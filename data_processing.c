/*
 * Copyright (c) 2024 Marzia Rivi
 *
 * This file is part of RadioLensfit.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, see <http://www.gnu.org/licenses/>.
 */


//  data_processing.c

#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "datatype.h"
#include "default_params.h"
#include "utils.h"
#include "galaxy_fitting.h"
#include "galaxy_visibilities.h"
#include "distributions.h"
#include "evaluate_uv_grid.h"
#include "source_extraction.h"
#include "data_processing.h"


#ifdef __cplusplus
extern "C" {
#endif

#ifdef FACET
void data_processing(bool re_fitting, unsigned int *bad_list, int nprocs, int rank, int nsources, double len, unsigned int num_coords,
                     FILE *pFile, likelihood_params *par, double *l, double *m, double *gflux, double *gscale, double *ge1, double *ge2, double *SNR_vis,
                     double *sum_w, complexd *visGal, complexd *visSkyMod, complexd *visData,
                     float *sigma2_vis, bool *flag, double *uu_metres, double *vv_metres, double *ww_metres, complexd *temp_facet_visData,
                     double *temp_facet_sigma2, double *temp_sum, double *com_time, double *fitting_time, int *bad)
#else
void data_processing(bool re_fitting, unsigned int *bad_list, int nprocs, int rank, int nsources, double len, unsigned int num_coords,
                     FILE *pFile, likelihood_params *par, double *l, double *m, double *gflux, double *gscale, double *ge1, double *ge2, double *SNR_vis,
                     double *sum_w, complexd *visGal, complexd *visSkyMod, complexd *visData,
                     float *sigma2_vis, bool *flag, double *uu_metres, double *vv_metres, double *ww_metres, double *fitting_time, int *bad)
#endif
{
#ifdef USE_MPI
    MPI_Status stat;
#endif
    unsigned int gal;
    double l0,m0;
    int k, nbad;
    long int my_g, ind; 
    double flux, mu, R_mu[nprocs];

    nbad = 0;
    unsigned int g = 0; // source global index
    while (g < nsources)
    { 
      k = 0; // source local index
      my_g = -1;  // if sources are finished, current task has my_g = -1 and has no source to fit

      // source extraction of nproc consecutive sources that will be fitted each one by a different task -----------------------
      for(int src = 0; src < nprocs && g < nsources; src++)
      {
        if (re_fitting) gal = bad_list[g];
        else gal = g;

        l0 = l[gal];  m0 = m[gal];
        flux = gflux[gal];
        mu = scale_mean(flux);
#ifndef SCALELENGTH_ON 
        R_mu[src] = exp(mu);
#else
        R_mu[src] = gscale[gal];
#endif
#ifdef FACET
        unsigned int facet = facet_size(R_mu[src],len);
        unsigned long int size = (unsigned long int) facet*facet;
#endif
        if (rank == src)  // proc src will fit the current source:
        {
          my_g = gal;  
          for (int nRo=1; nRo<par->numr; nRo++)   
            (par->rprior)[nRo] = rfunc(mu,R_STD,par->ro[nRo]);  // set log(prior) for scalelength of source gal
#ifdef FACET
           // extract visibilities, sigma2 and weights from my MS (already summed in the facet) for source gal  
           par->facet = facet;
           source_extraction(rank,facet,par,par->data,par->sigma2,sum_w,l0, m0, flux, R_mu[src], 0., 0., visSkyMod, visData, visGal, sigma2_vis, flag, num_coords, uu_metres, vv_metres, ww_metres, len);
#ifdef USE_MPI
           int n = 1;
           while (n<nprocs)
           {
             // collect and reduce facet vis, sigma2 and sum_w of the current source from the other procs (for their MS contribution)
             *com_time -= MPI_Wtime();
             MPI_Recv(temp_sum,size,MPI_DOUBLE,MPI_ANY_SOURCE,src,MPI_COMM_WORLD,&stat);
             MPI_Recv(temp_facet_sigma2,size,MPI_DOUBLE,MPI_ANY_SOURCE,nprocs+src,MPI_COMM_WORLD,&stat);
             MPI_Recv(temp_facet_visData,2*size,MPI_DOUBLE,MPI_ANY_SOURCE,2*nprocs+src,MPI_COMM_WORLD,&stat);
             *com_time += MPI_Wtime();

             for (unsigned long int i = 0; i<size; i++) sum_w[i] += temp_sum[i];
             for (unsigned long int i = 0; i<size; i++) (par->sigma2)[i] += temp_facet_sigma2[i];
             for (unsigned long int i = 0; i<size; i++) 
             {
                (par->data)[i].real += temp_facet_visData[i].real;
                (par->data)[i].imag += temp_facet_visData[i].imag;
             }
             n++;
           }
        }
        else 
        {
           // extract visibilities, sigma2 and sum_w from my MS (already summed in the facet) for source gal and send them to proc src
           source_extraction(rank,facet,par,temp_facet_visData, temp_facet_sigma2,temp_sum,l0, m0, flux, R_mu[src], 0., 0., visSkyMod, visData, visGal, sigma2_vis, flag, num_coords, uu_metres, vv_metres, ww_metres, len);

           *com_time -= MPI_Wtime();
           MPI_Send(temp_sum,size,MPI_DOUBLE,src,src,MPI_COMM_WORLD);
           MPI_Send(temp_facet_sigma2,size,MPI_DOUBLE,src,nprocs+src,MPI_COMM_WORLD);
           MPI_Send(temp_facet_visData,2*size,MPI_DOUBLE,src,2*nprocs+src,MPI_COMM_WORLD);
           *com_time += MPI_Wtime();
#endif
#else
           // extract source visibilities without faceting
           source_extraction(l0, m0, flux, R_mu[k], 0., 0., par, visSkyMod, visData, visGal, sigma2_vis, num_coords, uu_metres, vv_metres, ww_metres);
#endif
        }
        g++;
        k++;
      }
      
#ifdef FACET
      if (my_g >= 0)
      {
#ifdef USE_MPI
         // average facet visibilities (already summed within the cell)
         average_facets(par->facet*par->facet, par->data, par->sigma2, sum_w);
#endif
         // compute facet coordinates their number
         par->ncoords = evaluate_facet_coords(par->uu, par->vv, len, par->facet, sum_w);
      }
#endif

      // source fitting of this task --------------------------------------------------------------------------------------------------
#ifdef USE_MPI
      double start_fitting = MPI_Wtime();
#else
      long long start_fitting = current_timestamp();
#endif
      double mes_e1, mes_e2, maxL;
      double var_e1, var_e2, oneDimvar;
      if (my_g >= 0) 
      {
        source_fitting(rank, par, &mes_e1, &mes_e2, &var_e1, &var_e2, &oneDimvar, &maxL);
        printf("rank %d: n. %d flux = %f: measured e = %f , %f \n",rank,my_g,gflux[my_g],mes_e1,mes_e2);
      }
#ifdef USE_MPI
      *fitting_time += MPI_Wtime() - start_fitting;
#else
      *fitting_time += (double) (current_timestamp() -start_fitting)/1000.;
#endif

      // removal of the fitted sources visibilities from the MS data and sky model ----------------------------------------------------------
      double res[6]; // shape fitting results to be sent to the other procs
      ind = g - k;  // source global index
      k = 0;  // source local index
      for (int src = 0 ; src < nprocs && ind < nsources; src++)
      {
         if (re_fitting) gal = bad_list[ind];  //source global index
         else gal = ind;
         l0 = l[gal];  m0 = m[gal];
         flux = gflux[gal];

         if (rank == k)
         {
            res[0] = mes_e1; res[1] = mes_e2;
            res[2] = var_e1; res[3] = var_e2;
            res[4] = oneDimvar; res[5] = maxL;
         }
#ifdef USE_MPI 
         //MPI_Barrier(MPI_COMM_WORLD);
         // Bcast shape results from proc k to the others
         *com_time -= MPI_Wtime();
         MPI_Bcast(res,6,MPI_DOUBLE,k,MPI_COMM_WORLD);          
         *com_time += MPI_Wtime();
#endif
         if (res[5] <= -1e+10 || res[2] < VAR || res[3] < VAR || res[4] < VAR)  // bad measurement
         {
            if (re_fitting) 
            {
              if (rank == 0) fprintf(pFile, "%f | %f | %f | %e | %f | %f | %e | %e | %f | %f | %f \n",flux,ge1[gal],res[0],sqrt(res[2]),ge2[gal],res[1],sqrt(res[3]),res[4],SNR_vis[gal],l0/(ARCS2RAD),m0/(ARCS2RAD));   
              res[0] = 0.; res[1] = 0.; // remove round source model from original data
            } 
            else bad_list[nbad] = ind;  // each proc store the index of bad measurements to be fit again at the end
            nbad++;
         }
         else     
         {
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (unsigned int ch = 0; ch < par->nchannels; ch++)
            {
               unsigned long int ch_vis = (unsigned long int) ch*num_coords;
         
               // remove current round source model from the sky model
               data_galaxy_visibilities((par->spec)[ch], (par->wavenumbers)[ch], par->band_factor, par->acc_time, 0., 0., R_mu[k],
                                           flux, l0, m0, num_coords, uu_metres, vv_metres, ww_metres, &(visGal[ch_vis]));
                  
               for (unsigned long int i = ch_vis; i<ch_vis+num_coords; i++)
               {
                  visSkyMod[i].real -= visGal[i].real;
                  visSkyMod[i].imag -= visGal[i].imag;
               }
            
               // remove current source model fit from original data
               data_galaxy_visibilities((par->spec)[ch], (par->wavenumbers)[ch], par->band_factor, par->acc_time, res[0], res[1], R_mu[k],
                                           flux, l0, m0, num_coords, uu_metres, vv_metres, ww_metres, &(visGal[ch_vis]));
                  
               for (unsigned long int i = ch_vis; i<ch_vis+num_coords; i++)
               {
                  visData[i].real -= visGal[i].real;
                  visData[i].imag -= visGal[i].imag;
               }
            }
            if (rank == 0) fprintf(pFile, "%f | %f | %f | %e | %f | %f | %e | %e | %f | %f | %f \n",flux,ge1[gal],res[0],sqrt(res[2]),ge2[gal],res[1],sqrt(res[3]),res[4],SNR_vis[gal],l0/(ARCS2RAD),m0/(ARCS2RAD));
         }
         k++;
         ind++;
      }
    }

    *bad = nbad;
    return;
}


#ifdef __cplusplus
}
#endif

