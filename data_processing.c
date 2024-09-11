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

void data_processing(bool re_fitting, unsigned int *bad_list, int nsources, double len, unsigned int num_coords,
                     FILE *pFile, likelihood_params *par, double *l, double *m, double *gflux, double *gscale, double *ge1, double *ge2, double *SNR_vis,
                     double *sum_w, complexd *visGal, complexd *visSkyMod, complexd *visData,
                     float *sigma2_vis, bool *flag, double *uu_metres, double *vv_metres, double *ww_metres, double *fitting_time, int *bad)
{
    unsigned int gal;
    double mu, R_mu;
    double mes_e1, mes_e2, maxL;
    double var_e1, var_e2, oneDimvar; 
    double l0,m0,flux;
    int nbad = 0;
  
    for (unsigned int g = 0; g<nsources; g++)
    { 
        if (re_fitting) gal = bad_list[g];
        else gal = g;

        l0 = l[gal];  m0 = m[gal];
        flux = gflux[gal];
        mu = scale_mean(flux);
#ifndef SCALELENGTH_ON 
        R_mu = exp(mu);
#else
        R_mu = gscale[gal];
#endif
        // set log(prior) for scalelength of source gal
        for (int nRo=1; nRo<par->numr; nRo++)   
            (par->rprior)[nRo] = rfunc(mu,R_STD,par->ro[nRo]);   

#ifdef FACET
        unsigned int facet = facet_size(R_mu,len);
        // extract visibilities, sigma2 and weights from my MS (already summed in the facet) for source gal  
        par->facet = facet;
        source_extraction(0,facet,par,par->data,par->sigma2,sum_w,l0, m0, flux, R_mu, 0., 0., visSkyMod, visData, visGal, sigma2_vis, flag, num_coords, uu_metres, vv_metres, ww_metres, len);
        // compute facet coordinates their number
        par->ncoords = evaluate_facet_coords(par->uu, par->vv, len, par->facet, sum_w);
#else
        // extract source visibilities without faceting
        source_extraction(l0, m0, flux, R_mu, 0., 0., par, visSkyMod, visData, visGal, sigma2_vis, num_coords, uu_metres, vv_metres, ww_metres);
#endif

        // source fitting of this source --------------------------------------------------------------------------------------------------
        long long start_fitting = current_timestamp();
        int error = source_fitting(0, par, &mes_e1, &mes_e2, &var_e1, &var_e2, &oneDimvar, &maxL);
        printf("n. %d flux = %f: measured e = %f , %f \n",gal,flux,mes_e1,mes_e2); fflush(stdout);
        *fitting_time += (double) (current_timestamp() -start_fitting)/1000.;

        // removal of the fitted sources visibilities from the MS data and sky model ----------------------------------------------------------
        if (error) // bad measurement
        {  
            if (re_fitting) 
            {
              mes_e1 = 0.; mes_e2 = 0.;
              fprintf(pFile, "%f | %f | %f | %f | %f | %f | %f | %f | %f | %f | %f \n",flux,ge1[gal],mes_e1,sqrt(var_e1),ge2[gal],mes_e2,sqrt(var_e2),oneDimvar,SNR_vis[gal],l0/(ARCS2RAD),m0/(ARCS2RAD));
            } 
            else bad_list[nbad] = gal;  // store the index of bad measurements to be fit again at the end
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
               data_galaxy_visibilities((par->spec)[ch], (par->wavenumbers)[ch], par->band_factor, par->acc_time, 0., 0., R_mu,
                                           flux, l0, m0, num_coords, uu_metres, vv_metres, ww_metres, &(visGal[ch_vis]));
                  
               for (unsigned long int i = ch_vis; i<ch_vis+num_coords; i++)
               {
                  visSkyMod[i].real -= visGal[i].real;
                  visSkyMod[i].imag -= visGal[i].imag;
               }
            
               // remove current source model fit from original data
               data_galaxy_visibilities((par->spec)[ch], (par->wavenumbers)[ch], par->band_factor, par->acc_time, mes_e1, mes_e2, R_mu,
                                           flux, l0, m0, num_coords, uu_metres, vv_metres, ww_metres, &(visGal[ch_vis]));
                  
               for (unsigned long int i = ch_vis; i<ch_vis+num_coords; i++)
               {
                  visData[i].real -= visGal[i].real;
                  visData[i].imag -= visGal[i].imag;
               }
            }
            fprintf(pFile, "%f | %f | %f | %f | %f | %f | %f | %f | %f | %f | %f \n",flux,ge1[gal],mes_e1,sqrt(var_e1),ge2[gal],mes_e2,sqrt(var_e2),oneDimvar,SNR_vis[gal],l0/(ARCS2RAD),m0/(ARCS2RAD)); 
      }
    }

    *bad = nbad;
    return;
}

#ifdef __cplusplus
}
#endif

