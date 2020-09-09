/*
 * Copyright (c) 2018 Marzia Rivi
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

#include <new>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_rng.h>

#include "default_params.h"
#include "galaxy_visibilities.h"
#include "likelihood.h"
#include "marginalise_r.h"
#include "evaluate_uv_grid.h"
#include "galaxy_fitting.h"


// Source extraction --------------------------------------------------------
// Notice: instructions execution order is important as arrays are reused!!!!

#ifdef FACET
void source_extraction(int rank, unsigned int my_start_index, complexd *facet_vis, double *facet_sigma2, double l0, double m0, double flux, double mu, double e1, double e2, 
                       likelihood_params *par, complexd *visSkyMod, complexd *visData, complexd *visGal, double *sigma2_vis, unsigned int nchannels, 
                       unsigned long int num_coords, double *uu_metres, double *vv_metres, double *ww_metres, double len)
#else
void source_extraction(double l0, double m0, double flux, double mu, double e1, double e2, likelihood_params *par, complexd *visSkyMod, complexd *visData, complexd *visGal, 
                       unsigned long int num_coords, double *uu_metres, double *vv_metres, double *ww_metres)
#endif
{
   int facet = facet_size(mu,len);
   par->ncoords = evaluate_uv_grid(num_coords, par->uu, par->vv, uu_metres, vv_metres, len, facet, par->count);

#ifdef USE_MPI
   unsigned int my_freq_index = rank*nchannels; 
#else
   unsigned int my_freq_index = 0;
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
   for (unsigned int ch = 0; ch < nchannels; ch++)
   {
     unsigned long int ch_vis = ch*num_coords;
     unsigned int ch_ind = my_freq_index + ch;
     // get small round source at the current position and flux
     data_galaxy_visibilities((par->spec)[ch_ind], (par->wavenumbers)[ch_ind], par->band_factor, par->acc_time, e1, e2, mu,
                              flux, l0, m0, num_coords, uu_metres, vv_metres, ww_metres, &(visGal[ch_vis]));
    
     for (unsigned long int i = ch_vis; i<ch_vis+num_coords; i++)
     {
#ifdef USE_MPI
        // extract source data visibilities without modifying the sky model
        visGal[i].real = visData[i].real - visSkyMod[i].real + visGal[i].real;
        visGal[i].imag = visData[i].imag - visSkyMod[i].imag + visGal[i].imag;
#else
        // remove source model from the sky model
        visSkyMod[i].real -= visGal[i].real;
        visSkyMod[i].imag -= visGal[i].imag;
        
        // remove sky model from the original data (i.e. all other sources approximation) in order to have only the current galaxy data
        visGal[i].real = visData[i].real - visSkyMod[i].real;
        visGal[i].imag = visData[i].imag - visSkyMod[i].imag;
#endif 
     }
    
#ifdef FACET
     // Phase shift data visibilities (to be done after gridding because real data will be gridded)
     data_visibilities_phase_shift((par->wavenumbers)[ch_ind], l0, m0, num_coords, uu_metres, vv_metres, ww_metres, &(visGal[ch_vis]));
       
     // gridding visibilities
     unsigned int ch_visfacet = (my_start_index + ch)*par->ncoords;
     gridding_visibilities(num_coords,uu_metres,vv_metres,&(visGal[ch_vis]),&(sigma2_vis[ch_vis]),len,facet,&(facet_vis[ch_visfacet]),&(facet_sigma2[ch_visfacet]),par->count);
#else
     par->l0 = l0;
     par->m0 = m0;
#endif
   }
}


// Model fitting  -----------------------------------------------------
// Search for the maximum posterior to find starting ellipticity points

int source_fitting(int rank, likelihood_params *par, double *mes_e1, double *mes_e2, double *var_e1, double *var_e2, double *oneDimvar, double *maxL)
{
    int np_max = NP_MAX;  // min number of sampling points with likelihood above 5%ML
    
    gsl_multimin_function minex_func;
    minex_func.n = 2;
    minex_func.f = f_likelihood;
    minex_func.params = par;
    
    // use Simplex algorithm of Nelder and Mead provided by the GLS library to minimize -log(likelihood)
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    
    gsl_multimin_fminimizer *s = 0;
    gsl_vector *sx, *x;
    x = gsl_vector_alloc (2);
    sx = gsl_vector_alloc (2);
    s = gsl_multimin_fminimizer_alloc (T, 2);
    
    double start_e1 = 0.;
    double start_e2 = 0.;
    
    // Search for the maximum likelihood
    gsl_vector_set (x, 0, start_e1);
    gsl_vector_set (x, 1, start_e2);
    gsl_vector_set_all (sx, 0.1);
    
    gsl_multimin_fminimizer_set (s, &minex_func, x, sx);
    int iter = 0;
    int status;
    double size;
    
    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        
        if (status) break;
        
        size = gsl_multimin_fminimizer_size(s);
        status = gsl_multimin_test_size (size, TOL);
    }
    while (status == GSL_CONTINUE && iter < 50 && s->fval < 0.);
    
    *mes_e1 = gsl_vector_get(s->x, 0);
    *mes_e2 = gsl_vector_get(s->x, 1);
    *maxL= -s->fval;
    
    printf("rank %d:  Maximum log likelihood = %f for e = %f,%f\n",rank,*maxL,*mes_e1,*mes_e2);
    
    // Likelihood sampling to compute mean and variance
    int error = 0;
    *var_e1 = 0.;
    *var_e2 = 0.;
    *oneDimvar = 0.;
    if (*maxL > -1e+10)
    {
        likelihood_sampling(mes_e1, mes_e2, *maxL, par, np_max, var_e1, var_e2, oneDimvar);
        printf("rank %d: average: %f,%f variance: %e,%e 1Dvar: %e\n",rank,*mes_e1,*mes_e2,*var_e1,*var_e2,*oneDimvar);
        if (*var_e1 < 1e-4 || *var_e2 < 1e-4 || *oneDimvar < 1e-4)
        {
           printf("rank %d: ERROR likelihood sampling!\n",rank);
           error = 1;
        }
    }
    else error = 1;
    
    gsl_vector_free(x);
    gsl_vector_free(sx);
    gsl_multimin_fminimizer_free(s);
    
    return error;

}
