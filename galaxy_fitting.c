/*
 * Copyright (c) 2021 Marzia Rivi
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
#include "likelihood.h"
#include "marginalise_r.h"
#include "galaxy_fitting.h"


#ifdef __cplusplus
extern "C" {
#endif

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
   
    //printf("rank %d:  Maximum log likelihood = %f for e = %f,%f \n",rank,*maxL,*mes_e1,*mes_e2);
    if ((*mes_e1)*(*mes_e1)+(*mes_e2)*(*mes_e2) > 0.65) *maxL = -1e+10;

    // Likelihood sampling to compute mean and variance
    int error = 0;
    *var_e1 = 0.;
    *var_e2 = 0.;
    *oneDimvar = 0.;
    if (*maxL > -1.e+10)
    {
        likelihood_sampling(mes_e1, mes_e2, *maxL, par, np_max, var_e1, var_e2, oneDimvar);
        printf("rank %d: average: %f,%f variance: %e,%e 1Dvar: %e\n",rank,*mes_e1,*mes_e2,*var_e1,*var_e2,*oneDimvar);
        if (*var_e1 < VAR || *var_e2 < VAR || *oneDimvar < VAR)
        {
           printf("rank %d: ERROR likelihood sampling!\n",rank);
           error = 1;
        }
    }
    else 
    {
       printf("rank %d: BAD - NO likelihood sampling!\n",rank);
       error = 1;
    }
    
    gsl_vector_free(x);
    gsl_vector_free(sx);
    gsl_multimin_fminimizer_free(s);
    
    return error;
}

#ifdef __cplusplus
}
#endif

