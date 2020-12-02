/*
 * Copyright (c) 2020 Marzia Rivi
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


#ifndef _likelihood_h
#define _likelihood_h

#include "datatype.h"

#ifdef __cplusplus
extern "C" {
#endif
 

//double f_posterior (const gsl_vector *v, void *params);
double f_likelihood (const gsl_vector *v, void *params);
double loglikelihood(void *params, double ee1, double ee2, int *error);
    
#ifdef FACET
double loglikelihood_r(unsigned int nchannels, double band_factor, double acc_time, double* spec,
                           double* wavenumbers, double ee1, double ee2, double l, double m, double radius,
                           double scale,
                           unsigned long int n_uv_coords, unsigned long int* count, const double *variance,
                           double* uu_metres, double* vv_metres, double* weights, complexd* visData, double* visM);
    
int cross_correlation(unsigned int nchannels, double* wavenumbers, unsigned long int n_uv_coords,
                           const double* variance, double* uu_metres, double* vv_metres, double* weights, complexd* visData,
                           double* visMod, double* ho, double* det_sigma);
#else
double loglikelihood_r(unsigned int nchannels, double band_factor, double acc_time, double* spec,
                       double* wavenumbers, double ee1, double ee2, double l, double m, double radius,
                       double scale, unsigned long int n_uv_coords, unsigned long int* count,
                       const double *variance, double* uu_metres, double* vv_metres, double* ww_metres, complexd* visData, complexd* visM);
    
int cross_correlation(unsigned int nchannels, double* wavenumbers, unsigned long int n_uv_coords,
                       const double* variance, double* uu_metres, double* vv_metres, complexd* visData,
                       complexd* visMod, double* ho, double* det_sigma);
#endif

double marginalise_over_position_shift(double x);
    
void likelihood_sampling(double *mes_e1, double *mes_e2, double maxL, void *params, int np_max,
                         double *var_e1, double *var_e2, double *oneDimvar);
    
#ifdef __cplusplus
}
#endif

#endif
