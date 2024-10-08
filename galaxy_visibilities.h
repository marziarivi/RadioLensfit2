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
 

#ifndef _GalaxyVisibilities_h
#define _GalaxyVisibilities_h

#include <gsl/gsl_rng.h>
#include "datatype.h"

#ifdef __cplusplus
extern "C" {
#endif

double gaussian_shape(double k1, double k2, double scale_factor, double wave_factor);
double sersic_shape(double k1, double k2, double scale_factor, double wave_factor);    

void model_galaxy_visibilities_at_zero(double e1, double e2, double scale, unsigned int num_coords, 
                                       double* grid_u, double* grid_v,const double* sigma2, double* Modvis);

void model_galaxy_visibilities(unsigned int nchannels, double* spec, double* wavenumbers, double band_factor,
                               double acc_time, double e1, double e2, double scale, double l,
                               double m, unsigned int num_coords, double* uu_metres,
                               double* vv_metres, double* ww_metres, const double* sigma2, complexd* Modvis);
    
void data_galaxy_visibilities(double spectra, double wavenumber, double band_factor, double acc_time,
                              double e1, double e2, double scale, double flux, double l, double m,
                              unsigned int num_coords, double* uu_metres, double* vv_metres, double* ww_metres,
                              complexd* vis);

//double fq_smear(double band_factor, double phase);
//double t_smear(double acc_time, double phase);
//double primary_beam_profile(double ch_freq, double source_pos);

void add_system_noise(gsl_rng * gen, unsigned int num_coords, complexd* vis, double* sigma);

    
void data_visibilities_phase_shift(double wavenumber, double l, double m, unsigned int num_coords,
                                       double* uu_metres, double* vv_metres, double* ww_metres, complexd* vis);
   
#ifdef __cplusplus
}
#endif

#endif
