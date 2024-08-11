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

#ifndef _evaluate_uv_grid_h
#define _evaluate_uv_grid_h

#ifdef __cplusplus
extern "C" {
#endif
    
#include "datatype.h"
    
unsigned int facet_size(double theta_med, double len);
unsigned int evaluate_uv_grid_size(int rank, int nprocs, double len, double *wavenumbers, unsigned int num_channels, unsigned  int ncoords, double* u, double* v, unsigned int sizeg, bool* flag);
unsigned int evaluate_facet_coords(double* grid_u, double* grid_v, double len, unsigned int sizeg, double *count_w);
void gridding_visibilities(double *wavenumbers, unsigned int num_channels, unsigned int ncoords, double *u, double *v, complexd *vis, float *sigma2, double len, unsigned int sizeg, complexd *new_vis, double *new_sigma2, bool* flag, double* sum_w);
void average_facets(unsigned long int size, complexd* grid_vis, double* grid_sigma2, double *sum_w);
    
#ifdef __cplusplus
}
#endif

#endif
