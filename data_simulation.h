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

#ifndef ____data_simulation__
#define ____data_simulation__

 
#include "datatype.h"

void data_simulation(double *wavenumbers, double *spec, double channel_bandwidth_hz, int time_acc, unsigned int num_channels,
                     unsigned int num_baselines, double sigma, unsigned long int n_gal, double g1, double g2,
                     double *ge1, double *ge2, double *gflux, double *gscale, double *l, double *m,
                     double *SNR_vis, unsigned long int num_coords, double *uu_metres, double *vv_metres,
                     double *ww_metres, complexd *visGal, complexd* visData);

void sky_model(int rank, double *wavenumbers, double *spec,
                     double channel_bandwidth_hz, int time_acc, unsigned int num_channels, unsigned int num_baselines,
                     unsigned long int n_gal, double *gflux, double *gscale, double *l, double *m, unsigned long int num_coords,
                     double *uu_metres, double *vv_metres, double *ww_metres, complexd *visGal, complexd* visMod);



#endif /* defined(____data_simulation__) */
