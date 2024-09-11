/*
 * Copyright (c) 2024 Marzia Rivi
 *
 * This file is part of RadioLensfit2.
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

#ifndef DATA_PROCESSING_H_
#define DATA_PROCESSING_H_

/**
 * @file data_processing.h
*/

#ifdef __cplusplus
extern "C" {
#endif

#ifdef USE_MPI
void data_processing_mpi(bool re_fitting, unsigned int *bad_list, int nprocs, int rank, int nsources, double len, unsigned int num_coords,
                     FILE *pFile, likelihood_params *par, double *l, double *m, double *gflux, double *gscale, double *ge1, double *ge2, double *SNR_vis,
                     double *sum_w, complexd *visGal, complexd *visSkyMod, complexd *visData,
                     float *sigma2_vis, bool *flag, double *uu_metres, double *vv_metres, double *ww_metres, complexd *temp_facet_visData,
                     double *temp_facet_sigma2, double *temp_sum, double *com_time, double *fitting_time, int *bad);
#else
void data_processing(bool re_fitting, unsigned int *bad_list, int nsources, double len, unsigned int num_coords,
                     FILE *pFile, likelihood_params *par, double *l, double *m, double *gflux, double *gscale, double *ge1, double *ge2, double *SNR_vis,
                     double *sum_w, complexd *visGal, complexd *visSkyMod, complexd *visData,
                     float *sigma2_vis, bool *flag, double *uu_metres, double *vv_metres, double *ww_metres,double *fitting_time, int *bad);
#endif
 
#ifdef __cplusplus
}
#endif

#endif
