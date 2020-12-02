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
 

#ifndef ____source_extraction__
#define ____source_extraction__

#include "datatype.h"

#ifdef FACET
void source_extraction(int rank, unsigned int my_start_index, int facet, unsigned long int facet_ncoords, complexd *facet_vis, double *facet_sigma2, double l0, double m0, double flux, double mu, double e1, double e2, 
                       likelihood_params *par, complexd *visSkyMod, complexd *visData, complexd *visGal, double *sigma2_vis, unsigned int nchannels, 
                       unsigned long int num_coords, double *uu_metres, double *vv_metres, double *ww_metres, double len);
#else
void source_extraction(double l0, double m0, double flux, double mu, double e1, double e2, likelihood_params *par, complexd *visSkyMod, complexd *visData, complexd *visGal, 
                       unsigned long int num_coords, double *uu_metres, double *vv_metres, double *ww_metres);
#endif

#endif /* defined(____source_extraction__) */
