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
 

#ifndef ____galaxy_fitting__
#define ____galaxy_fitting__

#include "datatype.h"

#ifdef FACET
void source_extraction(double l0, double m0, double flux, double mu, double e1, double e2, likelihood_params *par, complexd *visSkyMod, complexd *visData, complexd *visGal, unsigned long int num_coords, double *uu_metres, double *vv_metres, int facet_size, double len);
#else
void source_extraction(double l0, double m0, double flux, double mu, double e1, double e2, likelihood_params *par, complexd *visSkyMod, complexd *visData, complexd *visGal, unsigned long int num_coords, double *uu_metres, double *vv_metres);
#endif

int source_fitting(int rank, likelihood_params *par, double *mes_e1, double *mes_e2, double *var_e1, double *var_e2, double *oneDimvar, double *maxL);

#endif /* defined(____galaxy_fitting__) */
