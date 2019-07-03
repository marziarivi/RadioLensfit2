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

#ifndef _generate_values_h
#define _generate_values_h

#include <gsl/gsl_rng.h>

#ifdef __cplusplus
extern "C" {
#endif

void generate_random_data(gsl_rng * gen, int nr, double *data, double min_value, double max_value, double (*CDFunc)(double,double), double param);
void generate_ellipticity(gsl_rng * gen, int ne, int NP, double *e1, double *e2);

#ifdef __cplusplus
}
#endif

#endif
