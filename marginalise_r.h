/*
 * Copyright (c) 2018 Lance Miller, Marzia Rivi
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


#ifndef _marginalise_r_h
#define _marginalise_r_h

//#ifdef __cplusplus
//extern "C" {
//#endif

void set_posterior_values(int numR, double* L_r, double* rprior, double* Ro, double* xmarvals, double* ymarvals, int* numvals);
double marginalise_posterior_r(int num_marvals, double *xmarvals, double *ymarvals);
    
double marf(double x, void *params);

//#ifdef __cplusplus
//}
//#endif

#endif
