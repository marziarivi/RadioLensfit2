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

#ifndef ____scalelength_dist__
#define ____scalelength_dist__

#ifdef __cplusplus
extern "C" {
#endif

const double scale_std = 0.3136;  // scalelength standard deviation
const double beta = -1.34; // flux prior: S^beta

double e_pdf(double e);
double flux_CDF(double alpha, double flux);
double scale_mean(double flux);
double r_CDF(double mean, double r);
double rfunc (const double mu, const double sigma, double r);
 
    
#ifdef __cplusplus
}
#endif

#endif /* defined(____scalelength_dist__) */
