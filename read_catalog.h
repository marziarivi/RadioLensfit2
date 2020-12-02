
/*
 * Copyright (c) 2020 Marzia Rivi
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

#ifndef ____read_catalog__
#define ____read_catalog__

#ifdef __cplusplus
extern "C" {
#endif

int read_catalog(unsigned long int nge, char *filename, double *gflux, double *gscale, double *ge1, double *ge2, double *l, double *m, double *SNR_vis);

#ifdef __cplusplus
}
#endif

#endif /* defined(____read_catalog__) */
