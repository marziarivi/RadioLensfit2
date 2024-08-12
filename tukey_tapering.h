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

#ifndef ____tukey_tapering__
#define ____tukey_tapering__


#ifdef __cplusplus
extern "C" {
#endif

double tukey_weight(double x, double n);
double tukey_tapering(double uvdist, double minuv_l, double maxuv_l, double inner_tukey, double outer_tukey);
    
#ifdef __cplusplus
}
#endif

#endif /* defined(____tukey_tapering__) */
