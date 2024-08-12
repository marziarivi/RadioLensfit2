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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "datatype.h"
#include "tukey_tapering.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Tukey tapering function 
 * It returns a value 0 <= w <= 1, in particular 0 when x = 0 and 1 when x=n) 
 * https://wsclean.readthedocs.io/en/latest/tapering.html#tukey-taper
*/
double tukey_weight(double x, double n)
{
  return 0.5*(1.0+cos((PI/n)*(x-n)));
}

double tukey_tapering(double uvdist, double minuv_l, double maxuv_l, double inner_tukey, double outer_tukey)
{
  double x,n;
  double w = 0;

  if (uvdist <= minuv_l) w = 0.;
  else if (uvdist >= maxuv_l) w = 0.;
  else if ((minuv_l + inner_tukey) < uvdist && uvdist < (maxuv_l - outer_tukey)) w= 1.;
  else if (minuv_l < uvdist && uvdist <= (minuv_l + inner_tukey)) 
    {
       x = uvdist - minuv_l;
       n = inner_tukey;
       w = tukey_weight(x,n);
    }
  else if ((maxuv_l - outer_tukey) <= uvdist && uvdist < maxuv_l) 
    {
       x = maxuv_l - uvdist;
       n = outer_tukey;
       w = tukey_weight(x,n);
    }
    return w;
} 


#ifdef __cplusplus
}
#endif
