/*
 * Copyright (c) 2024 Marzia Rivi, Lance Miller
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
#include <sys/time.h>
#include <string.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>

#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif

long long current_timestamp()
{
    struct timeval te;
    gettimeofday(&te, NULL); // get current time
    long long milliseconds = te.tv_sec*1000LL + te.tv_usec/1000; // calculate milliseconds
 
    return milliseconds;
}

/*
   Generates a random seed, either from the unix
   /dev/random hardware stream or from the clock if the
   hardware stream is not available unsigned long int random_seed()
*/
unsigned long int random_seed()
{
    unsigned int seed = 1;
    struct timeval tv;
    FILE *devrandom;
    if ((devrandom = fopen("/dev/random","r")) == NULL)
    {
        // use the clock
        gettimeofday(&tv,0);
        seed = tv.tv_sec + tv.tv_usec;
        printf("Got seed %u from gettimeofday()\n",seed);
    }
    else
    {
        // use the hardware random stream
        fread(&seed,sizeof(seed),1,devrandom);
        printf("Got seed %u from /dev/random\n",seed);
        fclose(devrandom);
    }
    return seed;
}

/*
   This function returns a random number from a Gaussian distribution with
   zero mean and unit variance.
   The random number generator may be seeded by calling srand() prior to this function.

   param another, if not NULL, is used to return a second random number.
*/
double random_gaussian(double* another, gsl_rng * gen)
{
    double x, y, r2, fac;

    unsigned long int maxG = gsl_rng_max(gen);
    unsigned long int minG = gsl_rng_min(gen);

    do
    {
        /* Choose x and y in a uniform square (-1, -1) to (+1, +1). */
        x = 2.0 * gsl_rng_uniform(gen) - 1.0;
        y = 2.0 * gsl_rng_uniform(gen) - 1.0;

        /* Check if this is in the unit circle. */
        r2 = x*x + y*y;
    } while (r2 >= 1.0 || r2 == 0.0);

    /* Box-Muller transform (polar form). */
    fac = sqrt(-2.0 * log(r2) / r2);
    x *= fac;
    if (another) *another = y * fac;

    /* Return the first random number. */
    return x;
}

#ifdef __cplusplus
}
#endif
