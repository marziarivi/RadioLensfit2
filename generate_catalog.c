
/*
 * Copyright (c) 2018 Marzia Rivi
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_rng.h>

#include "random_gaussian.h"
#include "distributions.h"
#include "generate_random_values.h"
#include "generate_catalog.h"

// Generate galaxy catalog ordered by source flux
// 2NP = number of sampled orientations (points on the circle of radius |e|) for each ellipticity module

unsigned long int galaxy_catalog(unsigned long int nge, int NP, double fov_eff, double Rmin, double Rmax, double Fmin, double Fmax, double* gflux, double* gscale, double* ge1, double *ge2, double *l, double *m)
{
    //setup random number generator
    const gsl_rng_type * G;
    gsl_rng * gen;
    G = gsl_rng_mt19937;  // Mersenne Twister
    gen = gsl_rng_alloc(G);
    
    unsigned long int seed = random_seed();
    gsl_rng_set(gen,seed);
    
#ifdef USE_MPI
    unsigned long int my_gal = nge/nprocs;
    unsigned long int rem = nge%nprocs;
    if (rem)
    if (rank < rem) my_gal++;
#else
    unsigned long int my_gal = nge;
#endif
    unsigned long int diffgal = my_gal/2*NP;
    my_gal = diffgal*2*NP;
    
    // generate flux values
    double *gflux2 = new double[diffgal];
    generate_random_data(gen,diffgal,gflux2,Fmin,Fmax,flux_CDF,beta);
    
    // sort flux values, so that to generate a population ordered by flux and therefore fitting sources by decreasing flux order
    gsl_sort(gflux2,1,diffgal); // sorting ascending order
    unsigned long int gal = 0;
    for (unsigned long int i=0; i<diffgal; i++)
        for (int k = 0; k < 2*NP; k++)
        {
            gflux[gal] = gflux2[diffgal-i-1];
            gal++;
        }
    delete[] gflux2;
    
    // generate scalelength
    double mu, scalelength;
    for (unsigned long int g = 0; g < my_gal; g++)
    {
        if (g%(2*NP) == 0)
        {
            mu = scale_mean(gflux[g]); //power law relation between flux and scalelength
            generate_random_data(gen,1,&scalelength,Rmin,Rmax,r_CDF,mu);
        }
        gscale[g] = scalelength;
    }
    
    // generate ellipticities
    generate_ellipticity(gen,diffgal,NP,ge1,ge2);
    
    // generate positions
    // uniformly random positions in RAD in a disk area
    // http://mathworld.wolfram.com/DiskPointPicking.html
    double radius,orient;
    
    for (unsigned long int gal=0; gal<my_gal; gal++)
    {
        radius = sqrt(gsl_rng_uniform(gen))*0.5*fov_eff;
        orient = gsl_rng_uniform(gen)*2*PI;
        
        l[gal] = radius*cos(orient);
        m[gal] = radius*sin(orient);
    }
    
    gsl_rng_free(gen);
    return my_gal;
}
