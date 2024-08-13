
/*
 * Copyright (c) 2024 Marzia Rivi
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

#include "datatype.h"
#include "utils.h"
#include "default_params.h"
#include "distributions.h"
#include "generate_random_values.h"

// Generate galaxy catalog ordered by source flux
// 2NP = number of sampled orientations (points on the circle of radius |e|) for each ellipticity module

//    Command line input parameters:
//    argv[1]  number of sources 
//    argv[2]  FoV eff (arcmin)
//    argv[3]  Fmin (uJy)

using namespace std;

int main(int argc, char *argv[])
//unsigned int galaxy_catalog(unsigned int nge, int NP, double fov_eff, double Rmin, double Rmax, double Fmin, double Fmax, double* gflux, double* gscale, double* ge1, double *ge2, double *l, double *m)
{
    unsigned int nge = atoi(argv[1]);
    double fov_eff = atof(argv[2])*60.*ARCS2RAD; //1.22*C0/(freq_start_hz*diameter);  // field of view in RAD
    //unsigned long int nge = ceil((flux_CDF(M_EXP, FMAX) - flux_CDF(M_EXP, Fmin))*fov_eff_arcmin*fov_eff_arcmin);
   
    //setup random number generator
    const gsl_rng_type * G;
    gsl_rng * gen;
    G = gsl_rng_mt19937;  // Mersenne Twister
    gen = gsl_rng_alloc(G);
    
    unsigned long int seed = random_seed();
    gsl_rng_set(gen,seed);

    int NP = NUM_ORIENT;
    unsigned int my_gal = nge;
    unsigned int diffgal = my_gal/(2*NP);
    my_gal = diffgal*2*NP;
    
    // generate flux values
    double *gflux = new double[nge];
    double *gflux2 = new double[diffgal];
    double Fmin = atof(argv[3]);
    generate_random_data(gen,diffgal,gflux2,Fmin,FMAX,flux_CDF,M_EXP);
    
    // sort flux values, so that to generate a population ordered by flux and therefore fitting sources by decreasing flux order
    gsl_sort(gflux2,1,diffgal); // sorting ascending order
    unsigned int gal = 0;
    for (unsigned int i=0; i<diffgal; i++)
        for (int k = 0; k < 2*NP; k++)
        {
            gflux[gal] = gflux2[diffgal-i-1];
            gal++;
        }
    delete[] gflux2;
    
    // generate scalelength
    double *gscale = new double[nge];
    double mu, scalelength;
    for (unsigned int g = 0; g < my_gal; g++)
    {
        if (g%(2*NP) == 0)
        {
            mu = scale_mean(gflux[g]); //power law relation between flux and scalelength
            generate_random_data(gen,1,&scalelength,RMIN,RMAX,r_CDF,mu);
        }
        gscale[g] = scalelength;
    }
    
    // generate ellipticities
    double *ge1 = new double[nge];
    double *ge2 = new double[nge];
    generate_ellipticity(gen,diffgal,NP,ge1,ge2);

    FILE *pFile = 0;
    char output[100];
    sprintf(output,"Catalog_%d.txt",nge);
    pFile = fopen(output,"w");
    fprintf(pFile, "l | m | flux | scale | e1 |  e2  \n");
    
    // generate positions
    // uniformly random positions in RAD in a disk area
    // http://mathworld.wolfram.com/DiskPointPicking.html
    double radius,orient,l,m;
    for (unsigned int gal=0; gal<my_gal; gal++)
    {
        radius = sqrt(gsl_rng_uniform(gen))*0.5*fov_eff;
        orient = gsl_rng_uniform(gen)*2*PI;
        
        l = radius*cos(orient);
        m = radius*sin(orient);
        fprintf(pFile, "%f | %f | %f  | %f | %f | %e | %e | %f | %f | %f \n",l,m,gflux[gal],gscale[gal],ge1[gal],ge2[gal]);
    }

    fclose(pFile);
    gsl_rng_free(gen);

    delete[] gflux;
    delete[] gscale;
    delete[] ge1;
    delete[] ge2;
    
    return nge;
}
