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

#ifdef _OPENMP
#include <omp.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "evaluate_uv_grid.h"

#ifdef __cplusplus
extern "C" {
#endif
    
// gaussian function with variances given by the ellipticity orthogonal to the uv coverage's one
 
double weight_func(double u, double v)
{
    double sigma1 = 1.;  
    double sigma2 = 1.;  
    double w = exp(-u*u/(2*sigma1) - v*v/(2*sigma2));
    return w;
}

// Compute grid of u,v coordinates (coordinates are put in the center of the cell).
unsigned long int evaluate_uv_grid(double len, unsigned long int ncoords, double* u, double* v, int sizeg, double** grid_u, double** grid_v, unsigned long int* count)
{
    unsigned int i,j;
    unsigned long int p,n;
    unsigned long int size = sizeg*sizeg;
    memset(count, 0, size*sizeof(unsigned long int));
    
    //double* weights = (double *) malloc(size*sizeof(double));
    //memset(weights, 0, size*sizeof(double));
    
    double* temp_grid = (double *) malloc(sizeg*sizeof(double));
    
    double inc = 2*len/sizeg;
    for (i = 0; i < sizeg; ++i) temp_grid[i] = (i+0.5)*inc - len;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (unsigned long int k = 0; k < ncoords; k++)
    {
        unsigned int pu = (unsigned int) ((u[k] + len) / inc);
        unsigned int pv = (unsigned int) ((v[k] + len) / inc);
        unsigned long int pc = (unsigned long int) pv * sizeg + pu;
#ifdef _OPENMP
#pragma omp critical
#endif
      {
        //weights[pc] += weight_func(u[k],v[k]);
        count[pc]++;
      }
    }
 
    n = 0;
  //  for (p = 0; p < size; p++)  if (weights[p]) n++;
    for (p = 0; p < size; p++)  if (count[p]) n++;
    
    if (!(*grid_u))
    {
        *grid_u = (double *) malloc(n*sizeof(double));
        *grid_v = (double *) malloc(n*sizeof(double));
    }
    
    n=0;
    for (p = 0; p < size; p++)
    {
        //if (weights[p])
        if (count[p])
        {
            j = p/sizeg;
            i = p%sizeg;
            (*grid_u)[n] = temp_grid[i];
            (*grid_v)[n] = temp_grid[j];
            //count[n] = weights[p];
            count[n] = count[p];
            n++;
        }
    }
 
    free(temp_grid);
    //free(weights);
    
    return n;
}

    
/*
 Compute visibilities at the pg = (ug,vg) grid points
 by adding all the original visibilities at the (u,v)
 points falling in the box centered in pg
 This is a gridding by convolution with the pillbox function 
 (Synthesis Imaging in Radio Astronomy II, p.143) with uniform weighting.
 */
void gridding_visibilities(unsigned long int ncoords, double *u, double *v, complexd *vis, double len, int sizeg, complexd *new_vis, unsigned long int *count)
{
    unsigned long int p,c;
    double  inc = 2*len/sizeg;
    unsigned long int size = sizeg*sizeg;
    
    complexd* temp_grid_vis = (complexd *) calloc(size,sizeof(complexd));
    
    for (unsigned long int k = 0; k < ncoords; k++)
    {
        unsigned int pu = (unsigned int) ((u[k] + len) / inc);
        unsigned int pv = (unsigned int) ((v[k] + len) / inc);
        unsigned long int pc = (unsigned long int) pv * sizeg + pu;
        
        double wp = 1.; //weight_func(u[k],v[k]);
 
        temp_grid_vis[pc].real += wp*vis[k].real;
        temp_grid_vis[pc].imag += wp*vis[k].imag;
    }
 
    c=0;
    for (p = 0; p < size; p++)
    {
        if (temp_grid_vis[p].real)
        {
            new_vis[c].real = temp_grid_vis[p].real/count[c];
            new_vis[c].imag = temp_grid_vis[p].imag/count[c];

            c++;
        }
    }
    free(temp_grid_vis);
}

    
    
    
/*
 Uniform gridding by convolution with the pillbox*sinc function, which is the FT of the rectangular function,
 (Synthesis Imaging in Radio Astronomy II, p.143) with uniform weighting.
 */
/*
void gridding_visibilities_sinc(unsigned long int ncoords, double *u, double *v, complexd *vis, double len, int sizeg, complexd *new_vis, unsigned long int *count)
    {
        unsigned long int p,c;
        double x,y,sinc;
        double  inc = 2*len/sizeg;
        unsigned long int size = sizeg*sizeg;
        
        complexd* temp_grid_vis = (complexd *) calloc(size,sizeof(complexd));
        
        for (unsigned long int k = 0; k < ncoords; k++)
        {
            unsigned int pu = (unsigned int) ((u[k] + len) / inc);
            unsigned int pv = (unsigned int) ((v[k] + len) / inc);
            unsigned long int pc = (unsigned long int) pv * sizeg + pu;
            
            x = PI*(pu-u[k])/inc;
            y = PI*(pv-v[k])/inc;
            sinc = sin(x)*sin(y)/(x*y);
            
            temp_grid_vis[pc].real += vis[k].real*sinc;
            temp_grid_vis[pc].imag += vis[k].imag*sinc;
        }
        
        c=0;
        for (p = 0; p < size; p++)
        {
            if (temp_grid_vis[p].real)
            {
                new_vis[c].real = temp_grid_vis[p].real/count[c];
                new_vis[c].imag = temp_grid_vis[p].imag/count[c];
                
                c++;
            }
        }
        free(temp_grid_vis);
    }


void convolve_with_sinc(unsigned long int ncoords, double *u, double *v, complexd *vis, double fov, double wavelength, complexd *new_vis)
{
    double x,y,sinc;
    double du = wavelength/fov*ARCS2RAD;
    
    for (unsigned long int i = 0; i < ncoords; i++)
    {
        new_vis[i].real = 0.;
        new_vis[i].imag = 0.;
        for (unsigned long int k = 0; k < ncoords; k++)
        {
            x = PI*(u[i]-u[k])/du;
            y = PI*(v[i]-v[k])/du;
            sinc = sin(x)*sin(y)/(x*y);
        
            new_vis[i].real += vis[k].real*sinc;
            new_vis[i].imag += vis[k].imag*sinc;
        }
    }
}
*/
#ifdef __cplusplus
}
#endif
