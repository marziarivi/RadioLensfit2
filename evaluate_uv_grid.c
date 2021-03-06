/*
 * Copyright (c) 2020 Marzia Rivi
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

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "evaluate_uv_grid.h"
#include "default_params.h"

#ifdef __cplusplus
extern "C" {
#endif
    
// gaussian function with variances given by the ellipticity orthogonal to the uv coverage's one
/* 
double weight_func(double u, double v)
{
    double sigma1 = 1.;  
    double sigma2 = 1.;  
    double w = exp(-u*u/(2*sigma1) - v*v/(2*sigma2));
    return w;
}
*/
// Compute facet size dependent on source flux
// use scalelength-flux relation: mu=log(theta_med[arcsec]) = ADD + ESP*log(flux[uJy]) to estimate galaxy scalelength
// multiply by a PSF factor to estimate facet FoV for that flux range and reference frequency
// compute facet cell size by the relation: Du[wavelength units]=1/FoV[rad] (Briggs 1999)
int facet_size(double theta_med, double len)
{
    double facet_du = 1./(PSF_NAT*theta_med*ARCS2RAD);
    int facet_size = 2*ceil(len/facet_du);
    return facet_size;
}

// Compute number of uv grid coordinates (coordinates are put in the center of the cell, only non-empty cells are considered)
unsigned long int evaluate_uv_grid_size(int rank, int nprocs, double len, double *wavenumbers, unsigned int num_channels, unsigned long int ncoords, double* u, double* v, int sizeg, bool *flag, unsigned long int* count)
{
    unsigned long int p,n;
    unsigned long int size = sizeg*sizeg;
    unsigned long int *temp_count = &(count[rank*size]);
    memset(temp_count, 0, size*sizeof(unsigned long int));
    
    double inc = 2*len/sizeg;

    n=0;
    for (unsigned int ch = 0; ch < num_channels; ch++)
    {
      double inv_lambda = wavenumbers[ch]/(2.*PI);
      for (unsigned long int k = 0; k < ncoords; k++)
      {
        if (!flag[n])
        {
          unsigned int pu = (unsigned int) ((u[k]*inv_lambda + len) / inc);
          unsigned int pv = (unsigned int) ((v[k]*inv_lambda + len) / inc);
          unsigned long int pc = (unsigned long int) pv * sizeg + pu;
          temp_count[pc]++;
        }
        n++;
      }
    }

#ifdef USE_MPI
    MPI_Allgather(temp_count,size,MPI_UNSIGNED_LONG,count,size,MPI_UNSIGNED_LONG,MPI_COMM_WORLD);
    n = size;
    for (int k = 1; k < nprocs; k++)
      for (p = 0; p < size; p++)
      {
        count[p] += count[n];
        n++;
      }
#endif
 
    n = 0;
    for (p = 0; p < size; p++)  if (count[p]) n++;
    return n;
}


// Compute uv facet coordinates in wavelength units and their number
unsigned long int evaluate_facet_coords(double* grid_u, double* grid_v, double len, int sizeg, unsigned long int *count)
{
    unsigned int i,j;
    unsigned long int p,n;
    double  inc = 2*len/sizeg;
    unsigned long int size = sizeg*sizeg;

    n=0;
    for (p = 0; p < size; p++)
    {
        if (count[p])
        {
            j = p/sizeg;
            i = p%sizeg;
            grid_u[n] = (i+0.5)*inc - len;
            grid_v[n] = (j+0.5)*inc - len;
            n++;
         }
    }
    return n;
}



/*
 Compute gridded visibilities by adding all the original visibilities at the (u,v)
 points falling in the same cell
 This is a gridding by convolution with the pillbox function
 (Synthesis Imaging in Radio Astronomy II, p.143) with uniform weighting.
 */
void gridding_visibilities(double *wavenumbers, unsigned int num_channels, unsigned long int ncoords, double *u, double *v, complexd *vis, double *sigma2, double len, int sizeg, complexd *new_vis, double *new_sigma2, bool *flag, unsigned long int *count)
{
    unsigned int i,j;
    unsigned long int p,n;
    double  inc = 2*len/sizeg;
    unsigned long int size = sizeg*sizeg;
    memset(count, 0, size*sizeof(unsigned long int));

#ifdef USE_MPI
    memset(new_vis, 0, size*sizeof(complexd));
    memset(new_sigma2, 0, size*sizeof(double));
#else
    complexd* temp_grid_vis = (complexd *) calloc(size,sizeof(complexd));
    double* temp_sigma2 = (double *) calloc(size,sizeof(double));
#endif

    n = 0;
    for (unsigned int ch = 0; ch < num_channels; ch++)
    {
      double inv_lambda = wavenumbers[ch]/(2.*PI);
      for (unsigned long int k = 0; k < ncoords; k++)
      {
         if (!flag[n])
         {
            unsigned int pu = (unsigned int) ((u[k]*inv_lambda + len) / inc);
            unsigned int pv = (unsigned int) ((v[k]*inv_lambda + len) / inc);
            unsigned long int pc = (unsigned long int) pv * sizeg + pu;
#ifdef USE_MPI
            new_vis[pc].real += vis[n].real;
            new_vis[pc].imag += vis[n].imag;
            new_sigma2[pc] += sigma2[n];
#else
            temp_grid_vis[pc].real += vis[n].real;
            temp_grid_vis[pc].imag += vis[n].imag;
            temp_sigma2[pc] += sigma2[n];
#endif
            count[pc]++;
         }
         n++;
      }
    }

// in the MPI version the average and removal of empty cells are performed after the collection of gridded visibilities from the other spectral windows
#ifndef USE_MPI
    n=0;
    for (p = 0; p < size; p++)
    {
        if (count[p])
        {
            new_vis[n].real = temp_grid_vis[p].real/count[p];
            new_vis[n].imag = temp_grid_vis[p].imag/count[p];
            new_sigma2[n] = temp_sigma2[p]/(count[p]*count[p]);
            n++;
        }
    }
    free(temp_grid_vis);
    free(temp_sigma2);
#endif
}


// average partial summed visibilities contained in the IFs facet
void average_facets(int nIFs, unsigned int size, complexd* facet_vis, double* facet_sigma2 ,unsigned long int *count)
{
   unsigned long int p;
   unsigned long int n = size;

   for (int k = 1; k < nIFs; k++)
     for (p = 0; p < size; p++)
     {
        facet_vis[p].real += facet_vis[n].real;
        facet_vis[p].imag += facet_vis[n].imag;
        n++;
     }

   n = size;
   for (int k = 1; k < nIFs; k++)
     for (p = 0; p < size; p++)
     {
        facet_sigma2[p] += facet_sigma2[n];
        n++;
     }

   n = size;
   for (int k = 1; k < nIFs; k++)
     for (p = 0; p < size; p++)
     {
        count[p] += count[n];
        n++;
     }

   n = 0;
   for (p = 0; p < size; p++)
   {
      if (count[p])
      {
         facet_vis[n].real = facet_vis[p].real/count[p];
         facet_vis[n].imag = facet_vis[p].imag/count[p];
         facet_sigma2[n] = facet_sigma2[p]/(count[p]*count[p]);
         n++;
      }
   }
}


// Uniform gridding by convolution with the pillbox*sinc function, which is the FT of the rectangular function,
// (Synthesis Imaging in Radio Astronomy II, p.143) with uniform weighting.
/*
void gridding_visibilities_sinc(unsigned long int ncoords, double *u, double *v, complexd *vis, double len, int sizeg, complexd *new_vis, unsigned long int *count)
 i   {
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
