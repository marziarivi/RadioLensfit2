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

#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
 
#include "utils.h"
#include "galaxy_visibilities.h"
#include "distributions.h"
#include "data_simulation.h"

#ifdef __cplusplus
extern "C" {
#endif

// Simulate data visibilities
void data_simulation(double *wavenumbers, double *spec, double channel_bandwidth_hz, 
                     int time_acc, unsigned int num_channels, unsigned int num_baselines, float sigma,
                     unsigned int n_gal, double g1, double g2, double *ge1, double *ge2, double *gflux, double *gscale,
                     double *l, double *m, double *SNR_vis, unsigned int num_coords, double *uu_metres, double *vv_metres,
                     double *ww_metres, complexd *visGal, complexd* visData)
{
    
    //setup random number generator
    const gsl_rng_type * G;
    gsl_rng * gen;
    G = gsl_rng_mt19937;  // Mersenne Twister
    gen = gsl_rng_alloc(G);
    
    unsigned long int seed = random_seed();
    gsl_rng_set(gen,seed);
   
    double ee1,ee2,l0,m0,den;
    complexd z1,z2; 
    double band_factor = channel_bandwidth_hz*PI/C0;

    for (unsigned int g=0; g<n_gal; g++)
    {
        l0 = l[g];
        m0 = m[g];
            
        ee1 = ge1[g];
        ee2 = ge2[g];
            
        // apply shear g
        z1.real = ee1+g1;            z1.imag = ee2+g2;         // z1 = e+g
        z2.real = 1.+ee1*g1+ee2*g2;  z2.imag = ee2*g1-ee1*g2;  // z2 = 1+conj(g)*e
        den = z2.real*z2.real+z2.imag*z2.imag;
        ee1 = (z1.real*z2.real + z1.imag*z2.imag)/den;
        ee2 = (z1.imag*z2.real - z1.real*z2.imag)/den;        // e = z1/z2
            
        ge1[g] = ee1;
        ge2[g] = ee2;
            
        SNR_vis[g] = 0.;

#ifdef _OPENMP
#pragma omp parallel for
#endif            
        for (unsigned int ch = 0; ch < num_channels; ch++)
        {
            // generate galaxy visibilities
            unsigned long int ch_vis = (unsigned long int) ch*num_coords;
            data_galaxy_visibilities(spec[ch], wavenumbers[ch], band_factor, time_acc, ee1, ee2, gscale[g],
                                         gflux[g], l0, m0, num_coords, uu_metres, vv_metres, ww_metres, &(visGal[ch_vis]));
                
            // compute signal-to-noise
            double SNR_ch = 0.;
            for (unsigned long int vs = ch_vis; vs < ch_vis+num_coords; vs++)
                SNR_ch += visGal[vs].real*visGal[vs].real + visGal[vs].imag*visGal[vs].imag;
#ifdef _OPENMP
#pragma omp critical
#endif 
            SNR_vis[g] += SNR_ch;
                
            // add visibilities of current galaxy to data vis
            for (unsigned long int i= ch_vis; i<ch_vis+num_coords; i++)
            {
                visData[i].real += visGal[i].real;
                visData[i].imag += visGal[i].imag;
            }
        }
    //    SNR_vis[g] /= sigma;
    //    SNR_vis[g] = sqrt(SNR_vis[g]);
    //    printf("n. %d flux = %f SNR = %f\n",g,gflux[g],SNR_vis[g]);
    }
    
    // Add a random Gaussian noise component to the visibilities.
    // same sigma to all baselines
    
    double* sigmab = (double *) malloc(sizeof(double)*num_baselines);
    for (unsigned int b=0; b<num_baselines; b++) sigmab[b] = sqrt(sigma);
   
    for (unsigned int ch = 0; ch < num_channels; ch++)
    {
        unsigned long int ch_vis = (unsigned long int) ch*num_coords;
        add_system_noise(gen, num_coords, &(visData[ch_vis]), sigmab);
    }
    
    gsl_rng_free(gen);
    free(sigmab);
    
}


// Simulate sky model visibilities
void sky_model(double *wavenumbers, double *spec, double channel_bandwidth_hz, int time_acc, unsigned int num_channels,
                     unsigned int n_gal, double *gflux, double *gscale, double *l, double *m, unsigned int num_coords, 
                     double *uu_metres, double *vv_metres, double *ww_metres, complexd *visGal, complexd* visMod)
{
    double R_mu, l0, m0;
    double band_factor = channel_bandwidth_hz*PI/C0;

    for (unsigned int g=0; g<n_gal; g++)
    {
#ifndef SCALELENGTH_ON
        R_mu = exp(scale_mean(gflux[g]));
#else
        R_mu = gscale[g];
#endif
        l0 = l[g];
        m0 = m[g];

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (unsigned int ch = 0; ch < num_channels; ch++)
        {
            // generate galaxy visibilities
            unsigned long int ch_vis = (unsigned long int) ch*num_coords;
            data_galaxy_visibilities(spec[ch], wavenumbers[ch], band_factor, time_acc, 0., 0., R_mu,
                                         gflux[g], l0, m0, num_coords, uu_metres, vv_metres, ww_metres, &(visGal[ch_vis]));
            
            for (unsigned long int i= ch_vis; i<ch_vis+num_coords; i++)
            {
                visMod[i].real += visGal[i].real;
                visMod[i].imag += visGal[i].imag;
            }
        }
    }
}

#ifdef __cplusplus
}
#endif
