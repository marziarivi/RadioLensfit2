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

#ifndef DEFAULT_H_
#define DEFAULT_H_

/**
 * @file default_params.h
 */

// Antennae noise parameters
#define EFFICIENCY 0.9     // system efficiency
#define SEFD_SKA 400e+6    // System Equivalent Flux Density (in micro-Jy) of each SKA1 antenna
#define SEFD_MKT 551e+6    // SEFD of each MeerKat antenna (in micro-Jy)

// Frequency reference of the source flux
#define REF_FREQ 1.4e+9

// Conversion of the beam frequency measure unit in GHz
#define BEAM_const 65e-9   

// Parameters for facet size computation
#define PSF_NAT 38  // PSF factor for natural weighting

// Parameter for catalog simulation
#define FMAX 200     // flux range maximum [micro-Jy]
#define NUM_ORIENT 5 // 2*NUM_ORIENT is the number of sampled orientations (points on the circle of radius |e|) for each ellipticity module
#define N_POINTS 1000 // number of points of the cumulative distribution function for random data generation    

// Parameters for numerical likelihood marginalisation of scalelength 
// scalelength range
#define RMIN 0.3   // scalelength range minimum [arcsec]
#define RMAX 3.5   // scalelength range maximum [arcsec]
#define NUM_R 29   // number of scalelength bins 

// Parameters of source scalelength-flux relation: log(alpha[arcsec]) = ADD + ESP*log(flux[uJy])
#define ADD -0.93
#define ESP 0.33

// Parameter for catalog simulation
#define FMAX 200     // flux range maximum [micro-Jy]
#define NUM_ORIENT 5 // 2*NUM_ORIENT is the number of sampled orientations (points on the circle of radius |e|) for each ellipticity module
#define N_POINTS 1000 // number of points of the cumulative distribution function for random data generation

// Prior distributions parameters
#define R_STD 0.3136  // scalelength prior - lognormal standard deviation

#define M_EXP -1.34   // flux prior - exponent of the flux power law: p(S) \propto S^(-M_EXP)    
#define NORM_S 3.0825 // flux prior - normalisation factor over a square degree

#define E_MAX 0.804   // ellipticity modulus prior - cutoff
#define E_0 0.0732    // ellipticity modulus prior - circularity parameter
#define DISP 0.2298   // ellipticity modulus prior - dispersion
#define NORM_E 2.595  // ellipticity modulus prior - normalisation factor per square degree 

// Parameters for the computation of the cumulative distribution (numerical integration)
#define JMAX 30
#define EPS 1.0e-5

// Parameter for likelihood sampling
#define NP_MAX 30  // min number of sampling points with likelihood above 5%ML

#endif
