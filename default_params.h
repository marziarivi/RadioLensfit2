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

// Parameters for numerical likelihood marginalisation of scalelength
// scalelength range
#define RMIN 0.3   // scalelength range minimum
#define RMAX 3.5   // scalelength range maximum
#define NUM_R 29   // number of scalelength bins 

// Parameters of source scalelength-flux relation: alpha[arcsec] = K_FAC*flux[uJy]^ESP
#define K_FAC 0.39
#define ESP 0.33

// Parameters for facet size computation
#define NUM_S 10      // number of flux bins
#define PSF_NAT 38.5  // PSF factor for natural weighting 

#endif
