/*
 * Copyright (c) 2018 Lance Miller, Marzia Rivi
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

#ifndef _random_gaussian_h
#define _random_gaussian_h

/**
 * @file random_gaussian.h
 */
#include <gsl/gsl_rng.h>

#ifdef __cplusplus
extern "C" {
#endif
    
    /**
     * @brief
     * Generates a random number from a Gaussian distribution with zero mean
     * and unit variance.
     *
     * @details
     * This function returns a random number from a Gaussian distribution with
     * zero mean and unit variance.
     *
     * The random number generator may be seeded by calling srand() prior to this
     * function.
     *
     * @param[in] another If not NULL, then this is used to return a second random number.
     */
    double random_gaussian(double* another, gsl_rng * gen);
    
   /**
     * @brief
     * this function generates a random seed, either from the unix
     * /dev/random hardware stream or from the clock if the
     * hardware stream is not available
     */
    unsigned long int random_seed();
    
#ifdef __cplusplus
}
#endif


#endif
