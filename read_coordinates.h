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

#ifndef _read_coordinates_h
#define _read_coordinates_h

#ifdef __cplusplus
extern "C" {
#endif
    
/**
 * @brief
 * Reads a coordinate file  
 *
 * @details
 * A coordinate file is an ASCII text file containing three
 * columns of comma- or space-separated values that represent the array index and the (x,y) coordinates. 
 * Each line corresponds to one point.
 */

int read_coords_oskar(const char* filename, unsigned long int ncoords, double* x, double* y, double *lenu, double *lenv);
    
double read_coord_ska(const char* filename1, const char* filename2, unsigned int ntimes, unsigned int *nbaselines,
                      double* x, double* y, double threshold, double *len);
    
unsigned long int read_uv_txt(const char* filename1, const char* filename2, unsigned int num_coords, double* x, double* y, double *len);
    
double read_uvw_txt(const char* filename1, const char* filename2, const char* filename3, unsigned int ntimes, unsigned int nbaselines, double* x, double* y, double* w, double *len);
    
#ifdef __cplusplus
}
#endif

#endif
