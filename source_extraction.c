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

#include "galaxy_visibilities.h"
#include "evaluate_uv_grid.h"
#include "source_extraction.h"


// Source extraction --------------------------------------------------------
// Notice: instructions execution order is important as arrays are reused!!!!

#ifdef FACET
void source_extraction(int rank, int facet, likelihood_params *par, complexd *facet_vis, double *facet_sigma2,
                       unsigned long int *count, double l0, double m0, double flux, double mu, double e1, double e2,  
                       complexd *visSkyMod, complexd *visData, complexd *visGal, double *sigma2_vis, unsigned int nchannels, 
                       unsigned long int num_coords, double *uu_metres, double *vv_metres, double *ww_metres, double len)
#else
void source_extraction(double l0, double m0, double flux, double mu, double e1, double e2, likelihood_params *par, complexd *visSkyMod, 
                       complexd *visData, complexd *visGal,double *sigma2_vis, unsigned int nchannels, unsigned long int num_coords, 
                       double *uu_metres, double *vv_metres, double *ww_metres)
#endif
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
   for (unsigned int ch = 0; ch < nchannels; ch++)
   {
     unsigned long int ch_vis = ch*num_coords;
     // get small round source at the current position and flux
     data_galaxy_visibilities((par->spec)[ch], (par->wavenumbers)[ch], par->band_factor, par->acc_time, e1, e2, mu,
                              flux, l0, m0, num_coords, uu_metres, vv_metres, ww_metres, &(visGal[ch_vis]));
    
     for (unsigned long int i = ch_vis; i<ch_vis+num_coords; i++)
     {
#ifdef USE_MPI
        // extract source data visibilities without modifying the sky model
        visGal[i].real += visData[i].real - visSkyMod[i].real;
        visGal[i].imag += visData[i].imag - visSkyMod[i].imag;
#else
        // remove source model from the sky model
        visSkyMod[i].real -= visGal[i].real;
        visSkyMod[i].imag -= visGal[i].imag;
        
        // remove sky model from the original data (i.e. all other sources approximation) in order to have only the current galaxy data
        visGal[i].real = visData[i].real - visSkyMod[i].real;
        visGal[i].imag = visData[i].imag - visSkyMod[i].imag;
#endif  
     }
    
#ifdef FACET
     // Phase shift data visibilities (to be done after gridding because real data will be gridded)
     data_visibilities_phase_shift((par->wavenumbers)[ch], l0, m0, num_coords, uu_metres, vv_metres, ww_metres, &(visGal[ch_vis]));
   }    
   // gridding visibilities
   gridding_visibilities(par->wavenumbers,nchannels,num_coords,uu_metres,vv_metres,visGal,sigma2_vis,len,facet,facet_vis,facet_sigma2,count);
#else
   }
  par->l0 = l0;
  par->m0 = m0;
#endif
  
}

