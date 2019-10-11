/*
 * Copyright (c) 2019 Marzia Rivi
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


//  Simulate.cpp
/*  
    Simulate a radio weak lensing observation with a given reduced cosmic shear and a source flux threshold.
    Radio telescope configuration and observing time sampling must be provided in a Measurement Set.
  
    The simulated visibilities are stored in the DATA column of the same Measurement Set.
    The text file called RWL_galaxy_catalog_<FoV>_<Fmin>.txt containing the simulated sourse catalog is created, 
    where <FoV> is the field of view provided in arcmin and <Fmin> is the flux threshold in uJy.
    The number of sources is estimated from the flux prior, the size of the field of view and Fmin.

    Sources are simulated according to the ring test to avoid shape noise: same source flux and size for 2*NUM_ORIENT 
    ellipiticity values simmetrically distributed along the same ring. 
  
    Command line input parameters:  

    argv[1]  filename Measurement Set
    argv[2]  effective field of view [arcmin]
    argv[3]  minimum galaxy flux [muJy]
    argv[4]  applied shear 1st component
    argv[5]  applied shear 2nd component
*/

#ifdef _OPENMP
#include <omp.h>
#endif

#include <iostream>
#include <new>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "datatype.h"
#include "utils.h"
#include "measurement_set.h"
#include "data_simulation.h"
#include "generate_catalog.h"
#include "distributions.h"

#define NUM_ORIENT 5  // 2*NUM_ORIENT is the number of sampled orientations (points on the circle of radius |e|) for each ellipticity module 
#define FMAX 300      // maximum source flux [uJy]
#define RMIN 0.3      // minimum source scalelength [arcsec]
#define RMAX 3.5      // maximun source scalelength [arcsec]

using namespace std;

int main(int argc, char *argv[])
{
    long long start_tot, end_tot;
    start_tot = current_timestamp();

#ifdef _OPENMP
    int num_threads;
#pragma omp parallel
    num_threads = omp_get_num_threads();
    cout << "Number of OpenMP threads = " << num_threads << endl;
#endif
    
    if (argc < 5)
    {
        cout << "ERROR: parameter missing!" << endl;
        cout << "usage: Simulate <filename MS> <effective field of view [arcmin]> <min flux [muJy]> <shear1> <shear2> " << endl;
        exit(EXIT_FAILURE);
    }

    // Read Measurement Set --------------------------------------------------------------------------------------------------------------------------------------------------------------------
    RL_MeasurementSet* ms = ms_open(argv[1]);

    //double RA = ms_phase_centre_ra_rad(ms);                 // Phase Centre coordinates
    //double Dec = ms_phase_centre_dec_rad(ms);   
    unsigned int num_stations = ms_num_stations(ms);        // Number of stations
    unsigned int num_channels = ms_num_channels(ms);        // Number of frequency channels
    unsigned int num_rows = ms_num_rows(ms);                // Number of rows 
    double freq_start_hz = ms_freq_start_hz(ms);            // Start Frequency, in Hz
    double channel_bandwidth_hz = ms_freq_inc_hz(ms);       // Frequency channel bandwidth, in Hz
    double full_bandwidth_hz = channel_bandwidth_hz * num_channels;  // Frequency total bandwidth, in Hz
    int time_acc = ms_time_inc_sec(ms);                     // accumulation time (sec)

    double efficiency = 0.9;     // system efficiency
    double SEFD_SKA = 400e+6;    // System Equivalent Flux Density (in micro-Jy) of each SKA1 antenna
    double SEFD_MKT = 551e+6;    // SEFD of each MeerKat antenna (in micro-Jy)

    double ref_frequency_hz = 1.4e+9;  //Reference frequency in Hz at which fluxes are measured
    
    unsigned int num_baselines = num_stations * (num_stations - 1) / 2;

    cout << "Number baselines: " << num_baselines << endl;
    cout << "Number of channels: " << num_channels << endl;
    cout << "Channels bandwidth (Hz): " << channel_bandwidth_hz << endl;
    cout << "Reference frequency (Hz): " << ref_frequency_hz << endl;
    cout << "Starting frequency (Hz): " << freq_start_hz << endl;
    cout << "Accumulation time (sec): " << time_acc << endl;
    
    double sizeGbytes, totGbytes = 0.;
    double fov_eff_arcmin = atof(argv[2]);
    double fov_eff = fov_eff_arcmin*60.*ARCS2RAD; //1.22*C0/(freq_start_hz*diameter);  // field of view in RAD
    printf("field of view: %e [rad] %f [arcsec] \n",fov_eff,fov_eff/(ARCS2RAD));
    
    // Allocate and read uv coordinates 
    unsigned long int num_coords = ms_num_rows(ms);
    double* uu_metres = new double[num_coords];
    double* vv_metres = new double[num_coords];
    double* ww_metres = new double[num_coords];
    sizeGbytes = 3*num_coords*sizeof(double)/((double)(1024*1024*1024));
    cout << "allocated original coordinates: " << sizeGbytes  << " GB" << endl;
    totGbytes += sizeGbytes;
    
    int status;
    double len = ms_read_coords(ms,0,num_coords,uu_metres,vv_metres,ww_metres,&status);
    
    // Allocate Galaxy and Sky Visibilities -----------------------------------------------------------------------------------------------------------------------
    unsigned long int num_rawvis  = (unsigned long int) num_channels * num_coords;
    complexd *visGal, *visData;
    try
    {
        visGal = new complexd[num_rawvis];
        sizeGbytes = num_rawvis*sizeof(complexd)/((double)(1024*1024*1024));
        cout << "allocated galaxy visibilities: " << num_rawvis << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "bad_alloc caught: " << ba.what() << '\n';
    }
    
    try
    {
        visData = new complexd[num_rawvis];
        cout << "allocated original data visibilities: " << num_rawvis << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "bad_alloc caught: " << ba.what() << '\n';
    }
    memset(visData, 0, num_rawvis*sizeof(complexd));
    
    long long start_data,end_data;
    start_data = current_timestamp();

    //---------------------------------------------------------------------------------------------------------------------------------------------------
    // Generate galaxy catalogue --------------------------------------------------------------------------------------------------------------------------

    double Fmin = atof(argv[3]);
    unsigned long int nge = ceil((flux_CDF(beta, FMAX) - flux_CDF(beta, Fmin))*fov_eff_arcmin*fov_eff_arcmin);
   
    double *gflux = new double[nge];
    double *gscale = new double[nge];
    double *ge1 = new double[nge];
    double *ge2 = new double[nge];
    double *l = new double[nge];
    double *m = new double[nge];
 
    unsigned long int mygalaxies = galaxy_catalog(nge, NUM_ORIENT, fov_eff, RMIN, RMAX, Fmin, FMAX, gflux, gscale,ge1,ge2,l,m);
    cout << "num gal: " << mygalaxies << endl;
    
    // Visibilities Simulation --------------------------------------------------------------------------------------------------------------------------
    double g1 = atof(argv[4]);  // shear to be applied
    double g2 = atof(argv[5]);
    
    double sigma = (SEFD_SKA*SEFD_SKA)/(2.*time_acc*channel_bandwidth_hz*efficiency*efficiency); // visibility noise variance
    cout << "sigma_vis  = " << sqrt(sigma) << " muJy" << endl;
    
    double *SNR_vis = new double[mygalaxies];

    // Pre-compute wavenumber and spectral factor for each channel 
    // They corresponds to the central frequency of each channel
    double *wavenumbers = new double[num_channels];
    double ch_freq = freq_start_hz + 0.5*channel_bandwidth_hz;
    double *spec = new double[num_channels];

    for (unsigned int ch = 0; ch < num_channels; ch++)
    {
        wavenumbers[ch] = 2.0 * PI * ch_freq / C0;
        spec[ch] = pow(ch_freq/ref_frequency_hz,-0.7);
        ch_freq += channel_bandwidth_hz;
    }
    
    data_simulation(ref_frequency_hz, wavenumbers, spec, channel_bandwidth_hz, time_acc, num_channels, 
                    num_baselines, sigma, mygalaxies, g1, g2, ge1, ge2, gflux, gscale, l, m, SNR_vis, num_coords, 
                    uu_metres, vv_metres, ww_metres, visGal, visData);
    
    end_data = current_timestamp();
    double data_time = (double)(end_data - start_data)/1000.;

    // Write visibilities on the DATA column of the Measurement Set
    ms_write_vis(ms,0,0,num_channels,num_coords,visData);

    // Write catalog: positions [rad], flux [uJy], scalelength [arcsec]
    FILE *pf;
    char filename[200];
    sprintf(filename,"RWL_galaxy_catalog_%.1f_%.1f.txt",fov_eff_arcmin,Fmin);
    pf = fopen(filename,"w");
    // fprintf(pf, "SNR | l | m | flux | scale | e1 | e2 \n");

    for (unsigned long int g = 0; g<mygalaxies; g++)
         fprintf(pf, "%f %f %f %f %f %f %f \n",SNR_vis[g],l[g],m[g],gflux[g],gscale[g],ge1[g],ge2[g]);
    fclose(pf);
 
    end_tot = current_timestamp();
    double total_time = (double)(end_tot - start_tot)/1000.;
   
    cout << "Data generation time (sec): " << data_time << endl;
    
    ms_close(ms);
    
    cout << "I/O time (sec): " << total_time - data_time << endl;
    cout << "Total time (sec): " << total_time << endl;
    
    // free memory ----------------------------------------------------------------------------------------------------------------
    delete[] visData;
    delete[] visGal;
    delete[] ge1;
    delete[] ge2;
    delete[] gflux;
    delete[] gscale;
    delete[] l;
    delete[] m;
    delete[] SNR_vis;
    delete[] uu_metres;
    delete[] vv_metres;
    delete[] ww_metres;
    delete[] wavenumbers;
    delete[] spec;

    return 0;
}
