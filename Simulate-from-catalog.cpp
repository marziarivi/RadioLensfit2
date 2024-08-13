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


//  Simulate-from-catalog.cpp
/*  
    Simulate a radio weak lensing observation of a given number of sources from a given galaxy catalog.
    Radio telescope configuration and observing time sampling must be provided in a Measurement Set.
    An input shear is applied to the sources shape.
  
    Command line input parameters:  

    argv[1]  filename Measurement Set
    argv[2]  galaxy catalog filename
    argv[3]  number of galaxies
    argv[4]  g1 shear component
    argv[5]  g2 shear component
*/

#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

#include <iostream>
#include <new>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "default_params.h"
#include "datatype.h"
#include "utils.h"
#include "measurement_set.h"
#include "data_simulation.h"
#include "read_catalog.h"
#include "distributions.h"

using namespace std;

int main(int argc, char *argv[])
{
    int nprocs, rank, num_threads=1;
#ifdef USE_MPI
    MPI_Init(&argc, &argv) ;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs) ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double start_tot = MPI_Wtime();
#else
    nprocs=1;
    rank=0;
    long long start_tot;
    start_tot = current_timestamp();
#endif

#ifdef _OPENMP
#pragma omp parallel
    num_threads = omp_get_num_threads();
    if (rank == 0) cout << "Number of OpenMP threads = " << num_threads << endl;
#endif
    
    if (argc < 5)
    {
      if (rank == 0)
      {
        cout << "ERROR: parameter missing!" << endl;
        cout << "usage: Simulate <filename MS> <galaxy catalog filename> <number of sources> <shear1> <shear2> " << endl;
      }
#ifdef USE_MPI
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
#else
      exit(EXIT_FAILURE);
#endif
    }

    // Read Measurement Set --------------------------------------------------------------------------------------------------------------------------------------------------------------------
    char filename[100];
    sprintf(filename,"%s%d.MS",argv[1],rank);
    RL_MeasurementSet* ms = ms_open(filename);
    cout << "rank " << rank << ": reading " << filename << "... " << endl;

    //double RA = ms_phase_centre_ra_rad(ms);                 // Phase Centre coordinates
    //double Dec = ms_phase_centre_dec_rad(ms);   
    unsigned int num_stations = ms_num_stations(ms);        // Number of stations
    unsigned int num_channels = ms_num_channels(ms);        // Number of frequency channels
    unsigned int num_rows = ms_num_rows(ms);                // Number of rows 
    double freq_start_hz = ms_freq_start_hz(ms);            // Start Frequency, in Hz
    double channel_bandwidth_hz = ms_freq_inc_hz(ms);       // Frequency channel bandwidth, in Hz
    double full_bandwidth_hz = channel_bandwidth_hz * num_channels;  // Frequency total bandwidth, in Hz
    int time_acc = ms_time_inc_sec(ms);                     // accumulation time (sec)
    unsigned int num_pols = ms_num_pols(ms);                          // number of polarizations

    double efficiency = EFFICIENCY;     // system efficiency
    double SEFD = SEFD_SKA;    // System Equivalent Flux Density (in micro-Jy) of each SKA antenna
    

    double ref_frequency_hz = 1.4e+9;  //Reference frequency in Hz at which fluxes are measured
    
    unsigned int num_baselines = num_stations * (num_stations - 1) / 2;

    if (num_pols != 1)
    {
      cout << "MS ERROR: number of polarizations is " << num_pols << endl;
      cout << "A single polarization is required!" << endl;
#ifdef USE_MPI
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
#else
      exit(EXIT_FAILURE);
#endif
    }

    if (rank == 0)
    {
      cout << "Number stations: " << num_stations << endl;
      cout << "Number of channels: " << num_channels << endl;
      cout << "Channels bandwidth (Hz): " << channel_bandwidth_hz << endl;
      cout << "Reference frequency (Hz): " << ref_frequency_hz << endl;
      cout << "Starting frequency (Hz): " << freq_start_hz << endl;
      cout << "Accumulation time (sec): " << time_acc << endl;
    }
    double sizeGbytes, totGbytes = 0.;
    
    // Allocate and read uv coordinates 
    unsigned int num_coords = ms_num_rows(ms);
    cout << "rank " << rank << ": number of rows: " << num_coords << endl;
    double* uu_metres = new double[num_coords];
    double* vv_metres = new double[num_coords];
    double* ww_metres = new double[num_coords];
    sizeGbytes = 3*num_coords*sizeof(double)/((double)(1024*1024*1024));
    cout << "rank " << rank << ": allocated original coordinates: " << sizeGbytes  << " GB" << endl;
    totGbytes += sizeGbytes;
    
    int status;
    double len = ms_read_coords(ms,0,num_coords,uu_metres,vv_metres,ww_metres,&status);
    if (status) 
    {
        cout << "rank " << rank << ": ERROR reading MS - uvw points: " << status << endl;
#ifdef USE_MPI
        MPI_Abort(MPI_COMM_WORLD,1);
        MPI_Finalize();
#else
        exit(EXIT_FAILURE);
#endif
    }

    // Read galaxy catalogue --------------------------------------------------------------------------------------------------------------------------
    // For catalog generation use the code generate_catalog.c  

    unsigned int nge = atoi(argv[3]);

    double *gflux = new double[nge];
    double *gscale = new double[nge];
    double *ge1 = new double[nge];
    double *ge2 = new double[nge];
    double *l = new double[nge];
    double *m = new double[nge];
    double *temp_SNR = new double[nge];
    bool readSNR = false;

    unsigned int ngalaxies = read_catalog(nge, argv[2],gflux,gscale,ge1,ge2,l,m,temp_SNR,readSNR);
    if (rank == 0) cout << "Read catalog. Number of sources: " << ngalaxies << endl;

    // Allocate Galaxy and Sky Visibilities -----------------------------------------------------------------------------------------------------------------------
    unsigned long int num_rawvis  = (unsigned long int) num_channels * num_coords;
    complexd *visGal, *visData;
    try
    {
        visGal = new complexd[num_rawvis];
        sizeGbytes = num_rawvis*sizeof(complexd)/((double)(1024*1024*1024));
        cout << "rank " << rank << ": allocated galaxy visibilities: " << num_rawvis << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "rank " << rank << ": bad_alloc caught: " << ba.what() << '\n';
    }
    
    try
    {
        visData = new complexd[num_rawvis];
        cout << "rank " << rank << ": allocated original data visibilities: " << num_rawvis << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "rank " << rank << ": bad_alloc caught: " << ba.what() << '\n';
    }
    memset(visData, 0, num_rawvis*sizeof(complexd));
    
    // Visibilities Simulation --------------------------------------------------------------------------------------------------------------------------
    // shear to be applied
    double g1 = atof(argv[4]); 
    double g2 = atof(argv[5]);
 
#ifdef USE_MPI
    double start_sim = MPI_Wtime();
#else
    long long start_sim = current_timestamp();
#endif
    float sigma = (SEFD*SEFD)/(2.*time_acc*channel_bandwidth_hz*efficiency*efficiency); // visibility noise variance
    if (rank == 0) cout << "sigma_vis  = " << sqrt(sigma) << " muJy" << endl;
    
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
    
    data_simulation(wavenumbers, spec, channel_bandwidth_hz, time_acc, num_channels, 
                    num_baselines, sigma, ngalaxies, g1, g2, ge1, ge2, gflux, gscale, l, m, temp_SNR, num_coords, 
                    uu_metres, vv_metres, ww_metres, visGal, visData);
   
#ifdef USE_MPI
    double *SNR_vis = new double[ngalaxies];
    MPI_Reduce(temp_SNR,SNR_vis,ngalaxies,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    double simulation_time = MPI_Wtime() - start_sim;
#else
    double *SNR_vis = temp_SNR;
    long long end_sim = current_timestamp();
    double simulation_time = (double)(end_sim - start_sim)/1000.;
#endif

    // Write visibilities on the DATA column of the Measurement Set
    ms_write_vis(ms,0,0,num_channels,num_coords,visData);
    ms_close(ms);
   
    // Write catalog: SNR, positions [rad], flux [uJy], scalelength [arcsec] 
    if (rank == 0)
    {
      FILE *pf;
      sprintf(filename,"%s_SNR",argv[2]);
      pf = fopen(filename,"w");
      // fprintf(pf, "SNR | l | m | flux | scale | e1 | e2 \n"); 
      for (unsigned int g = 0; g < ngalaxies; g++)
      {
        SNR_vis[g] /= sigma;
        SNR_vis[g] = sqrt(SNR_vis[g]);
  
        fprintf(pf, "%f %f %f %f %f %f %f \n",SNR_vis[g],l[g],m[g],gflux[g],gscale[g],ge1[g],ge2[g]);
#ifdef GAUSSIAN
        cout << "n. " << g << " flux = " << gflux[g] << ", sigma = " << gscale[g] << ", SNR = " << SNR_vis[g] << endl;
#else
        cout << "n. " << g << " flux = " << gflux[g] << ", scale-length = " << gscale[g] << ", SNR = " << SNR_vis[g] << endl;
#endif
      }
      fclose(pf);
    }
#ifdef USE_MPI
    double total_time = MPI_Wtime() - start_tot;
#else
    long long end_tot = current_timestamp();
    double total_time = (double)(end_tot - start_tot)/1000.;
#endif
   
    if (rank == 0)
    {
      cout << "Simulation time (sec): " << simulation_time << endl;
      cout << "I/O time (sec): " << total_time - simulation_time << endl;
      cout << "Total time (sec): " << total_time << endl;
    }

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
#ifdef USE_MPI
    delete[] temp_SNR;
#endif  
    delete[] uu_metres;
    delete[] vv_metres;
    delete[] ww_metres;
    delete[] wavenumbers;
    delete[] spec;

#ifdef USE_MPI
     MPI_Finalize() ;
#endif
    return 0;
}
