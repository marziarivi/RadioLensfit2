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


//  RadioLensfit2-mpi.cpp
/*
    Measure star forming galaxy ellipticies from a radio weak lensing observation.
    A model fitting approach is adopted where the likelihood is marginalised over position, 
    flux and scalelength source parameters.

    Data visibilities and observation configuration must be provided in a Measurement Set.
    The number of galaxies and the corresponding source catalog (ordered by decreasing flux) 
    containing source SNR, position and flux must be provided.  

    A text file containing the list of the galaxies with the measured ellipticities will be generated.
 
    Command line input parameters:
    argv[1]  source catalog filename 
    argv[2]  number of sources
    argv[3]  MS1 filename
    ....   
*/

// MPI version enabled only for the FACET case!
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
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

#include "datatype.h"
#include "default_params.h"
#include "utils.h"
#include "measurement_set.h"
#include "data_simulation.h"
#include "read_catalog.h"
#include "galaxy_fitting-mpi.h"
#include "galaxy_visibilities.h"
#include "distributions.h"
#include "evaluate_uv_grid.h"

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
    if (rank==0) cout << "Number of OpenMP threads = " << num_threads << endl;
#endif
    
    if (argc != nprocs+3)
    {
      if (rank == 0)
      {
        cout << "ERROR: bad number of parameters!" << endl;
        cout << "usage: RadioLensfit2-MS <source catalog filename> <num_sources> <filename MS1> <filename MS2> .... " << endl;
        cout << "number of MS must be equal to the number of MPI tasks" << endl;
      }  
      exit(EXIT_FAILURE);
    }

    double data_time = 0.;
    double extraction_time = 0.;
    double fitting_time = 0.;
#ifdef USE_MPI
    double start_data,end_data,start_fitting,end_fitting,start_extraction,end_extraction;
    start_data = MPI_Wtime();
#else
    long long start_data,end_data,start_fitting,end_fitting,start_extraction,end_extraction;
    start_data = current_timestamp();
#endif

    // Read Measurement Set --------------------------------------------------------------------------------------------------------------------------------------------------------------------
    RL_MeasurementSet* ms = ms_open(argv[3+rank]);

    //double RA = ms_phase_centre_ra_rad(ms);                 // Phase Centre coordinates
    //double Dec = ms_phase_centre_dec_rad(ms);   
    const unsigned int num_stations = ms_num_stations(ms);        // Number of stations
    const unsigned int num_channels = ms_num_channels(ms);        // Number of frequency channels
    const unsigned long int num_rows = ms_num_rows(ms);                // Number of rows 
    const double freq_start_hz = ms_freq_start_hz(ms); //1280e+6;          // Start Frequency, in Hz
    const double channel_bandwidth_hz = ms_freq_inc_hz(ms); //240e+6      // Frequency channel bandwidth, in Hz
    const double full_bandwidth_hz = channel_bandwidth_hz * num_channels;  // Frequency total bandwidth, in Hz
    const int time_acc = ms_time_inc_sec(ms);                     // accumulation time (sec)
    const unsigned int num_baselines = num_stations * (num_stations - 1) / 2;
    const unsigned int tot_nchannels = num_channels*nprocs;

    const double efficiency = EFFICIENCY;     // system efficiency
    const double SEFD = SEFD_SKA;    // System Equivalent Flux Density (in micro-Jy) of each SKA1 antenna

    const double ref_frequency_hz = REF_FREQ;  // Reference frequency in Hz at which fluxes are measured    
    double freq0 = freq_start_hz;        // Initial frequency for all dataset
#ifdef USE_MPI
    double com_time = MPI_Wtime(); 
    // Bcast from 0 to other procs the starting frequency of the all dataset
    MPI_Bcast(&freq0,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    com_time -= MPI_Wtime();
#endif

    if (rank == 0)
    { 
      cout << "Reference frequency (Hz): " << ref_frequency_hz << endl;
      cout << "rank " << rank << ": Number baselines: " << num_baselines << endl;
      cout << "rank " << rank << ": Number of channels: " << tot_nchannels << endl;
      cout << "rank " << rank << ": Channels bandwidth (Hz): " << channel_bandwidth_hz << endl;
      cout << "rank " << rank << ": Inital frequency (Hz): " << freq0 << endl;
      cout << "rank " << rank << ": Accumulation time (sec): " << time_acc << endl;
      cout << "rank " << rank << ": Number of rows: " << num_rows << endl;
    }
    cout << "rank " << rank << ": my starting frequency (Hz): " << freq_start_hz << endl;
    double sizeGbytes, totGbytes = 0.;
    
    // Allocate and read uv coordinates
    unsigned long int num_coords = num_rows; 
    double* uu_metres = new double[num_coords];
    double* vv_metres = new double[num_coords];
    double* ww_metres = new double[num_coords];
    sizeGbytes = 3*num_coords*sizeof(double)/((double)(1024*1024*1024));
    cout << "rank " << rank << ": allocated original coordinates: " << sizeGbytes  << " GB" << endl;
    totGbytes += sizeGbytes;
    
    int status = 0;
    double len = ms_read_coords(ms,0,num_coords,uu_metres,vv_metres,ww_metres,&status);
    if (status) 
    {
        cout << "rank " << rank << ": ERROR reading MS - uvw points: " << status << endl;
        exit(EXIT_FAILURE);
    }
    
    // Allocate and read Data visibilities of the current MS 
    unsigned long int num_vis  = (unsigned long int) num_channels * num_coords;
    complexd *visData;
    try
    {
        visData = new complexd[num_vis];
        sizeGbytes = num_vis*sizeof(complexd)/((double)(1024*1024*1024));
        cout << "rank " << rank << ": allocated original data visibilities: " << num_vis << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "rank " << rank << ": bad_alloc caught: " << ba.what() << '\n';
    }

    ms_read_vis(ms, 0, 0, num_channels, num_rows, "DATA", visData, &status);
    if (status) 
    {
        cout << "rank " << rank << ": ERROR reading MS - DATA column: " << status << endl;
        exit(EXIT_FAILURE);
    }

    double *sigma2_vis;
    try
    {
        sigma2_vis = new double[num_vis];
        sizeGbytes = num_vis*sizeof(double)/((double)(1024*1024*1024));
        cout << "rank " << rank << ": allocated sigma2 visibilities: " << num_vis << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "rank " << rank << ": bad_alloc caught: " << ba.what() << '\n';
    }

    //ms_read_sigma(ms, 0, num_coords, sigma2_vis, &status);
    if (status)
    {
        cout << "rank " << rank << ": ERROR reading MS - sigma: " << status << endl;
        exit(EXIT_FAILURE);
    }
    for (unsigned long int i = 0; i<num_vis; i++)
        sigma2_vis[i] = (SEFD*SEFD)/(2.*time_acc*channel_bandwidth_hz*efficiency*efficiency); // visibility noise variance
 
    ms_close(ms); 

    // Read galaxy catalogue --------------------------------------------------------------------------------------------------------------------------
    unsigned long int nge = atof(argv[2]);
    
    double *gflux = new double[nge];
    double *l = new double[nge];
    double *m = new double[nge];
    double *ge1 = new double[nge];
    double *ge2 = new double[nge];
    double *gscale = new double[nge];
    double *SNR_vis = new double[nge];
 
    unsigned long int ngalaxies = read_catalog(nge, argv[1],gflux,gscale,ge1,ge2,l,m,SNR_vis);
    if (rank == 0) cout << "Number of sources: " << ngalaxies << endl;
   
#ifdef USE_MPI
    end_data = MPI_Wtime();
    data_time = end_data - start_data;
    double start_model = MPI_Wtime();
#else
    end_data = current_timestamp();
    data_time = (double)(end_data - start_data)/1000.;
    long long start_model = current_timestamp();
#endif
    
    // Sky model visibilities computation --------------------------------------------------------------------------------------------------------------------------
    // Allocate Galaxy and Sky Model Visibilities for this  MS
    complexd *visGal, *visSkyMod;
    try
    {
        visGal = new complexd[num_vis];
        cout << "rank " << rank << ": allocated galaxy visibilities: " << num_vis << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "rank " << rank << ": bad_alloc caught: " << ba.what() << '\n';
    }

    try
    {
        visSkyMod = new complexd[num_vis];
        cout << "rank " << rank << ": allocated sky model visibilities: " << num_vis << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "rank " << rank << ": bad_alloc caught: " << ba.what() << '\n';
    }

    memset(visSkyMod, 0, num_vis*sizeof(complexd));    

    // Pre-compute wavenumber and spectral factor for each channel of all MS 
    // They corresponds to the central frequency of each channel
    double *wavenumbers = new double[tot_nchannels];
    double ch_freq = freq0 + 0.5*channel_bandwidth_hz;
    double *spec = new double[tot_nchannels];

    for (unsigned int ch = 0; ch < tot_nchannels; ch++)
    {
        wavenumbers[ch] = 2.0 * PI * ch_freq / C0;
        spec[ch] = pow(ch_freq/ref_frequency_hz,-0.7);
        ch_freq += channel_bandwidth_hz;
    }

    // compute sky model for this MS
    const unsigned int my_freq_index = rank*num_channels;
    sky_model(&(wavenumbers[my_freq_index]), &(spec[my_freq_index]), channel_bandwidth_hz, time_acc, num_channels, 
               ngalaxies, gflux, gscale, l, m, num_coords, uu_metres, vv_metres, ww_metres, visGal, visSkyMod);
    
#ifdef USE_MPI
    double model_time = MPI_Wtime() - start_model;
#else
    long long end_model = current_timestamp();
    double model_time = (double)(end_model - start_model)/1000.;
#endif

    // Setup Model Fitting ----------------------------------------------------------------------------------------------------------------------------------------
    // define steps in galaxy scalelength (Ro in ARCSEC)
    const double Rmin = RMIN;
    double Rmax = RMAX;
    int numR = NUM_R;
    double* Ro = new double[numR];
    double* rprior = new double[numR];
    Ro[0] = 0.;
    Ro[1] = Rmin;

    rprior[0]= 0.;
    int nRo=2;
    while (nRo<numR && Ro[nRo-1] < Rmax)
    {
        // quadratic spacing of samples
        double Rinterval = 0.08 + 0.4*pow( ((Ro[nRo-1]-Rmin)/(Rmax-Rmin)), 2);
        Ro[nRo] = Ro[nRo-1] + Rinterval;
        nRo++;
    }
    numR = nRo;
    if (Ro[nRo-1]>Rmax) Rmax=Ro[nRo-1];
 
    int num_models = numR-1;
    
    // Set likelihood computation parameters
    likelihood_params par;
    par.numr = numR;
    par.ro = Ro;
    par.rprior = rprior;
    par.nchannels = tot_nchannels;  // all channels for source fitting
    par.band_factor = channel_bandwidth_hz*PI/C0;
    par.acc_time = time_acc;
    par.spec = spec;
    par.wavenumbers = wavenumbers; // wavenumbers for the model
       
#ifdef FACET
    // Faceting uv coordinates ----------------------------------------------------------------------------------------
    int facet = facet_size(RMAX,len);  
    unsigned long int ncells = facet*facet;
    unsigned long int* count = new unsigned long int[ncells];
    unsigned long int facet_ncoords = evaluate_max_uv_grid_size(len,num_coords, uu_metres, vv_metres, facet, count);

    double* facet_u = new double[facet_ncoords];
    double* facet_v = new double[facet_ncoords];
    sizeGbytes = (2*facet_ncoords*sizeof(double)+ncells*sizeof(unsigned long int))/((double)(1024*1024*1024));
    cout << "rank " << rank << ": allocated grid coordinates and array counter: " << sizeGbytes  << " GB" << endl;
    totGbytes += sizeGbytes;
  
#ifdef USE_MPI
    // allocate temporary facet vis and sigma2 for the gridded visibilities of this MS to be sent to the task that will process the current source 
    unsigned long int temp_facet_nvis = facet_ncoords*num_channels;
    complexd* temp_facet_visData;
    double* temp_facet_sigma2;
    try
    {
        temp_facet_visData = new complexd[temp_facet_nvis];
        temp_facet_sigma2 = new double[temp_facet_nvis];
        sizeGbytes = (temp_facet_nvis*(sizeof(complexd)+sizeof(double)))/((double)(1024*1024*1024));
        cout << "rank " << rank << ": allocated temporary gridded visibilities and variances: " << temp_facet_nvis << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "rank " << rank << ": bad_alloc caught: " << ba.what() << '\n';
    }
#endif

    // allocate facet data and model visibilities for all MS (same ncoords and num channels)
    // for the fitting of all visibilities of the current source by this task
    unsigned long int facet_nvis = facet_ncoords*tot_nchannels;
    complexd* facet_visData;
    double* facet_sigma2;
    try
    {
        facet_visData = new complexd[facet_nvis];
        facet_sigma2 = new double[facet_nvis];
        sizeGbytes = (facet_nvis*sizeof(complexd)+facet_ncoords*sizeof(double))/((double)(1024*1024*1024));
        cout << "rank " << rank << ": allocated gridded visibilities and variances: " << facet_nvis << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "rank " << rank << ": bad_alloc caught: " << ba.what() << '\n';
    }

    double* visMod;
    try
    {
        unsigned long int model_ncoords = facet_ncoords;
        visMod = new double[num_models*model_ncoords*tot_nchannels];
        sizeGbytes = num_models*model_ncoords*tot_nchannels*sizeof(double)/((double)(1024*1024*1024));
        cout << "rank " << rank << ": allocated models: num_models= " << num_models << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "rank " << rank << ": bad_alloc caught: " << ba.what() << '\n';
    }

    par.ncoords = facet_ncoords;
    par.uu = facet_u;
    par.vv = facet_v;
    //par.weights = weights;
    par.data = facet_visData;
    par.count = count;
    par.mod = visMod;
    par.sigma2 = facet_sigma2;   
 
#else
    par_sigma2 = sigma2_vis;

    complexd* visMod;
    try
    {
        unsigned long int model_ncoords = num_coords;
        visMod = new complexd[num_models*model_ncoords*tot_nchannels];
        sizeGbytes = num_models*model_ncoords*tot_nchannels*sizeof(complexd)/((double)(1024*1024*1024));
        cout << "rank " << rank << ": allocated models: num_models= " << num_models << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "rank " << rank << ": bad_alloc caught: " << ba.what() << '\n';
    }

    par.ncoords = num_coords;
    par.uu = uu_metres;
    par.vv = vv_metres;
    par.ww = ww_metres;
    par.data = visGal;
    par.count = 0;
    par.mod = visMod;
#endif
    
    if (rank==0) 
    {
      cout << "Total Visibilities GBytes per rank: " << totGbytes << endl;
      cout << num_models << " samples in galaxy scale-length, " << Rmin << " < r0 < " << Rmax << " arcsec" << endl;
    }

    // Data Processing -----------------------------------------------------------------------------------------------------------------------------------
    // Only rank 0 writes the output file 
    FILE *pFile = 0;
    if (rank == 0)
    {
      pFile = fopen("ellipticities.txt","w");
      fprintf(pFile, "flux | e1 | m_e1 | err1 | e2 | m_e2 | err2 | 1D var | SNR |   l  |  m  | \n");
    }
#ifdef USE_MPI
    start_extraction = MPI_Wtime();
    double *recv_facet_buffer, *recv_sigma2_buffer;
    double *send_facet_buffer, *send_sigma2_buffer;
#else
    start_extraction = current_timestamp();
#endif

    double l0,m0;
    unsigned long int bad_list[ngalaxies];  // all tasks update the bad_list to avoid communication
    int k, bad = 0;
    long int my_g, ind; 
    double mu, R_mu[nprocs];

    long int g = 0; // source global index
    while (g < ngalaxies)
    { 
      k = 0; // source local index
      my_g = -1;  // if sources are finished, current task has my_g = -1 and has no source to fit

      // extract facet visibilities for nproc consecutive sources that will be fit each one by a different task -----------------------
      for(int src = 0; src < nprocs && g < ngalaxies; src++)
      {
        l0 = l[g];  m0 = m[g];
        mu = scale_mean(gflux[g]);
        R_mu[k] = exp(mu);
        //R_mu[k] = gscale[g];

        if (rank == k)  // proc k will fit the current source
        {
          my_g = g;     
          for (int nRo=1; nRo<numR; nRo++)   
            rprior[nRo] = rfunc(mu,R_STD,Ro[nRo]);  // set log(prior) for scalelength of source g
#ifdef FACET
           // extract my averaged visibilities and sigma2 contribution (my MS) for source g and store them in the corresponding section of my source facet array 
           source_extraction(rank, my_freq_index, facet_visData, facet_sigma2,l0, m0, gflux[g], R_mu[k], 0., 0., &par, visSkyMod, visData, visGal, sigma2_vis, num_channels, num_coords, uu_metres, vv_metres, ww_metres, len);
#ifdef USE_MPI
           // Buffer facet vis and sigma2 for collection of source g from the other procs (for their MS contribution) 
           recv_facet_buffer = (double *) facet_visData;
           send_facet_buffer = (double *) &(facet_visData[my_freq_index*facet_ncoords]);
           recv_sigma2_buffer = facet_sigma2;
           send_sigma2_buffer = &(facet_sigma2[my_freq_index*facet_ncoords]);
        }
        else 
        {
           // extract my averaged visibilities and sigma2 contribution (my MS) for source g and store them in the temporary IF facet array
           source_extraction(rank, 0, temp_facet_visData, temp_facet_sigma2,l0, m0, gflux[g], R_mu[k], 0., 0., &par, visSkyMod, visData, visGal, sigma2_vis, num_channels, num_coords, uu_metres, vv_metres, ww_metres, len);
           // Buffer temp facet vis and sigma2 (my MS) of source g to send to proc = k 
           recv_facet_buffer = 0;
           send_facet_buffer = (double *) temp_facet_visData;
           recv_sigma2_buffer = 0;
           send_sigma2_buffer = temp_facet_sigma2;
#endif
#else
           source_extraction(rank, my_freq_index, l0, m0, gflux[g], R_mu[k], 0., 0., &par, visSkyMod, visData, visGal, sigma2_vis, num_channels, num_coords, uu_metres, vv_metres, ww_metres);
#endif
        }
#ifdef USE_MPI
        // Proc k collects facet vis and sigma2 of the current source from the other procs (for their MS contribution)
        com_time -= MPI_Wtime();
        MPI_Gather(send_facet_buffer,2*temp_facet_nvis,MPI_DOUBLE,recv_facet_buffer,2*facet_nvis,MPI_DOUBLE,k,MPI_COMM_WORLD);
        MPI_Gather(send_sigma2_buffer,temp_facet_nvis,MPI_DOUBLE,recv_sigma2_buffer,facet_nvis,MPI_DOUBLE,k,MPI_COMM_WORLD);
        com_time += MPI_Wtime();
#endif
        g++;
        k++;
      }
      
      // source fitting of this task --------------------------------------------------------------------------------------------------
#ifdef USE_MPI
      start_fitting = MPI_Wtime();
#else
      start_fitting = current_timestamp();
#endif
      double mes_e1, mes_e2, maxL;
      double var_e1, var_e2, oneDimvar;
      if (my_g >= 0)  //this task will fit the current source
      {
        source_fitting(rank, &par, &mes_e1, &mes_e2, &var_e1, &var_e2, &oneDimvar, &maxL);
        cout << "rank " << rank << ": n. " << my_g << " flux = " << gflux[my_g] << ": measured e = " << mes_e1 << "," << mes_e2 << endl;
      }
#ifdef USE_MPI
      fitting_time += MPI_Wtime() - start_fitting;
#else
      fitting_time += (double) (current_timestamp() -start_fitting)/1000.;
#endif

      // removal of the fit sources visibilities from the MS data and sky model ----------------------------------------------------------
      double res[6]; // shape fitting results to be sent to the other procs
      ind = g - k;  // source global index
      k = 0;  // source local index
      for (int src = 0 ; src < nprocs && ind < ngalaxies; src++)
      {
         l0 = l[ind];  m0 = m[ind];
#ifdef USE_MPI 
         // Bcast shape results from proc k to the others
         if (rank == k)
         {
            res[0] = mes_e1; res[1] = mes_e2;
            res[2] = var_e1; res[3] = var_e2;
            res[4] = oneDimvar; res[5] = maxL;
         }

         com_time -= MPI_Wtime();
         MPI_Bcast(res,6,MPI_DOUBLE,k,MPI_COMM_WORLD);          
         com_time += MPI_Wtime();

         if (rank != k)
         {
            mes_e1 = res[0]; mes_e2 = res[1];
            var_e1 = res[2]; var_e2 = res[3];
            oneDimvar = res[4]; maxL = res[5];
         }
#endif
         if (maxL < -1e+10 || var_e1 < 1e-4 || var_e2 < 1e-4 || oneDimvar < 1e-4)  // bad measurement
         {
            bad_list[bad] = ind;  // each proc store the index of bad measurements to be fit again at the end
            bad++;
         }
         else     
         {
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (unsigned int ch = 0; ch < num_channels; ch++)
            {
               unsigned int ch_ind = my_freq_index + ch;
               unsigned long int ch_vis = ch*num_coords;
         
               // remove current round source model from the sky model
               data_galaxy_visibilities(spec[ch_ind], wavenumbers[ch_ind], par.band_factor, time_acc, 0., 0., R_mu[k],
                                           gflux[ind], l0, m0, num_coords, uu_metres, vv_metres, ww_metres, &(visGal[ch_vis]));
                  
               for (unsigned int i = ch_vis; i<ch_vis+num_coords; i++)
               {
                  visSkyMod[i].real -= visGal[i].real;
                  visSkyMod[i].imag -= visGal[i].imag;
               }
            
               // remove current source model fit from original data
               data_galaxy_visibilities(spec[ch_ind], wavenumbers[ch_ind], par.band_factor, time_acc, mes_e1, mes_e2, R_mu[k],
                                           gflux[ind], l0, m0, num_coords, uu_metres, vv_metres, ww_metres, &(visGal[ch_vis]));
                  
               for (unsigned int i = ch_vis; i<ch_vis+num_coords; i++)
               {
                  visData[i].real -= visGal[i].real;
                  visData[i].imag -= visGal[i].imag;
               }
            }
            if (rank == 0) fprintf(pFile, "%f | %f | %f | %f | %f | %f | %f | %f | %f | %f | %f \n",gflux[ind],ge1[ind],mes_e1,sqrt(var_e1),ge2[ind],mes_e2,sqrt(var_e2),oneDimvar,SNR_vis[ind],l0/(ARCS2RAD),m0/(ARCS2RAD));
         }
         k++;
         ind++;
      }
    }

    // Re-fitting bad sources again -------------------------------------------------------------------------------------------------------------------------------------
    if (rank == 0) cout << "Re-fitting " << bad << " bad sources" << endl;
    
    int nsources = bad;
    bad = 0;
    unsigned long int b = 0; 
    while (b < nsources)
    {
      k = 0;  // source local index
      my_g = -1;
      
      for (int src = 0; src < nprocs && b < nsources; src++)
      {
        unsigned long int gal = bad_list[b]; // source global index
        double flux = gflux[gal];
        mu = scale_mean(flux);
        R_mu[k] = exp(mu);
        //R_mu[k] = gscale[gal];
        l0 = l[gal];  m0 = m[gal];

        if (rank == k) // proc k will perform the fitting of the current source
        {
          my_g = gal;
          for (int nRo=1; nRo<numR; nRo++)
                rprior[nRo] = rfunc(mu,R_STD,Ro[nRo]);
#ifdef FACET
          source_extraction(rank, my_freq_index, facet_visData, facet_sigma2,l0, m0, flux, R_mu[k], 0., 0., &par, visSkyMod, visData, visGal, sigma2_vis, num_channels, num_coords, uu_metres, vv_metres, ww_metres, len);
#ifdef USE_MPI
          // Buffer facet vis and sigma2 of source g for collection from the other procs of their MS contribution
          recv_facet_buffer = (double *) facet_visData;
          send_facet_buffer = (double *) &(facet_visData[my_freq_index*facet_ncoords]);
          recv_sigma2_buffer = facet_sigma2;
          send_sigma2_buffer = &(facet_sigma2[my_freq_index*facet_ncoords]);
        }
        else
        {
          source_extraction(rank, 0., temp_facet_visData, temp_facet_sigma2,l0, m0, flux, R_mu[k], 0., 0., &par, visSkyMod, visData, visGal, sigma2_vis, num_channels, num_coords, uu_metres, vv_metres, ww_metres, len);
          // Buffer temp facet vis and sigma2 (my MS section) of source g to send to proc = k
          recv_facet_buffer = 0;
          send_facet_buffer = (double *) temp_facet_visData;
          recv_sigma2_buffer = 0;
          send_sigma2_buffer = temp_facet_sigma2;
#endif
#else
          source_extraction(rank, my_freq_index, l0, m0, flux, R_mu[k], 0., 0., &par, visSkyMod, visData, visGal, sigma2_vis, num_channels, num_coords, uu_metres, vv_metres, ww_metres);
#endif
        }
#ifdef USE_MPI
        // Proc k collects facet vis and sigma2 of the current source from the other procs (for their MS contribution)
        com_time -= MPI_Wtime();
        MPI_Gather(send_facet_buffer,2*temp_facet_nvis,MPI_DOUBLE,recv_facet_buffer,2*facet_nvis,MPI_DOUBLE,k,MPI_COMM_WORLD);
        MPI_Gather(send_sigma2_buffer,temp_facet_nvis,MPI_DOUBLE,recv_sigma2_buffer,facet_nvis,MPI_DOUBLE,k,MPI_COMM_WORLD);
        com_time += MPI_Wtime();
#endif
        b++;
        k++;
      }
        
#ifdef USE_MPI
      start_fitting = MPI_Wtime();
#else
      start_fitting = current_timestamp();
#endif
      double mes_e1, mes_e2, maxL;
      double var_e1, var_e2, oneDimvar;
      if (my_g >= 0) //this task will fit the current source
      {
        source_fitting(rank, &par, &mes_e1, &mes_e2, &var_e1, &var_e2, &oneDimvar, &maxL);
        cout << "rank " << rank << ": n. " << my_g << " flux = " << gflux[my_g] << "): measured e = " << mes_e1 << "," << mes_e2 << endl;
      }
#ifdef USE_MPI
      fitting_time += MPI_Wtime() - start_fitting;
#else
      fitting_time += (double) (current_timestamp() - start_fitting)/1000.;
#endif

      // Removal of the fit sources visibilities from the MS data and sky model
      double res[6];
      ind = b - k;  
      k = 0; // source local index
      for (int proc = 0 ; proc < nprocs && ind < nsources; proc++)
      {
          unsigned long int gal = bad_list[ind];  //source global index
          l0 = l[gal];  m0 = m[gal];

#ifdef USE_MPI   // Bcast shape results from proc k to the other
          if (rank == k)
          {
            res[0] = mes_e1; res[1] = mes_e2;
            res[2] = var_e1; res[3] = var_e2;
            res[4] = oneDimvar; res[5] = maxL;
          }

          com_time -= MPI_Wtime();
          MPI_Bcast(res,6,MPI_DOUBLE,k,MPI_COMM_WORLD);
          com_time += MPI_Wtime();

          if (rank != k)
          {
           mes_e1 = res[0]; mes_e2 = res[1];
           var_e1 = res[2]; var_e2 = res[3];
           oneDimvar = res[4]; maxL = res[5];
          }
#endif    
          if (rank == 0)   // write shape measurement (always!)
             fprintf(pFile, "%f | %f | %f | %f | %f | %f | %f | %f | %f | %f | %f  \n",gflux[gal],ge1[gal],mes_e1,sqrt(var_e1),ge2[gal],mes_e2,sqrt(var_e2),oneDimvar,SNR_vis[gal],l0/(ARCS2RAD),m0/(ARCS2RAD));
          
          if (maxL <= -1e+10 || var_e1 < 1e-4 || var_e2 < 1e-4 || oneDimvar < 1e-4) // bad measurement
          {
             bad++; 
             mes_e1 = 0.; mes_e2 = 0.; // remove round source model from original data
          }
#ifdef _OPENMP
#pragma omp parallel for
#endif
          for (unsigned int ch = 0; ch < num_channels; ch++)
          {
            unsigned int ch_ind = my_freq_index + ch;
            unsigned long int ch_vis = ch*num_coords;

            // remove current round source model from the sky model
            data_galaxy_visibilities(spec[ch_ind], wavenumbers[ch_ind], par.band_factor, time_acc, 0., 0., R_mu[k],
                                       gflux[gal], l0, m0, num_coords, uu_metres, vv_metres, ww_metres, &(visGal[ch_vis]));

            for (unsigned int i = ch_vis; i<ch_vis+num_coords; i++)
            {
               visSkyMod[i].real -= visGal[i].real;
               visSkyMod[i].imag -= visGal[i].imag;
            }
 
            // remove current source model from original data
            data_galaxy_visibilities(spec[ch_ind], wavenumbers[ch_ind], par.band_factor, time_acc, mes_e1, mes_e2, R_mu[k],
                                     gflux[gal], l0, m0, num_coords, uu_metres, vv_metres, ww_metres, &(visGal[ch_vis]));
 
            for (unsigned int i = ch_vis; i<ch_vis+num_coords; i++)
            {
              visData[i].real -= visGal[i].real;
              visData[i].imag -= visGal[i].imag;
            }
          }
          k++;
          ind++;
      }
    }
    
    if (rank == 0 && pFile != 0) fclose(pFile);
    
#ifdef USE_MPI
    end_extraction = MPI_Wtime();
    extraction_time = end_extraction - start_extraction - fitting_time - com_time;
    double end_tot = MPI_Wtime();
    double total_time = end_tot - start_tot;
#else
    end_extraction = current_timestamp();
    extraction_time = (double)(end_extraction - start_extraction)/1000. - fitting_time;
    long long end_tot = current_timestamp();
    double total_time = (double)(end_tot - start_tot)/1000.;
#endif

    if (rank == 0)
    {
       cout << "Removed " << bad << " bad data galaxies" << endl << endl;
       cout << "Set up time (sec): " << total_time - data_time - model_time - extraction_time - fitting_time << endl;
       cout << "Data reading time (sec): " << data_time << endl;
       cout << "Sky model visibilities time (sec): " << model_time << endl;
       cout << "Source extraction computation time (sec): " << extraction_time << endl;
       cout << "Data fitting computation time (sec): " << fitting_time << endl;
#ifdef USE_MPI
       cout << "Communication time (sec): " << com_time << endl;
#endif
       cout << "Total time (sec): " << total_time << endl;
    }
    // free memory ----------------------------------------------------------------------------------------------------------------
    delete[] visMod;
    delete[] visGal;
    delete[] visSkyMod;
    delete[] visData;
    delete[] sigma2_vis;
    delete[] Ro;
    delete[] rprior;
    delete[] gflux;
    delete[] gscale;
    delete[] ge1;
    delete[] ge2;
    delete[] l;
    delete[] m;
    delete[] SNR_vis;
    delete[] uu_metres;
    delete[] vv_metres;
    delete[] ww_metres;
    delete[] wavenumbers;
    delete[] spec;
#ifdef FACET
    delete[] facet_u;
    delete[] facet_v;
    delete[] facet_visData;
    delete[] facet_sigma2;
#ifdef USE_MPI
    delete[] temp_facet_visData;
    delete[] temp_facet_sigma2;
#endif
#endif

#ifdef USE_MPI
     MPI_Finalize() ;
#endif
    return 0;
}
