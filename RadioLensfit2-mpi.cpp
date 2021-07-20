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
    containing source (SNR), position and flux must be provided.  

    A text file containing the list of the galaxies with the measured ellipticities will be generated.
 
    This version allows MS with spectral windows having different number of rows and frequency channels. 

    Command line input parameters:
    argv[1]  source catalog filename 
    argv[2]  number of sources
    argv[3]  MS filename prefix

    MPI version is implemented only for the faceting option!
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
#include "galaxy_fitting.h"
#include "galaxy_visibilities.h"
#include "distributions.h"
#include "evaluate_uv_grid.h"
#include "source_extraction.h"

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
    
    if (argc != 4)
    {
      if (rank == 0)
      {
        cout << "ERROR: bad number of parameters!" << endl;
        cout << "usage: RadioLensfit2-mpi <source catalog filename> <num_sources> <MS filename prefix>" << endl;
        cout << "number of MS must be equal to the number of MPI tasks" << endl;
        cout << "filename of input MS should be <prefix><IF number>.MS" << endl;
      }  
#ifdef USE_MPI
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
#else
      exit(EXIT_FAILURE);
#endif
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
    char filename[100];
    sprintf(filename,"%s%d.MS",argv[3],rank);
    RL_MeasurementSet* ms = ms_open(filename);
    cout << "rank " << rank << ": reading " << filename << "... " << endl;

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

    const double efficiency = EFFICIENCY;     // system efficiency
    const double SEFD = SEFD_SKA;    // System Equivalent Flux Density (in micro-Jy) of each SKA1 antenna

    const double ref_frequency_hz = REF_FREQ;  // Reference frequency in Hz at which fluxes are measured    

    if (rank == 0)
    { 
      cout << "Number baselines: " << num_baselines << endl;
      cout << "Accumulation time (sec): " << time_acc << endl;
      cout << "Number of channels per IF: " << num_channels << endl;
      cout << "Channels bandwidth (Hz): " << channel_bandwidth_hz << endl;
      cout << "Inital frequency (Hz): " << freq_start_hz << endl;
      cout << "Reference frequency (Hz): " << ref_frequency_hz << endl;
      cout << "Number of polarizations: " << ms_num_pols(ms) << endl;
    }
    cout << "rank " << rank << ": starting frequency (Hz): " << freq_start_hz << endl;
    cout << "rank " << rank << ": number of rows: " << num_rows << endl;
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
#ifdef USE_MPI
        MPI_Abort(MPI_COMM_WORLD,1);
        MPI_Finalize();
#else
        exit(EXIT_FAILURE);
#endif
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
#ifdef USE_MPI
        MPI_Abort(MPI_COMM_WORLD,1);
        MPI_Finalize();
#else
        exit(EXIT_FAILURE);
#endif
    }
  
    // allocate and read FLAG column
    bool *flag;
    try
    {
        flag = new bool[num_vis];
        sizeGbytes = num_vis*sizeof(bool)/((double)(1024*1024*1024));
        cout << "rank " << rank << ": allocated flag column: " << num_vis << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "rank " << rank << ": bad_alloc caught: " << ba.what() << '\n';
    }

    unsigned long int nF = ms_read_Flag(ms, 0, 0, num_channels, num_rows, "FLAG",flag, &status);
    if (status)
    {
      cout << "rank " << rank << ": ERROR reading MS - flag: " << status << endl;
#ifdef USE_MPI
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
#else
      exit(EXIT_FAILURE);
#endif
    }
    else
      cout << "rank " << rank << ": percentage of flagged visibilities: " << round(nF*100./num_vis) << "%" << endl;

    // allocate and read SIGMA column
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
#ifdef USE_MPI
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
#else
        exit(EXIT_FAILURE);
#endif
    }
    double sigma2 = (SEFD*SEFD)/(2.*time_acc*channel_bandwidth_hz*efficiency*efficiency);
    for (unsigned long int i = 0; i<num_vis; i++)
        sigma2_vis[i] = sigma2; // visibility noise variance
 
    cout << "rank " << rank << ": MS data GBytes: " << totGbytes << endl;
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
    sizeGbytes = nge*sizeof(double)/((double)(1024*1024*1024));
    totGbytes += sizeGbytes;

    bool readSNR = true;
    unsigned long int ngalaxies = read_catalog(nge, argv[1],gflux,gscale,ge1,ge2,l,m,SNR_vis,readSNR);
    if (rank == 0)  cout << "Read catalog. Number of sources: " << ngalaxies << endl;

#ifdef USE_MPI
    end_data = MPI_Wtime();
    data_time = end_data - start_data;
#else
    end_data = current_timestamp();
    data_time = (double)(end_data - start_data)/1000.;
#endif
    
    // Sky model visibilities computation --------------------------------------------------------------------------------------------------------------------------
    // Allocate Galaxy and Sky Model Visibilities for this  MS
    complexd *visGal, *visSkyMod;
    try
    {
        visGal = new complexd[num_vis];
        cout << "rank " << rank << ": allocated galaxy visibilities: " << num_vis << ", size = " << sizeGbytes  << " GB" << endl;
        sizeGbytes = num_vis*sizeof(complexd)/((double)(1024*1024*1024));
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

    // Pre-compute wavenumber and spectral factor for each channel of the current MS 
    // They corresponds to the central frequency of each channel
    double *wavenumbers = new double[num_channels];
    double *spec = new double[num_channels];
    double ch_freq = freq_start_hz + 0.5*channel_bandwidth_hz;

    for (unsigned int i = 0; i < num_channels; i++)
    {
       wavenumbers[i] = 2.0 * PI * ch_freq / C0;
       spec[i] = pow(ch_freq/ref_frequency_hz,-0.7);
       ch_freq += channel_bandwidth_hz;
      }
    
#ifdef USE_MPI
    double start_model = MPI_Wtime();
#else
    long long start_model = current_timestamp();
#endif
    // compute sky model for this MS
    sky_model(wavenumbers,spec, channel_bandwidth_hz, time_acc, num_channels,
               ngalaxies, gflux, gscale, l, m, num_coords, uu_metres, vv_metres, ww_metres, visGal, visSkyMod);
    cout << "rank " << rank << ": computed sky model " << endl;
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
    if (rank==0) cout << num_models << " samples in galaxy scale-length, " << Rmin << " < r0 < " << Rmax << " arcsec" << endl;
    
    // Set likelihood computation parameters
    likelihood_params par;
    par.numr = numR;
    par.ro = Ro;
    par.rprior = rprior;
    par.nchannels = num_channels;  
    par.band_factor = channel_bandwidth_hz*PI/C0;
    par.acc_time = time_acc;
    par.spec = spec;
    par.wavenumbers = wavenumbers; 
    par.count = 0;

#ifdef FACET
    // Faceting uv coordinates ----------------------------------------------------------------------------------------
    // Compute max uv distance in units of wavelenghts (for all spectral windows)
    len = len*wavenumbers[num_channels-1]/(2*PI);
#ifdef USE_MPI
    double send_len[nprocs],recv_len[nprocs];          // Will contain the maximum uv distance (in units of wavelengths) for each spectral window
    for (int spw = 0; spw < nprocs; spw++) send_len[spw] = len;
    MPI_Barrier(MPI_COMM_WORLD);
    double com_time = -MPI_Wtime();
    // Send and receive uv distance for each spectral window from each task 
    MPI_Alltoall(send_len,1,MPI_DOUBLE,recv_len,1,MPI_DOUBLE,MPI_COMM_WORLD);
    com_time += MPI_Wtime();
    for (int spw = 0; spw < nprocs; spw++) len = fmax(len,recv_len[spw]); 
#endif
    int facet = facet_size(RMAX,len);
    unsigned long int ncells = facet*facet;
    unsigned long int* count;
    try
    {
        count = new unsigned long int[ncells*nprocs];
        sizeGbytes = ncells*nprocs*sizeof(unsigned long int)/((double)(1024*1024*1024));
        totGbytes += sizeGbytes;
    }
     catch (bad_alloc& ba)
    {
        cerr << "rank " << rank << ": bad_alloc caught: " << ba.what() << '\n';
    }
    cout << "rank " << rank << ": allocated  array counter: " << sizeGbytes  << " GB" << endl;

    unsigned long int facet_ncoords = evaluate_uv_grid_size(rank,nprocs,len,wavenumbers,num_channels,num_coords, uu_metres, vv_metres, facet, flag, count);

    double* facet_u;
    double* facet_v;
    try
    {
        facet_u = new double[facet_ncoords];
        facet_v = new double[facet_ncoords];
        sizeGbytes = 2*facet_ncoords*sizeof(double)/((double)(1024*1024*1024));
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "rank " << rank << ": bad_alloc caught: " << ba.what() << '\n';
    }
    cout << "rank " << rank << ": allocated facet coordinates: " << sizeGbytes  << " GB" << endl;
    totGbytes += sizeGbytes;

    // allocate facet model visibilities 
    double* visMod;
    try
    {
      visMod = new double[num_models*facet_ncoords];
      sizeGbytes = num_models*facet_ncoords*sizeof(double)/((double)(1024*1024*1024));
      cout << "rank " << rank << ": allocated facet model visibilities, size = " << sizeGbytes  << " GB" << endl;
      totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
      cerr << "rank " << rank << ": bad_alloc caught: " << ba.what() << '\n';
    }
  
#ifdef USE_MPI
    // allocate facet vis and sigma2 for the gridded visibilities of this MS to be sent to the task that will process the current source 
    unsigned long int temp_facet_nvis = ncells;
    unsigned long int* temp_count;
    complexd* temp_facet_visData;
    double* temp_facet_sigma2;
    try
    {
        temp_count = new unsigned long int[temp_facet_nvis];
        temp_facet_visData = new complexd[temp_facet_nvis];
        temp_facet_sigma2 = new double[temp_facet_nvis];
        sizeGbytes = temp_facet_nvis*(sizeof(complexd)+sizeof(double)+sizeof(unsigned long int))/((double)(1024*1024*1024));
        cout << "rank " << rank << ": allocated my facet visibilities, variances and count: " << temp_facet_nvis << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "rank " << rank << ": bad_alloc caught: " << ba.what() << '\n';
    }
#endif

    // allocate facet data visibilities 
    // for the fitting of all visibilities of the current source by this task
    unsigned long int facet_nvis = ncells*nprocs; // and for reduction of partial facets from the other tasks
    complexd* facet_visData;
    double* facet_sigma2;
    try
    {
        facet_visData = new complexd[facet_nvis];
        facet_sigma2 = new double[facet_nvis];
        sizeGbytes = (facet_nvis*(sizeof(complexd)+sizeof(double)))/((double)(1024*1024*1024));
        cout << "rank " << rank << ": allocated gridded visibilities and variances: " << facet_nvis << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "rank " << rank << ": bad_alloc caught: " << ba.what() << '\n';
    }

    par.mod = visMod;
    par.uu = facet_u;
    par.vv = facet_v;
    //par.weights = weights;
    par.data = facet_visData;
    par.sigma2 = facet_sigma2;   
#else
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
    par.mod = visMod;
    par_sigma2 = sigma2_vis;
#endif
    cout << "rank " << rank << ": Total Visibilities GBytes: " << totGbytes << endl;
      
    // Data Processing -----------------------------------------------------------------------------------------------------------------------------------
    // Only rank 0 writes the output file 
    FILE *pFile = 0;
    if (rank == 0)
    {
      char output[100];
      sprintf(output,"ellipticities_%dIF-%dch.txt",nprocs,num_channels);
      pFile = fopen(output,"w");
      fprintf(pFile, "flux | e1 | m_e1 | err1 | e2 | m_e2 | err2 | 1D var | SNR |   l  |  m  | \n");
    }
#ifdef USE_MPI
    start_extraction = MPI_Wtime();
    double *recv_facet_buffer, *recv_sigma2_buffer;
    double *send_facet_buffer, *send_sigma2_buffer;
    unsigned long int *send_count_buffer, *recv_count_buffer;
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
        R_mu[src] = exp(mu);
        //R_mu[src] = gscale[g];
#ifdef FACET
        int facet = facet_size(R_mu[src],len);
        unsigned long int size = facet*facet;
#endif
        if (rank == src)  // proc k will fit the current source
        {
          my_g = g;  
          for (int nRo=1; nRo<numR; nRo++)   
            rprior[nRo] = rfunc(mu,R_STD,Ro[nRo]);  // set log(prior) for scalelength of source g
#ifdef FACET
           // extract my averaged visibilities, sigma2 and count contribution (my MS) for source g and store them in the corresponding section of my source facet array 
           par.facet = facet;
           source_extraction(rank,facet,&par,&(facet_visData[rank*size]),&(facet_sigma2[rank*size]),&(count[rank*size]),l0, m0, gflux[g], R_mu[src], 0., 0., visSkyMod, visData, visGal, sigma2_vis, flag, num_channels, num_coords, uu_metres, vv_metres, ww_metres, len);
#ifdef USE_MPI
           // Buffer facet vis, sigma2 and count for collection of source g from the other procs (for their MS contribution) 
           recv_facet_buffer = (double *) facet_visData;
           send_facet_buffer = (double *) &(facet_visData[rank*size]);
           recv_sigma2_buffer = facet_sigma2;
           send_sigma2_buffer = &(facet_sigma2[rank*size]);
           recv_count_buffer = count;
           send_count_buffer = &(count[rank*size]);
        }
        else 
        {
           // extract my averaged visibilities, sigma2 and count contribution (my MS) for source g and store them in the temporary IF facet array
           source_extraction(rank,facet,&par,temp_facet_visData, temp_facet_sigma2,temp_count,l0, m0, gflux[g], R_mu[src], 0., 0., visSkyMod, visData, visGal, sigma2_vis, flag, num_channels, num_coords, uu_metres, vv_metres, ww_metres, len);
           // Buffer facet vis, sigma2 and count (my MS) of source g to send to proc = k 
           recv_facet_buffer = 0;
           send_facet_buffer = (double *) temp_facet_visData;
           recv_sigma2_buffer = 0;
           send_sigma2_buffer = temp_facet_sigma2;
           recv_count_buffer = 0;
           send_count_buffer = temp_count;
#endif
#else
           source_extraction(l0, m0, gflux[g], R_mu[k], 0., 0., &par, visSkyMod, visData, visGal, sigma2_vis, num_channels, num_coords, uu_metres, vv_metres, ww_metres);
#endif
        }
#ifdef USE_MPI
        //MPI_Barrier(MPI_COMM_WORLD);
        // Proc k collects facet vis and sigma2 (and their counts) of the current source from the other procs (for their MS contribution)
        com_time -= MPI_Wtime();
        double time_Gather = MPI_Wtime();
        MPI_Gather(send_facet_buffer,2*size,MPI_DOUBLE,recv_facet_buffer,2*size,MPI_DOUBLE,k,MPI_COMM_WORLD);
        MPI_Gather(send_sigma2_buffer,size,MPI_DOUBLE,recv_sigma2_buffer,size,MPI_DOUBLE,k,MPI_COMM_WORLD);
        MPI_Gather(send_count_buffer,size,MPI_UNSIGNED_LONG,recv_count_buffer,size,MPI_UNSIGNED_LONG,k,MPI_COMM_WORLD);
        com_time += MPI_Wtime();
#endif
        g++;
        k++;
      }
      
#ifdef FACET
      if (my_g >= 0)
      {
#ifdef USE_MPI
         // average partial sum of visibilities
         average_facets(nprocs, par.facet*par.facet, facet_visData, facet_sigma2, count);
#endif
         // compute facet coordinates their number
         par.ncoords = evaluate_facet_coords(par.uu, par.vv, len, par.facet, count);
      }
#endif

      // source fitting of this task --------------------------------------------------------------------------------------------------
#ifdef USE_MPI
      start_fitting = MPI_Wtime();
#else
      start_fitting = current_timestamp();
#endif
      double mes_e1, mes_e2, maxL;
      double var_e1, var_e2, oneDimvar;
      if (my_g >= 0) 
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

         if (rank == k)
         {
            res[0] = mes_e1; res[1] = mes_e2;
            res[2] = var_e1; res[3] = var_e2;
            res[4] = oneDimvar; res[5] = maxL;
         }
#ifdef USE_MPI 
         //MPI_Barrier(MPI_COMM_WORLD);
         // Bcast shape results from proc k to the others
         com_time -= MPI_Wtime();
         double time_Bcast = MPI_Wtime();
         MPI_Bcast(res,6,MPI_DOUBLE,k,MPI_COMM_WORLD);          
         com_time += MPI_Wtime();
#endif
         if (res[5] < -1e+10 || res[2] < VAR || res[3] < VAR || res[4] < VAR)  // bad measurement
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
               unsigned long int ch_vis = ch*num_coords;
         
               // remove current round source model from the sky model
               data_galaxy_visibilities(spec[ch], wavenumbers[ch], par.band_factor, time_acc, 0., 0., R_mu[k],
                                           gflux[ind], l0, m0, num_coords, uu_metres, vv_metres, ww_metres, &(visGal[ch_vis]));
                  
               for (unsigned int i = ch_vis; i<ch_vis+num_coords; i++)
               {
                  visSkyMod[i].real -= visGal[i].real;
                  visSkyMod[i].imag -= visGal[i].imag;
               }
            
               // remove current source model fit from original data
               data_galaxy_visibilities(spec[ch], wavenumbers[ch], par.band_factor, time_acc, res[0], res[1], R_mu[k],
                                           gflux[ind], l0, m0, num_coords, uu_metres, vv_metres, ww_metres, &(visGal[ch_vis]));
                  
               for (unsigned int i = ch_vis; i<ch_vis+num_coords; i++)
               {
                  visData[i].real -= visGal[i].real;
                  visData[i].imag -= visGal[i].imag;
               }
            }
            if (rank == 0) fprintf(pFile, "%f | %f | %f | %e | %f | %f | %e | %e | %f | %f | %f \n",gflux[ind],ge1[ind],res[0],sqrt(res[2]),ge2[ind],res[1],sqrt(res[3]),res[4],SNR_vis[ind],l0/(ARCS2RAD),m0/(ARCS2RAD));
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

#ifdef FACET
        int facet = facet_size(R_mu[k],len);
        unsigned long int size = facet*facet;
#endif
        if (rank == src) // proc k will perform the fitting of the current source
        {
          my_g = gal;
          for (int nRo=1; nRo<numR; nRo++)
                rprior[nRo] = rfunc(mu,R_STD,Ro[nRo]);
#ifdef FACET
          par.facet = facet;
          source_extraction(rank,facet, &par, &(facet_visData[rank*size]), &(facet_sigma2[rank*size]), &(count[rank*size]),l0, m0, flux, R_mu[src], 0., 0., visSkyMod, visData, visGal, sigma2_vis, flag, num_channels, num_coords, uu_metres, vv_metres, ww_metres, len);
#ifdef USE_MPI
          // Buffer facet vis and sigma2 of source g for collection from the other procs of their MS contribution
          recv_facet_buffer = (double *) facet_visData;
          send_facet_buffer = (double *) &(facet_visData[rank*size]);
          recv_sigma2_buffer = facet_sigma2;
          send_sigma2_buffer = &(facet_sigma2[rank*size]);
          recv_count_buffer = count;
          send_count_buffer = &(count[rank*size]);
        }
        else
        {
          source_extraction(rank,facet,&par,temp_facet_visData, temp_facet_sigma2,temp_count,l0, m0, flux, R_mu[src], 0., 0., visSkyMod, visData, visGal, sigma2_vis, flag, num_channels, num_coords, uu_metres, vv_metres, ww_metres, len);
          // Buffer temp facet vis and sigma2 (my MS section) of source g to send to proc = k
          recv_facet_buffer = 0;
          send_facet_buffer = (double *) temp_facet_visData;
          recv_sigma2_buffer = 0;
          send_sigma2_buffer = temp_facet_sigma2;
          recv_count_buffer = 0;
          send_count_buffer = temp_count;
#endif
#else
          source_extraction(l0, m0, flux, R_mu[k], 0., 0., &par, visSkyMod, visData, visGal, sigma2_vis, num_channels, num_coords, uu_metres, vv_metres, ww_metres);
#endif
        }
#ifdef USE_MPI
        //MPI_Barrier(MPI_COMM_WORLD);
        // Proc k collects facet vis and sigma2 of the current source from the other procs (for their MS contribution)
        com_time -= MPI_Wtime();
        double time_gather = MPI_Wtime();
        MPI_Gather(send_facet_buffer,2*size,MPI_DOUBLE,recv_facet_buffer,2*size,MPI_DOUBLE,k,MPI_COMM_WORLD);
        MPI_Gather(send_sigma2_buffer,size,MPI_DOUBLE,recv_sigma2_buffer,size,MPI_DOUBLE,k,MPI_COMM_WORLD);
        MPI_Gather(send_count_buffer,size,MPI_UNSIGNED_LONG,recv_count_buffer,size,MPI_UNSIGNED_LONG,k,MPI_COMM_WORLD);
        com_time += MPI_Wtime();
#endif
        b++;
        k++;
      }
        
      if (my_g >= 0)
      {
#ifdef USE_MPI
        // average partial sum of visibilities and count coordinates
        average_facets(nprocs, par.facet*par.facet, facet_visData, facet_sigma2, count);
#endif
        par.ncoords = evaluate_facet_coords(par.uu, par.vv, len, par.facet, count);
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
        cout << "rank " << rank << ": n. " << my_g << " flux = " << gflux[my_g] << ": measured e = " << mes_e1 << "," << mes_e2 << endl;
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

          if (rank == k)
          {
            res[0] = mes_e1; res[1] = mes_e2;
            res[2] = var_e1; res[3] = var_e2;
            res[4] = oneDimvar; res[5] = maxL;
          }
#ifdef USE_MPI  // Bcast shape results from proc k to the other
          //MPI_Barrier(MPI_COMM_WORLD);
          com_time -= MPI_Wtime();
          double time_Bcast = MPI_Wtime();
          MPI_Bcast(res,6,MPI_DOUBLE,k,MPI_COMM_WORLD);
          com_time += MPI_Wtime();
#endif    
          if (rank == 0)   // write shape measurement (always!)
             fprintf(pFile, "%f | %f | %f | %e | %f | %f | %e | %e | %f | %f | %f  \n",gflux[gal],ge1[gal],res[0],sqrt(res[2]),ge2[gal],res[1],sqrt(res[3]),res[4],SNR_vis[gal],l0/(ARCS2RAD),m0/(ARCS2RAD));
          if (res[5] < -1e+10 || res[2] < VAR || res[3] < VAR || res[4] < VAR) // bad measurement 
          {
             bad++; 
             res[0] = 0.; res[1] = 0.; // remove round source model from original data
          }
#ifdef _OPENMP
#pragma omp parallel for
#endif
          for (unsigned int ch = 0; ch < num_channels; ch++)
          {
            unsigned long int ch_vis = ch*num_coords;

            // remove current round source model from the sky model
            data_galaxy_visibilities(spec[ch], wavenumbers[ch], par.band_factor, time_acc, 0., 0., R_mu[k],
                                       gflux[gal], l0, m0, num_coords, uu_metres, vv_metres, ww_metres, &(visGal[ch_vis]));

            for (unsigned int i = ch_vis; i<ch_vis+num_coords; i++)
            {
               visSkyMod[i].real -= visGal[i].real;
               visSkyMod[i].imag -= visGal[i].imag;
            }
 
            // remove current source model from original data
            data_galaxy_visibilities(spec[ch], wavenumbers[ch], par.band_factor, time_acc, res[0], res[1], R_mu[k],
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
       cout << "Total time (sec): " << total_time << endl;
       cout << "Data reading time (sec): " << data_time << endl;
#ifdef USE_MPI 
       cout << "Set up time (sec): " << total_time - data_time - model_time - extraction_time - fitting_time - com_time << endl;
    }
    cout << "rank " << rank << ": Communication time (sec): " << com_time << endl;
#else   
       cout << "Set up time (sec): " << total_time - data_time - model_time - extraction_time - fitting_time << endl;
    }
#endif
    cout << "rank " << rank << ": Sky model visibilities computation time (sec): " << model_time << endl;
    cout << "rank " << rank << ": Source extraction computation time (sec): " << extraction_time << endl;
    cout << "rank " << rank << ": Data fitting computation time (sec): " << fitting_time << endl;

#ifdef USE_MPI
     MPI_Finalize() ;
#endif

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
    delete[] count;
    delete[] facet_u;
    delete[] facet_v;
    delete[] facet_visData;
    delete[] facet_sigma2;
#ifdef USE_MPI
    delete[] temp_count;
    delete[] temp_facet_visData;
    delete[] temp_facet_sigma2;
#endif
#endif

    return 0;
}
