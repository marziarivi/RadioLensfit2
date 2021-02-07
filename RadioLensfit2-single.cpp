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


//  RadioLensfit2-single.cpp
//
//  Simulation of SKA1-MID observation and fitting of single galaxies 
//  in the field of view
//
//  argv[1]  Measurement Set filename
//  argv[2]  Galaxy catalog filename
//  argv[3]  number of galaxies
//  argv[4]  applied shear 1st component
//  argv[5]  applied shear 2nd component

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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>

#include "datatype.h"
#include "generate_random_values.h"
#include "utils.h"
#include "measurement_set.h"
#include "galaxy_visibilities.h"
#include "likelihood.h"
#include "marginalise_r.h"
#include "distributions.h"
#include "evaluate_uv_grid.h"
#include "utils.h"
#include "default_params.h"
#include "read_catalog.h"

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

    if (argc != 6)
    {
        cout << "ERROR: bad number of parameters!" << endl;
        cout << "usage: RadioLensfit.x <MS filename> <source catalog filename> <nge> <shear g1> <shear g2>" << endl;
        exit(EXIT_FAILURE);
    }

    // Read Measurement Set --------------------------------------------------------------------------------------------------------------------------------------------------------------------
    RL_MeasurementSet* ms = ms_open(argv[1]);

    //double RA = ms_phase_centre_ra_rad(ms);                 // Phase Centre coordinates
    //double Dec = ms_phase_centre_dec_rad(ms);   
    const unsigned int num_stations = ms_num_stations(ms);        // Number of stations
    const unsigned int num_channels = ms_num_channels(ms);        // Number of frequency channels
    const unsigned int num_rows = ms_num_rows(ms);                // Number of rows 
    const double freq_start_hz = ms_freq_start_hz(ms);          // Start Frequency, in Hz
    const double channel_bandwidth_hz = ms_freq_inc_hz(ms);      // Frequency channel bandwidth, in Hz
    const double full_bandwidth_hz = channel_bandwidth_hz * num_channels;  // Frequency total bandwidth, in Hz
    const int time_acc = ms_time_inc_sec(ms);                     // accumulation time (sec)

    const double ref_frequency_hz = REF_FREQ;  //Reference frequency in Hz at which fluxes are measured
    
    double efficiency = EFFICIENCY;     // system efficiency
    double SEFD = SEFD_JVLA;      // System Equivalent Flux Density (in micro-Jy) of each SKA1 antenna
    
    unsigned int num_baselines = num_stations * (num_stations - 1) / 2;
    if (rank==0)
    {
        cout << "Number baselines: " << num_baselines << endl;
        cout << "Number of channels: " << num_channels << endl;
        cout << "Channels bandwidth (Hz): " << channel_bandwidth_hz << endl;
        cout << "Reference frequency (Hz): " << ref_frequency_hz << endl;
        cout << "Starting frequency (Hz): " << freq_start_hz << endl;
        cout << "Accumulation time (sec): " << time_acc << endl;
        cout << "Number of rows: " << num_rows << endl;
    }
    
    double sizeGbytes, totGbytes = 0.;
    
    // Allocate and read uv coordinates ------------------------------------------------------------------------------
    unsigned long int num_coords = ms_num_rows(ms);
    double* uu_metres = new double[num_coords];
    double* vv_metres = new double[num_coords];
    double* ww_metres = new double[num_coords];
    sizeGbytes = 3*num_coords*sizeof(double)/((double)(1024*1024*1024));
    cout << "rank " << rank << ": allocated original coordinates: " << sizeGbytes  << " GB" << endl;
    totGbytes += sizeGbytes;
    
    int status = 0;
    double len = ms_read_coords(ms,0,num_coords,uu_metres,vv_metres,ww_metres,&status);
    ms_close(ms);
    
    // Allocate and read Data visibilities
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
    memset(visData, 0, num_vis*sizeof(complexd));

    double *sigma2_vis;
    try
    {
        sigma2_vis = new double[num_vis];
        sizeGbytes = num_vis*sizeof(double)/((double)(1024*1024*1024));
        cout << "allocated sigma2 visibilities: " << num_vis << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << " bad_alloc caught: " << ba.what() << '\n';
    }
    double sigma2 = (SEFD*SEFD)/(2.*time_acc*channel_bandwidth_hz*efficiency*efficiency);
    for (unsigned long int i = 0; i<num_vis; i++)
        sigma2_vis[i] = sigma2; // visibility noise variance

    if (rank==0) cout << "sigma_vis  = " << sqrt(sigma2) << " muJy" << endl;

    // Pre-compute wavenumber and spectral factor for each channel ---------------------------------------------------------------------
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
 
    // define steps in galaxy scalelength (Ro in ARCSEC) ------------------------------
    double Rmin = RMIN;
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
    if (rank==0) cout << numR-1 << " samples in galaxy scale-length, " << Rmin << " < r0 < " << Rmax << " arcsec" << endl;


    // Read galaxy catalogue --------------------------------------------------------------------------------------------------------------------------
    unsigned long int my_gal = atof(argv[3]);

    double *gflux = new double[my_gal];
    double *l = new double[my_gal];
    double *m = new double[my_gal];
    double *ge1 = new double[my_gal];
    double *ge2 = new double[my_gal];
    double *gscale = new double[my_gal];
    double *SNR_vis = new double[my_gal];

    read_catalog(my_gal, argv[2],gflux,gscale,ge1,ge2,l,m,SNR_vis);

    //setup random number generator
    const gsl_rng_type * G;
    gsl_rng * gen;
    G = gsl_rng_mt19937;  // Mersenne Twister
    gen = gsl_rng_alloc(G);
    
    unsigned long int seed = random_seed();
    gsl_rng_set(gen,seed);

    // Faceting uv coordinates ----------------------------------------------------------------------------------------
#ifdef FACET
    len = len*wavenumbers[num_channels-1]/(2*PI);
    int facet = facet_size(RMAX,len);
    unsigned long int ncells = facet*facet;
    unsigned long int* count = new unsigned long int[ncells];

    unsigned long int facet_ncoords = evaluate_uv_grid_size(len,wavenumbers, num_channels,num_coords, uu_metres, vv_metres, facet, count);
    double* facet_u = new double[facet_ncoords];
    double* facet_v = new double[facet_ncoords];
    sizeGbytes = (2*facet_ncoords*sizeof(double)+ncells*sizeof(unsigned long int))/((double)(1024*1024*1024));
    cout << "rank " << rank << ": allocated grid coordinates and array counter: " << sizeGbytes  << " GB" << endl;
    totGbytes += sizeGbytes;
   
    unsigned long int facet_nvis = facet_ncoords;
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
    if (rank==0) cout << "grid length = " << 2*len << ", grid size = " << facet << endl;  

#endif
    // Allocate Model Visibilities ---------------------------------------------------------------------------------------------
    int num_models = numR-1;

#ifdef FACET
    double* visMod;
    try
    {
        unsigned long int model_ncoords = facet_ncoords;
        visMod = new double[num_models*model_ncoords];
        sizeGbytes = num_models*model_ncoords*sizeof(complexd)/((double)(1024*1024*1024));
#else
        complexd* visMod;
        try
        {
            unsigned long int model_ncoords = num_coords;
            visMod = new complexd[num_models*model_ncoords*num_channels];
            sizeGbytes = num_models*model_ncoords*num_channels*sizeof(complexd)/((double)(1024*1024*1024));
#endif
        cout << "rank " << rank << ": allocated models: num_models= " << num_models << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "rank " << rank << ": bad_alloc caught: " << ba.what() << '\n';
    }
    if (rank==0) cout << "Total Visibilities GBytes: " << totGbytes << endl;    


    // shear to be applied
    double g1 = atof(argv[4]);
    double g2 = atof(argv[5]);
    double* sigmab = new double[num_baselines];
    
    // Set function to be minimized
    likelihood_params par;
    par.numr = numR;
    par.ro = Ro;
    par.rprior = rprior;
    par.mod = visMod;
#ifdef FACET
    par.ncoords = facet_ncoords;
    par.uu = facet_u;
    par.vv = facet_v;
    par.data = facet_visData;
    par.sigma2 = facet_sigma2;
#else
    par.ncoords = num_coords;
    par.uu = uu_metres;
    par.vv = vv_metres;
    par.ww = ww_metres;
    par.data = visData;
    par.count = 0;
    par.sigma2 = sigma2_vis; // visibility noise variance
    par.nchannels = num_channels;
    par.band_factor = channel_bandwidth_hz*PI/C0;
    par.acc_time = time_acc;
    par.spec = spec;
    par.wavenumbers = wavenumbers; // wavenumbers for the model
#endif

    FILE *pFile;
    char filename[100];
    sprintf(filename,"ellipticities-%dch.txt",num_channels);
    pFile = fopen(filename,"w");
    fprintf(pFile, "flux | e1 | m_e1 | err1 | e2 | m_e2 | err2 | 1D var | SNR |   l  |  m  | \n");

    int bad = 0;
    int np_max = NP_MAX;  // min number of sampling points with likelihood above 5%ML
    
    gsl_multimin_function minex_func;
    minex_func.n = 2;
    minex_func.f = f_likelihood;
    minex_func.params = &par;
    
    // use Simplex algorithm of Nelder and Mead provided by the GLS library to minimize -log(likelihood)
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    
    gsl_multimin_fminimizer *s = 0;
    gsl_vector *sx, *x;
    x = gsl_vector_alloc (2);
    sx = gsl_vector_alloc (2);
    s = gsl_multimin_fminimizer_alloc (T, 2);
    
#ifdef USE_MPI
    double data_time = 0.;
    double fitting_time = 0.;
    double start_data,end_data,start_fitting,end_fitting;
#else
    double data_time = 0;
    double fitting_time = 0;
    long long start_data,end_data,start_fitting,end_fitting;
#endif
    
    double l0,m0,ee1,ee2,den;
    complexd z1,z2;
  
    for (unsigned long int g=0; g<my_gal; g++)
    {
       // set log(prior) for scalelength
          double mu = scale_mean(gflux[g]);
          for (int nRo=1; nRo<numR; nRo++)
              rprior[nRo] = rfunc(mu,R_STD,Ro[nRo]);
          double R_mu = exp(mu);
 
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
    
#ifdef USE_MPI
        start_data = MPI_Wtime();
#else
        start_data = current_timestamp();
#endif
    
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (unsigned int ch = 0; ch < num_channels; ch++)
        {
           // generate galaxy visibilities
           unsigned long int ch_vis = ch*num_coords;
           data_galaxy_visibilities(spec[ch], wavenumbers[ch], par.band_factor, time_acc, ee1, ee2, gscale[g],
                                      gflux[g], l[g], m[g], num_coords, uu_metres, vv_metres, ww_metres, &(visData[ch_vis]));
 
           // compute signal-to-noise
           double SNR_ch = 0.;
           for (unsigned long int vs = ch_vis; vs < ch_vis+num_coords; vs++)
               SNR_ch += visData[vs].real*visData[vs].real + visData[vs].imag*visData[vs].imag;
           SNR_ch /= sigma2;
#ifdef _OPENMP
#pragma omp critical
#endif
          {
             SNR_vis[g] += SNR_ch;
             add_system_noise(gen, num_coords, &(visData[ch_vis]), sigmab);
          }
            
#ifdef FACET
            // Phase shift data visibilities (to be done after gridding because real data will be gridded)
            data_visibilities_phase_shift(wavenumbers[ch], l0, m0, num_coords, uu_metres, vv_metres, ww_metres, &(visData[ch_vis]));
        }
        // gridding visibilities
        facet = facet_size(R_mu,len);
        gridding_visibilities(wavenumbers,num_channels,num_coords,uu_metres,vv_metres,visData,sigma2_vis,len,facet,facet_visData,facet_sigma2,count);
        par.ncoords = evaluate_facet_coords(par.uu, par.vv, len, facet, count);
#else
        }
        par.l0 = l0;
        par.m0 = m0;
#endif
        cout << "SNR: " << sqrt(SNR_vis[g]) << endl;
          
#ifdef USE_MPI
        end_data = MPI_Wtime();
        data_time += end_data - start_data;
        start_fitting = MPI_Wtime();
#else
        end_data = current_timestamp();
        data_time += (double)(end_data - start_data)/1000.;
        start_fitting = current_timestamp();
#endif
        
        // Model fitting  -----------------------------------------------------
        // Search for the maximum posterior to find starting ellipticity points
        int iter = 0;
        int status;
        double size;
 
        double start_e1 = 0.;
        double start_e2 = 0.;
        
        // Search for the maximum likelihood
        gsl_vector_set (x, 0, start_e1);
        gsl_vector_set (x, 1, start_e2);
        gsl_vector_set_all (sx, 0.1);
        
        gsl_multimin_fminimizer_set (s, &minex_func, x, sx);
        iter = 0;

        do
        {
                iter++;
                status = gsl_multimin_fminimizer_iterate(s);
            
                if (status) break;
            
                size = gsl_multimin_fminimizer_size(s);
                status = gsl_multimin_test_size (size, 1e-3);
        }
        while (status == GSL_CONTINUE && iter < 50 && s->fval < 0.);
        
        double mes_e1, mes_e2, maxL;
        mes_e1 = gsl_vector_get(s->x, 0);
        mes_e2 = gsl_vector_get(s->x, 1);
        maxL= -s->fval;
        cout << "rank:" << rank << " n. " << g << " flux = " << gflux[g] << " scalelength = " << gscale[g] << " position [arcsec] (" << l0/(ARCS2RAD) << "," << m0/(ARCS2RAD) << "): Maximum log likelihood = " << maxL << " n.iter = " << iter << " for e = " << mes_e1 << "," << mes_e2 <<  "  original e = " << ge1[g] << "," << ge2[g] << endl;
            
        // Likelihood sampling to compute mean and variance
        double var_e1, var_e2, oneDimvar;
        likelihood_sampling(&mes_e1, &mes_e2, maxL, &par, np_max, &var_e1, &var_e2, &oneDimvar);

#ifdef USE_MPI
        end_fitting = MPI_Wtime();
        fitting_time += end_fitting - start_fitting;
#else
        end_fitting = current_timestamp();
        fitting_time += (double)(end_fitting - start_fitting)/1000.;
#endif

        fprintf(pFile, "%f | %f | %f | %f | %f | %f | %f | %f | %f | %f | %f \n",gflux[g],ge1[g],mes_e1,sqrt(var_e1), ge2[g],mes_e2,sqrt(var_e2),oneDimvar,sqrt(SNR_vis[g]),l0/(ARCS2RAD),m0/(ARCS2RAD));
        if (var_e1 < 1e-4 || var_e2 < 1e-4 || oneDimvar < 1e-4)
        {
              cout << "ERROR likelihood sampling!" << endl;
              bad++;
        }
    }
    
    gsl_vector_free(x);
    gsl_vector_free(sx);
    gsl_multimin_fminimizer_free(s);
    
    gsl_rng_free(gen);
 
#ifdef USE_MPI
    double end_tot = MPI_Wtime();
#else
    long long end_tot = current_timestamp();
#endif

    cout << "rank : " << rank << "removed " << bad << " bad data galaxies" << endl << endl;
#ifdef USE_MPI
    double total_time = end_tot - start_tot;
#else
    double total_time = (double)(end_tot - start_tot)/1000.;
#endif
    cout << "rank: " << rank << " set up time (sec): " << total_time - data_time - fitting_time << endl;
    cout << "rank: " << rank << " data generation time (sec): " << data_time << endl;
    cout << "rank: " << rank << " data fitting computation time (sec): " << fitting_time << endl;
    cout << "rank: " << rank << " Total time (sec): " << total_time << endl;
        
    if (pFile != 0) fclose(pFile);

    
    // free memory ----------------------------------------------------------------------------------------------------------------
    delete[] visMod;
    delete[] visData;
    delete[] Ro;
    delete[] rprior;
    delete[] sigmab;
    delete[] ge1;
    delete[] ge2;
    delete[] gflux;
    delete[] gscale;
    delete[] l;
    delete[] m;
    delete[] SNR_vis;
    delete[] uu_metres;
    delete[] vv_metres;
    delete[] wavenumbers;
    delete[] spec;
    delete[] sigma2_vis;
#ifdef FACET
    delete[] facet_u;
    delete[] facet_v;
    delete[] facet_visData;
    delete[] facet_sigma2;
    delete[] count;
#endif

#ifdef USE_MPI
      MPI_Finalize() ;
#endif
    return 0;
}
