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
//  argv[1]  name of the file contaning u coordinates
//  argv[2]  name of the file contaning v coordinates
//  argv[3]  number of galaxies
//  argv[4]  applied shear 1st component
//  argv[5]  applied shear 2nd component
//  argv[6]  minimum galaxy flux in muJy
//  argv[7]  maximum galaxy flux in muJy

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
#include "read_coordinates.h"
#include "utils.h"
#include "galaxy_visibilities.h"
#include "likelihood.h"
#include "marginalise_r.h"
#include "distributions.h"
#include "evaluate_uv_grid.h"
#include "utils.h"
#include "default_params.h"

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

    if (argc < 7)
    {
        cout << "ERROR: parameter missing!" << endl;
        cout << "usage: RadioLensfit.x <filename u-coord> <filename v-coord> <nge> <shear1> <shear2> <min flux> <max flux>" << endl;
        cout << "corresponding to ng = nge*10 galaxies, g1 = shear1, g2 = shear2 and minimum galaxy flux [muJy]" << endl;
        exit(EXIT_FAILURE);
    }

    // Initialise input data
    unsigned int num_stations = 197;      // Number of stations
    unsigned int num_channels = 1;        // Number of frequency channels
    unsigned int num_times = 480; //1920; // Number of time samples
    double freq_start_hz = 950e+6;        // Start Frequency, in Hz
    //double freq_inc_hz = 1e+6;          // Frequency increment, in Hz
    double full_bandwidth_hz = 240e+6;    // Frequency total bandwidth, in Hz
    double ref_frequency_hz = 1.4e+9;  //Reference frequency in Hz at which fluxes are measured
    int time_acc = 60; //15;     // accumulation time
    double efficiency = 0.9;     // system efficiency
    double SEFD = SEFD_SKA;      // System Equivalent Flux Density (in micro-Jy) of each SKA1 antenna
    
    unsigned int num_baselines = num_stations * (num_stations - 1) / 2;
    double channel_bandwidth_hz = full_bandwidth_hz/num_channels; // Frequency channel bandwidth, in Hz
    if (rank==0)
    {
        cout << "Number baselines: " << num_baselines << endl;
        cout << "Number of time snapshots: " << num_times << endl;
        cout << "Number of channels: " << num_channels << endl;
        cout << "Channels bandwidth (Hz): " << channel_bandwidth_hz << endl;
        cout << "Reference frequency (Hz): " << ref_frequency_hz << endl;
        cout << "Starting frequency (Hz): " << freq_start_hz << endl;
        cout << "Accumulation time (sec): " << time_acc << endl;
    }
    
    double sizeGbytes, totGbytes = 0.;
    double fov = 3600*ARCS2RAD; //1.22*C0/(freq_start_hz*diameter);  // 1 degree field of view in RAD
    printf("field of view: %e [rad] %f [arcsec] \n",fov,fov/(ARCS2RAD));
    
    // Allocate and read uv coordinates ------------------------------------------------------------------------------
    // coordinates in the file are ordered as nbaselines x ntimes
    // coordinates in the array will be ordered as ntimes x nbaselines
    unsigned long int num_coords = num_times * num_baselines;
    double* uu_metres = new double[num_coords];
    double* vv_metres = new double[num_coords];
    sizeGbytes = 2*num_coords*sizeof(double)/((double)(1024*1024*1024));
    cout << "rank " << rank << ": allocated original coordinates: " << sizeGbytes  << " GB" << endl;
    totGbytes += sizeGbytes;
    
    double len, threshold = 0.; //only uv-points above this threshold [metres] will be used
    // read only the baselines above the threshold and update their number
    double maxB = read_coord_ska(argv[1], argv[2], num_times, &num_baselines, uu_metres, vv_metres, threshold, &len);
    if (rank==0)
    {
        cout << "max baseline: " << maxB << endl;
        cout << "New number baselines: " << num_baselines << endl;
        cout << "resolution: " << (1.22*C0/(freq_start_hz*maxB))/(ARCS2RAD) << " arcsec" << endl;
    }
    num_coords = num_times * num_baselines;
    
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
 
    
    // Allocate Data Visibilities ------------------------------------------------------------------------------------------
    unsigned long int num_rawvis  = (unsigned long int) num_channels * num_coords;
    if (rank==0)
        cout << "Num visibilities: " << num_rawvis << endl;
    
    complexd* visData;
    try
    {
        visData = new complexd[num_rawvis];
        sizeGbytes = num_rawvis*sizeof(complexd)/((double)(1024*1024*1024));
        cout << "rank " << rank << ": allocated original visibilities: " << num_rawvis << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "rank " << rank << ": bad_alloc caught: " << ba.what() << '\n';
    }
    
    memset(visData, 0, num_rawvis*sizeof(complexd));
    
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

    
    // Faceting uv coordinates ----------------------------------------------------------------------------------------
    double Fmin = atof(argv[6]);
    double Fmax = atof(argv[7]);

#ifdef FACET
    double Rmu_max = scale_mean(Fmax);
    Rmu_max = exp(Rmu_max); 
    int facet = facet_size(Rmu_max,len);
    unsigned long int ncells = facet*facet;
    unsigned long int* count = new unsigned long int[ncells];
    
    unsigned long int facet_ncoords = evaluate_max_uv_grid(len, num_coords, uu_metres, vv_metres, facet, count);
    //unsigned long int facet_ncoords = evaluate_max_uv_circular_grid(len,num_coords, uu_metres, vv_metres, facet, count);
    double* facet_u = new double[facet_ncoords];
    double* facet_v = new double[facet_ncoords];
    sizeGbytes = (2*facet_ncoords*sizeof(double)+ncells*sizeof(unsigned long int))/((double)(1024*1024*1024));
    cout << "rank " << rank << ": allocated grid coordinates and array counter: " << sizeGbytes  << " GB" << endl;
    totGbytes += sizeGbytes;
   
    unsigned long int facet_nvis = num_channels*facet_ncoords;
    complexd* facet_visData;
    try
    {
       facet_visData = new complexd[facet_nvis];
       sizeGbytes = facet_nvis*sizeof(complexd)/((double)(1024*1024*1024));
       cout << "rank " << rank << ": allocated gridded visibilities: " << facet_nvis << ", size = " << sizeGbytes  << " GB" << endl;
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
        visMod = new double[num_models*model_ncoords*num_channels];
        sizeGbytes = num_models*model_ncoords*num_channels*sizeof(complexd)/((double)(1024*1024*1024));
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


    // Generate fake galaxies ---------------------------------------------------------------------------------------------------
    // shear to be applied
    double g1 = atof(argv[4]);
    double g2 = atof(argv[5]);
    double* sigmab = new double[num_baselines];
    
    int NP = 1;    // 2NP = number of sampled orientations (points on the circle of radius |e|) for each ellipticity module
    int nge = atoi(argv[3]);
    
    //setup random number generator
    const gsl_rng_type * G;
    gsl_rng * gen;
    G = gsl_rng_mt19937;  // Mersenne Twister
    gen = gsl_rng_alloc(G);
    
    unsigned long int seed = 4151155891; //random_seed();
    gsl_rng_set(gen,seed);

#ifdef USE_MPI
    int my_gal = nge/nprocs;
    int rem = nge%nprocs;
    if (rem)
        if (rank < rem) my_gal++;
#else
    int my_gal = nge;
#endif
    int mygalaxies = my_gal*2*NP;

    // generate flux values
    const double beta = -1.34; // flux prior: S^beta
    double *gflux = new double[my_gal];
    double *gflux2 = new double[my_gal];
    generate_random_data(gen,my_gal,gflux2,Fmin,Fmax,flux_CDF,beta);
    
    // sort flux values, so that to generate a population ordered by flux and therefore fitting sources by decreasing flux order
    gsl_sort(gflux2,1,my_gal); // sorting ascending order
    for (unsigned long int i=0; i<my_gal; i++)
        gflux[i] = gflux2[my_gal-i-1];
    
    delete[] gflux2;
    
    // generate scalelength
    double *gscale = new double[my_gal];
    for (unsigned long int g=0; g<my_gal; g++)
    {
        double mu = scale_mean(gflux[g]); //power law relation between flux and scalelength
        generate_random_data(gen,1,&(gscale[g]),Rmin,Rmax,r_CDF,mu);
    }
    
    // generate ellipticities
    double *ge1 = new double[mygalaxies];
    double *ge2 = new double[mygalaxies];
    generate_ellipticity(gen,my_gal,NP,ge1,ge2);
   
    // uniformly random positions in RAD in a disk area
    //  http://mathworld.wolfram.com/DiskPointPicking.html
    double *l = new double[mygalaxies];
    double *m = new double[mygalaxies];
    //double radius,orient;

    for (unsigned long int gal=0; gal<mygalaxies; gal++)
    {
       //  radius = sqrt(gsl_rng_uniform(gen))*0.5*fov;
       //  orient = gsl_rng_uniform(gen)*2*PI;

       l[gal] = gsl_rng_uniform(gen)*fov - fov/2.; //radius*cos(orient);
       m[gal] = gsl_rng_uniform(gen)*fov - fov/2.; //radius*sin(orient);
    }
  
    double *SNR_vis = new double[mygalaxies];
    
    // Set function to be minimized
    likelihood_params par;
    par.numr = numR;
    par.ro = Ro;
    par.rprior = rprior;
#ifdef FACET
    par.ncoords = facet_ncoords;
    par.uu = facet_u;
    par.vv = facet_v;
    par.data = facet_visData;
    par.count = count;
#else
    par.ncoords = num_coords;
    par.uu = uu_metres;
    par.vv = vv_metres;
    par.data = visData;
    par.count = 0;
#endif
    par.nchannels = num_channels;
    par.nbaselines = num_baselines;
    par.band_factor = channel_bandwidth_hz*PI/C0;
    par.acc_time = time_acc;
    par.spec = spec;
    par.wavenumbers = wavenumbers; // wavenumbers for the model
    par.mod = visMod;
    par.sigma = (SEFD*SEFD)/(2.*time_acc*channel_bandwidth_hz*efficiency*efficiency); // visibility noise variance

    if (rank==0) cout << "sigma_vis  = " << sqrt(par.sigma) << " muJy" << endl;
    for (unsigned int b=0; b<num_baselines; b++) sigmab[b] = sqrt(par.sigma);
    
    
    FILE *pFile;
    char filename[100];
    sprintf(filename,"ellipticities%d.txt",rank);
    pFile = fopen(filename,"w");
    fprintf(pFile, "flux | e1 | m_e1 | err1 | e2 | m_e2 | err2 | 1D var | SNR |   l  |  m  | \n");

    int bad = 0;
    unsigned long int gal = 0;
    int np_max = 30;  // min number of sampling points with likelihood above 5%ML
    
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
 
      for (int ang=0; ang<2*NP; ang++)
      {
        l0 = l[gal];
        m0 = m[gal];
          
        ee1 = ge1[gal];
        ee2 = ge2[gal];

        // apply shear g
        z1.real = ee1+g1;            z1.imag = ee2+g2;         // z1 = e+g
        z2.real = 1.+ee1*g1+ee2*g2;  z2.imag = ee2*g1-ee1*g2;  // z2 = 1+conj(g)*e
        den = z2.real*z2.real+z2.imag*z2.imag;
        ee1 = (z1.real*z2.real + z1.imag*z2.imag)/den;
        ee2 = (z1.imag*z2.real - z1.real*z2.imag)/den;        // e = z1/z2
        
        ge1[gal] = ee1;
        ge2[gal] = ee2;
        
        SNR_vis[gal] = 0.;
    
#ifdef USE_MPI
        start_data = MPI_Wtime();
#else
        start_data = current_timestamp();
#endif
    
#ifdef FACET
       facet = facet_size(R_mu,len);
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (unsigned int ch = 0; ch < num_channels; ch++)
        {
           // generate galaxy visibilities
           unsigned long int ch_vis = ch*num_coords;
           data_galaxy_visibilities2D(spec[ch], wavenumbers[ch], par.band_factor, time_acc, ee1, ee2, gscale[g],
                                      gflux[g], l[gal], m[gal], num_coords, uu_metres, vv_metres, &(visData[ch_vis]));
 
           // compute signal-to-noise
           double SNR_ch = 0.;
           for (unsigned long int vs = ch_vis; vs < ch_vis+num_coords; vs++)
               SNR_ch += visData[vs].real*visData[vs].real + visData[vs].imag*visData[vs].imag;
           SNR_ch /= par.sigma;
#ifdef _OPENMP
#pragma omp critical
#endif
          {
             SNR_vis[gal] += SNR_ch;
             add_system_noise(gen, num_coords, &(visData[ch_vis]), sigmab);
          }
            
#ifdef FACET
            // Phase shift data visibilities (to be done after gridding because real data will be gridded)
            data_visibilities_phase_shift2D(wavenumbers[ch], l0, m0, num_coords, uu_metres, vv_metres, &(visData[ch_vis]));

            // gridding visibilities
            unsigned int ch_visfacet = ch*facet_ncoords;
            par.ncoords = gridding_visibilities(num_coords,facet_u,facet_v,uu_metres,vv_metres,&(visData[ch_vis]),len,facet,&(facet_visData[ch_visfacet]),count,ch);
            //par.ncoords = circular_gridding_visibilities(num_coords,facet_u,facet_v,uu_metres,vv_metres,&(visData[ch_vis]),len,facet,&(facet_visData[ch_visfacet]),count,ch);
#else
            par.l0 = l0;
            par.m0 = m0;
#endif
        }
        cout << "SNR: " << sqrt(SNR_vis[gal]) << endl;
          
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
        cout << "rank:" << rank << " n. " << gal << " flux = " << gflux[g] << " scalelength = " << gscale[g] << " position [arcsec] (" << l0/(ARCS2RAD) << "," << m0/(ARCS2RAD) << "): Maximum log likelihood = " << maxL << " n.iter = " << iter << " for e = " << mes_e1 << "," << mes_e2 <<  "  original e = " << ge1[gal] << "," << ge2[gal] << endl;
            
        // Likelihood sampling to compute mean and variance
        double var_e1, var_e2, cov_e;
        int error = likelihood_sampling(rank,&mes_e1, &mes_e2, maxL, &par, np_max, &var_e1, &var_e2, &cov_e);
        double oneDimvar = sqrt(var_e1*var_e2-cov_e*cov_e);

#ifdef USE_MPI
        end_fitting = MPI_Wtime();
        fitting_time += end_fitting - start_fitting;
#else
        end_fitting = current_timestamp();
        fitting_time += (double)(end_fitting - start_fitting)/1000.;
#endif

        fprintf(pFile, "%f | %f | %f | %f | %f | %f | %f | %f | %f | %f | %f \n",gflux[g],ge1[gal],mes_e1,sqrt(var_e1), ge2[gal],mes_e2,sqrt(var_e2),oneDimvar,sqrt(SNR_vis[gal]),l0/(ARCS2RAD),m0/(ARCS2RAD));
        if (error)
        {
              cout << "ERROR likelihood sampling!" << endl;
              bad++;
        }
        gal++;
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
#ifdef FACET
    delete[] facet_u;
    delete[] facet_v;
    delete[] facet_visData;
#endif

#ifdef USE_MPI
      MPI_Finalize() ;
#endif
    return 0;
}
