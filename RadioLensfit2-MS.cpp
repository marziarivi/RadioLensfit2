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


//  RadioLensfit.cpp
//
//  argv[1]  filename Measurement Set
//  argv[2]  filename source catalog containing l,m,flux
//  argv[3]  number of sources
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

#include "datatype.h"
#include "utils.h"
#include "measurement_set.h"
#include "data_simulation.h"
#include "read_catalog.h"
#include "galaxy_fitting.h"
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
    
    if (argc < 5)
    {
        cout << "ERROR: parameter missing!" << endl;
        cout << "usage: RadioLensfit.x <filename MS> <min flux> <shear1> <shear2> " << endl;
        exit(EXIT_FAILURE);
    }

    // Read Measurement Set --------------------------------------------------------------------------------------------------------------------------------------------------------------------
    RL_MeasurementSet* ms = ms_open(argv[1]);

    double RA = ms_phase_centre_ra_rad(ms);                 // Phase Centre coordinates
    double Dec = ms_phase_centre_dec_rad(ms);   
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
    if (rank==0)
    {
        cout << "Number baselines: " << num_baselines << endl;
        cout << "Number of channels: " << num_channels << endl;
        cout << "Channels bandwidth (Hz): " << channel_bandwidth_hz << endl;
        cout << "Reference frequency (Hz): " << ref_frequency_hz << endl;
        cout << "Starting frequency (Hz): " << freq_start_hz << endl;
        cout << "Accumulation time (sec): " << time_acc << endl;
    }
    
    double sizeGbytes, totGbytes = 0.;
    double fov_eff_arcmin = 60.;
    double fov_eff = fov_eff_arcmin*60.*ARCS2RAD; //1.22*C0/(freq_start_hz*diameter);  // 1 degree field of view in RAD
    printf("field of view: %e [rad] %f [arcsec] \n",fov_eff,fov_eff/(ARCS2RAD));
    
    // Allocate and read uv coordinates 
    unsigned long int num_coords = ms_num_rows(ms);
    double* uu_metres = new double[num_coords];
    double* vv_metres = new double[num_coords];
    double* ww_metres = new double[num_coords];
    sizeGbytes = 2*num_coords*sizeof(double)/((double)(1024*1024*1024));
    cout << "rank " << rank << ": allocated original coordinates: " << sizeGbytes  << " GB" << endl;
    totGbytes += sizeGbytes;
    
    int status;
    double len = ms_read_coords(ms,0,num_coords,uu_metres,vv_metres,ww_metres,&status);
    
    ms_close(ms);
    delete[] ww_metres;
  
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
 
    
    // Allocate Galaxy and Raw Data Visibilities -----------------------------------------------------------------------------------------------------------------------
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
    
    // Allocate Sky Model Visibilities ------------------------------------------------------------------------------------------
    unsigned long int num_vis  = (unsigned long int) num_channels * num_coords;
    complexd *visSkyMod;
    try
    {
        visSkyMod = new complexd[num_vis];
        sizeGbytes = num_vis*sizeof(complexd)/((double)(1024*1024*1024));
        cout << "rank " << rank << ": allocated sky model visibilities: " << num_vis << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "rank " << rank << ": bad_alloc caught: " << ba.what() << '\n';
    }
    
    memset(visSkyMod, 0, num_vis*sizeof(complexd));

    // define steps in galaxy scalelength (Ro in ARCSEC) --------------------------------------------------------------------------
    double Rmin = 0.3;
    double Rmax = 3.5;
    int numR = 29;
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
    
    //---------------------------------------------------------------------------------------------------------------------------------------------------
    // Read galaxy catalogue --------------------------------------------------------------------------------------------------------------------------

    unsigned long int nge = atof(argv[3]);
    
    double *gflux = new double[nge];
    double *gscale = new double[nge];
    double *ge1 = new double[nge];
    double *ge2 = new double[nge];
    double *l = new double[nge];
    double *m = new double[nge];
    double *SNR_vis = new double[nge];
 
    int NP = 1;    // 2NP = number of sampled orientations (points on the circle of radius |e|) for each ellipticity module
    unsigned long int mygalaxies = read_catalog(nge, argv[2], gflux, gscale,ge1,ge2,l,m);
    cout << "num gal: " << mygalaxies << endl;
    
#ifdef USE_MPI
    double data_time = 0.;
    double fitting_time = 0.;
    double start_data,end_data,start_fitting,end_fitting;
    start_data = MPI_Wtime();
#else
    double data_time = 0;
    double fitting_time = 0;
    long long start_data,end_data,start_fitting,end_fitting;
    start_data = current_timestamp();
#endif
    
    // Visibilities Simulation --------------------------------------------------------------------------------------------------------------------------
    double g1 = atof(argv[4]);  // shear to be applied
    double g2 = atof(argv[5]);
    
    double sigma = (SEFD_SKA*SEFD_SKA)/(2.*time_acc*channel_bandwidth_hz*efficiency*efficiency); // visibility noise variance
    if (rank==0) cout << "sigma_vis  = " << sqrt(sigma) << " muJy" << endl;
    
    data_simulation(freq_start_hz,ref_frequency_hz, wavenumbers, spec, channel_bandwidth_hz, time_acc, num_channels, num_baselines,
                    sigma, mygalaxies, g1, g2, ge1, ge2, gflux, gscale, l, m, SNR_vis, num_coords, uu_metres, vv_metres,
                    visGal, visData);

    sky_model(freq_start_hz,ref_frequency_hz, wavenumbers, spec, channel_bandwidth_hz, time_acc, num_channels, num_baselines,
                    mygalaxies, gflux, l, m, num_coords, uu_metres, vv_metres, visGal, visSkyMod);
    
#ifdef USE_MPI
    end_data = MPI_Wtime();
    data_time = end_data - start_data;
#else
    end_data = current_timestamp();
    data_time = (double)(end_data - start_data)/1000.;
#endif
    
   //---------------------------------------------------------------------------------------------------------------------------------------------------
   //---------------------------------------------------------------------------------------------------------------------------------------------------
    
    int num_models = numR-1;
    
    // Set likelihood computation parameters
    likelihood_params par;
    par.numr = numR;
    par.ro = Ro;
    par.rprior = rprior;
    par.nchannels = num_channels;
    par.nbaselines = num_baselines;
    par.band_factor = channel_bandwidth_hz*PI/C0;
    par.acc_time = time_acc;
    par.spec = spec;
    par.wavenumbers = wavenumbers; // wavenumbers for the model
    par.sigma = sigma; // visibility noise variance
    
#ifdef FACET
    // lower limit flux
    double threshold_flux[10];
    threshold_flux[0] = 150.;
    threshold_flux[1] = 100.;
    threshold_flux[2] = 80.;
    threshold_flux[3] = 60.;
    threshold_flux[4] = 40.;
    threshold_flux[5] = 20.;
    threshold_flux[6] = 10.;
    
    // corresponding facet size
    int facet_size[10];
    facet_size[0] = 600;
    facet_size[1] = 550;
    facet_size[2] = 500;
    facet_size[3] = 460;
    facet_size[4] = 420;
    facet_size[5] = 350;
    facet_size[6] = 280;
    
    // Faceting uv coordinates ----------------------------------------------------------------------------------------
    int ind = 0;
    while (gflux[0]<threshold_flux[ind]) ind++;
    int facet = facet_size[ind];
    double* facet_u = 0;
    double* facet_v = 0;
    unsigned long int ncells = facet*facet;
    unsigned long int* count = new unsigned long int[ncells];
    
    unsigned long int facet_ncoords = evaluate_uv_grid(len, num_coords, uu_metres, vv_metres, facet, &facet_u, &facet_v, count);
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
    if (rank==0)
    {
        cout << "grid length = " << 2*len << ", grid size = " << facet << endl;
        cout << num_models << " samples in galaxy scale-length, " << Rmin << " < r0 < " << Rmax << " arcsec" << endl;
    }
    
    // Allocate Facet Model Visibilities ------------------------------------------------------------------------------------------
    double* visMod;
    try
    {
        unsigned long int model_ncoords = facet_ncoords;
        visMod = new double[num_models*model_ncoords*num_channels];
        sizeGbytes = num_models*model_ncoords*num_channels*sizeof(double)/((double)(1024*1024*1024));
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
    par.data = facet_visData;
    par.count = count;
    par.mod = visMod;
    
#else
    complexd* visMod;
    try
    {
        unsigned long int model_ncoords = num_coords;
        visMod = new complexd[num_models*model_ncoords*num_channels];
        sizeGbytes = num_models*model_ncoords*num_channels*sizeof(complexd)/((double)(1024*1024*1024));
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
    par.data = visGal;
    par.count = 0.;
    par.mod = visMod;
#endif
    
#ifdef USE_MPI
    start_fitting = MPI_Wtime();
#else
    start_fitting = current_timestamp();
#endif
    
    if (rank==0) cout << "Total Visibilities GBytes per rank: " << totGbytes << endl;

    // Data Fitting -----------------------------------------------------------------------------------------------------------------------------------
    FILE *pFile;
    char filename[100];
    sprintf(filename,"ellipticities%d.txt",rank);
    pFile = fopen(filename,"w");
    fprintf(pFile, "flux | scale | e1 | m_e1 | err1 | e2 | m_e2 | err2 | 1D var | SNR |   l  |  m  | \n");
    
    double l0,m0;
    unsigned long int bad_list[mygalaxies];
    int bad = 0;
    
    for (unsigned long int g=0; g<mygalaxies; g++)
    {

      if (SNR_vis[g] >= 10.)
      {
#ifdef FACET
        if (gflux[g] < threshold_flux[ind])
        {
            ind++; facet = facet_size[ind];
            par.ncoords = evaluate_uv_grid(len, num_coords, uu_metres, vv_metres, facet, &facet_u, &facet_v, count);
            if (rank==0) cout << " new facet size: " << facet << endl;
        }
#endif
        // set log(prior) for scalelength
        double mu = scale_mean(gflux[g]);
        for (int nRo=1; nRo<numR; nRo++)
           rprior[nRo] = rfunc(mu,scale_std,Ro[nRo]);
        double R_mu = exp(mu);
        
          l0 = l[g];  m0 = m[g];
          
#ifdef FACET
          source_extraction(l0, m0, gflux[g], exp(mu), 0., 0., &par, visSkyMod, visData, visGal, num_coords, uu_metres, vv_metres, facet, len);
#else
          source_extraction(l0, m0, gflux[g], exp(mu), 0., 0., &par, visSkyMod, visData, visGal, num_coords, uu_metres, vv_metres);
#endif
 
          double mes_e1, mes_e2, maxL;
          double var_e1, var_e2, oneDimvar;
          int error = source_fitting(rank, &par, &mes_e1, &mes_e2, &var_e1, &var_e2, &oneDimvar, &maxL);
 
          cout << "rank:" << rank << " n. " << g << " flux = " << gflux[g] << " scalelength = " << gscale[g] << " position [arcsec] (" << l0/(ARCS2RAD) << "," << m0/(ARCS2RAD) << "): measured e = " << mes_e1 << "," << mes_e2 <<  ",  original e = " << ge1[g] << "," << ge2[g] << endl;
          
          if (error)
          {
              mes_e1 = 0.; mes_e2 = 0.;
              bad_list[bad] = g;  // store index bad sources to try to fit again at the end
              bad++;
          }
          else fprintf(pFile, "%f | %f | %f | %f | %f | %f | %f | %f | %f | %f | %f | %f \n",gflux[g],gscale[g],ge1[g],mes_e1,sqrt(var_e1), ge2[g],mes_e2,sqrt(var_e2),oneDimvar,SNR_vis[g],l0/(ARCS2RAD),m0/(ARCS2RAD));
              
          
#ifdef _OPENMP
#pragma omp parallel for
#endif
          for (unsigned int ch = 0; ch < num_channels; ch++)
          {
              unsigned long int ch_vis = ch*num_coords;
              if (error)
              {
                  // Try fitting at the end (when almost all galaxies are fitted and removed from the data):
                  // add current source model back to the sky model
                  data_galaxy_visibilities(spec[ch], wavenumbers[ch], par.band_factor, time_acc, 0., 0., R_mu,
                                           gflux[g], l0, m0, num_coords, uu_metres, vv_metres, &(visGal[ch_vis]));
                  
                  for (unsigned int i = ch_vis; i<ch_vis+num_coords; i++)
                  {
                      visSkyMod[i].real += visGal[i].real;
                      visSkyMod[i].imag += visGal[i].imag;
                  }
            
              }
              else
              {
                  // get current source model fit
                  data_galaxy_visibilities(spec[ch], wavenumbers[ch], par.band_factor, time_acc, mes_e1, mes_e2, R_mu,
                                           gflux[g], l0, m0, num_coords, uu_metres, vv_metres, &(visGal[ch_vis]));
                  
                  for (unsigned int i = ch_vis; i<ch_vis+num_coords; i++)
                  {
                      // remove it from original data
                      visData[i].real -= visGal[i].real;
                      visData[i].imag -= visGal[i].imag;
                  }
              }
          }
      }
    }
    
    // Re-fitting bad sources again -------------------------------------------------------------------------------------------------------------------------------------
    cout << "bad sources: " << bad << endl << endl;
    cout << "Re-fitting bad sources" << endl;
    
    int nsources = bad;
    bad = 0;
    for (unsigned long int b = 0; b < nsources; b++)
    {
        unsigned long int gal = bad_list[b];
        double flux = gflux[gal];
        
        // set log(prior) for scalelength
        double mu = scale_mean(flux);
        for (int nRo=1; nRo<numR; nRo++)
            rprior[nRo] = rfunc(mu,scale_std,Ro[nRo]);
        double R_mu = exp(mu);
        
        l0 = l[gal];  m0 = m[gal];
#ifdef FACET
        ind = 0;
        while (flux < threshold_flux[ind]) ind++;
        facet = facet_size[ind];
        par.ncoords = evaluate_uv_grid(len, num_coords, uu_metres, vv_metres, facet, &facet_u, &facet_v, count);
        
        source_extraction(l0, m0, flux, exp(mu), 0., 0., &par, visSkyMod, visData, visGal, num_coords, uu_metres, vv_metres, facet, len);
#else
        source_extraction(l0, m0, flux, exp(mu), 0., 0., &par, visSkyMod, visData, visGal, num_coords, uu_metres, vv_metres);
#endif
        
        double mes_e1, mes_e2, maxL;
        double var_e1, var_e2, oneDimvar;
        int error = source_fitting(rank, &par, &mes_e1, &mes_e2, &var_e1, &var_e2, &oneDimvar, &maxL);
        
        cout << "rank:" << rank << " n. " << gal << " flux = " << flux << " scalelength = " << gscale[gal] << " position [arcsec] (" << l0/(ARCS2RAD) << "," << m0/(ARCS2RAD) << "): measured e = " << mes_e1 << "," << mes_e2 <<  ",  original e = " << ge1[gal] << "," << ge2[gal] << endl;
        
        fprintf(pFile, "%f | %f | %f | %f | %f | %f | %f | %f | %f | %f | %f | %f \n",flux,gscale[gal],ge1[gal],mes_e1,sqrt(var_e1), ge2[gal],mes_e2,sqrt(var_e2),oneDimvar,SNR_vis[gal],l0/(ARCS2RAD),m0/(ARCS2RAD));
        
        if (error)
        {
            mes_e1 = 0.; mes_e2 = 0.;
            bad++;
        }
        
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (unsigned int ch = 0; ch < num_channels; ch++)
        {
            unsigned long int ch_vis = ch*num_coords;
            // get current source model fit
            data_galaxy_visibilities(spec[ch], wavenumbers[ch], par.band_factor, time_acc, mes_e1, mes_e2, R_mu,
                                     flux, l0, m0, num_coords, uu_metres, vv_metres, &(visGal[ch_vis]));
 
            for (unsigned int i = ch_vis; i<ch_vis+num_coords; i++)
            {
                // remove it from original data
                visData[i].real -= visGal[i].real;
                visData[i].imag -= visGal[i].imag;
                
            }
        }
    }
    
    cout << "rank: " << rank << " removed " << bad << " bad data galaxies" << endl << endl;
    if (pFile != 0) fclose(pFile);
    
#ifdef USE_MPI
    end_fitting = MPI_Wtime();
    fitting_time = end_fitting - start_fitting;
    double end_tot = MPI_Wtime();
    double total_time = end_tot - start_tot;
#else
    end_fitting = current_timestamp();
    fitting_time = (double)(end_fitting - start_fitting)/1000.;
    long long end_tot = current_timestamp();
    double total_time = (double)(end_tot - start_tot)/1000.;
#endif

    cout << "rank: " << rank << " set up time (sec): " << total_time - data_time - fitting_time << endl;
    cout << "rank: " << rank << " data generation time (sec): " << data_time << endl;
    cout << "rank: " << rank << " data fitting computation time (sec): " << fitting_time << endl;
    cout << "rank: " << rank << " Total time (sec): " << total_time << endl;
    
    // free memory ----------------------------------------------------------------------------------------------------------------
    delete[] visMod;
    delete[] visGal;
    delete[] visSkyMod;
    delete[] visData;
    delete[] Ro;
    delete[] rprior;
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
