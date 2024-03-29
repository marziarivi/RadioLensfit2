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


//  RadioLensfit2.cpp (serial version)
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
    argv[3]  MS filename
*/

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
    int num_threads, rank = 0;
    long long start_tot;
    start_tot = current_timestamp();
#ifdef _OPENMP
#pragma omp parallel
    num_threads = omp_get_num_threads();
    cout << "Number of OpenMP threads = " << num_threads << endl;
#endif
    
    if (argc != 4)
    {
      cout << "ERROR: bad number of parameters!" << endl;
      cout << "usage: RadioLensfit2 <source catalog filename> <num_sources> <filename MS> " << endl;
      exit(EXIT_FAILURE);
    }

    double data_time = 0.;
    double extraction_time = 0.;
    double fitting_time = 0.;
    long long start_data,end_data,start_fitting,end_fitting,start_extraction,end_extraction;
    start_data = current_timestamp();

    // Read Measurement Set --------------------------------------------------------------------------------------------------------------------------------------------------------------------
    RL_MeasurementSet* ms = ms_open(argv[3]);

    //double RA = ms_phase_centre_ra_rad(ms);                 // Phase Centre coordinates
    //double Dec = ms_phase_centre_dec_rad(ms);   
    const unsigned int num_stations = ms_num_stations(ms);        // Number of stations
    const unsigned int num_channels = ms_num_channels(ms);        // Number of frequency channels
    const unsigned int num_rows = ms_num_rows(ms);                // Number of rows 
    const double freq_start_hz = ms_freq_start_hz(ms); //1280e+6;          // Start Frequency, in Hz
    const double channel_bandwidth_hz = ms_freq_inc_hz(ms); //240e+6;      // Frequency channel bandwidth, in Hz
    const double full_bandwidth_hz = channel_bandwidth_hz * num_channels;  // Frequency total bandwidth, in Hz
    const int time_acc = ms_time_inc_sec(ms);                     // accumulation time (sec)

    const double efficiency = EFFICIENCY;     // system efficiency
    const double SEFD = SEFD_SKA;    // System Equivalent Flux Density (in micro-Jy) of each SKA1 antenna

    const double ref_frequency_hz = REF_FREQ;  //Reference frequency in Hz at which fluxes are measured    
    unsigned int num_baselines = num_stations * (num_stations - 1) / 2;

    cout << "Reference frequency (Hz): " << ref_frequency_hz << endl;
    cout << "Number baselines: " << num_baselines << endl;
    cout << "Number of channels: " << num_channels << endl;
    cout << "Channels bandwidth (Hz): " << channel_bandwidth_hz << endl;
    cout << "Starting frequency (Hz): " << freq_start_hz << endl;
    cout << "Accumulation time (sec): " << time_acc << endl;
    cout << "Number of rows: " << num_rows << endl;
    cout << "Number of polarizations: " << ms_num_pols(ms) << endl;   
 
    double sizeGbytes, totGbytes = 0.;
    
    // Allocate and read uv coordinates
    unsigned long int num_coords = num_rows; 
    double* uu_metres = new double[num_coords];
    double* vv_metres = new double[num_coords];
    double* ww_metres = new double[num_coords];
    sizeGbytes = 3*num_coords*sizeof(double)/((double)(1024*1024*1024));
    cout <<  "allocated original coordinates: " << sizeGbytes  << " GB" << endl;
    totGbytes += sizeGbytes;
    
    int status = 0;
    double len = ms_read_coords(ms,0,num_coords,uu_metres,vv_metres,ww_metres,&status);
    if (status) 
    {
        cout <<"ERROR reading MS - uvw points: " << status << endl;
        exit(EXIT_FAILURE);
    }
 
    // Allocate and read Data visibilities
    unsigned long int num_vis  = (unsigned long int) num_channels * num_coords;
    complexd *visData;
    try
    {
        visData = new complexd[num_vis];
        sizeGbytes = num_vis*sizeof(complexd)/((double)(1024*1024*1024));
        cout << "allocated original data visibilities: " << num_vis << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "bad_alloc caught: " << ba.what() << '\n';
    }

    ms_read_vis(ms, 0, 0, num_channels, num_rows, "DATA", visData, &status);
    if (status) 
    {
        cout << "ERROR reading MS - DATA column: " << status << endl;
        exit(EXIT_FAILURE);
    }

    // Allocate and read FLAG column
    bool *flag;
    try
    {
        flag = new bool[num_vis];
        sizeGbytes = num_vis*sizeof(bool)/((double)(1024*1024*1024));
        cout << "allocated flag column: " << num_vis << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "bad_alloc caught: " << ba.what() << '\n';
    }

    unsigned long int nF = ms_read_Flag(ms, 0, 0, num_channels, num_rows, "FLAG",flag, &status);
    if (status)
    {
      cout << "ERROR reading MS - flag: " << status << endl;
      exit(EXIT_FAILURE);
    }
    else
      cout << "Number of flagged visibilities: " << nF << endl;

    // Allocate and read SIGMA column
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

    //ms_read_sigma(ms, 0, num_coords, sigma2_vis, &status);
    if (status)
    {
        cout << "ERROR reading MS - sigma: " << status << endl;
        exit(EXIT_FAILURE);
    }
    double sigma2 = (SEFD*SEFD)/(2.*time_acc*channel_bandwidth_hz*efficiency*efficiency);
    for (unsigned long int i = 0; i<num_vis; i++)
        sigma2_vis[i] = sigma2; // visibility noise variance
 
    ms_close(ms); 

    // Read galaxy catalogue --------------------------------------------------------------------------------------------------------------------------
    unsigned int nge = atof(argv[2]);
    
    double *gflux = new double[nge];
    double *l = new double[nge];
    double *m = new double[nge];
    double *ge1 = new double[nge];
    double *ge2 = new double[nge];
    double *gscale = new double[nge];
    double *SNR_vis = new double[nge];
    bool readSNR = true;
 
    unsigned int mygalaxies = read_catalog(nge, argv[1],gflux,gscale,ge1,ge2,l,m,SNR_vis,readSNR);
    cout << "Read catalog. Number of sources: " << mygalaxies << endl;
   
    end_data = current_timestamp();
    data_time = (double)(end_data - start_data)/1000.;
    
    // Sky model visibilities computation --------------------------------------------------------------------------------------------------------------------------
    // Allocate Galaxy and Sky Model Visibilities
    complexd *visGal, *visSkyMod;
    try
    {
        visGal = new complexd[num_vis];
        cout << "allocated galaxy visibilities: " << num_vis << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "bad_alloc caught: " << ba.what() << '\n';
    }

    try
    {
        visSkyMod = new complexd[num_vis];
        cout << "allocated sky model visibilities: " << num_vis << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "bad_alloc caught: " << ba.what() << '\n';
    }

    memset(visSkyMod, 0, num_vis*sizeof(complexd));    

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

    long long start_model = current_timestamp();
    sky_model(wavenumbers, spec, channel_bandwidth_hz, time_acc, num_channels,
               mygalaxies, gflux, gscale, l, m, num_coords, uu_metres, vv_metres, ww_metres, visGal, visSkyMod);
    long long end_model = current_timestamp();
    double model_time = (double)(end_model - start_model)/1000.;

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
    par.nchannels = num_channels;
    par.band_factor = channel_bandwidth_hz*PI/C0;
    par.acc_time = time_acc;
    par.spec = spec;
    par.wavenumbers = wavenumbers; // wavenumbers for the model
       
#ifdef FACET
    // Faceting uv coordinates ----------------------------------------------------------------------------------------
    len = len*wavenumbers[num_channels-1]/(2*PI);    // max len in wavelength units
    int facet = facet_size(RMAX,len);
    unsigned long int ncells = facet*facet;
    unsigned long int* count = new unsigned long int[ncells];
    unsigned int facet_ncoords = evaluate_uv_grid_size(0,1,len,wavenumbers,num_channels,num_coords, uu_metres, vv_metres, facet, flag, count);

    double *facet_u, *facet_v;
    try
    {
      facet_u = new double[facet_ncoords];
      facet_v = new double[facet_ncoords];
      sizeGbytes = (2*facet_ncoords*sizeof(double)+ncells*sizeof(unsigned long int))/((double)(1024*1024*1024));
      cout << "allocated grid coordinates and array counter: " << sizeGbytes  << " GB" << endl;
      totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "bad_alloc caught: " << ba.what() << '\n';
    }    

    unsigned int facet_nvis = facet_ncoords;
    complexd* facet_visData;
    double* facet_sigma2;
    try
    {
        facet_visData = new complexd[facet_nvis];
        facet_sigma2 = new double[facet_nvis];
        sizeGbytes = (facet_nvis*sizeof(complexd)+facet_ncoords*sizeof(double))/((double)(1024*1024*1024));
        cout << "allocated gridded visibilities and variances: " << facet_nvis << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "bad_alloc caught: " << ba.what() << '\n';
    }

    // Allocate Facet Model Visibilities ------------------------------------------------------------------------------------------
    double* visMod;
    try
    {
        unsigned int model_ncoords = facet_ncoords;
        visMod = new double[num_models*model_ncoords];
        sizeGbytes = num_models*model_ncoords*sizeof(double)/((double)(1024*1024*1024));
        cout << "allocated models: num_models= " << num_models << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "bad_alloc caught: " << ba.what() << '\n';
    }

    par.uu = facet_u;
    par.vv = facet_v;
    //par.weights = weights;
    par.data = facet_visData;
    par.count = count;
    par.mod = visMod;
    par.sigma2 = facet_sigma2;   
 
#else
    complexd* visMod;
    try
    {
        unsigned int model_ncoords = num_coords;
        visMod = new complexd[num_models*model_ncoords*num_channels];
        sizeGbytes = num_models*model_ncoords*num_channels*sizeof(complexd)/((double)(1024*1024*1024));
        cout << "allocated models: num_models= " << num_models << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "bad_alloc caught: " << ba.what() << '\n';
    }

    par.ncoords = num_coords;
    par.uu = uu_metres;
    par.vv = vv_metres;
    par.ww = ww_metres;
    par.data = visGal;
    par.count = 0;
    par.mod = visMod;
    par.sigma2 = sigma2_vis;
#endif
    
    cout << "Total Visibilities GBytes per rank: " << totGbytes << endl;
    cout << num_models << " samples in galaxy scale-length, " << Rmin << " < r0 < " << Rmax << " arcsec" << endl;

    start_extraction = current_timestamp();
 
    // Data Fitting -----------------------------------------------------------------------------------------------------------------------------------
    FILE *pFile;
    char output[100];
    sprintf(output,"ellipticities-%dch.txt",num_channels);
    pFile = fopen(output,"w");
    fprintf(pFile, "flux | e1 | m_e1 | err1 | e2 | m_e2 | err2 | 1D var | SNR |   l  |  m  | \n");
    
    double l0,m0;
    double mes_e1, mes_e2, maxL;
    double var_e1, var_e2, oneDimvar;
    unsigned int bad_list[mygalaxies];
    int bad = 0;
    
    for (unsigned int g=0; g<mygalaxies; g++)
    {
        double mu = scale_mean(gflux[g]);
        double R_mu = exp(mu);
     //   double R_mu = gscale[g];

        // set log(prior) for scalelength
        for (int nRo=1; nRo<numR; nRo++)
           rprior[nRo] = rfunc(mu,R_STD,Ro[nRo]);
        
        l0 = l[g];  m0 = m[g];
#ifdef FACET
        unsigned int facet = facet_size(R_mu,len);
        source_extraction(0,facet,&par,par.data,par.sigma2,count,l0, m0, gflux[g], R_mu, 0.,0., visSkyMod, visData, visGal, sigma2_vis, flag, num_channels, num_coords, uu_metres, vv_metres, ww_metres, len);
        par.ncoords = evaluate_facet_coords(par.uu, par.vv, len, facet, count);
#else
        source_extraction(l0, m0, gflux[g], R_mu, 0., 0., &par, visSkyMod, visData, visGal, sigma2_vis, num_channels, num_coords, uu_metres, vv_metres, ww_metres);
#endif

        start_fitting = current_timestamp();
        int error = source_fitting(rank, &par, &mes_e1, &mes_e2, &var_e1, &var_e2, &oneDimvar, &maxL);
        cout << "n. " << g << " flux = " << gflux[g] << ": measured e = " << mes_e1 << "," << mes_e2 << endl;
        fitting_time += (double) (current_timestamp() -start_fitting)/1000.;
 
       if (error)
        {
              bad_list[bad] = g;  // store index bad sources to try to fit again at the end
              bad++;

              // Try fitting at the end (when almost all galaxies are fitted and removed from the data):
              // add current source round model back to the sky model
#ifdef _OPENMP
#pragma omp parallel for
#endif
              for (unsigned int ch = 0; ch < num_channels; ch++)
              {             
                 unsigned long int ch_vis = (unsigned long int) ch*num_coords;
                 data_galaxy_visibilities(spec[ch], wavenumbers[ch], par.band_factor, time_acc, 0., 0., R_mu,
                                             gflux[g], l0, m0, num_coords, uu_metres, vv_metres, ww_metres, &(visGal[ch_vis]));
               
                 for (unsigned int i = ch_vis; i<ch_vis+num_coords; i++)
                 {
                      visSkyMod[i].real += visGal[i].real;
                      visSkyMod[i].imag += visGal[i].imag;
                 }

              }
        }
        else
        {
             fprintf(pFile, "%f | %f | %f | %f | %f | %f | %f | %f | %f | %f | %f \n",gflux[g],ge1[g],mes_e1,sqrt(var_e1),ge2[g],mes_e2,sqrt(var_e2),oneDimvar,SNR_vis[g],l0/(ARCS2RAD),m0/(ARCS2RAD));
             
#ifdef _OPENMP
#pragma omp parallel for
#endif
             // get current source model fit and remove it from the original data
             for (unsigned int ch = 0; ch < num_channels; ch++)
             {
                unsigned long int ch_vis = (unsigned long int) ch*num_coords;
                data_galaxy_visibilities(spec[ch], wavenumbers[ch], par.band_factor, time_acc, mes_e1, mes_e2, R_mu,
                                           gflux[g], l0, m0, num_coords, uu_metres, vv_metres, ww_metres, &(visGal[ch_vis]));
                  
                for (unsigned int i = ch_vis; i<ch_vis+num_coords; i++)
                {
                      visData[i].real -= visGal[i].real;
                      visData[i].imag -= visGal[i].imag;
                }
              }
        }
    }
    
    // Re-fitting bad sources again -------------------------------------------------------------------------------------------------------------------------------------
    cout << "Re-fitting " << bad << " bad sources" << endl;
    
    int nsources = bad;
    bad = 0;
    for (unsigned int b = 0; b < nsources; b++)
    {
        unsigned int gal = bad_list[b];
        double flux = gflux[gal];
        double mu = scale_mean(flux);
        double R_mu = exp(mu);
        //double R_mu = gscale[gal];

        // set log(prior) for scalelength
        for (int nRo=1; nRo<numR; nRo++)
            rprior[nRo] = rfunc(mu,R_STD,Ro[nRo]);
        
        l0 = l[gal];  m0 = m[gal];
#ifdef FACET
        unsigned int facet = facet_size(R_mu,len);
        //par.ncoords = evaluate_uv_grid_size(len,wavenumbers, num_channels,num_coords, uu_metres, vv_metres, facet, count);
        source_extraction(0,facet,&par,par.data,par.sigma2,count,l0, m0, flux, R_mu, 0., 0., visSkyMod, visData, visGal, sigma2_vis, flag, num_channels, num_coords, uu_metres, vv_metres, ww_metres, len);
        par.ncoords = evaluate_facet_coords(par.uu, par.vv, len, facet, count);
#else
        source_extraction(l0, m0, flux, R_mu, 0., 0., &par, visSkyMod, visData, visGal, sigma2_vis, num_channels, num_coords, uu_metres, vv_metres, ww_metres);
#endif
        start_fitting = current_timestamp();
        double mes_e1, mes_e2, maxL;
        double var_e1, var_e2, oneDimvar;
        int error = source_fitting(rank, &par, &mes_e1, &mes_e2, &var_e1, &var_e2, &oneDimvar, &maxL);
        
        cout << "n. " << gal << " flux = " << flux << "): measured e = " << mes_e1 << "," << mes_e2 << endl;
        fprintf(pFile, "%f | %f | %f | %f | %f | %f | %f | %f | %f | %f | %f  \n",flux,ge1[gal],mes_e1,sqrt(var_e1),ge2[gal],mes_e2,sqrt(var_e2),oneDimvar,SNR_vis[gal],l0/(ARCS2RAD),m0/(ARCS2RAD));
      
        fitting_time += (double) (current_timestamp() - start_fitting)/1000.;
  
        if (error)
        {
            mes_e1 = 0.; mes_e2 = 0.;  // remove round source model from original data
            bad++;
        }

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (unsigned int ch = 0; ch < num_channels; ch++)
        {
            unsigned long int ch_vis = (unsigned int ) ch*num_coords;
            // get current source model fit
            data_galaxy_visibilities(spec[ch], wavenumbers[ch], par.band_factor, time_acc, mes_e1, mes_e2, R_mu,
                                     flux, l0, m0, num_coords, uu_metres, vv_metres, ww_metres, &(visGal[ch_vis]));
 
            for (unsigned int i = ch_vis; i<ch_vis+num_coords; i++)
            {
                // remove it from original data
                visData[i].real -= visGal[i].real;
                visData[i].imag -= visGal[i].imag;
            }
        }
    }
    
    cout << "removed " << bad << " bad data galaxies" << endl << endl;
    if (pFile != 0) fclose(pFile);
    
    end_extraction = current_timestamp();
    extraction_time = (double)(end_extraction - start_extraction)/1000. - fitting_time;
    long long end_tot = current_timestamp();
    double total_time = (double)(end_tot - start_tot)/1000.;

    cout << "set up time (sec): " << total_time - data_time - model_time - extraction_time - fitting_time << endl;
    cout << "data reading time (sec): " << data_time << endl;
    cout << "sky model visibilities time (sec): " << model_time << endl;
    cout << "source extraction computation time (sec): " << extraction_time << endl;
    cout << "data fitting computation time (sec): " << fitting_time << endl;
    cout << " Total time (sec): " << total_time << endl;
    
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
#endif
}
