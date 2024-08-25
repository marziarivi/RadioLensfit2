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
#include "evaluate_uv_grid.h"
#include "data_processing.h"

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
    const double freq_start_hz = ms_freq_start_hz(ms);            // Start Frequency, in Hz
    const double channel_bandwidth_hz = ms_freq_inc_hz(ms);       // Frequency channel bandwidth, in Hz
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
    float *sigma2_vis;
    try
    {
        sigma2_vis = new float[num_vis];
        sizeGbytes = num_vis*sizeof(float)/((double)(1024*1024*1024));
        cout << "allocated sigma2 visibilities: " << num_vis << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << " bad_alloc caught: " << ba.what() << '\n';
    }
    float sigma2 = (SEFD*SEFD)/(2.*time_acc*channel_bandwidth_hz*efficiency*efficiency);
    ms_read_sigma(ms, 0, num_coords, sigma2_vis, &status);
    if (status)
    {   
        cout << "ERROR reading MS - sigma: " << status << endl;
        for (unsigned long int i = 0; i<num_vis; i++)
             sigma2_vis[i] = sigma2; // visibility noise variance
        cout << "Use theoretical rms: " << sqrt(sigma2) << " uJy" << endl;
    }
    else
    {
      for (unsigned long int i = 0; i<num_vis; i++)
         if (!sigma2_vis[i])  // if sigma is 0 the use theoretical value
              sigma2_vis[i] = sigma2; // visibility noise variance
    }    

    cout << "rank " << rank << ": MS data total GBytes: " << totGbytes << endl;
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
    unsigned int facet_ncoords = evaluate_uv_grid_size(0,1,len,wavenumbers,num_channels,num_coords, uu_metres, vv_metres, facet, flag);

    // allocate partial weights sum (per cell for weighted average)
    double *sum_w;
    try
    {
        sum_w = new double[ncells];
        sizeGbytes = ncells*sizeof(double)/((double)(1024*1024*1024));
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "rank " << rank << ": bad_alloc caught: " << ba.what() << '\n';
    }
    cout << "rank " << rank << ": allocated  array weights sum: " << ncells << ", size = " << sizeGbytes  << " GB" << endl;
    
    // allocate facet arrays
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
        cout << "allocated models: num_models = " << num_models << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "bad_alloc caught: " << ba.what() << '\n';
    }

    par.uu = facet_u;
    par.vv = facet_v;
    par.data = facet_visData;
    par.mod = visMod;
    par.sigma2 = facet_sigma2;   
#else
    complexd* visMod;
    try
    {
        unsigned int model_ncoords = num_coords;
        visMod = new complexd[num_models*model_ncoords*num_channels];
        sizeGbytes = num_models*model_ncoords*num_channels*sizeof(complexd)/((double)(1024*1024*1024));
        cout << "allocated models: num_models = " << num_models << ", size = " << sizeGbytes  << " GB" << endl;
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
    
    cout << "Total Visibilities GBytes: " << totGbytes << endl;
    cout << num_models << " samples in galaxy scale-length, " << Rmin << " < r0 < " << Rmax << " arcsec" << endl;
 
    // Data processing -----------------------------------------------------------------------------------------------------------------------------------
    FILE *pFile;
    char output[100];
    sprintf(output,"ellipticities-%dch.txt",num_channels);
    pFile = fopen(output,"w");
    fprintf(pFile, "flux | e1 | m_e1 | err1 | e2 | m_e2 | err2 | 1D var | SNR |   l  |  m  | \n");
    start_extraction = current_timestamp();

    unsigned int bad_list[mygalaxies];
    int bad = 0;
    data_processing(false, bad_list, 1, 0, mygalaxies, len, num_coords, pFile, &par, l, m, gflux, gscale, ge1, ge2, SNR_vis,
                    sum_w, visGal, visSkyMod, visData, sigma2_vis, flag, uu_metres, vv_metres, ww_metres, &fitting_time, &bad);
    
    // Re-fitting bad sources again -------------------------------------------------------------------------------------------------------------------------------------
    if (bad > 0)
    {    
        cout << "Re-fitting " << bad << " bad sources" << endl;
        data_processing(true, bad_list, 1, 0, bad, len, num_coords, pFile, &par, l, m, gflux, gscale, ge1, ge2, SNR_vis,
                      sum_w, visGal, visSkyMod, visData, sigma2_vis, flag, uu_metres, vv_metres, ww_metres, &fitting_time, &bad);

    }
    end_extraction = current_timestamp();
    extraction_time = (double)(end_extraction - start_extraction)/1000. - fitting_time;
    long long end_tot = current_timestamp();
    double total_time = (double)(end_tot - start_tot)/1000.;

    if (pFile != 0) fclose(pFile);
    cout << "Removed " << bad << " bad data galaxies" << endl << endl;
    cout << "Total time (sec): " << total_time << endl;
    cout << "Data reading time (sec): " << data_time << endl;
    cout << "Set up time (sec): " << total_time - data_time - model_time - extraction_time - fitting_time << endl;
    cout << "Sky model visibilities time (sec): " << model_time << endl;
    cout << "Source extraction computation time (sec): " << extraction_time << endl;
    cout << "Data fitting computation time (sec): " << fitting_time << endl;
    
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
    delete[] sum_w;
    delete[] facet_u;
    delete[] facet_v;
    delete[] facet_visData;
    delete[] facet_sigma2;
#endif
}
