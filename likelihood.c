/*
 * Copyright (c) 2018 Marzia Rivi
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

//
//  Computation of the likelihood function through galaxy model-fitting
//  It depends on 5 parameters: ellipcticy (e1,e2), scalelength (r), position (x), flux (S).
//  Marginalisation over S, x is adopted to reduce it to a function of the ellipticy and scalelength only.

#ifdef _OPENMP
#include <omp.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_vector.h>

#include "datatype.h"
#include "galaxy_visibilities.h"
#include "likelihood.h"
#include "marginalise_r.h"
#include "distributions.h"

#ifdef __cplusplus
extern "C" {
#endif
    
    
/*
double f_posterior (const gsl_vector *v, void *params)
{
        double log_post,ee1, ee2, scale;
        int error = 0;
        likelihood_params *par = (likelihood_params *)params;
        ee1 = gsl_vector_get(v, 0);
        ee2 = gsl_vector_get(v, 1);
        scale = gsl_vector_get(v, 2);
    
        double e_mod = sqrt(ee1*ee1+ee2*ee2);
    
        if(e_mod <= 0.8 && e_mod > 0.)
        {
            double L_r = loglikelihood_r(par->nchannels, par->band_factor, par->acc_time, par->spec, par->wavenumbers, ee1, ee2, par->l0, par->m0, par->radius, scale, par->ncoords, par->count,par->sigma,par->uu,par->vv,par->data,par->mod);
            log_post = L_r+log(e_pdf(e_mod));
        }
        else log_post = -1.e10;
    
        return -log_post;
}
*/
    
/*
 *  Computation of the likelihood for a fixed value of the ellipticity (ee1,ee2)
 *
 *  First computate the likelihood for a fixed value of the ellipticity (e1,e2) and scalelength r:
 *  
 *  chi_square and likelihood L defined as in Paper I (Miller et al., 2007 p. 319)
 *  after integrating analitically over the flux using a uniform prior as in Paper I, we get
 *   
 *     logL(e1, e2, r, X) = const + h^2(X)/(2*sigma_pix^2)
 *  
 *   where h(X) is the cross-correlation and sigma is the uncertainty of the data value.
 *   L must be then numerically marginalised over the shift X of the galaxy position.
 *
 *   Then marginalise over the scalelenght.
 */
   
    
double f_likelihood (const gsl_vector *v, void *params)
{
    double L_e,ee1, ee2;
    int error = 0;
    ee1 = gsl_vector_get(v, 0);
    ee2 = gsl_vector_get(v, 1);
    
    if(ee1*ee1+ee2*ee2 <= 0.64)
        L_e = loglikelihood(params, ee1, ee2, &error);
    
    else L_e = -1.e10;
    
    return -L_e;
}

    
/*
 *  Likelihood computation as function of ellipticity (marginalise over position and scalelength)
 */

double loglikelihood(void *params, double ee1, double ee2, int *error)
{
    likelihood_params *par = (likelihood_params *)params;
    int numR = par->numr;
    double L_e;
    int numvals;
    double* xmarvals = (double *) malloc(sizeof(double)*(numR+10));
    double* ymarvals = (double *) malloc(sizeof(double)*(numR+10));
    double* L_r = (double *) malloc(sizeof(double)*numR);
    *error = 0;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int nRo = 1; nRo < numR; nRo++)
    {
       unsigned long int index = (nRo-1)*(par->ncoords)*(par->nchannels);

#ifdef FACET
       L_r[nRo] = loglikelihood_r(par->nchannels, par->band_factor, par->acc_time, par->spec, par->wavenumbers, ee1, ee2, par->l0, par->m0, par->radius,(par->ro)[nRo], par->ncoords, par->count, par->sigma, par->uu, par->vv, par->weights, par->data, &((par->mod)[index]));
#else
       L_r[nRo] = loglikelihood_r(par->nchannels, par->band_factor, par->acc_time, par->spec, par->wavenumbers, ee1, ee2, par->l0, par->m0, par->radius,(par->ro)[nRo], par->ncoords, par->count, par->sigma, par->uu, par->vv, par->ww, par->data, &((par->mod)[index]));
#endif
    }
    
    // marginalisation over scalelength
    set_posterior_values(numR, L_r, par->rprior, par->ro, xmarvals, ymarvals, &numvals);
    if (numvals > 4)  L_e = marginalise_posterior_r(numvals, xmarvals, ymarvals);
    else
    {
        printf("likelihood error for e1 = %f, e2 = %f: too few points (n=%d) for marginalisation over the scalelength\n",ee1,ee2,numvals);
        *error = 1; L_e = -1.e10;
    }
    if (L_e == -1.e10) *error = 1;
    
    free(L_r);
    free(xmarvals);
    free(ymarvals);
  
    return L_e;
}
   
    
/*
 *  Likelihood computation as function of ellipticity and scalelength (marginalise over position)
 */
    
#ifdef FACET
double loglikelihood_r(unsigned int nchannels, double band_factor, double acc_time, double* spec,
                       double* wavenumbers, double ee1, double ee2, double l, double m, double radius,
                       double scale,
                       unsigned long int n_uv_coords, unsigned long int* count, const double variance,
                       double* uu_metres, double* vv_metres, double* weights, complexd* visData, double* visM)
#else
double loglikelihood_r(unsigned int nchannels, double band_factor, double acc_time, double* spec,
                double* wavenumbers, double ee1, double ee2, double l, double m, double radius,
                double scale, unsigned long int n_uv_coords, unsigned long int* count,
                const double variance, double* uu_metres, double* vv_metres, double* ww_metres, complexd* visData, complexd* visM)
#endif
{
    // generate model
#ifdef FACET
    model_galaxy_visibilities_at_zero(nchannels, spec, wavenumbers, ee1, ee2, scale, radius, n_uv_coords, uu_metres, vv_metres, count, visM);
#else
    model_galaxy_visibilities(nchannels, spec, wavenumbers, band_factor, acc_time, ee1, ee2, scale, l,m, radius, n_uv_coords, uu_metres, vv_metres, ww_metres, count, visM);
#endif
        
    // Compute log(likelihood) dependend only on ellipticity and scale-length
    double L_er,ho, det_sigma;
    int error = cross_correlation(nchannels, wavenumbers, n_uv_coords, count, uu_metres, vv_metres, weights, visData, visM, &ho, &det_sigma);
 //   if (det_sigma <= 0) L_er = -1.e10; //printf("det_sigma = %e\n",det_sigma);fflush(stdout);}
    if (error) L_er = -1.e10;
    else
    {
        double B = (ho*ho)/(2.*variance);
        L_er = marginalise_over_position_shift(B);
        if (isnan(L_er)) L_er = -1.e10;
        else L_er -= 0.5*log(det_sigma); // = log(sqrt(1/det_sigma))
    }
        
    return L_er;
}

 
/*
 *  the cross-correlation (take only the real part) is approximated as a 2D gaussian
 *  centered in its maximum ho = h(xo,yo) and det(covariance matrix) = 1/det(Hessian(xo,yo)).
 */
#ifdef FACET
int cross_correlation(unsigned int nchannels, double* wavenumbers, unsigned long int n_uv_coords,
                           unsigned long int* count, double* uu_metres, double* vv_metres, double* weights, complexd* visData,
                           double* visMod, double* ho, double* det_sigma)
#else
int cross_correlation(unsigned int nchannels, double* wavenumbers, unsigned long int n_uv_coords,
                       unsigned long int* count, double* uu_metres, double* vv_metres, complexd* visData,
                       complexd* visMod, double* ho, double* det_sigma)
#endif
{
    double a, real_part, imag_part, res_arc, wavenumber, wavenumber2, u,v;
    double det, d2h_dx2_xo, d2h_dy2_yo, dh_dx_xo, dh_dy_yo, d2h_dxdy_o, incx, incy;
    unsigned long int i,k;
    int error = 0;
    /* cross-correlation h in the Fourier domain: conj(visData)*visMod */
    complexd* h_F = (complexd*) malloc(nchannels*n_uv_coords*sizeof(complexd));
    
    d2h_dx2_xo = 0.;
    d2h_dy2_yo = 0.;
    dh_dx_xo = 0.;
    dh_dy_yo = 0.;
    d2h_dxdy_o = 0;
 
    k=0;
    double value = 0;
    for (unsigned int ch=0; ch<nchannels; ch++)
    {
        wavenumber = wavenumbers[ch];
        wavenumber2 = wavenumber*wavenumber;
        
        for (i=0; i < n_uv_coords; i++)
        {
#ifdef FACET
            h_F[k].real = visData[k].real*visMod[k];
            h_F[k].imag = - visData[k].imag*visMod[k];
            
            h_F[k].real *= count[i];
            h_F[k].imag *= count[i];
#else
            h_F[k].real = (visData[k].real*visMod[k].real + visData[k].imag*visMod[k].imag);
            h_F[k].imag = (visData[k].real*visMod[k].imag - visData[k].imag*visMod[k].real);
#endif
            value += h_F[k].real;
            
            u = uu_metres[i];
            v = vv_metres[i];
            
            dh_dx_xo -= h_F[k].imag * wavenumber * u;
            dh_dy_yo -= h_F[k].imag * wavenumber * v;
            
            d2h_dx2_xo -= h_F[k].real * wavenumber2 * u*u;
            d2h_dy2_yo -= h_F[k].real * wavenumber2 * v*v;
            d2h_dxdy_o -= h_F[k].real * wavenumber2 * u*v;
            k++;
        }
    }
    
    double h_0 = value;
    det  = d2h_dx2_xo*d2h_dy2_yo - d2h_dxdy_o*d2h_dxdy_o;
 
    /*
     * we compute (xo,yo) as the critical point of the quadratic obtained by the Taylor
     * expansion of h(X)=h_F*e^{ikX} near the origin up to the second order (Newton's method).
     */
    double xo=0., yo=0.;
    double tolerance=1.e-7;
    double norm = 1.;
    int iter = 0;
    
    while (det > 0 && (norm > tolerance) && iter < 50)
    {
        iter++;
        
        // inc = - Hess(h)^{-1}*grad(h)
        incx = (dh_dx_xo*d2h_dy2_yo - dh_dy_yo*d2h_dxdy_o)/det;
        incy = (dh_dy_yo*d2h_dx2_xo - dh_dx_xo*d2h_dxdy_o)/det;
        
        xo -= incx;
        yo -= incy;
        
        d2h_dx2_xo = 0.;
        d2h_dy2_yo = 0.;
        dh_dx_xo = 0.;
        dh_dy_yo = 0.;
        d2h_dxdy_o = 0;
        value = 0;
        
        for (unsigned int ch=0; ch<nchannels; ch++)
        {
          wavenumber = wavenumbers[ch];
          wavenumber2 = wavenumber*wavenumber;
            
          for (i = 0; i < n_uv_coords; ++i)
          {
             u = uu_metres[i];
             v = vv_metres[i];
             a = wavenumber * (xo*u + yo*v);
        
             // complex multiply, take only the real part
             k = ch*n_uv_coords + i;
             real_part = (h_F[k].real * cos(a) - h_F[k].imag * sin(a));
             imag_part = (h_F[k].real * sin(a) + h_F[k].imag * cos(a));
             value += real_part;
        
             dh_dx_xo -= imag_part * wavenumber * u;
             dh_dy_yo -= imag_part * wavenumber * v;
        
             d2h_dx2_xo -= real_part * wavenumber2 * u*u;
             d2h_dy2_yo -= real_part * wavenumber2 * v*v;
             d2h_dxdy_o -= real_part * wavenumber2 * u*v;
          }
        }
        det  = d2h_dx2_xo*d2h_dy2_yo - d2h_dxdy_o*d2h_dxdy_o;
        
        norm = sqrt(incx*incx+incy*incy);
    }
 
    *ho = value;
    *det_sigma = det;
 
    free(h_F);
    if (det <= 0) return 1;
    return 0;
}
    
    
/*
 * Marginalisation of L over the position shift:
 * we assume a uniform prior;
 * after moving to polar coordinates, it corresponds to the
 * exact integration of \int_0^r_max r * exp[B*exp(-r^2)]dr.
 * By the substitution t=-B*exp(-r^2), this integral is equal to Ei(B)-Ei(B*exp(-r_max^2)),
 * where Ei is the exponential integral function.
 * The logL(0) - logL(r_max) = DeltaChiSquared/2 ~ 3 condition assumes we want to integrate out to
 * a radius r_max where the positional confidence interval is 95 percent.
 * Uses gnu scientific library routine for ExpIntegralEi.
 *
 * This function corresponding to likelfunc(double x) in lensfit
 * the parameter x is B = logL(0).
 * It returns log(L(e1, e2, r))
 */
double marginalise_over_position_shift(double x)
{
    double y, Delta, DeltaChiSquared;
    gsl_sf_result yc, ys;
    
    DeltaChiSquared = 5.991;
    
    // Note that Delta = logL(0) - logL(r_max) = B - B*exp(r^2), i.e. the integration interval
    Delta = 0.5*DeltaChiSquared;
    
    /*
     * if chi-squared max below 5.99 the 95 percent confidence region
     * on the galaxy's position is unconstrained so the shape measurement
     * carries no signal
     */
    
    if (x <= Delta)
    {
        y = -1.e10;
        return y;
    }
    
    // deal with very large values using asymptotic approximation
    if (x > 100.)
    {
        y = x - log(Delta) - exp(-Delta);
        return y;
    }
  
    // do not let x be too close to Delta
    if ( (x-Delta)<0.0001 ) x = Delta+0.0001;
    
    // evaluate exponential integral and normalise by prior area (=\pi*r_max^2)
    if ( gsl_sf_expint_Ei_e(x, &yc) == 0 &&
        gsl_sf_expint_Ei_e((x-Delta), &ys) == 0 )
    {
        y = log( (yc.val - ys.val)/log(x/(x-Delta)) );
    }
    else
    {
        // error from gsl_sf: don't halt, just return "bad" value
        y =  -1.e10;
    }
    return y;
    
}
    
    
/*
 * Evaluation of the likelihood in a neighborhood of its maximum
 * Take the mean and variance of the values, for which L[e] >= 0.05 maxL, as measure and uncertainty of the ellipticity
 */
int likelihood_sampling(int rank, double *mes_e1, double *mes_e2, double maxL, void *params, int np_max, double *var_e1, double *var_e2, double *cov_e)
{
    int error = 0;
    int ie1,ie2,np;
    double totL, L_e, average_1,average_2, cov, var1, var2,x_e1,x_e2;
    double xL1, xL2, mod2 = 0;
    likelihood_params *par = (likelihood_params *)params;
    
    double L = 1.;
    double threshold = 0.01; //0.05;
    double max_e1 = *mes_e1;
    double max_e2 = *mes_e2;
    average_1 = max_e1;
    average_2 = max_e2;
    var1 = max_e1*max_e1;
    var2 = max_e2*max_e2;
    cov = max_e1*max_e2;
    totL = 1.; //L=exp(maxL-maxL)
    np = 1;
    
    //search for the range where sampling the likelihood
    float oldsampling = 1.;     //dummy value to make sure loop is entered
    float sampling = 0.05;
    
      //search threshold along the first component
      int k1=0;
      x_e1 = max_e1;
      while (L>threshold && mod2<=1. && !error)
      {
        x_e1 += sampling;
        mod2 = x_e1*x_e1+max_e2*max_e2;
        if(mod2 <= 1.)
        {
            L_e = loglikelihood(par, x_e1, max_e2, &error);
            if (!error)
            {
                
                // dividing each likelihood by exp(Lmax) the mean and variance doesn't change
                L=exp(L_e-maxL);
                k1++;
            }
        }
      }
 
      //search threshold along the second component
      L = 1.;
      int k2=0;
      mod2=0.;
      x_e2 = max_e2;
      while (L>threshold && mod2<=1. && !error)
      {
        x_e2 += sampling;
        mod2 = x_e2*x_e2+max_e1*max_e1;
        if(mod2 <= 1.)
        {
            L_e = loglikelihood(par, max_e1, x_e2, &error);
            if (!error)
            {
                
                // dividing each likelihood by exp(Lmax) the mean and variance doesn't change
                L=exp(L_e-maxL);
                k2++;
            }
        }
      }
   
    float edim1 = k1*sampling;
    float edim2 = k2*sampling;
    float start1 = edim1;
    float start2 = edim2;
    
    int k=0;
    while (oldsampling > sampling && sampling > 0.003)
    {
        for(ie1 = -start1/sampling; ie1 <= start1/sampling; ie1 +=2)
        {
            for(ie2 = -start2/sampling; ie2 <= start2/sampling; ie2 += 2)
            {
               if(ie1!=0 || ie2!=0)
               {
                 x_e1 = max_e1+ie1 * sampling;
                 x_e2 = max_e2+ie2 * sampling;
                 mod2 = x_e1*x_e1+x_e2*x_e2;
                 if(mod2 <= 1.)
                 {
                    L_e = loglikelihood(par, x_e1, x_e2, &error);
                    if (!error)
                    {
                        // dividing each likelihood by exp(Lmax) the mean and variance doesn't change
                        L=exp(L_e-maxL);
                        if (L >= threshold)
                        {
                            np++;
                            
                            xL1 = x_e1*L;
                            xL2 = x_e2*L;
                            average_1 += xL1;
                            average_2 += xL2;
                            var1 += x_e1*xL1;
                            var2 += x_e2*xL2;
                            cov += x_e1*xL2;
                            totL += L;
                        }
                        k++;
                    }
                 }
                
               }
            }
        }
        oldsampling = sampling;
        if (np < np_max)
        {
            sampling /= 2;
            start1 = edim1 - sampling;
            start2 = edim2 - sampling;
        }
    }
    
    average_1 /= totL;
    average_2 /= totL;
    var1 /= totL; var1 -= average_1*average_1;
    var2 /= totL; var2 -= average_2*average_2;
    cov /= totL; cov -= average_1*average_2;
    
    *mes_e1 = average_1;
    *mes_e2 = average_2;
    
    *var_e1 = var1;
    *var_e2 = var2;
    *cov_e = cov;
    
    printf("rank %d: n. points: %d, %d average: %f,%f variance: %e,%e cov: %e\n",rank,np,k,average_1,average_2,var1,var2,cov);
    
    if (var1 < 1e-4 || var2 <1e-4) return 1;
    else return 0;
                
}


#ifdef __cplusplus
}
#endif
