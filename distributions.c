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


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_sf_erf.h>

#include "distributions.h"
#include "default_params.h"

#ifdef __cplusplus
extern "C" {
#endif
    
// probability distribution function of the ellipticity modulus
// parameters fitted from VLA-COSMOS
double e_pdf(double e)
{
    const double e_max = E_MAX;   // ellipticity cutoff
    const double e_0 = E_0;       // circularity parameter
    const double a = DISP;        // dispersion
        
    const double A = NORM_E;      // normalization factor 
        
    double p_e = A*e*(1.-exp((e-e_max)/a))/((1.+e)*sqrt(e*e+e_0*e_0));
    return p_e;
}

// Cumulative Distribution Function of a power function
// flux prior is a power law whose CDF is well know: x^{alpha + 1}/(alpha + 1)
double flux_CDF(double alpha, double flux)
{
    const double a = alpha+1.;
    const double norm = NORM_S;  // normalisation factor over 1 square degree
    // double norm = 1.37;  // norm = 1/729878 normalization factor referred to an area of 1600 arcmin^2
    // norm *= 2.25;  // = 3600/1600, to obtain % of N(<flux) per square degree
    return norm*pow(flux,a)/a;
}

//power law relation between flux and scalelength
double scale_mean(double flux)
{
    return ADD+ESP*log(flux);
}


/* prior function for scalelength as a function of flux, derived from Deep SWIRE 20cm continuum catalog
 
 NB returns ln(prior)
 
 parameters assume r is measured in arcsec
 */

double rfunc (const double mu, const double sigma, double r)
{
    double logprior;
    if (r<=0.)
    {
        fflush(stdout);
        fprintf(stderr," illegal value given to rfunc: r=%f \n",r);
        exit(EXIT_FAILURE);
    }
    
    //lognormal distribution
    double x = (log(r) - mu)/sigma;
    logprior = -(x*x)*0.5 - log(r*sigma);
    
    return logprior;
}


/* Cumulative Distribution Function of a Lognormal
   scalelength prior is a lognormal distribution whose CDF is well know: 1/2+erf((ln(x)-mu)/(sqrt(2)*sigma))/2
*/

double r_CDF(double mean, double r)
{
    //lognormal distribution
    double x = (log(r) - mean)/(sqrt(2.)*R_STD);
    
    double F_r = (1.+gsl_sf_erf(x))*0.5;
    return F_r;
}


// Cumulative Distribution Function of a given probability density function: \int_0^b pdf(x)
// apply an extended trapezoidal rule (Numerical Recipes p.137)
double CDF(double (*pdf)(double), double b)
{
    double x, tnm, del, sum;
    double s, olds = 0.;
    int it,j,k;

    s=0.5*b*(*pdf)(b);
    for (k=2; k<=JMAX; k++)
    {
        for (it=1,j=1;j<k-1;j++) it <<= 1;
        tnm = (double) it;
        del = b/tnm;      //This is the spacing of the points to be added.
        x = 0.5*del;
        for (sum=0.,j=1; j<it; j++,x+=del) sum += (*pdf)(x);
        s = 0.5*(s+del*sum);  //This replaces s by its refined value. NB. pdf(0)=0

        if (k > 5)          //Avoid spurious early convergence.
            if (fabs(s-olds) < EPS*fabs(olds) ||
                (s == 0.0 && olds == 0.0)) return s;
        olds=s;
    }
    printf("\n CDF: Too many steps\n");
    return 0.0;
}
    
#ifdef __cplusplus
}
#endif


