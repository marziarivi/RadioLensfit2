/*
 * Copyright (c) 2020 Lance Miller, Marzia Rivi
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

#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "datatype.h"
#include "marginalise_r.h"


void set_posterior_values(int numR, double* L_r, double* rprior, double* Ro, double* xmarvals, double* ymarvals, int* numvals)
{
  int num_marvals = 0;    

  // find first valid measurement
  int nRo = 1;
  while (nRo < numR && ( L_r[nRo] <= -1.e+10 ) ) nRo++;
    
  /*
   * if found OK, make linearly interpolated array at intervals R/5 from  R/10 to
   * the first measured value R using zerolikel as the r=0 value.  As the prior
   * has not been defined on this finer grid, and as the log(posterior) value
   * at r=0 is -infinity, interpolate the likelihood and then apply the prior.
   */

  if (nRo < numR)
  {
    /*
     * go through the measured points, in order of increasing R, adding into the
     * interpolation array
     */
    double likelval;
    for(nRo = 1; nRo < numR; nRo++)
    {
      likelval = L_r[nRo];
      // if likelihood previously measured successfully...
      if (likelval > -1.e+10)
      {
          xmarvals[num_marvals] = Ro[nRo];
          ymarvals[num_marvals] = likelval + rprior[nRo];
          num_marvals++;
      }
    }
  }
  *numvals = num_marvals;
}
 
     
/*
 * Marginalise likelihood over scalelength
 * 1. find interval of integration
 * 2. interpolate log(likelihood) points with a cubic spline
 * 3. integrate using QAG adaptive integration
 */


// GSL interpolation functions and integration workspace
gsl_interp_accel *facc;
gsl_spline *fspline;
int gsl_work_size={1000};
gsl_integration_workspace *gsl_work;

double marginalise_posterior_r(int num_marvals, double *xmarvals, double *ymarvals)
{
    double sump, error;
    int gsl_status;
    
    gsl_set_error_handler_off();
    
    //numerical integration function kernels
    gsl_function marF;
    marF.function = &marf;
    //marF.params = &nt;
    
    // find the position of the maximum
    double ymax = -10.0;
    int ymaxpos = 0;
    int ii;
    
    for(ii = 0; ii < num_marvals; ii++)
    {
        if (ymarvals[ii] > ymax)
        {
            ymax = ymarvals[ii];
            ymaxpos = ii;
        }
    }
    // set integration limits
    double xmin, xmax;
    xmin = xmax = xmarvals[ymaxpos];
    if (ymaxpos < num_marvals-1) xmax = xmarvals[ymaxpos+1];
    if (ymaxpos > 0) xmin = xmarvals[ymaxpos-1];
    for (ii=0; ii < num_marvals; ii++)
    {
        {
            if (xmarvals[ii] < xmin) xmin = xmarvals[ii];
            if (xmarvals[ii] > xmax) xmax = xmarvals[ii];
        }
        ymarvals[ii] -= ymax;  // trick to integrate a very large values of the likelihood
    }
    
  // check number of interpolation points and integration limits
  if (xmax <= xmin)
  {
    fflush(stdout);
    fprintf(stderr, " error in function \n");
    fprintf(stderr, " num = %d xmin,max = %f %f \n", num_marvals, xmin, xmax);
    exit(EXIT_FAILURE);
  }

  // allocate the interpolation objects (in the gsl routine this object
  // must be equal in dimension to the data array so we cannot simply keep
  // reusing a single object, sadly)
    
  // allocate spline interpolation for marginalisation
  facc = gsl_interp_accel_alloc();
  if (facc == NULL)
  {
     fflush(stdout);
     fprintf(stderr, " error from gsl_interp_accel_alloc, facc \n");
     exit(EXIT_FAILURE);
  }
    
  // allocate the workspace for the integrator
  gsl_work = gsl_integration_workspace_alloc(gsl_work_size);
  if (gsl_work == NULL)
  {
     fflush(stdout);
     fprintf(stderr, " error from gsl_integration_workspace, work \n");
     exit(EXIT_FAILURE);
  }
    
  // cspline for log(posterior) values
  fspline = gsl_spline_alloc(gsl_interp_cspline,num_marvals);
  if (fspline == NULL)
  {
    fflush(stdout);
    fprintf(stderr," error from gsl_interp_accel_alloc, fspline \n");
    exit(EXIT_FAILURE);
  }
 
  // load R and log(posterior) values into interpolation arrays
  gsl_status = gsl_spline_init(fspline,xmarvals,ymarvals,num_marvals);
  if (gsl_status)
  {
    fflush(stdout);
    fprintf(stderr," error gsl_spline_init fspline \n");
    fprintf(stderr," %s \n",gsl_strerror(gsl_status));
    exit(EXIT_FAILURE);
  }

  // integrate the exponential of the interpolated log(posterior) array
  gsl_status = gsl_integration_qag(&marF, xmin, xmax,0., 1.e-3, gsl_work_size,
                                 GSL_INTEG_GAUSS41,gsl_work, &sump, &error);
 
  if (gsl_status)
  {
    // occasionally cannot integrate, try again with linear interpolation
    gsl_spline_free(fspline);
    fspline = gsl_spline_alloc(gsl_interp_linear,num_marvals);
    if (fspline == NULL)
    {
        fflush(stdout);
        fprintf(stderr," error from linear gsl_interp_accel_alloc, fspline \n");
        exit(EXIT_FAILURE);
    }
    // reset the accelerator
    gsl_status = gsl_interp_accel_reset(facc);
    if (gsl_status)
    {
        fflush(stdout);
        fprintf(stderr," error linear gsl_interp_accel_reset facc \n");
        fprintf(stderr," %s \n",gsl_strerror(gsl_status));
        exit(EXIT_FAILURE);
    }
    // load R and log(posterior) values into interpolation arrays
    gsl_status = gsl_spline_init(fspline,xmarvals,ymarvals,num_marvals);
    if (gsl_status)
    {
        fflush(stdout);
        fprintf(stderr," error linear gsl_spline_init fspline \n");
        fprintf(stderr," %s \n",gsl_strerror(gsl_status));
        exit(EXIT_FAILURE);
    }
    // integrate the exponential of the interpolated log(posterior) array
    gsl_status = gsl_integration_qag(&marF, xmin, xmax,0., 1.e-3, gsl_work_size,
                                         GSL_INTEG_GAUSS41,gsl_work, &sump, &error);
  }
    
  // free the gsl objects
  gsl_spline_free(fspline);
  gsl_interp_accel_free(facc);
  gsl_integration_workspace_free(gsl_work);
    
  if (gsl_status)
  {
      // return a null value, no marginalised values will be computed at this point
      fflush(stdout);
      fprintf(stderr," cannot integrate interpolated posterior \n");
      fprintf(stderr," %s \n",gsl_strerror(gsl_status));
      return -1.e10;
  }

  return log(sump)+ymax;   // return log(L), trick to deal with a very large values of the likelihood
}
    
    
/* integration kernel used by marginalisation step for integrating posterior */
double marf(double x, void *params)
{
    // interpolate log(posterior)
    double yval;
    int gsl_status = gsl_spline_eval_e(fspline, x, facc, &yval);
    if (gsl_status)
        yval = 0.;
    else
        // return exponentiated function
        yval = exp(yval);
    return yval;
}
