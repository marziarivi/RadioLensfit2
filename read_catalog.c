
/*
 * Copyright (c) 2020 Marzia Rivi
 *
 * This file is part of RadioLensfit2.
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
#include "read_catalog.h"

#ifdef __cplusplus
extern "C" {
#endif
    
int read_catalog(unsigned long int nge, char *filename, double *gflux, double *gscale, double *ge1, double *ge2, double *l, double *m, double *SNR_vis, bool readSNR)
{
    FILE *fp;
    char *token;
    char line[1000];
    
    fp = fopen(filename,"r");
    if (!fp)
    {
        printf("ERROR: Unable to open the file %s\n", filename);
        exit(EXIT_FAILURE);
        return 1;
    }
    // SNR | l | m | flux | scale | e1 | e2 (last 3 params are provided for simulations)
    double SNR,ll,mm,flux;
    double scale,e1,e2;
    
    unsigned long int g = 0;
    while (!feof(fp) && g<nge)
    {
        fgets(line, 1000, fp);
        if (readSNR)  
        { 
          sscanf(line, "%lf %lf %lf %lf %lf %lf %lf",&SNR,&ll,&mm,&flux,&scale, &e1, &e2);
          SNR_vis[g] = SNR;
        }
        else sscanf(line, "%lf %lf %lf %lf %lf %lf",&ll,&mm,&flux,&scale, &e1, &e2);

        l[g] = ll;
        m[g] = mm;
        gflux[g] = flux;
        gscale[g] = scale;
        ge1[g] = e1;
        ge2[g] = e2;
        
        g++;
        
    }
    fclose(fp);
    return g;
}

#ifdef __cplusplus
}
#endif
