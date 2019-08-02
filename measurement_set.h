#ifndef MEASUREMENT_SET_H_
#define MEASUREMENT_SET_H_

#include <ms/MeasurementSets.h>
#include <stddef.h>
#include "datatype.h"

enum MS_ERROR_CODES
{
        ERR_MS_COLUMN_NOT_FOUND        = -1,
        ERR_MS_OUT_OF_RANGE            = -2,
        ERR_MS_NO_DATA                 = -3
};

struct RL_MeasurementSet
{
    casacore::MeasurementSet* ms;   // Pointer to the Measurement Set.
    casacore::MSColumns* msc;       // Pointer to the sub-tables.
    casacore::MSMainColumns* msmc;  // Pointer to the main columns.
    unsigned int *a1, *a2;
    unsigned int num_pols, num_channels, num_stations;
    double freq_start_hz, freq_inc_hz;
    double phase_centre_ra, phase_centre_dec;
    double start_time, end_time, time_inc_sec;
};

typedef struct RL_MeasurementSet RL_MeasurementSet;

RL_MeasurementSet* ms_open(const char* filename);
void ms_close(RL_MeasurementSet* p);

void ms_ensure_num_rows(RL_MeasurementSet* p, unsigned int num);
double ms_freq_inc_hz(const RL_MeasurementSet* p);
double ms_freq_start_hz(const RL_MeasurementSet* p);
unsigned int ms_num_channels(const RL_MeasurementSet* p);
unsigned int ms_num_pols(const RL_MeasurementSet* p);
unsigned int ms_num_rows(const RL_MeasurementSet* p);
unsigned int ms_num_stations(const RL_MeasurementSet* p);
double ms_phase_centre_ra_rad(const RL_MeasurementSet* p);
double ms_phase_centre_dec_rad(const RL_MeasurementSet* p);
double ms_time_inc_sec(const RL_MeasurementSet* p);

double ms_read_coords(RL_MeasurementSet* p,
                      unsigned int start_row, unsigned int num_coords,
                      double* uu, double* vv, double* ww, int* status);

void ms_read_sigma(RL_MeasurementSet* p,
                      unsigned int start_row, unsigned int num_coords,
                      float* sigma2, int* status);

void ms_read_Flag(RL_MeasurementSet* p,
		  unsigned int start_row, unsigned int start_channel,
                  unsigned int num_channels, unsigned int num_coords,
		  const char* column, bool* flag, int* status);

void ms_read_vis(RL_MeasurementSet* p,
                 unsigned int start_row, unsigned int start_channel,
                 unsigned int num_channels, unsigned int num_coords,
                 const char* column, complexd* vis, int* status);
void ms_write_vis(RL_MeasurementSet* p,
                  unsigned int start_row, unsigned int start_channel,
                  unsigned int num_channels, unsigned int num_coords,
                  const complexd* vis);


#endif /* MEASUREMENT_SET_H_ */
