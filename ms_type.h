#include <ms/MeasurementSets.h>

#ifndef MEASUREMENT_SET_TYPEDEF_
#define MEASUREMENT_SET_TYPEDEF_

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

#endif /* MEASUREMENT_SET_TYPEDEF_ */
