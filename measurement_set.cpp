#include "measurement_set.h"

#include <tables/Tables.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>

#include <cstdlib>
#include <cstdio>
#include <cmath>


using namespace casacore;

RL_MeasurementSet* ms_open(const char* filename)
{
    RL_MeasurementSet* p = (RL_MeasurementSet*) calloc(1, sizeof(RL_MeasurementSet));

    try
    {
        // Create the MeasurementSet. Storage managers are recreated as needed.
        p->ms = new MeasurementSet(filename,
                                   TableLock(TableLock::NoLocking), Table::Update);

        // Create the MSMainColumns and MSColumns objects for accessing data
        // in the main table and subtables.
        p->msc = new MSColumns(*(p->ms));
        p->msmc = new MSMainColumns(*(p->ms));
    }
    catch (AipsError& e)
    {
        fprintf(stderr, "Caught AipsError: %s\n", e.what());
        fflush(stderr);
        ms_close(p);
        return 0;
    }

    // Refuse to open if there is more than one spectral window.
    if (p->ms->spectralWindow().nrow() != 1)
    {
        fprintf(stderr, "RadioLensfit can read Measurement Sets with one spectral "
                "window only. Use 'split' or 'mstransform' in CASA to select "
                "the spectral window first.\n");
        fflush(stderr);
        ms_close(p);
        return 0;
    }

    // Refuse to open if there is more than one field.
    if (p->ms->field().nrow() != 1)
    {
        fprintf(stderr, "RadioLensfit can read Measurement Sets with one target "
                "field only. Use 'split' or 'mstransform' in CASA to select "
                "the target field first.\n");
        fflush(stderr);
        ms_close(p);
        return 0;
    }

    // Get the data dimensions.
    p->num_pols = 0;
    p->num_channels = 0;
    if (p->ms->polarization().nrow() > 0)
        p->num_pols = p->msc->polarization().numCorr().get(0);
    if (p->num_pols != 1)
    {
        fprintf(stderr, "RadioLensfit can only read Measurement Sets with a sigle polarization (I stokes)\n");
        fflush(stderr);
        ms_close(p);
        return 0;
    }

    if (p->ms->spectralWindow().nrow() > 0)
    {
        p->num_channels = p->msc->spectralWindow().numChan().get(0);
        p->freq_start_hz = p->msc->spectralWindow().refFrequency().get(0);
        p->freq_inc_hz = (p->msc->spectralWindow().chanWidth().get(0))(IPosition(1, 0));
    }
    p->num_stations = p->ms->antenna().nrow();
    if (p->ms->nrow() > 0)
        p->time_inc_sec = p->msc->interval().get(0);

    // Get the phase centre.
    p->phase_centre_ra = 0.0;
    p->phase_centre_dec = 0.0;
    if (p->ms->field().nrow() > 0)
    {
        Vector<MDirection> dir;
        p->msc->field().phaseDirMeasCol().get(0, dir, true);
        if (dir.size() > 0)
        {
            Vector<Double> v = dir(0).getAngle().getValue();
            p->phase_centre_ra = v(0);
            p->phase_centre_dec = v(1);
        }
    }

    // Get the time range.
    Vector<Double> range(2, 0.0);
    if (p->msc->observation().nrow() > 0)
        p->msc->observation().timeRange().get(0, range);
    p->start_time = range[0];
    p->end_time = range[1];

    return p;
}


void ms_close(RL_MeasurementSet* p)
{
    if (!p) return;
    if (p->msmc)
        delete p->msmc;
    if (p->msc)
        delete p->msc;
    if (p->ms)
        delete p->ms;
    free(p->a1);
    free(p->a2);
    free(p);
}


// Ensures the specified number of rows exist in the Measurement Set, adding extra ones if necessary.
void ms_ensure_num_rows(RL_MeasurementSet* p, unsigned int num)
{
    if (!p->ms) return;
    int rows_to_add = (int)num - (int)(p->ms->nrow());
    if (rows_to_add > 0)
        p->ms->addRow((unsigned int)rows_to_add);
}

// Return the channel bandwidth in Hz 
double ms_freq_inc_hz(const RL_MeasurementSet* p)
{
    return p->freq_inc_hz;
}

double ms_freq_start_hz(const RL_MeasurementSet* p)
{
    return p->freq_start_hz;
}

unsigned int ms_num_channels(const RL_MeasurementSet* p)
{
    return p->num_channels;
}

unsigned int ms_num_pols(const RL_MeasurementSet* p)
{
    return p->num_pols;
}

// Returns number of rows in the main table
unsigned int ms_num_rows(const RL_MeasurementSet* p)
{
    if (!p->ms) return 0;
    return p->ms->nrow();
}

unsigned int ms_num_stations(const RL_MeasurementSet* p)
{
    return p->num_stations;
}

double ms_phase_centre_ra_rad(const RL_MeasurementSet* p)
{
    return p->phase_centre_ra;
}

double ms_phase_centre_dec_rad(const RL_MeasurementSet* p)
{
    return p->phase_centre_dec;
}

double ms_time_inc_sec(const RL_MeasurementSet* p)
{
    return p->time_inc_sec;
}

/*
   Reads a list of baseline coordinates from the main table of the Measurement Set. 
   The coordinate arrays must be allocated to the correct size on entry.
   Return the max coordinate for gridding
*/
double ms_read_coords(RL_MeasurementSet* p,
                          unsigned int start_row, unsigned int num_coords,
                          double* uu, double* vv, double* ww, int* status)
{
    if (!p->ms || !p->msmc || num_coords == 0) return 0.;

    // Check that the row is within the table bounds.
    unsigned int total_rows = p->ms->nrow();
    if (start_row >= total_rows)
    {
        *status = ERR_MS_OUT_OF_RANGE;
        return 0.;
    }
    if (start_row + num_coords > total_rows)
        num_coords = total_rows - start_row;

    // Read the coordinate data and copy it into the supplied arrays.
    Slice slice(start_row, num_coords, 1);
    Array<Double> column_range = p->msmc->uvw().getColumnRange(slice);
    Matrix<Double> matrix;
    matrix.reference(column_range);

    double umax = 0.;
    double vmax = 0.;

    for (unsigned int i = 0; i < num_coords; ++i)
    {
        uu[i] = matrix(0, i);
        vv[i] = matrix(1, i);
        ww[i] = matrix(2, i);
        if (fabs(uu[i]) > umax) umax = uu[i];
        if (fabs(vv[i]) > vmax) vmax = vv[i];
    }

    return ceil(fmax(umax,vmax));
}


/*
   Reads noise sigtma from the main table of the Measurement Set.
   The sigma array must be allocated to the correct size on entry.
*/

void ms_read_sigma(RL_MeasurementSet* p,
                      unsigned int start_row, unsigned int num_coords,
                      float* sigma2, int* status)
{
    if (!p->ms || !p->msmc || num_coords == 0) return;

    // Check that the row is within the table bounds.
    unsigned int total_rows = p->ms->nrow();
    if (start_row >= total_rows)
    {
        *status = ERR_MS_OUT_OF_RANGE;
        return;
    }
    if (start_row + num_coords > total_rows)
        num_coords = total_rows - start_row;

    // Read the coordinate data and copy it into the supplied arrays.
    Slice slice(start_row, num_coords, 1);
    Array<Float> column_range = p->msmc->sigma().getColumnRange(slice);

    const float* in = (const float*) column_range.data();
    for (unsigned int i = 0; i < num_coords; ++i)
    {
        sigma2[i] = in[i]*in[i];
    }
}


/*
 *    Read a block of flags from the specified column of the main table of t
 *    he Measurement Set.
 *       
 *    It is assumed the MS contains a single polarization component: the flux I.
 *    The dimensionality of the boolean flag block is: num_channels * num_coords 
 *    with num_coords the fastest varying dimension, and num_channels the slowest.
 *    The flag array must be allocated to the correct size on entry. 
 */
void ms_read_Flag(RL_MeasurementSet* p,
		  unsigned int start_row, unsigned int start_channel,
	          unsigned int num_channels, unsigned int num_coords,
	          const char* column, bool* flag, int* status)
{
  if (!p->ms || !p->msmc || num_coords == 0 || num_channels == 0) return;
  
  // Check that the column exists.
  if (!p->ms->tableDesc().isColumn(column))
  {
    *status = ERR_MS_COLUMN_NOT_FOUND;
    return;
  }

  // Check that the row is within the table bounds.
  unsigned int total_rows = p->ms->nrow();
  if (start_row >= total_rows)
  {
    *status = ERR_MS_OUT_OF_RANGE;
    return;
  }
  if (start_row + num_coords > total_rows)
    num_coords = total_rows - start_row;

    // Create the slicers for the column.
    unsigned int num_pols = 1;
    IPosition start1(1, start_row);
    IPosition length1(1, num_coords);
    Slicer row_range(start1, length1);
    IPosition start2(2, 0, start_channel);
    IPosition length2(2, num_pols, num_channels);
    Slicer array_section(start2, length2);

    // Read the data.
    ArrayColumn<Bool> ac(*(p->ms), column);
    Array<Bool> column_range = ac.getColumnRange(row_range, array_section);

    // Copy the visibility data into the supplied array,
    // swapping coords and channel dimensions.
    const bool* in = (const bool*) column_range.data();
    for (unsigned int c = 0; c < num_channels; ++c)
    {
       for (unsigned int b = 0; b < num_coords; ++b)
       {
         unsigned int i = (b * num_channels + c);
         unsigned int j = (c * num_coords + b);
         flag[j] = in[i];
       }
    }
}


/*
   Read a block of visibilities from the specified column of the main table of the Measurement Set.
   
   It is assumed the MS contains a single stokes component: the flux I.
   The dimensionality of the complex vis data block is: num_channels * num_coords 
   with num_coords the fastest varying dimension, and num_channels the slowest.
   The vis array must be allocated to the correct size on entry. 
*/
void ms_read_vis(RL_MeasurementSet* p,
        unsigned int start_row, unsigned int start_channel,
        unsigned int num_channels, unsigned int num_coords,
        const char* column, complexd* vis, int* status)
{
    if (!p->ms || !p->msmc || num_coords == 0 || num_channels == 0) return;

    // Check that the column exists.
    if (!p->ms->tableDesc().isColumn(column))
    {
        *status = ERR_MS_COLUMN_NOT_FOUND;
        return;
    }

    // Check that the row is within the table bounds.
    unsigned int total_rows = p->ms->nrow();
    if (start_row >= total_rows)
    {
        *status = ERR_MS_OUT_OF_RANGE;
        return;
    }
    if (start_row + num_coords > total_rows)
        num_coords = total_rows - start_row;

    // Create the slicers for the column.
    unsigned int num_pols = 1;
    IPosition start1(1, start_row);
    IPosition length1(1, num_coords);
    Slicer row_range(start1, length1);
    IPosition start2(2, 0, start_channel);
    IPosition length2(2, num_pols, num_channels);
    Slicer array_section(start2, length2);

    // Read the data.
    ArrayColumn<Complex> ac(*(p->ms), column);
    Array<Complex> column_range = ac.getColumnRange(row_range, array_section);

    // Copy the visibility data into the supplied array,
    // swapping coords and channel dimensions.
    const float* in = (const float*) column_range.data();
    for (unsigned int c = 0; c < num_channels; ++c)
    {
        for (unsigned int b = 0; b < num_coords; ++b)
        {
                unsigned int i = (b * num_channels + c) << 1;
                unsigned int j = (c * num_coords + b);
                vis[j].real = in[i];
                vis[j].imag = in[i + 1];
        }
    }
}


/* 
   Writes the given block of visibility data to the data column of the Measurement Set, 
   extending it if necessary.
  
   It is assumed to write on a MS with a single polarization column.
   The dimensionality of the complex vis data block is: num_channels * num_coords,
   with num_coords the fastest varying dimension and num_channels the slowest.
*/
void ms_write_vis(RL_MeasurementSet* p,
        unsigned int start_row, unsigned int start_channel,
        unsigned int num_channels, unsigned int num_coords, const complexd* vis)
{
    MSMainColumns* msmc = p->msmc;
    if (!msmc) return;

    // Allocate storage for the block of visibility data.
    unsigned int num_pols = 1;
    IPosition shape(3, num_pols, num_channels, num_coords);
    Array<Complex> vis_data(shape);

    // Copy visibility data into the array,
    // swapping baseline and channel dimensions.
    float* out = (float*) vis_data.data();
    for (unsigned int c = 0; c < num_channels; ++c)
    {
        for (unsigned int b = 0; b < num_coords; ++b)
        {
                unsigned int i = (c * num_coords + b);
                unsigned int j = (b * num_channels + c) << 1;
                out[j]     = vis[i].real;
                out[j + 1] = vis[i].imag;
        }
    }

    // Add new rows if required.
    ms_ensure_num_rows(p, start_row + num_coords);

    // Create the slicers for the column.
    IPosition start1(1, start_row);
    IPosition length1(1, num_coords);
    Slicer row_range(start1, length1);
    IPosition start2(2, 0, start_channel);
    IPosition length2(2, num_pols, num_channels);
    Slicer array_section(start2, length2);

    // Write visibilities to DATA column.
    ArrayColumn<Complex>& col_data = msmc->data();
    col_data.putColumnRange(row_range, array_section, vis_data);
}

