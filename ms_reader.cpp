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
                                   TableLock(TableLock::PermanentLocking), Table::Update);
        
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


template<typename T>
void copy_array(const RL_MeasurementSet* p, const char* column,
                unsigned int start_row, unsigned int num_rows,
                size_t data_size_bytes, void* data, size_t* required_size,
                int* status)
{
    try
    {
        Slice slice(start_row, num_rows, 1);
        ArrayColumn<T> ac(*(p->ms), column);
        Array<T> a = ac.getColumnRange(slice);
        *required_size = a.size() * sizeof(T);
        if (data_size_bytes >= *required_size)
            memcpy(data, a.data(), *required_size);
        else
            *status = ERR_MS_OUT_OF_RANGE;
    }
    catch (...)
    {
        *status = ERR_MS_NO_DATA;
    }
}

template<typename T>
void copy_scalar(const RL_MeasurementSet* p, const char* column,
                 unsigned int start_row, unsigned int num_rows,
                 size_t data_size_bytes, void* data, size_t* required_size,
                 int* status)
{
    try
    {
        Slice slice(start_row, num_rows, 1);
        ScalarColumn<T> ac(*(p->ms), column);
        Array<T> a = ac.getColumnRange(slice);
        *required_size = a.size() * sizeof(T);
        if (data_size_bytes >= *required_size)
            memcpy(data, a.data(), *required_size);
        else
            *status = ERR_MS_OUT_OF_RANGE;
    }
    catch (...)
    {
        *status = ERR_MS_NO_DATA;
    }
}

void ms_read_column(const RL_MeasurementSet* p, const char* column,
                          unsigned int start_row, unsigned int num_rows,
                          size_t data_size_bytes, void* data, size_t* required_size_bytes,
                          int* status)
{
    if (*status || !p->ms) return;
    
    // Check that the column exists.
    if (!p->ms->tableDesc().isColumn(column))
    {
        *status = ERR_MS_COLUMN_NOT_FOUND;
        return;
    }
    
    // Check that some data are selected.
    if (num_rows == 0) return;
    
    // Check that the row is within the table bounds.
    unsigned int total_rows = p->ms->nrow();
    if (start_row >= total_rows)
    {
        *status = ERR_MS_OUT_OF_RANGE;
        return;
    }
    if (start_row + num_rows > total_rows)
        num_rows = total_rows - start_row;
    
    // Get column description and data type.
    const ColumnDesc& cdesc = p->ms->tableDesc().columnDesc(column);
    DataType dtype = cdesc.dataType();
    
    if (cdesc.isScalar())
    {
        switch (dtype)
        {
            case TpBool:
                copy_scalar<Bool>(p, column, start_row, num_rows,
                                  data_size_bytes, data, required_size_bytes, status); break;
            case TpChar:
                copy_scalar<Char>(p, column, start_row, num_rows,
                                  data_size_bytes, data, required_size_bytes, status); break;
            case TpUChar:
                copy_scalar<uChar>(p, column, start_row, num_rows,
                                   data_size_bytes, data, required_size_bytes, status); break;
            case TpShort:
                copy_scalar<Short>(p, column, start_row, num_rows,
                                   data_size_bytes, data, required_size_bytes, status); break;
            case TpUShort:
                copy_scalar<uShort>(p, column, start_row, num_rows,
                                    data_size_bytes, data, required_size_bytes, status); break;
            case TpInt:
                copy_scalar<Int>(p, column, start_row, num_rows,
                                 data_size_bytes, data, required_size_bytes, status); break;
            case TpUInt:
                copy_scalar<uInt>(p, column, start_row, num_rows,
                                  data_size_bytes, data, required_size_bytes, status); break;
            case TpFloat:
                copy_scalar<Float>(p, column, start_row, num_rows,
                                   data_size_bytes, data, required_size_bytes, status); break;
            case TpDouble:
                copy_scalar<Double>(p, column, start_row, num_rows,
                                    data_size_bytes, data, required_size_bytes, status); break;
            case TpComplex:
                copy_scalar<Complex>(p, column, start_row, num_rows,
                                     data_size_bytes, data, required_size_bytes, status); break;
            case TpDComplex:
                copy_scalar<DComplex>(p, column, start_row, num_rows,
                                      data_size_bytes, data, required_size_bytes, status); break;
            default:
                *status = ERR_MS_UNKNOWN_DATA_TYPE; break;
        }
    }
    else
    {
        switch (dtype)
        {
            case TpBool:
                copy_array<Bool>(p, column, start_row, num_rows,
                                 data_size_bytes, data, required_size_bytes, status); break;
            case TpChar:
                copy_array<Char>(p, column, start_row, num_rows,
                                 data_size_bytes, data, required_size_bytes, status); break;
            case TpUChar:
                copy_array<uChar>(p, column, start_row, num_rows,
                                  data_size_bytes, data, required_size_bytes, status); break;
            case TpShort:
                copy_array<Short>(p, column, start_row, num_rows,
                                  data_size_bytes, data, required_size_bytes, status); break;
            case TpUShort:
                copy_array<uShort>(p, column, start_row, num_rows,
                                   data_size_bytes, data, required_size_bytes, status); break;
            case TpInt:
                copy_array<Int>(p, column, start_row, num_rows,
                                data_size_bytes, data, required_size_bytes, status); break;
            case TpUInt:
                copy_array<uInt>(p, column, start_row, num_rows,
                                 data_size_bytes, data, required_size_bytes, status); break;
            case TpFloat:
                copy_array<Float>(p, column, start_row, num_rows,
                                  data_size_bytes, data, required_size_bytes, status); break;
            case TpDouble:
                copy_array<Double>(p, column, start_row, num_rows,
                                   data_size_bytes, data, required_size_bytes, status); break;
            case TpComplex:
                copy_array<Complex>(p, column, start_row, num_rows,
                                    data_size_bytes, data, required_size_bytes, status); break;
            case TpDComplex:
                copy_array<DComplex>(p, column, start_row, num_rows,
                                     data_size_bytes, data, required_size_bytes, status); break;
            default:
                *status = ERR_MS_UNKNOWN_DATA_TYPE; break;
        }
    }
}

double ms_read_coords(RL_MeasurementSet* p,
                          unsigned int start_row, unsigned int num_baselines,
                          double* uu, double* vv, double* ww, int* status)
{
    if (!p->ms || !p->msmc || num_baselines == 0) return 0.;
    
    // Check that the row is within the table bounds.
    unsigned int total_rows = p->ms->nrow();
    if (start_row >= total_rows)
    {
        *status = ERR_MS_OUT_OF_RANGE;
        return 0.;
    }
    if (start_row + num_baselines > total_rows)
        num_baselines = total_rows - start_row;
    
    // Read the coordinate data and copy it into the supplied arrays.
    Slice slice(start_row, num_baselines, 1);
    Array<Double> column_range = p->msmc->uvw().getColumnRange(slice);
    Matrix<Double> matrix;
    matrix.reference(column_range);

    double umax = 0.;
    double vmax = 0.;

    for (unsigned int i = 0; i < num_baselines; ++i)
    {
        uu[i] = matrix(0, i);
        vv[i] = matrix(1, i);
        ww[i] = matrix(2, i);
        if (fabs(uu[i]) > umax) umax = uu[i];
        if (fabs(vv[i]) > vmax) vmax = vv[i];
    }
    
    return ceil(fmax(umax,vmax));
}


void ms_read_vis(RL_MeasurementSet* p,
                       unsigned int start_row, unsigned int start_channel,
                       unsigned int num_channels, unsigned int num_baselines,
                       const char* column, double* vis, int* status)
{
    if (!p->ms || !p->msmc || num_baselines == 0 || num_channels == 0) return;
    
    // Check that the column exists.
    if (!p->ms->tableDesc().isColumn(column))
    {
        *status = ERR_MS_COLUMN_NOT_FOUND;
        return;
    }
    if (strcmp(column, "DATA") && strcmp(column, "CORRECTED_DATA") &&
        strcmp(column, "MODEL_DATA"))
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
    if (start_row + num_baselines > total_rows)
        num_baselines = total_rows - start_row;
    
    // Create the slicers for the column.
    unsigned int num_pols = p->num_pols;
    IPosition start1(1, start_row);
    IPosition length1(1, num_baselines);
    Slicer row_range(start1, length1);
    IPosition start2(2, 0, start_channel);
    IPosition length2(2, num_pols, num_channels);
    Slicer array_section(start2, length2);
    
    // Read the data.
    ArrayColumn<Complex> ac(*(p->ms), column);
    Array<Complex> column_range = ac.getColumnRange(row_range, array_section);
    
    // Copy the visibility data into the supplied array,
    // swapping baseline and channel dimensions.
    const float* in = (const float*) column_range.data();
    for (unsigned int c = 0; c < num_channels; ++c)
    {
        for (unsigned int b = 0; b < num_baselines; ++b)
        {
            for (unsigned int p = 0; p < num_pols; ++p)
            {
                unsigned int i = (num_pols * (b * num_channels + c) + p) << 1;
                unsigned int j = (num_pols * (c * num_baselines + b) + p) << 1;
                vis[j]     = in[i];
                vis[j + 1] = in[i + 1];
            }
        }
    }
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
