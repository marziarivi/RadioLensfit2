#include "measurement_set.h"

#include <tables/Tables.h>
#include <casa/Arrays/Vector.h>

using namespace casacore;

size_t ms_column_element_size(const RL_MeasurementSet* p, const char* column)
{
    if (!p->ms || !p->ms->tableDesc().isColumn(column)) return 0;
    
    switch (p->ms->tableDesc().columnDesc(column).dataType())
    {
        case TpBool:     return sizeof(Bool);
        case TpChar:     return sizeof(Char);
        case TpUChar:    return sizeof(uChar);
        case TpShort:    return sizeof(Short);
        case TpUShort:   return sizeof(uShort);
        case TpInt:      return sizeof(Int);
        case TpUInt:     return sizeof(uInt);
        case TpFloat:    return sizeof(Float);
        case TpDouble:   return sizeof(Double);
        case TpComplex:  return sizeof(Complex);
        case TpDComplex: return sizeof(DComplex);
        default:         return 0;
    }
    return 0;
}

int ms_column_element_type(const RL_MeasurementSet* p, const char* column)
{
    if (!p->ms || !p->ms->tableDesc().isColumn(column))
        return MS_UNKNOWN_TYPE;
    
    switch (p->ms->tableDesc().columnDesc(column).dataType())
    {
        case TpBool:     return MS_BOOL;
        case TpChar:     return MS_CHAR;
        case TpUChar:    return MS_UCHAR;
        case TpShort:    return MS_SHORT;
        case TpUShort:   return MS_USHORT;
        case TpInt:      return MS_INT;
        case TpUInt:     return MS_UINT;
        case TpFloat:    return MS_FLOAT;
        case TpDouble:   return MS_DOUBLE;
        case TpComplex:  return MS_COMPLEX;
        case TpDComplex: return MS_DCOMPLEX;
        default:         return MS_UNKNOWN_TYPE;
    }
    return MS_UNKNOWN_TYPE;
}

size_t* ms_column_shape(const RL_MeasurementSet* p, const char* column, size_t* ndim)
{
    size_t i = 0, *t = 0;
    if (!p->ms || !p->ms->tableDesc().isColumn(column)) return 0;
    
    const ColumnDesc& cdesc = p->ms->tableDesc().columnDesc(column);
    const IPosition& shape = cdesc.shape();
    if (shape.size() > 0)
    {
        *ndim = (int) shape.size();
        t = (size_t*) calloc(*ndim, sizeof(size_t));
        for (i = 0; i < *ndim; ++i) t[i] = shape(i);
    }
    else if (p->ms->nrow() > 0)
    {
        // If shape is not fixed, return shape of first cell instead.
        TableColumn tc(*(p->ms), column);
        IPosition shape = tc.shape(0);
        if (shape.size() > 0)
        {
            *ndim = (int) shape.size();
            t = (size_t*) calloc(*ndim, sizeof(size_t));
            for (i = 0; i < *ndim; ++i) t[i] = shape(i);
        }
    }
    return t;
}

void ms_ensure_num_rows(RL_MeasurementSet* p, unsigned int num)
{
    if (!p->ms) return;
    int rows_to_add = (int)num - (int)(p->ms->nrow());
    if (rows_to_add > 0)
        p->ms->addRow((unsigned int)rows_to_add);
}

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

double ms_time_start_mjd_utc(const RL_MeasurementSet* p)
{
    return p->start_time / 86400.0;
}

