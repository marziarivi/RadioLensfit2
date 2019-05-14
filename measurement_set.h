#ifndef MEASUREMENT_SET_H_
#define MEASUREMENT_SET_H_

#ifdef __cplusplus
extern "C" {
#endif
    
//struct RL_MeasurementSet;
//#ifndef MEASUREMENT_SET_TYPEDEF_
//#define MEASUREMENT_SET_TYPEDEF_
//typedef struct RL_MeasurementSet RL_MeasurementSet;
//#endif /* MEASUREMENT_SET_TYPEDEF_ */
    

enum MS_ERROR_CODES
{
        ERR_MS_COLUMN_NOT_FOUND        = -200,
        ERR_MS_OUT_OF_RANGE            = -201,
        ERR_MS_UNKNOWN_DATA_TYPE       = -202,
        ERR_MS_NO_DATA                 = -203
};
    
enum MS_TYPES
{
        MS_UNKNOWN_TYPE = -1,
        MS_BOOL,
        MS_CHAR,
        MS_UCHAR,
        MS_SHORT,
        MS_USHORT,
        MS_INT,
        MS_UINT,
        MS_FLOAT,
        MS_DOUBLE,
        MS_COMPLEX,
        MS_DCOMPLEX
};
    
#ifdef __cplusplus
}
#endif

#include "ms_utils.h"
#include "ms_reader.h"

#endif /* MEASUREMENT_SET_H_ */
