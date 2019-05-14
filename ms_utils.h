#ifndef MS_UTILS_H_
#define MS_UTILS_H_

#include <stddef.h>
#include "ms_type.h"

#ifdef __cplusplus
extern "C" {
#endif
    
    /**
     * @brief
     * Returns the size in bytes of an element in a column of the Measurement Set.
     *
     * @details
     * Returns the size in bytes of an element in a column of the Measurement Set.
     *
     * @param[in] column     Column name.
     */
    size_t ms_column_element_size(const RL_MeasurementSet* p,
                                        const char* column);
    
    /**
     * @brief
     * Returns the data type of an element in a column of the Measurement Set.
     *
     * @details
     * Returns the data type of an element in a column of the Measurement Set.
     *
     * This is one of the values from the OSKAR_MS_TYPE enumerator.
     *
     * @param[in] column     Column name.
     */
    int ms_column_element_type(const RL_MeasurementSet* p,
                                     const char* column);
    
    /**
     * @brief
     * Returns the shape of an element in a column of the Measurement Set.
     *
     * @details
     * Returns the shape of an element in a column of the Measurement Set.
     *
     * Note that the returned array must be freed by the caller using free().
     *
     * @param[in] column     Column name.
     * @param[out] ndim      Number of dimensions.
     *
     * @return Array containing the size of each dimension.
     * Must be freed by the caller using free().
     */
    size_t* ms_column_shape(const RL_MeasurementSet* p, const char* column,
                                  size_t* ndim);
    
    /**
     * @brief
     * Ensures the specified number of rows exist in the Measurement Set.
     *
     * @details
     * Ensures the specified number of rows exist in the Measurement Set,
     * adding extra ones if necessary.
     *
     * @param[in] num    Total number of rows in the Measurement Set.
     */
    void ms_ensure_num_rows(RL_MeasurementSet* p, unsigned int num);
    
    /**
     * @brief
     * Returns the channel separation in the Measurement Set.
     *
     * @details
     * Returns the channel separation in the Measurement Set.
     */
    double ms_freq_inc_hz(const RL_MeasurementSet* p);
    
    /**
     * @brief
     * Returns the reference frequency in the Measurement Set.
     *
     * @details
     * Returns the reference frequency in the Measurement Set.
     */
    double ms_freq_start_hz(const RL_MeasurementSet* p);
    
    /**
     * @brief
     * Returns the number of channels in the Measurement Set.
     *
     * @details
     * Returns the number of channels in the Measurement Set.
     */
     unsigned int ms_num_channels(const RL_MeasurementSet* p);
    
    /**
     * @brief
     * Returns the number of polarisations in the Measurement Set.
     *
     * @details
     * Returns the number of polarisations in the Measurement Set.
     */
    unsigned int ms_num_pols(const RL_MeasurementSet* p);
    
    /**
     * @brief
     * Returns the number of rows in the main table.
     *
     * @details
     * Returns the number of rows in the main table.
     */
    unsigned int ms_num_rows(const RL_MeasurementSet* p);
    
    /**
     * @brief
     * Returns the number of stations in the Measurement Set.
     *
     * @details
     * Returns the number of stations in the Measurement Set.
     */
    unsigned int ms_num_stations(const RL_MeasurementSet* p);
    
    /**
     * @brief
     * Returns the phase centre RA in the Measurement Set.
     *
     * @details
     * Returns the phase centre RA in the Measurement Set.
     */
    double ms_phase_centre_ra_rad(const RL_MeasurementSet* p);
    
    /**
     * @brief
     * Returns the phase centre Dec in the Measurement Set.
     *
     * @details
     * Returns the phase centre Dec in the Measurement Set.
     */
    double ms_phase_centre_dec_rad(const RL_MeasurementSet* p);
    
    
    /**
     * @brief
     * Returns the time increment in the Measurement Set.
     *
     * @details
     * Returns the time increment in the Measurement Set.
     */
     double ms_time_inc_sec(const RL_MeasurementSet* p);
    
    /**
     * @brief
     * Returns the start time in the Measurement Set.
     *
     * @details
     * Returns the start time in the Measurement Set.
     */
    double ms_time_start_mjd_utc(const RL_MeasurementSet* p);
    
#ifdef __cplusplus
}
#endif

#endif /* MS_UTILS_H_ */
