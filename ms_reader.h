#ifndef MS_READER_H_
#define MS_READER_H_

#include "ms_type.h"

#ifdef __cplusplus
extern "C" {
#endif
    
    /**
     * @brief Opens an existing Measurement Set.
     *
     * @details
     * Opens an existing Measurement Set.
     */
    RL_MeasurementSet* ms_open(const char* filename);
    
    /**
     * @brief Gets data from one column in a Measurement Set.
     *
     * @details
     * Gets data from one column in a Measurement Set.
     *
     * @param[in] p                     Pointer to opened Measurement Set.
     * @param[in] column                Name of required column in main table.
     * @param[in] start_row             Start row.
     * @param[in] num_rows              Number of rows to return.
     * @param[in] data_size_bytes       Data size of allocated block, in bytes.
     * @param[in,out] data              Data block to fill.
     * @param[out] required_size_bytes  Required size of the data block, in bytes.
     * @param[in,out] status            Status return code.
     */
     void ms_read_column(const RL_MeasurementSet* p, const char* column,
                              unsigned int start_row, unsigned int num_rows,
                              size_t data_size_bytes, void* data, size_t* required_size_bytes,
                              int* status);
    
    /**
     * @details
     * Reads baseline coordinate data from the main table.
     *
     * @details
     * This function reads a list of baseline coordinates from
     * the main table of the Measurement Set. The coordinate arrays must be
     * allocated to the correct size on entry.
     *
     * @param[in] start_row     The start row from which to read (zero-based).
     * @param[in] num_baselines Number of baselines to read from the main table.
     * @param[in,out] uu        Baseline u-coordinates, in metres.
     * @param[in,out] vv        Baseline v-coordinates, in metres.
     * @param[in,out] ww        Baseline w-coordinates, in metres.
     * @param[in,out] status    Status return code.
     * return max coordinate for gridding
     */
    double ms_read_coords(RL_MeasurementSet* p,
                                unsigned int start_row, unsigned int num_baselines,
                                double* uu, double* vv, double* ww, int* status);
    
    
    /**
     * @details
     * Reads visibility data from the main table.
     *
     * @details
     * This function reads a block of visibility data from the specified column of
     * the main table of the Measurement Set. The \p vis array must be
     * allocated to the correct size on entry.
     *
     * The dimensionality of the complex \p vis data block is:
     * (num_channels * num_baselines * num_pols)
     * with num_pols the fastest varying dimension, then num_baselines,
     * and num_channels the slowest.
     *
     * @param[in] start_row     The start row from which to read (zero-based).
     * @param[in] start_channel The start channel index to read (zero-based).
     * @param[in] num_channels  Number of channels to read.
     * @param[in] num_baselines Number of baselines to read from the main table.
     * @param[in] column        Name of column (DATA, MODEL_DATA or CORRECTED_DATA).
     * @param[in,out] vis       Visibility data.
     * @param[in,out] status    Status return code.
     */
    void ms_read_vis(RL_MeasurementSet* p,
                             unsigned int start_row, unsigned int start_channel,
                             unsigned int num_channels, unsigned int num_baselines,
                             const char* column, double* vis, int* status);
    
    /**
     * @brief Closes the Measurement Set.
     *
     * @details
     * Closes the Measurement Set and flushes any pending write operations to disk.
     */
    void ms_close(RL_MeasurementSet* p);
    
#ifdef __cplusplus
}
#endif

#endif /* MS_READER_H_ */
