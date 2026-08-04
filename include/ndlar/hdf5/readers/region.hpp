#pragma once

#include <cstddef>
#include <vector>

#include <hdf5.h>

#include "ndlar/hdf5/types.hpp"

namespace ndlar::hdf5 {

// Reader for ref_region datasets containing [start, stop) spans.
struct RawRefRegionReader {
    hid_t dset = H5I_INVALID_HID;
    hid_t mem_type = H5I_INVALID_HID;
    hid_t filespace = H5I_INVALID_HID;
    size_t row_count = 0;

    RawRefRegionReader(hid_t file_id, const char* dataset_path);
    ~RawRefRegionReader();

    bool read_rows(size_t first_idx, size_t count, std::vector<RefRegion>& out) const;
};

}  // namespace ndlar::hdf5
