#pragma once

#include <cstddef>
#include <vector>

#include <hdf5.h>

#include "ndlar/hdf5/types.hpp"

namespace ndlar::hdf5 {

// Reader for mc_truth/segments/data.
struct RawTrueSegmentReader {
    hid_t dset = H5I_INVALID_HID;
    hid_t mem_type = H5I_INVALID_HID;
    hid_t filespace = H5I_INVALID_HID;
    size_t row_count = 0;

    RawTrueSegmentReader(hid_t file_id, const char* dataset_path);
    ~RawTrueSegmentReader();

    bool read_rows(size_t first_idx, size_t count, std::vector<TrueSegment>& out) const;
};

}  // namespace ndlar::hdf5
