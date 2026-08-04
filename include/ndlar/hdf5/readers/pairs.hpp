#pragma once

#include <cstddef>
#include <vector>

#include <hdf5.h>

#include "ndlar/hdf5/types.hpp"

namespace ndlar::hdf5 {

// Reader for ref datasets containing uint32 source/target pairs.
struct RawRefPairReader {
    hid_t dset = H5I_INVALID_HID;
    hid_t filespace = H5I_INVALID_HID;
    hid_t pair_array_type = H5I_INVALID_HID;
    bool is_2d = false;
    size_t row_count = 0;

    RawRefPairReader(hid_t file_id, const char* dataset_path);
    ~RawRefPairReader();

    bool read_rows(size_t first_idx, size_t count, std::vector<RefPair>& out) const;
};

}  // namespace ndlar::hdf5
