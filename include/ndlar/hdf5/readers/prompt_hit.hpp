#pragma once

#include <cstddef>
#include <vector>

#include <hdf5.h>

#include "ndlar/hdf5/types.hpp"

namespace ndlar::hdf5 {

// Reader for prompt hits.
struct RawPromptHitReader {
    hid_t dset = H5I_INVALID_HID;
    hid_t mem_type = H5I_INVALID_HID;
    hid_t filespace = H5I_INVALID_HID;
    size_t row_count = 0;

    explicit RawPromptHitReader(hid_t file_id);
    ~RawPromptHitReader();

    bool read_rows(size_t first_idx, size_t count, std::vector<PromptHit>& out) const;
};

}  // namespace ndlar::hdf5
