#pragma once

#include <cstddef>
#include <vector>

#include <hdf5.h>

#include "ndlar/hdf5/types.hpp"

namespace ndlar::hdf5 {

// Reader for mc_truth/trajectories/data.
struct RawTrajectoryReader {
    hid_t dset = H5I_INVALID_HID;
    hid_t vec3_type = H5I_INVALID_HID;
    hid_t mem_type = H5I_INVALID_HID;
    hid_t event_id_mem_type = H5I_INVALID_HID;
    hid_t filespace = H5I_INVALID_HID;
    size_t row_count = 0;

    explicit RawTrajectoryReader(hid_t file_id);
    ~RawTrajectoryReader();

    bool read_rows(size_t first_idx, size_t count, std::vector<Trajectory>& out) const;
    bool read_event_ids(size_t first_idx, size_t count, std::vector<int64_t>& out) const;
};

}  // namespace ndlar::hdf5
