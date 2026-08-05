#pragma once

#include <cstddef>
#include <vector>

#include <hdf5.h>

#include "ndlar/hdf5/readers/base.hpp"
#include "ndlar/hdf5/types.hpp"

namespace ndlar::hdf5
{

// Reader for mc_truth/trajectories/data.
class RawTrajectoryReader : public RawReaderBase
{
public:
    RawTrajectoryReader(hid_t file_id);
    bool read_rows(size_t first_idx, size_t count, std::vector<Trajectory> &out) const;
    bool read_event_ids(size_t first_idx, size_t count, std::vector<int64_t> &out) const;

private:
    UniqueHID event_id_mem_type_;
};

} // namespace ndlar::hdf5
