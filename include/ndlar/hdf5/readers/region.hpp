#pragma once

#include <cstddef>
#include <vector>

#include <hdf5.h>

#include "ndlar/hdf5/readers.hpp"
#include "ndlar/hdf5/types.hpp"

namespace ndlar::hdf5
{

// Reader for ref_region datasets containing [start, stop) spans.
class RawRefRegionReader : public RawReaderBase
{
public:
    RawRefRegionReader(hid_t file_id, const char *dataset_path);
    bool read_rows(size_t first_idx, size_t count, std::vector<RefRegion> &out) const;
};

} // namespace ndlar::hdf5
