#pragma once

#include <cstddef>
#include <vector>

#include <hdf5.h>

#include "ndlar/hdf5/readers/base.hpp"
#include "ndlar/hdf5/types.hpp"

namespace ndlar::hdf5
{

// Reader for mc_truth/segments/data.
class RawTrueSegmentReader : public RawReaderBase
{
public:
    RawTrueSegmentReader(hid_t file_id, const char *dataset_path);
    bool read_rows(size_t first_idx, size_t count, std::vector<TrueSegment> &out) const;
};

} // namespace ndlar::hdf5
