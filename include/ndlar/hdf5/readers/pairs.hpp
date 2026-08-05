#pragma once

#include <cstddef>
#include <vector>

#include <hdf5.h>

#include "ndlar/hdf5/readers/base.hpp"
#include "ndlar/hdf5/types.hpp"

namespace ndlar::hdf5
{

class RawRefPairReader : public RawReaderBase
{
public:
    RawRefPairReader(hid_t file_id, const char *dataset_path);
    bool read_rows(size_t first_idx, size_t count, std::vector<RefPair> &out) const;

private:
    bool is_2d = false;
    UniqueHID pair_array_type_;
};

} // namespace ndlar::hdf5
