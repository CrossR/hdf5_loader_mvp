#pragma once

#include <cstddef>
#include <vector>

#include <hdf5.h>

#include "ndlar/hdf5/readers.hpp"
#include "ndlar/hdf5/types.hpp"

namespace ndlar::hdf5
{

// Reader for prompt hits.
class RawPromptHitReader : public RawReaderBase
{
public:
    RawPromptHitReader(hid_t file_id);
    bool read_rows(size_t first_idx, size_t count, std::vector<PromptHit> &out) const;
};

} // namespace ndlar::hdf5
