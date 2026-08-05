#include "ndlar/hdf5/readers/prompt_hit.hpp"
#include "ndlar/hdf5/paths.hpp"

namespace ndlar::hdf5
{

RawPromptHitReader::RawPromptHitReader(hid_t file_id) :
    RawReaderBase(file_id, paths::dataset::kPromptHits, "prompt_hits dataset")
{
    mem_type_.reset(H5Tcreate(H5T_COMPOUND, sizeof(PromptHit)), H5Tclose);
    H5Tinsert(mem_type_.get(), "id", HOFFSET(PromptHit, id), H5T_NATIVE_UINT);
    H5Tinsert(mem_type_.get(), "x", HOFFSET(PromptHit, x), H5T_NATIVE_FLOAT);
    H5Tinsert(mem_type_.get(), "y", HOFFSET(PromptHit, y), H5T_NATIVE_FLOAT);
    H5Tinsert(mem_type_.get(), "z", HOFFSET(PromptHit, z), H5T_NATIVE_FLOAT);
    H5Tinsert(mem_type_.get(), "Q", HOFFSET(PromptHit, Q), H5T_NATIVE_FLOAT);
    H5Tinsert(mem_type_.get(), "E", HOFFSET(PromptHit, E), H5T_NATIVE_FLOAT);
    H5Tinsert(mem_type_.get(), "ts_pps", HOFFSET(PromptHit, ts_pps), H5T_NATIVE_FLOAT);
}

bool RawPromptHitReader::read_rows(size_t first_idx, size_t count, std::vector<PromptHit> &out) const
{
    return read_rows_1d_hyperslab(dset_.get(), mem_type_.get(), filespace_.get(), row_count, first_idx, count, out, "prompt_hits");
}

} // namespace ndlar::hdf5
