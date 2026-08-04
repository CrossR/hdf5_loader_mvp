
#include "ndlar/hdf5/readers/prompt_hit.hpp"
#include "ndlar/hdf5/paths.hpp"
#include "ndlar/hdf5/readers.hpp"

namespace ndlar::hdf5
{

RawPromptHitReader::RawPromptHitReader(hid_t file_id)
{
    dset = open_dataset_or_throw(file_id, paths::dataset::kPromptHits, "prompt_hits dataset");
    row_count = dataset_row_count_or_throw(dset, "prompt_hits dataset");
    filespace = get_filespace_or_throw(dset, "prompt_hits dataset");

    mem_type = H5Tcreate(H5T_COMPOUND, sizeof(PromptHit));
    H5Tinsert(mem_type, "id", HOFFSET(PromptHit, id), H5T_NATIVE_UINT);
    H5Tinsert(mem_type, "x", HOFFSET(PromptHit, x), H5T_NATIVE_FLOAT);
    H5Tinsert(mem_type, "y", HOFFSET(PromptHit, y), H5T_NATIVE_FLOAT);
    H5Tinsert(mem_type, "z", HOFFSET(PromptHit, z), H5T_NATIVE_FLOAT);
    H5Tinsert(mem_type, "Q", HOFFSET(PromptHit, Q), H5T_NATIVE_FLOAT);
    H5Tinsert(mem_type, "E", HOFFSET(PromptHit, E), H5T_NATIVE_FLOAT);
    H5Tinsert(mem_type, "ts_pps", HOFFSET(PromptHit, ts_pps), H5T_NATIVE_FLOAT);
}

RawPromptHitReader::~RawPromptHitReader()
{
    if (mem_type >= 0)
        H5Tclose(mem_type);
    if (dset >= 0)
        H5Dclose(dset);
    if (filespace >= 0)
        H5Sclose(filespace);
}

bool RawPromptHitReader::read_rows(size_t first_idx, size_t count, std::vector<PromptHit> &out) const
{
    return read_rows_1d_hyperslab(dset, mem_type, filespace, row_count, first_idx, count, out, "prompt_hits");
}

} // namespace ndlar::hdf5
