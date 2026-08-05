#include "ndlar/hdf5/readers/region.hpp"

namespace ndlar::hdf5
{

RawRefRegionReader::RawRefRegionReader(hid_t file_id, const char *dataset_path) :
    RawReaderBase(file_id, dataset_path, "ref_region dataset")
{
    mem_type_.reset(H5Tcreate(H5T_COMPOUND, sizeof(RefRegion)), H5Tclose);
    H5Tinsert(mem_type_.get(), "start", HOFFSET(RefRegion, start), H5T_NATIVE_INT);
    H5Tinsert(mem_type_.get(), "stop", HOFFSET(RefRegion, stop), H5T_NATIVE_INT);
}

bool RawRefRegionReader::read_rows(size_t first_idx, size_t count, std::vector<RefRegion> &out) const
{
    return read_rows_1d_hyperslab(dset_.get(), mem_type_.get(), filespace_.get(), row_count, first_idx, count, out, "ref_region");
}

} // namespace ndlar::hdf5
