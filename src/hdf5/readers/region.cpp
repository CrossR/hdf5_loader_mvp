
#include "ndlar/hdf5/readers.hpp"
#include "ndlar/hdf5/readers/region.hpp"

namespace ndlar::hdf5
{

RawRefRegionReader::RawRefRegionReader(hid_t file_id, const char *dataset_path)
{
    dset = open_dataset_or_throw(file_id, dataset_path, "ref_region dataset");
    filespace = get_filespace_or_throw(dset, dataset_path);
    row_count = dataset_row_count_or_throw(dset, dataset_path);

    mem_type = H5Tcreate(H5T_COMPOUND, sizeof(RefRegion));
    H5Tinsert(mem_type, "start", HOFFSET(RefRegion, start), H5T_NATIVE_INT);
    H5Tinsert(mem_type, "stop", HOFFSET(RefRegion, stop), H5T_NATIVE_INT);
}

RawRefRegionReader::~RawRefRegionReader()
{
    if (mem_type >= 0)
        H5Tclose(mem_type);
    if (dset >= 0)
        H5Dclose(dset);
    if (filespace >= 0)
        H5Sclose(filespace);
}

bool RawRefRegionReader::read_rows(size_t first_idx, size_t count, std::vector<RefRegion> &out) const
{
    return read_rows_1d_hyperslab(dset, mem_type, filespace, row_count, first_idx, count, out, "ref_region");
}

} // namespace ndlar::hdf5
