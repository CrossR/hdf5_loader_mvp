
#include "ndlar/hdf5/readers.hpp"
#include "ndlar/hdf5/readers/segment.hpp"

namespace ndlar::hdf5
{

RawTrueSegmentReader::RawTrueSegmentReader(hid_t file_id, const char *dataset_path)
{
    dset = open_dataset_or_throw(file_id, dataset_path, "true_segment dataset");
    row_count = dataset_row_count_or_throw(dset, dataset_path);
    filespace = get_filespace_or_throw(dset, dataset_path);

    mem_type = H5Tcreate(H5T_COMPOUND, sizeof(TrueSegment));
    H5Tinsert(mem_type, "segment_id", HOFFSET(TrueSegment, segment_id), H5T_NATIVE_UINT);
    H5Tinsert(mem_type, "pdg_id", HOFFSET(TrueSegment, pdg_id), H5T_NATIVE_INT);
    H5Tinsert(mem_type, "file_traj_id", HOFFSET(TrueSegment, file_traj_id), H5T_NATIVE_UINT);
    H5Tinsert(mem_type, "traj_id", HOFFSET(TrueSegment, traj_id), H5T_NATIVE_UINT);
    H5Tinsert(mem_type, "vertex_id", HOFFSET(TrueSegment, vertex_id), H5T_NATIVE_ULLONG);
    H5Tinsert(mem_type, "event_id", HOFFSET(TrueSegment, event_id), H5T_NATIVE_LLONG);
}

RawTrueSegmentReader::~RawTrueSegmentReader()
{
    if (mem_type >= 0)
        H5Tclose(mem_type);
    if (dset >= 0)
        H5Dclose(dset);
    if (filespace >= 0)
        H5Sclose(filespace);
}

bool RawTrueSegmentReader::read_rows(size_t first_idx, size_t count, std::vector<TrueSegment> &out) const
{
    return read_rows_1d_hyperslab(dset, mem_type, filespace, row_count, first_idx, count, out, "true_segment");
}

} // namespace ndlar::hdf5
