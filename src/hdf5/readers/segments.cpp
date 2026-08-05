#include "ndlar/hdf5/readers/segment.hpp"

namespace ndlar::hdf5
{

RawTrueSegmentReader::RawTrueSegmentReader(hid_t file_id, const char *dataset_path) :
    RawReaderBase(file_id, dataset_path, "true_segment dataset")
{
    mem_type_.reset(H5Tcreate(H5T_COMPOUND, sizeof(TrueSegment)), H5Tclose);
    H5Tinsert(mem_type_.get(), "segment_id", HOFFSET(TrueSegment, segment_id), H5T_NATIVE_UINT);
    H5Tinsert(mem_type_.get(), "pdg_id", HOFFSET(TrueSegment, pdg_id), H5T_NATIVE_INT);
    H5Tinsert(mem_type_.get(), "file_traj_id", HOFFSET(TrueSegment, file_traj_id), H5T_NATIVE_UINT);
    H5Tinsert(mem_type_.get(), "traj_id", HOFFSET(TrueSegment, traj_id), H5T_NATIVE_UINT);
    H5Tinsert(mem_type_.get(), "vertex_id", HOFFSET(TrueSegment, vertex_id), H5T_NATIVE_ULLONG);
    H5Tinsert(mem_type_.get(), "event_id", HOFFSET(TrueSegment, event_id), H5T_NATIVE_LLONG);
}

bool RawTrueSegmentReader::read_rows(size_t first_idx, size_t count, std::vector<TrueSegment> &out) const
{
    return read_rows_1d_hyperslab(dset_.get(), mem_type_.get(), filespace_.get(), row_count, first_idx, count, out, "true_segment");
}

} // namespace ndlar::hdf5
