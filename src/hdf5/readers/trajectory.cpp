#include "ndlar/hdf5/readers/trajectory.hpp"
#include "ndlar/hdf5/paths.hpp"

namespace ndlar::hdf5
{

RawTrajectoryReader::RawTrajectoryReader(hid_t file_id) :
    RawReaderBase(file_id, paths::dataset::kTrajectories, "trajectories dataset")
{
    // Fix magic number
    constexpr size_t kVec3Size = sizeof(Trajectory::xyz_start) / sizeof(Trajectory::xyz_start[0]);
    hsize_t vec3_dims[1] = {kVec3Size};

    UniqueHID vec3_type(H5Tarray_create2(H5T_NATIVE_FLOAT, 1, vec3_dims), H5Tclose);

    mem_type_.reset(H5Tcreate(H5T_COMPOUND, sizeof(Trajectory)), H5Tclose);
    H5Tinsert(mem_type_.get(), "event_id", HOFFSET(Trajectory, event_id), H5T_NATIVE_LLONG);
    H5Tinsert(mem_type_.get(), "xyz_start", HOFFSET(Trajectory, xyz_start), vec3_type.get());
    H5Tinsert(mem_type_.get(), "xyz_end", HOFFSET(Trajectory, xyz_end), vec3_type.get());
    H5Tinsert(mem_type_.get(), "file_traj_id", HOFFSET(Trajectory, file_traj_id), H5T_NATIVE_UINT);
    H5Tinsert(mem_type_.get(), "traj_id", HOFFSET(Trajectory, traj_id), H5T_NATIVE_UINT);
    H5Tinsert(mem_type_.get(), "pdg_id", HOFFSET(Trajectory, pdg_id), H5T_NATIVE_INT);
    H5Tinsert(mem_type_.get(), "E_start", HOFFSET(Trajectory, E_start), H5T_NATIVE_FLOAT);
    H5Tinsert(mem_type_.get(), "pxyz_start", HOFFSET(Trajectory, pxyz_start), vec3_type.get());
    H5Tinsert(mem_type_.get(), "vertex_id", HOFFSET(Trajectory, vertex_id), H5T_NATIVE_ULLONG);
    H5Tinsert(mem_type_.get(), "parent_id", HOFFSET(Trajectory, parent_id), H5T_NATIVE_LLONG);

    event_id_mem_type_.reset(H5Tcreate(H5T_COMPOUND, sizeof(int64_t)), H5Tclose);
    H5Tinsert(event_id_mem_type_.get(), "event_id", 0, H5T_NATIVE_LLONG);
}

bool RawTrajectoryReader::read_rows(size_t first_idx, size_t count, std::vector<Trajectory> &out) const
{
    return read_rows_1d_hyperslab(dset_.get(), mem_type_.get(), filespace_.get(), row_count, first_idx, count, out, "trajectories");
}

bool RawTrajectoryReader::read_event_ids(size_t first_idx, size_t count, std::vector<int64_t> &out) const
{
    return read_rows_1d_hyperslab(dset_.get(), event_id_mem_type_.get(), filespace_.get(), row_count, first_idx, count, out, "trajectory_ids");
}

} // namespace ndlar::hdf5
