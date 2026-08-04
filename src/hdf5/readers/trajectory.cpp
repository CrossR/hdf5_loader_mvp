
#include "ndlar/hdf5/paths.hpp"
#include "ndlar/hdf5/readers.hpp"

namespace ndlar::hdf5
{

RawTrajectoryReader::RawTrajectoryReader(hid_t file_id)
{
    dset = open_dataset_or_throw(file_id, paths::dataset::kTrajectories, "trajectories dataset");
    row_count = dataset_row_count_or_throw(dset, "trajectories dataset");
    filespace = get_filespace_or_throw(dset, "trajectories dataset");

    hsize_t vec3_dims[1] = {3};
    vec3_type = H5Tarray_create2(H5T_NATIVE_FLOAT, 1, vec3_dims);
    mem_type = H5Tcreate(H5T_COMPOUND, sizeof(Trajectory));
    H5Tinsert(mem_type, "event_id", HOFFSET(Trajectory, event_id), H5T_NATIVE_LLONG);
    H5Tinsert(mem_type, "xyz_start", HOFFSET(Trajectory, xyz_start), vec3_type);
    H5Tinsert(mem_type, "xyz_end", HOFFSET(Trajectory, xyz_end), vec3_type);
    H5Tinsert(mem_type, "file_traj_id", HOFFSET(Trajectory, file_traj_id), H5T_NATIVE_UINT);
    H5Tinsert(mem_type, "traj_id", HOFFSET(Trajectory, traj_id), H5T_NATIVE_UINT);
    H5Tinsert(mem_type, "pdg_id", HOFFSET(Trajectory, pdg_id), H5T_NATIVE_INT);
    H5Tinsert(mem_type, "E_start", HOFFSET(Trajectory, E_start), H5T_NATIVE_FLOAT);
    H5Tinsert(mem_type, "pxyz_start", HOFFSET(Trajectory, pxyz_start), vec3_type);
    H5Tinsert(mem_type, "vertex_id", HOFFSET(Trajectory, vertex_id), H5T_NATIVE_ULLONG);
    H5Tinsert(mem_type, "parent_id", HOFFSET(Trajectory, parent_id), H5T_NATIVE_LLONG);

    event_id_mem_type = H5Tcreate(H5T_COMPOUND, sizeof(int64_t));
    H5Tinsert(event_id_mem_type, "event_id", 0, H5T_NATIVE_LLONG);
}

RawTrajectoryReader::~RawTrajectoryReader()
{
    if (event_id_mem_type >= 0)
        H5Tclose(event_id_mem_type);
    if (mem_type >= 0)
        H5Tclose(mem_type);
    if (vec3_type >= 0)
        H5Tclose(vec3_type);
    if (dset >= 0)
        H5Dclose(dset);
    if (filespace >= 0)
        H5Sclose(filespace);
}

bool RawTrajectoryReader::read_rows(size_t first_idx, size_t count, std::vector<Trajectory> &out) const
{
    return read_rows_1d_hyperslab(dset, mem_type, filespace, row_count, first_idx, count, out, "trajectories");
}

bool RawTrajectoryReader::read_event_ids(size_t first_idx, size_t count, std::vector<int64_t> &out) const
{
    return read_rows_1d_hyperslab(dset, event_id_mem_type, filespace, row_count, first_idx, count, out, "trajectory_ids");
}

} // namespace ndlar::hdf5
