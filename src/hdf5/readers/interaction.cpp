#include "ndlar/hdf5/paths.hpp"

#include "ndlar/hdf5/readers.hpp"

namespace ndlar::hdf5
{

RawInteractionReader::RawInteractionReader(hid_t file_id)
{
    dset = open_dataset_or_throw(file_id, paths::dataset::kInteractions, "interactions dataset");
    row_count = dataset_row_count_or_throw(dset, "interactions dataset");
    filespace = get_filespace_or_throw(dset, "interactions dataset");

    hsize_t vec4_dims[1] = {4};
    vec4_type = H5Tarray_create2(H5T_NATIVE_FLOAT, 1, vec4_dims);
    mem_type = H5Tcreate(H5T_COMPOUND, sizeof(Interaction));
    H5Tinsert(mem_type, "event_id", HOFFSET(Interaction, event_id), H5T_NATIVE_LLONG);
    H5Tinsert(mem_type, "vertex_id", HOFFSET(Interaction, vertex_id), H5T_NATIVE_ULLONG);
    H5Tinsert(mem_type, "x_vert", HOFFSET(Interaction, x_vert), H5T_NATIVE_FLOAT);
    H5Tinsert(mem_type, "y_vert", HOFFSET(Interaction, y_vert), H5T_NATIVE_FLOAT);
    H5Tinsert(mem_type, "z_vert", HOFFSET(Interaction, z_vert), H5T_NATIVE_FLOAT);
    H5Tinsert(mem_type, "Enu", HOFFSET(Interaction, Enu), H5T_NATIVE_FLOAT);
    H5Tinsert(mem_type, "nu_pdg", HOFFSET(Interaction, nu_pdg), H5T_NATIVE_INT);
    H5Tinsert(mem_type, "nu_4mom", HOFFSET(Interaction, nu_4mom), vec4_type);
    H5Tinsert(mem_type, "isCC", HOFFSET(Interaction, isCC), H5T_NATIVE_UCHAR);
    H5Tinsert(mem_type, "isQES", HOFFSET(Interaction, isQES), H5T_NATIVE_UCHAR);
    H5Tinsert(mem_type, "isRES", HOFFSET(Interaction, isRES), H5T_NATIVE_UCHAR);
    H5Tinsert(mem_type, "isDIS", HOFFSET(Interaction, isDIS), H5T_NATIVE_UCHAR);
    H5Tinsert(mem_type, "isMEC", HOFFSET(Interaction, isMEC), H5T_NATIVE_UCHAR);
    H5Tinsert(mem_type, "isCOH", HOFFSET(Interaction, isCOH), H5T_NATIVE_UCHAR);

    event_id_mem_type = H5Tcreate(H5T_COMPOUND, sizeof(int64_t));
    H5Tinsert(event_id_mem_type, "event_id", 0, H5T_NATIVE_LLONG);
}

RawInteractionReader::~RawInteractionReader()
{
    if (event_id_mem_type >= 0)
        H5Tclose(event_id_mem_type);
    if (mem_type >= 0)
        H5Tclose(mem_type);
    if (vec4_type >= 0)
        H5Tclose(vec4_type);
    if (dset >= 0)
        H5Dclose(dset);
    if (filespace >= 0)
        H5Sclose(filespace);
}

bool RawInteractionReader::read_rows(size_t first_idx, size_t count, std::vector<Interaction> &out) const
{
    return read_rows_1d_hyperslab(dset, mem_type, filespace, row_count, first_idx, count, out, "interactions");
}

bool RawInteractionReader::read_event_ids(size_t first_idx, size_t count, std::vector<int64_t> &out) const
{
    return read_rows_1d_hyperslab(dset, event_id_mem_type, filespace, row_count, first_idx, count, out, "interaction_ids");
}

} // namespace ndlar::hdf5
