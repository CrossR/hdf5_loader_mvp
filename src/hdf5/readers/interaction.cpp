
#include "ndlar/hdf5/readers/interaction.hpp"
#include "ndlar/hdf5/paths.hpp"

namespace ndlar::hdf5
{

RawInteractionReader::RawInteractionReader(hid_t file_id) :
    RawReaderBase(file_id, paths::dataset::kInteractions, "interactions dataset")
{
    // Dynamically size based on the struct
    constexpr size_t kMomArraySize = sizeof(Interaction::nu_4mom) / sizeof(Interaction::nu_4mom[0]);
    hsize_t vec4_dims[1] = {kMomArraySize};

    UniqueHID vec4_type(H5Tarray_create2(H5T_NATIVE_FLOAT, 1, vec4_dims), H5Tclose);

    mem_type_.reset(H5Tcreate(H5T_COMPOUND, sizeof(Interaction)), H5Tclose);
    H5Tinsert(mem_type_.get(), "event_id", HOFFSET(Interaction, event_id), H5T_NATIVE_LLONG);
    H5Tinsert(mem_type_.get(), "vertex_id", HOFFSET(Interaction, vertex_id), H5T_NATIVE_ULLONG);
    H5Tinsert(mem_type_.get(), "x_vert", HOFFSET(Interaction, x_vert), H5T_NATIVE_FLOAT);
    H5Tinsert(mem_type_.get(), "y_vert", HOFFSET(Interaction, y_vert), H5T_NATIVE_FLOAT);
    H5Tinsert(mem_type_.get(), "z_vert", HOFFSET(Interaction, z_vert), H5T_NATIVE_FLOAT);
    H5Tinsert(mem_type_.get(), "Enu", HOFFSET(Interaction, Enu), H5T_NATIVE_FLOAT);
    H5Tinsert(mem_type_.get(), "nu_pdg", HOFFSET(Interaction, nu_pdg), H5T_NATIVE_INT);
    H5Tinsert(mem_type_.get(), "nu_4mom", HOFFSET(Interaction, nu_4mom), vec4_type.get());
    H5Tinsert(mem_type_.get(), "isCC", HOFFSET(Interaction, isCC), H5T_NATIVE_UCHAR);
    H5Tinsert(mem_type_.get(), "isQES", HOFFSET(Interaction, isQES), H5T_NATIVE_UCHAR);
    H5Tinsert(mem_type_.get(), "isRES", HOFFSET(Interaction, isRES), H5T_NATIVE_UCHAR);
    H5Tinsert(mem_type_.get(), "isDIS", HOFFSET(Interaction, isDIS), H5T_NATIVE_UCHAR);
    H5Tinsert(mem_type_.get(), "isMEC", HOFFSET(Interaction, isMEC), H5T_NATIVE_UCHAR);
    H5Tinsert(mem_type_.get(), "isCOH", HOFFSET(Interaction, isCOH), H5T_NATIVE_UCHAR);

    event_id_mem_type_.reset(H5Tcreate(H5T_COMPOUND, sizeof(int64_t)), H5Tclose);
    H5Tinsert(event_id_mem_type_.get(), "event_id", 0, H5T_NATIVE_LLONG);
}

bool RawInteractionReader::read_rows(size_t first_idx, size_t count, std::vector<Interaction> &out) const
{
    return read_rows_1d_hyperslab(dset_.get(), mem_type_.get(), filespace_.get(), row_count, first_idx, count, out, "interactions");
}

bool RawInteractionReader::read_event_ids(size_t first_idx, size_t count, std::vector<int64_t> &out) const
{
    return read_rows_1d_hyperslab(dset_.get(), event_id_mem_type_.get(), filespace_.get(), row_count, first_idx, count, out, "interaction_ids");
}

} // namespace ndlar::hdf5
