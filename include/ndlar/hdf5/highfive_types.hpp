#pragma once

#include <cstddef> // for offsetof
#include <hdf5.h>
#include <highfive/H5File.hpp>
#include <type_traits> // for std::is_same_v

#include "ndlar/hdf5/types.hpp"

namespace ndlar::hdf5
{

struct H5ArrayWrapper : public HighFive::DataType
{
    explicit H5ArrayWrapper(hid_t id) :
        HighFive::DataType(id)
    {
    }
};

// Direct mapping from C++ types to guaranteed HDF5 native constants
template <typename T>
inline hid_t get_native_hdf5_type()
{
    if constexpr (std::is_same_v<T, float>)
        return H5T_NATIVE_FLOAT;
    else if constexpr (std::is_same_v<T, double>)
        return H5T_NATIVE_DOUBLE;
    else if constexpr (std::is_same_v<T, int64_t>)
        return H5T_NATIVE_INT64;
    else if constexpr (std::is_same_v<T, uint64_t>)
        return H5T_NATIVE_UINT64;
    else if constexpr (std::is_same_v<T, int32_t>)
        return H5T_NATIVE_INT32;
    else if constexpr (std::is_same_v<T, uint32_t>)
        return H5T_NATIVE_UINT32;
    else if constexpr (std::is_same_v<T, uint8_t>)
        return H5T_NATIVE_UINT8;
    else if constexpr (std::is_same_v<T, int8_t>)
        return H5T_NATIVE_INT8;
}

template <typename T>
inline HighFive::DataType create_array_type(size_t size)
{
    hsize_t dims[1] = {size};
    hid_t base_type = get_native_hdf5_type<T>();
    hid_t array_id = H5Tarray_create2(base_type, 1, dims);
    return H5ArrayWrapper(array_id);
}

inline HighFive::CompoundType create_compound_EventIdOnly()
{
    return {{"event_id", HighFive::AtomicType<int64_t>{}}};
}

inline HighFive::CompoundType create_compound_PromptHit()
{
    return HighFive::CompoundType(
        {{"id", HighFive::AtomicType<uint32_t>{}, offsetof(PromptHit, id)}, {"x", HighFive::AtomicType<float>{}, offsetof(PromptHit, x)},
            {"y", HighFive::AtomicType<float>{}, offsetof(PromptHit, y)}, {"z", HighFive::AtomicType<float>{}, offsetof(PromptHit, z)},
            {"Q", HighFive::AtomicType<float>{}, offsetof(PromptHit, Q)}, {"E", HighFive::AtomicType<float>{}, offsetof(PromptHit, E)},
            {"ts_pps", HighFive::AtomicType<float>{}, offsetof(PromptHit, ts_pps)}},
        sizeof(PromptHit));
}

inline HighFive::CompoundType create_compound_EventRow()
{
    return HighFive::CompoundType({{"id", HighFive::AtomicType<int64_t>{}, offsetof(EventRow, id)},
                                      {"ts_start", HighFive::AtomicType<int64_t>{}, offsetof(EventRow, ts_start)},
                                      {"ts_end", HighFive::AtomicType<int64_t>{}, offsetof(EventRow, ts_end)},
                                      {"unix_ts", HighFive::AtomicType<int64_t>{}, offsetof(EventRow, unix_ts)},
                                      {"unix_ts_usec", HighFive::AtomicType<int64_t>{}, offsetof(EventRow, unix_ts_usec)}},
        sizeof(EventRow));
}

inline HighFive::CompoundType create_compound_ExtTrig()
{
    return {{"iogroup", HighFive::AtomicType<int32_t>{}}};
}

inline HighFive::CompoundType create_compound_Trajectory()
{
    return HighFive::CompoundType({{"event_id", HighFive::AtomicType<int64_t>{}, offsetof(Trajectory, event_id)},
                                      {"xyz_start", create_array_type<float>(3), offsetof(Trajectory, xyz_start)},
                                      {"xyz_end", create_array_type<float>(3), offsetof(Trajectory, xyz_end)},
                                      {"file_traj_id", HighFive::AtomicType<uint32_t>{}, offsetof(Trajectory, file_traj_id)},
                                      {"traj_id", HighFive::AtomicType<uint32_t>{}, offsetof(Trajectory, traj_id)},
                                      {"pdg_id", HighFive::AtomicType<int32_t>{}, offsetof(Trajectory, pdg_id)},
                                      {"E_start", HighFive::AtomicType<float>{}, offsetof(Trajectory, E_start)},
                                      {"pxyz_start", create_array_type<float>(3), offsetof(Trajectory, pxyz_start)},
                                      {"vertex_id", HighFive::AtomicType<uint64_t>{}, offsetof(Trajectory, vertex_id)},
                                      {"parent_id", HighFive::AtomicType<int64_t>{}, offsetof(Trajectory, parent_id)}},
        sizeof(Trajectory));
}

inline HighFive::CompoundType create_compound_Interaction()
{
    return HighFive::CompoundType({{"event_id", HighFive::AtomicType<int64_t>{}, offsetof(Interaction, event_id)},
                                      {"vertex_id", HighFive::AtomicType<uint64_t>{}, offsetof(Interaction, vertex_id)},
                                      {"x_vert", HighFive::AtomicType<float>{}, offsetof(Interaction, x_vert)},
                                      {"y_vert", HighFive::AtomicType<float>{}, offsetof(Interaction, y_vert)},
                                      {"z_vert", HighFive::AtomicType<float>{}, offsetof(Interaction, z_vert)},
                                      {"Enu", HighFive::AtomicType<float>{}, offsetof(Interaction, Enu)},
                                      {"nu_pdg", HighFive::AtomicType<int32_t>{}, offsetof(Interaction, nu_pdg)},
                                      {"nu_4mom", create_array_type<float>(4), offsetof(Interaction, nu_4mom)},
                                      {"isCC", HighFive::AtomicType<uint8_t>{}, offsetof(Interaction, isCC)},
                                      {"isQES", HighFive::AtomicType<uint8_t>{}, offsetof(Interaction, isQES)},
                                      {"isRES", HighFive::AtomicType<uint8_t>{}, offsetof(Interaction, isRES)},
                                      {"isDIS", HighFive::AtomicType<uint8_t>{}, offsetof(Interaction, isDIS)},
                                      {"isMEC", HighFive::AtomicType<uint8_t>{}, offsetof(Interaction, isMEC)},
                                      {"isCOH", HighFive::AtomicType<uint8_t>{}, offsetof(Interaction, isCOH)}},
        sizeof(Interaction));
}

inline HighFive::CompoundType create_compound_PacketFraction()
{
    return HighFive::CompoundType({{"segment_ids", create_array_type<int64_t>(20), offsetof(PacketFraction, segment_ids)},
                                      {"fraction", create_array_type<double>(20), offsetof(PacketFraction, fraction)}},
        sizeof(PacketFraction));
}

inline HighFive::CompoundType create_compound_RefRegion()
{
    return HighFive::CompoundType({{"start", HighFive::AtomicType<int32_t>{}, offsetof(RefRegion, start)},
                                      {"stop", HighFive::AtomicType<int32_t>{}, offsetof(RefRegion, stop)}},
        sizeof(RefRegion));
}

inline HighFive::CompoundType create_compound_TrueSegment()
{
    return HighFive::CompoundType({{"segment_id", HighFive::AtomicType<uint32_t>{}, offsetof(TrueSegment, segment_id)},
                                      {"pdg_id", HighFive::AtomicType<int32_t>{}, offsetof(TrueSegment, pdg_id)},
                                      {"file_traj_id", HighFive::AtomicType<uint32_t>{}, offsetof(TrueSegment, file_traj_id)},
                                      {"traj_id", HighFive::AtomicType<uint32_t>{}, offsetof(TrueSegment, traj_id)},
                                      {"vertex_id", HighFive::AtomicType<uint64_t>{}, offsetof(TrueSegment, vertex_id)},
                                      {"event_id", HighFive::AtomicType<int64_t>{}, offsetof(TrueSegment, event_id)}},
        sizeof(TrueSegment));
}

} // namespace ndlar::hdf5

// Register the compound types with HighFive for serialization/deserialization.
//
// Must be outside of the namespace, otherwise HighFive will not be able to find
// them.
HIGHFIVE_REGISTER_TYPE(ndlar::hdf5::PromptHit, ndlar::hdf5::create_compound_PromptHit)
HIGHFIVE_REGISTER_TYPE(ndlar::hdf5::EventRow, ndlar::hdf5::create_compound_EventRow)
HIGHFIVE_REGISTER_TYPE(ndlar::hdf5::ExtTrig, ndlar::hdf5::create_compound_ExtTrig)
HIGHFIVE_REGISTER_TYPE(ndlar::hdf5::Trajectory, ndlar::hdf5::create_compound_Trajectory)
HIGHFIVE_REGISTER_TYPE(ndlar::hdf5::Interaction, ndlar::hdf5::create_compound_Interaction)
HIGHFIVE_REGISTER_TYPE(ndlar::hdf5::PacketFraction, ndlar::hdf5::create_compound_PacketFraction)
HIGHFIVE_REGISTER_TYPE(ndlar::hdf5::RefRegion, ndlar::hdf5::create_compound_RefRegion)
HIGHFIVE_REGISTER_TYPE(ndlar::hdf5::TrueSegment, ndlar::hdf5::create_compound_TrueSegment)
HIGHFIVE_REGISTER_TYPE(ndlar::hdf5::EventIdOnly, ndlar::hdf5::create_compound_EventIdOnly)
