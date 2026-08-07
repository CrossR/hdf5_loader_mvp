#pragma once

#include <hdf5.h>
#include <highfive/H5File.hpp>

#include "ndlar/hdf5/types.hpp"

namespace ndlar::hdf5
{

// Helper to manually build HDF5 array types since HighFive
// lacks the built-in traits for std::array inside compounds.
struct H5ArrayWrapper : public HighFive::DataType {
    explicit H5ArrayWrapper(hid_t id) : HighFive::DataType(id) {}
};

template <typename T>
inline HighFive::DataType create_array_type(size_t size) {
    hsize_t dims[1] = { size };
    hid_t array_id = H5Tarray_create(HighFive::create_datatype<T>().getId(), 1, dims);
    return H5ArrayWrapper(array_id);
}

// These inline functions simply create a mapping between the C++ struct and the
// HDF5 compound type, so that HighFive can automatically serialize/deserialize
// the data for us.
//
// They are registered at the bottom of this file.

inline HighFive::CompoundType create_compound_EventIdOnly()
{
    return {{"event_id", HighFive::AtomicType<int64_t>{}}};
}

inline HighFive::CompoundType create_compound_PromptHit()
{
    return {{"id", HighFive::AtomicType<uint32_t>{}}, {"x", HighFive::AtomicType<float>{}}, {"y", HighFive::AtomicType<float>{}},
        {"z", HighFive::AtomicType<float>{}}, {"Q", HighFive::AtomicType<float>{}}, {"E", HighFive::AtomicType<float>{}},
        {"ts_pps", HighFive::AtomicType<float>{}}};
}

inline HighFive::CompoundType create_compound_EventRow()
{
    return {{"id", HighFive::AtomicType<int64_t>{}}, {"ts_start", HighFive::AtomicType<int64_t>{}}, {"ts_end", HighFive::AtomicType<int64_t>{}},
        {"unix_ts", HighFive::AtomicType<int64_t>{}}, {"unix_ts_usec", HighFive::AtomicType<int64_t>{}}};
}

inline HighFive::CompoundType create_compound_ExtTrig()
{
    return {{"iogroup", HighFive::AtomicType<int32_t>{}}};
}

inline HighFive::CompoundType create_compound_Trajectory()
{
    return {{"event_id", HighFive::AtomicType<int64_t>{}}, {"xyz_start", create_array_type<float>(3)}, {"xyz_end", create_array_type<float>(3)},
        {"file_traj_id", HighFive::AtomicType<uint32_t>{}}, {"traj_id", HighFive::AtomicType<uint32_t>{}},
        {"pdg_id", HighFive::AtomicType<int32_t>{}}, {"E_start", HighFive::AtomicType<float>{}}, {"pxyz_start", create_array_type<float>(3)},
        {"vertex_id", HighFive::AtomicType<uint64_t>{}}, {"parent_id", HighFive::AtomicType<int64_t>{}}};
}

inline HighFive::CompoundType create_compound_Interaction()
{
    return {{"event_id", HighFive::AtomicType<int64_t>{}}, {"vertex_id", HighFive::AtomicType<uint64_t>{}},
        {"x_vert", HighFive::AtomicType<float>{}}, {"y_vert", HighFive::AtomicType<float>{}}, {"z_vert", HighFive::AtomicType<float>{}},
        {"Enu", HighFive::AtomicType<float>{}}, {"nu_pdg", HighFive::AtomicType<int32_t>{}}, {"nu_4mom", create_array_type<float>(4)},
        {"isCC", HighFive::AtomicType<uint8_t>{}}, {"isQES", HighFive::AtomicType<uint8_t>{}}, {"isRES", HighFive::AtomicType<uint8_t>{}},
        {"isDIS", HighFive::AtomicType<uint8_t>{}}, {"isMEC", HighFive::AtomicType<uint8_t>{}}, {"isCOH", HighFive::AtomicType<uint8_t>{}}};
}

inline HighFive::CompoundType create_compound_PacketFraction()
{
    return {{"segment_ids", create_array_type<int64_t>(20)}, {"fraction", create_array_type<double>(20)}};
}

inline HighFive::CompoundType create_compound_RefRegion()
{
    return {{"start", HighFive::AtomicType<int32_t>{}}, {"stop", HighFive::AtomicType<int32_t>{}}};
}

inline HighFive::CompoundType create_compound_TrueSegment()
{
    return {{"segment_id", HighFive::AtomicType<uint32_t>{}}, {"pdg_id", HighFive::AtomicType<int32_t>{}},
        {"file_traj_id", HighFive::AtomicType<uint32_t>{}}, {"traj_id", HighFive::AtomicType<uint32_t>{}},
        {"vertex_id", HighFive::AtomicType<uint64_t>{}}, {"event_id", HighFive::AtomicType<int64_t>{}}};
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
