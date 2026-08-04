#pragma once

#include <highfive/H5File.hpp>

#include "ndlar/hdf5/types.hpp"

namespace ndlar::hdf5 {

inline HighFive::CompoundType create_compound_PromptHit() {
    return {{"id", HighFive::AtomicType<uint32_t>{}}, {"x", HighFive::AtomicType<float>{}}, {"y", HighFive::AtomicType<float>{}},
            {"z", HighFive::AtomicType<float>{}}, {"Q", HighFive::AtomicType<float>{}}, {"E", HighFive::AtomicType<float>{}},
            {"ts_pps", HighFive::AtomicType<float>{}}};
}

inline HighFive::CompoundType create_compound_EventRow() {
    return {{"id", HighFive::AtomicType<int64_t>{}},          {"ts_start", HighFive::AtomicType<int64_t>{}},
            {"ts_end", HighFive::AtomicType<int64_t>{}},      {"unix_ts", HighFive::AtomicType<int64_t>{}},
            {"unix_ts_usec", HighFive::AtomicType<int64_t>{}}};
}

inline HighFive::CompoundType create_compound_ExtTrig() {
    return {{"iogroup", HighFive::AtomicType<int32_t>{}}};
}

inline HighFive::CompoundType create_compound_RefRegion() {
    return {{"start", HighFive::AtomicType<int32_t>{}}, {"stop", HighFive::AtomicType<int32_t>{}}};
}

inline HighFive::CompoundType create_compound_TrueSegment() {
    return {{"segment_id", HighFive::AtomicType<uint32_t>{}}, {"pdg_id", HighFive::AtomicType<int32_t>{}},
            {"file_traj_id", HighFive::AtomicType<uint32_t>{}}, {"traj_id", HighFive::AtomicType<uint32_t>{}},
            {"vertex_id", HighFive::AtomicType<uint64_t>{}}, {"event_id", HighFive::AtomicType<int64_t>{}}};
}

}  // namespace ndlar::hdf5

HIGHFIVE_REGISTER_TYPE(ndlar::hdf5::PromptHit, ndlar::hdf5::create_compound_PromptHit)
HIGHFIVE_REGISTER_TYPE(ndlar::hdf5::EventRow, ndlar::hdf5::create_compound_EventRow)
HIGHFIVE_REGISTER_TYPE(ndlar::hdf5::ExtTrig, ndlar::hdf5::create_compound_ExtTrig)
HIGHFIVE_REGISTER_TYPE(ndlar::hdf5::RefRegion, ndlar::hdf5::create_compound_RefRegion)
HIGHFIVE_REGISTER_TYPE(ndlar::hdf5::TrueSegment, ndlar::hdf5::create_compound_TrueSegment)
