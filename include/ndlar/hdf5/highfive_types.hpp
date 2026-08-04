#pragma once

#include <highfive/H5File.hpp>

#include "ndlar/hdf5/types.hpp"

namespace ndlar::hdf5 {

/*
 * Creates a compound type for the PromptHit struct.
 *
 * @return A HighFive::CompoundType representing the PromptHit struct.
 */
inline HighFive::CompoundType create_compound_PromptHit() {
    return {{"id", HighFive::AtomicType<uint32_t>{}}, {"x", HighFive::AtomicType<float>{}}, {"y", HighFive::AtomicType<float>{}},
            {"z", HighFive::AtomicType<float>{}}, {"Q", HighFive::AtomicType<float>{}}, {"E", HighFive::AtomicType<float>{}},
            {"ts_pps", HighFive::AtomicType<float>{}}};
}

/*
 * Creates a compound type for the EventRow struct.
 *
 * @return A HighFive::CompoundType representing the EventRow struct.
 */
inline HighFive::CompoundType create_compound_EventRow() {
    return {{"id", HighFive::AtomicType<int64_t>{}},          {"ts_start", HighFive::AtomicType<int64_t>{}},
            {"ts_end", HighFive::AtomicType<int64_t>{}},      {"unix_ts", HighFive::AtomicType<int64_t>{}},
            {"unix_ts_usec", HighFive::AtomicType<int64_t>{}}};
}

/*
 * Creates a compound type for the ExtTrig struct.
 *
 * @return A HighFive::CompoundType representing the ExtTrig struct.
 */
inline HighFive::CompoundType create_compound_ExtTrig() {
    return {{"iogroup", HighFive::AtomicType<int32_t>{}}};
}

/*
 * Creates a compound type for the RefRegion struct.
 *
 * @return A HighFive::CompoundType representing the RefRegion struct.
 */
inline HighFive::CompoundType create_compound_RefRegion() {
    return {{"start", HighFive::AtomicType<int32_t>{}}, {"stop", HighFive::AtomicType<int32_t>{}}};
}

/*
 * Creates a compound type for the TrueSegment struct.
 *
 * @return A HighFive::CompoundType representing the TrueSegment struct.
 */
inline HighFive::CompoundType create_compound_TrueSegment() {
    return {{"segment_id", HighFive::AtomicType<uint32_t>{}}, {"pdg_id", HighFive::AtomicType<int32_t>{}},
            {"file_traj_id", HighFive::AtomicType<uint32_t>{}}, {"traj_id", HighFive::AtomicType<uint32_t>{}},
            {"vertex_id", HighFive::AtomicType<uint64_t>{}}, {"event_id", HighFive::AtomicType<int64_t>{}}};
}

}  // namespace ndlar::hdf5

// INFO: Register the compound types with HighFive for serialization/deserialization.
HIGHFIVE_REGISTER_TYPE(ndlar::hdf5::PromptHit, ndlar::hdf5::create_compound_PromptHit)
HIGHFIVE_REGISTER_TYPE(ndlar::hdf5::EventRow, ndlar::hdf5::create_compound_EventRow)
HIGHFIVE_REGISTER_TYPE(ndlar::hdf5::ExtTrig, ndlar::hdf5::create_compound_ExtTrig)
HIGHFIVE_REGISTER_TYPE(ndlar::hdf5::RefRegion, ndlar::hdf5::create_compound_RefRegion)
HIGHFIVE_REGISTER_TYPE(ndlar::hdf5::TrueSegment, ndlar::hdf5::create_compound_TrueSegment)
