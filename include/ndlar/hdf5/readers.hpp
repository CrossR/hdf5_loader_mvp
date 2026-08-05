#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <hdf5.h>

// Umbrella header for all HDF5 readers.
#include "ndlar/hdf5/readers/interaction.hpp"
#include "ndlar/hdf5/readers/packet_fraction.hpp"
#include "ndlar/hdf5/readers/pairs.hpp"
#include "ndlar/hdf5/readers/prompt_hit.hpp"
#include "ndlar/hdf5/readers/region.hpp"
#include "ndlar/hdf5/readers/segment.hpp"
#include "ndlar/hdf5/readers/trajectory.hpp"

namespace ndlar::hdf5
{

// Group sorted indices into [start, length] spans, allowing optional small gaps.
std::vector<std::array<size_t, 2>> contiguous_spans(const std::vector<size_t> &indices, size_t max_gap = 0);

// Build event_id -> row_indices map from trajectory rows.
std::unordered_map<int64_t, std::vector<size_t>> build_event_index_from_rows(const RawTrajectoryReader &reader);

// Build event_id -> row_indices map from interaction rows.
std::unordered_map<int64_t, std::vector<size_t>> build_event_index_from_rows(const RawInteractionReader &reader);

// Build an event index from rows of a dataset, using the provided Reader type.
template <typename Reader>
std::unordered_map<int64_t, std::vector<size_t>> build_event_index_from_rows_impl(const Reader &reader);

} // namespace ndlar::hdf5
