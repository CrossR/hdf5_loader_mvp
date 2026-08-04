#pragma once

#include <cstddef>
#include <memory>
#include <unordered_map>
#include <vector>

#include <hdf5.h>

#include "ndlar/hdf5/types.hpp"

namespace ndlar::hdf5 {

/**
 * A reader for fetching raw packet fraction data from an HDF5 file.
 */
struct RawPacketFractionReader {
    hid_t dset = H5I_INVALID_HID;
    hid_t seg_array = H5I_INVALID_HID;
    hid_t frac_array = H5I_INVALID_HID;
    hid_t mem_type = H5I_INVALID_HID;
    size_t row_count = 0;

    explicit RawPacketFractionReader(hid_t file_id);
    ~RawPacketFractionReader();

    bool read_rows(size_t first_idx, size_t count, std::vector<PacketFraction>& out) const;
};

/**
 * A reader for fetching raw reference region data from an HDF5 file.
 */
struct RawRefRegionReader {
    hid_t dset = H5I_INVALID_HID;
    hid_t mem_type = H5I_INVALID_HID;
    size_t row_count = 0;

    RawRefRegionReader(hid_t file_id, const char* dataset_path);
    ~RawRefRegionReader();

    bool read_rows(size_t first_idx, size_t count, std::vector<RefRegion>& out) const;
};

/**
 * A reader for fetching raw reference pair data from an HDF5 file.
 */
struct RawRefPairReader {
    hid_t dset = H5I_INVALID_HID;
    hid_t pair_array_type = H5I_INVALID_HID;
    bool is_2d = false;
    size_t row_count = 0;

    RawRefPairReader(hid_t file_id, const char* dataset_path);
    ~RawRefPairReader();

    bool read_rows(size_t first_idx, size_t count, std::vector<RefPair>& out) const;
};

/**
 * A reader for fetching raw true segment data from an HDF5 file.
 */
struct RawTrueSegmentReader {
    hid_t dset = H5I_INVALID_HID;
    hid_t mem_type = H5I_INVALID_HID;
    size_t row_count = 0;

    RawTrueSegmentReader(hid_t file_id, const char* dataset_path);
    ~RawTrueSegmentReader();

    bool read_rows(size_t first_idx, size_t count, std::vector<TrueSegment>& out) const;
};

/**
 * A reader for fetching raw trajectory data from an HDF5 file.
 */
struct RawTrajectoryReader {
    hid_t dset = H5I_INVALID_HID;
    hid_t vec3_type = H5I_INVALID_HID;
    hid_t mem_type = H5I_INVALID_HID;
    hid_t event_id_mem_type = H5I_INVALID_HID;
    size_t row_count = 0;

    explicit RawTrajectoryReader(hid_t file_id);
    ~RawTrajectoryReader();

    bool read_rows(size_t first_idx, size_t count, std::vector<Trajectory>& out) const;
    bool read_event_ids(size_t first_idx, size_t count, std::vector<int64_t>& out) const;
};

/**
 * A reader for fetching raw interaction data from an HDF5 file.
 */
struct RawInteractionReader {
    hid_t dset = H5I_INVALID_HID;
    hid_t vec4_type = H5I_INVALID_HID;
    hid_t mem_type = H5I_INVALID_HID;
    hid_t event_id_mem_type = H5I_INVALID_HID;
    size_t row_count = 0;

    explicit RawInteractionReader(hid_t file_id);
    ~RawInteractionReader();

    bool read_rows(size_t first_idx, size_t count, std::vector<Interaction>& out) const;
    bool read_event_ids(size_t first_idx, size_t count, std::vector<int64_t>& out) const;
};

/**
 * Groups consecutive indices into contiguous spans.
 *
 * @param indices The vector of indices to group.
 * @param max_gap The maximum gap between indices to consider them contiguous.
 * @return A vector of arrays representing the contiguous spans.
 */
std::vector<std::array<size_t, 2>> contiguous_spans(const std::vector<size_t>& indices, size_t max_gap = 0);

/**
 * Builds an event index from the rows of a trajectory reader.
 *
 * @param reader The trajectory reader to build the index from.
 * @return A map from event IDs to vectors of row indices.
 */
std::unordered_map<int64_t, std::vector<size_t>> build_event_index_from_rows(const RawTrajectoryReader& reader);

/**
 * Builds an event index from the rows of an interaction reader.
 *
 * @param reader The interaction reader to build the index from.
 * @return A map from event IDs to vectors of row indices.
 */
std::unordered_map<int64_t, std::vector<size_t>> build_event_index_from_rows(const RawInteractionReader& reader);

}  // namespace ndlar::hdf5
