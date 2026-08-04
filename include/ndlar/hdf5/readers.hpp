#pragma once

#include <cstddef>
#include <memory>
#include <unordered_map>
#include <vector>

#include <hdf5.h>

#include "ndlar/hdf5/types.hpp"

namespace ndlar::hdf5 {

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

struct RawRefRegionReader {
    hid_t dset = H5I_INVALID_HID;
    hid_t mem_type = H5I_INVALID_HID;
    size_t row_count = 0;

    RawRefRegionReader(hid_t file_id, const char* dataset_path);
    ~RawRefRegionReader();

    bool read_rows(size_t first_idx, size_t count, std::vector<RefRegion>& out) const;
};

struct RawRefPairReader {
    hid_t dset = H5I_INVALID_HID;
    hid_t pair_array_type = H5I_INVALID_HID;
    bool is_2d = false;
    size_t row_count = 0;

    RawRefPairReader(hid_t file_id, const char* dataset_path);
    ~RawRefPairReader();

    bool read_rows(size_t first_idx, size_t count, std::vector<RefPair>& out) const;
};

struct RawTrueSegmentReader {
    hid_t dset = H5I_INVALID_HID;
    hid_t mem_type = H5I_INVALID_HID;
    size_t row_count = 0;

    RawTrueSegmentReader(hid_t file_id, const char* dataset_path);
    ~RawTrueSegmentReader();

    bool read_rows(size_t first_idx, size_t count, std::vector<TrueSegment>& out) const;
};

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

std::vector<std::array<size_t, 2>> contiguous_spans(const std::vector<size_t>& indices, size_t max_gap = 0);

std::unordered_map<int64_t, std::vector<size_t>> build_event_index_from_rows(const RawTrajectoryReader& reader);
std::unordered_map<int64_t, std::vector<size_t>> build_event_index_from_rows(const RawInteractionReader& reader);

}  // namespace ndlar::hdf5
