#include <algorithm>
#include <array>
#include <climits>
#include <cstdint>
#include <chrono>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "hep_evd.h"
#include <hdf5.h>
#include <highfive/H5File.hpp>

namespace {

using RefPair = std::array<uint32_t, 2>;

// ---------------------------------------------------------
// HepEVD Geometry Constants
// ---------------------------------------------------------
static const std::vector<float> NDLAR_X_VALUES = {
    -323.7130126953125,
    -276.2869873046875,
    -223.71299743652344,
    -176.28700256347656,
    -123.71299743652344,
    -76.28700256347656,
    -23.71299934387207,
    23.71299934387207,
    76.28700256347656,
    123.71299743652344,
    176.28700256347656,
    223.71299743652344,
    276.2869873046875,
    323.7130126953125,
};
static const std::vector<float> NDLAR_Z_VALUES = {
    465.7558898925781,
    565.7559204101563,
    665.7559204101563,
    765.7559204101563,
    865.7559204101563,
};
static const float NDLAR_Y_CONST = -66.87129974365234;
static const float NDLAR_X_WIDTH = 46.79100036621094;
static const float NDLAR_Y_WIDTH = 299.5989990234375;
static const float NDLAR_Z_WIDTH = 95.6635971069336;
static constexpr int kDebugMatchPrintLimit = 3;
static constexpr int32_t kInvalidTrigger = INT32_MAX;
static constexpr float kMeVToGeV = 1.0e-3f;

using SteadyClock = std::chrono::steady_clock;

double elapsed_ms(const SteadyClock::time_point& start, const SteadyClock::time_point& end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

struct PromptHit {
    uint32_t id;
    float x, y, z, Q, E, ts_pps;
};

struct EventRow {
    int64_t id;
    int64_t ts_start;
    int64_t ts_end;
    int64_t unix_ts;
    int64_t unix_ts_usec;
};

struct ExtTrig {
    int32_t iogroup;
};

struct Trajectory {
    int64_t event_id;
    std::array<float, 3> xyz_start;
    std::array<float, 3> xyz_end;
    uint32_t file_traj_id;
    uint32_t traj_id;
    int32_t pdg_id;
    float E_start;
    std::array<float, 3> pxyz_start;
    uint64_t vertex_id;
    int64_t parent_id;
};

struct Interaction {
    int64_t event_id;
    uint64_t vertex_id;
    float x_vert;
    float y_vert;
    float z_vert;
    float Enu;
    int32_t nu_pdg;
    std::array<float, 4> nu_4mom;
    uint8_t isCC;
    uint8_t isQES;
    uint8_t isRES;
    uint8_t isDIS;
    uint8_t isMEC;
    uint8_t isCOH;
};

struct PacketFraction {
    std::array<int64_t, 20> segment_ids;
    std::array<double, 20> fraction;
};

struct RefRegion {
    int32_t start;
    int32_t stop;
};

struct TrueSegment {
    uint32_t segment_id;
    int32_t pdg_id;
    uint32_t file_traj_id;
    uint32_t traj_id;
    uint64_t vertex_id;
    int64_t event_id;
};

struct EventProducts {
    int32_t trigger_id = kInvalidTrigger;
    int64_t spill_id = -1;

    std::vector<float> hit_x;
    std::vector<float> hit_y;
    std::vector<float> hit_z;
    std::vector<float> hit_charge;
    std::vector<float> hit_E;
    std::vector<float> hit_ts;
    std::vector<uint16_t> hit_matches;

    std::vector<float> hit_packetFrac;
    std::vector<int32_t> hit_pdg;
    std::vector<int32_t> hit_segmentID;
    std::vector<int64_t> hit_particleID;
    std::vector<int64_t> hit_particleIDLocal;
    std::vector<int64_t> hit_vertexID;

    std::vector<float> mcp_startx;
    std::vector<float> mcp_starty;
    std::vector<float> mcp_startz;
    std::vector<float> mcp_endx;
    std::vector<float> mcp_endy;
    std::vector<float> mcp_endz;
    std::vector<int64_t> mcp_id;
    std::vector<int64_t> mcp_idLocal;
    std::vector<int32_t> mcp_pdg;
    std::vector<float> mcp_energy;
    std::vector<float> mcp_px;
    std::vector<float> mcp_py;
    std::vector<float> mcp_pz;
    std::vector<int64_t> mcp_nuid;
    std::vector<int64_t> mcp_mother;

    std::vector<int64_t> nuID;
    std::vector<float> nue;
    std::vector<int32_t> nuPDG;
    std::vector<float> nupx;
    std::vector<float> nupy;
    std::vector<float> nupz;
    std::vector<float> nuvtxx;
    std::vector<float> nuvtxy;
    std::vector<float> nuvtxz;
    std::vector<int32_t> ccnc;
    std::vector<int32_t> mode;
};

bool is_valid_region(const RefRegion& region) {
    return region.start >= 0 && region.stop > region.start;
}

int region_size(const RefRegion& region) {
    if (!is_valid_region(region)) {
        return 0;
    }
    return region.stop - region.start;
}

float resolve_packet_fraction(const PacketFraction& row, uint32_t segment_id) {
    for (size_t i = 0; i < row.segment_ids.size(); ++i) {
        if (row.segment_ids[i] == static_cast<int64_t>(segment_id)) {
            return static_cast<float>(row.fraction[i]);
        }
    }
    return 0.0f;
}

struct RawPacketFractionReader {
    hid_t dset = H5I_INVALID_HID;
    hid_t seg_array = H5I_INVALID_HID;
    hid_t frac_array = H5I_INVALID_HID;
    hid_t mem_type = H5I_INVALID_HID;
    size_t row_count = 0;

    explicit RawPacketFractionReader(hid_t file_id) {
        dset = H5Dopen2(file_id, "mc_truth/packet_fraction/data", H5P_DEFAULT);
        if (dset < 0) {
            throw std::runtime_error("Failed to open packet_fraction dataset");
        }

        hid_t space = H5Dget_space(dset);
        if (space < 0) {
            throw std::runtime_error("Failed to get packet_fraction dataspace");
        }
        const hssize_t nrows = H5Sget_simple_extent_npoints(space);
        H5Sclose(space);
        if (nrows < 0) {
            throw std::runtime_error("Failed to get packet_fraction row count");
        }
        row_count = static_cast<size_t>(nrows);

        hsize_t dims[1] = {20};
        seg_array = H5Tarray_create2(H5T_NATIVE_LLONG, 1, dims);
        frac_array = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, dims);
        mem_type = H5Tcreate(H5T_COMPOUND, sizeof(PacketFraction));
        H5Tinsert(mem_type, "segment_ids", HOFFSET(PacketFraction, segment_ids), seg_array);
        H5Tinsert(mem_type, "fraction", HOFFSET(PacketFraction, fraction), frac_array);
    }

    ~RawPacketFractionReader() {
        if (mem_type >= 0) H5Tclose(mem_type);
        if (frac_array >= 0) H5Tclose(frac_array);
        if (seg_array >= 0) H5Tclose(seg_array);
        if (dset >= 0) H5Dclose(dset);
    }

    bool read_rows(size_t first_idx, size_t count, std::vector<PacketFraction>& out) const {
        if (count == 0 || first_idx >= row_count || first_idx + count > row_count) {
            return false;
        }

        hid_t filespace = H5Dget_space(dset);
        if (filespace < 0) {
            throw std::runtime_error("Failed to get packet_fraction filespace");
        }

        const hsize_t start[1] = {first_idx};
        const hsize_t hcount[1] = {count};
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, nullptr, hcount, nullptr);
        hid_t memspace = H5Screate_simple(1, hcount, nullptr);

        out.resize(count);
        const herr_t status = H5Dread(dset, mem_type, memspace, filespace, H5P_DEFAULT, out.data());
        H5Sclose(memspace);
        H5Sclose(filespace);

        return status >= 0;
    }
};

struct RawTrajectoryReader {
    hid_t dset = H5I_INVALID_HID;
    hid_t vec3_type = H5I_INVALID_HID;
    hid_t mem_type = H5I_INVALID_HID;
    hid_t event_id_mem_type = H5I_INVALID_HID;
    size_t row_count = 0;

    explicit RawTrajectoryReader(hid_t file_id) {
        dset = H5Dopen2(file_id, "mc_truth/trajectories/data", H5P_DEFAULT);
        if (dset < 0) {
            throw std::runtime_error("Failed to open trajectories dataset");
        }

        hid_t space = H5Dget_space(dset);
        if (space < 0) {
            throw std::runtime_error("Failed to get trajectories dataspace");
        }
        const hssize_t nrows = H5Sget_simple_extent_npoints(space);
        H5Sclose(space);
        if (nrows < 0) {
            throw std::runtime_error("Failed to get trajectories row count");
        }
        row_count = static_cast<size_t>(nrows);

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

    ~RawTrajectoryReader() {
        if (event_id_mem_type >= 0) H5Tclose(event_id_mem_type);
        if (mem_type >= 0) H5Tclose(mem_type);
        if (vec3_type >= 0) H5Tclose(vec3_type);
        if (dset >= 0) H5Dclose(dset);
    }

    bool read_rows(size_t first_idx, size_t count, std::vector<Trajectory>& out) const {
        if (count == 0 || first_idx >= row_count || first_idx + count > row_count) {
            return false;
        }

        hid_t filespace = H5Dget_space(dset);
        if (filespace < 0) {
            throw std::runtime_error("Failed to get trajectories filespace");
        }

        const hsize_t start[1] = {first_idx};
        const hsize_t hcount[1] = {count};
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, nullptr, hcount, nullptr);
        hid_t memspace = H5Screate_simple(1, hcount, nullptr);

        out.resize(count);
        const herr_t status = H5Dread(dset, mem_type, memspace, filespace, H5P_DEFAULT, out.data());
        H5Sclose(memspace);
        H5Sclose(filespace);

        return status >= 0;
    }

    bool read_event_ids(size_t first_idx, size_t count, std::vector<int64_t>& out) const {
        if (count == 0 || first_idx >= row_count || first_idx + count > row_count) {
            return false;
        }

        hid_t filespace = H5Dget_space(dset);
        if (filespace < 0) {
            throw std::runtime_error("Failed to get trajectories filespace");
        }

        const hsize_t start[1] = {first_idx};
        const hsize_t hcount[1] = {count};
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, nullptr, hcount, nullptr);
        hid_t memspace = H5Screate_simple(1, hcount, nullptr);

        out.resize(count);
        const herr_t status = H5Dread(dset, event_id_mem_type, memspace, filespace, H5P_DEFAULT, out.data());
        H5Sclose(memspace);
        H5Sclose(filespace);

        return status >= 0;
    }
};

struct RawInteractionReader {
    hid_t dset = H5I_INVALID_HID;
    hid_t vec4_type = H5I_INVALID_HID;
    hid_t mem_type = H5I_INVALID_HID;
    hid_t event_id_mem_type = H5I_INVALID_HID;
    size_t row_count = 0;

    explicit RawInteractionReader(hid_t file_id) {
        dset = H5Dopen2(file_id, "mc_truth/interactions/data", H5P_DEFAULT);
        if (dset < 0) {
            throw std::runtime_error("Failed to open interactions dataset");
        }

        hid_t space = H5Dget_space(dset);
        if (space < 0) {
            throw std::runtime_error("Failed to get interactions dataspace");
        }
        const hssize_t nrows = H5Sget_simple_extent_npoints(space);
        H5Sclose(space);
        if (nrows < 0) {
            throw std::runtime_error("Failed to get interactions row count");
        }
        row_count = static_cast<size_t>(nrows);

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

    ~RawInteractionReader() {
        if (event_id_mem_type >= 0) H5Tclose(event_id_mem_type);
        if (mem_type >= 0) H5Tclose(mem_type);
        if (vec4_type >= 0) H5Tclose(vec4_type);
        if (dset >= 0) H5Dclose(dset);
    }

    bool read_rows(size_t first_idx, size_t count, std::vector<Interaction>& out) const {
        if (count == 0 || first_idx >= row_count || first_idx + count > row_count) {
            return false;
        }

        hid_t filespace = H5Dget_space(dset);
        if (filespace < 0) {
            throw std::runtime_error("Failed to get interactions filespace");
        }

        const hsize_t start[1] = {first_idx};
        const hsize_t hcount[1] = {count};
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, nullptr, hcount, nullptr);
        hid_t memspace = H5Screate_simple(1, hcount, nullptr);

        out.resize(count);
        const herr_t status = H5Dread(dset, mem_type, memspace, filespace, H5P_DEFAULT, out.data());
        H5Sclose(memspace);
        H5Sclose(filespace);

        return status >= 0;
    }

    bool read_event_ids(size_t first_idx, size_t count, std::vector<int64_t>& out) const {
        if (count == 0 || first_idx >= row_count || first_idx + count > row_count) {
            return false;
        }

        hid_t filespace = H5Dget_space(dset);
        if (filespace < 0) {
            throw std::runtime_error("Failed to get interactions filespace");
        }

        const hsize_t start[1] = {first_idx};
        const hsize_t hcount[1] = {count};
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, nullptr, hcount, nullptr);
        hid_t memspace = H5Screate_simple(1, hcount, nullptr);

        out.resize(count);
        const herr_t status = H5Dread(dset, event_id_mem_type, memspace, filespace, H5P_DEFAULT, out.data());
        H5Sclose(memspace);
        H5Sclose(filespace);

        return status >= 0;
    }
};

std::vector<std::array<size_t, 2>> contiguous_spans(const std::vector<size_t>& indices) {
    std::vector<std::array<size_t, 2>> spans;
    if (indices.empty()) {
        return spans;
    }

    size_t span_start = indices[0];
    size_t prev = indices[0];
    for (size_t i = 1; i < indices.size(); ++i) {
        const size_t cur = indices[i];
        if (cur == prev + 1) {
            prev = cur;
            continue;
        }

        spans.push_back({span_start, prev - span_start + 1});
        span_start = cur;
        prev = cur;
    }
    spans.push_back({span_start, prev - span_start + 1});

    return spans;
}

template <typename T>
bool read_rows_h5(const HighFive::DataSet& dset, size_t first_idx, size_t count, size_t nrows, std::vector<T>& out) {
    if (count == 0 || first_idx >= nrows || first_idx + count > nrows) {
        out.clear();
        return false;
    }

    out.clear();
    dset.select({first_idx}, {count}).read(out);
    return out.size() == count;
}

bool read_refpair_rows_h5(const HighFive::DataSet& dset, size_t first_idx, size_t count, size_t nrows, std::vector<RefPair>& out) {
    if (count == 0 || first_idx >= nrows || first_idx + count > nrows) {
        out.clear();
        return false;
    }

    const std::vector<size_t> dims = dset.getDimensions();
    out.clear();
    if (dims.size() == 2 && dims[1] == 2) {
        dset.select({first_idx, 0}, {count, 2}).read(out);
        return out.size() == count;
    }

    return read_rows_h5<RefPair>(dset, first_idx, count, nrows, out);
}

struct StreamingContext {
    std::vector<EventRow> events;
    std::vector<ExtTrig> ext_trigs;
    std::vector<RefRegion> hit_event_bounds;
    std::vector<RefRegion> event_to_exttrig_reg;
    std::vector<RefPair> event_to_exttrig_ref;

    HighFive::DataSet dset_hits;
    HighFive::DataSet dset_hit_to_pkt_reg;
    HighFive::DataSet dset_hit_to_pkt_ref;
    HighFive::DataSet dset_pkt_to_seg_reg;
    HighFive::DataSet dset_pkt_to_seg_ref;
    HighFive::DataSet dset_pkt_to_frac_reg;
    HighFive::DataSet dset_pkt_to_frac_ref;
    HighFive::DataSet dset_segments;

    size_t n_hit_to_pkt_reg = 0;
    size_t n_hit_to_pkt_ref = 0;
    size_t n_pkt_to_seg_reg = 0;
    size_t n_pkt_to_seg_ref = 0;
    size_t n_pkt_to_frac_reg = 0;
    size_t n_pkt_to_frac_ref = 0;
    size_t n_segments = 0;

    std::unordered_map<int64_t, std::vector<size_t>> traj_rows_by_event;
    std::unordered_map<int64_t, std::vector<size_t>> int_rows_by_event;

    std::unordered_map<size_t, TrueSegment> segment_cache;
    std::unordered_map<size_t, PacketFraction> fraction_cache;
};

std::unordered_map<int64_t, std::vector<size_t>> build_event_index_from_rows(const RawTrajectoryReader& reader) {
    std::unordered_map<int64_t, std::vector<size_t>> by_event;
    static constexpr size_t kChunkRows = 65536;
    std::vector<int64_t> ids;
    for (size_t base = 0; base < reader.row_count; base += kChunkRows) {
        const size_t n = std::min(kChunkRows, reader.row_count - base);
        if (!reader.read_event_ids(base, n, ids)) {
            continue;
        }
        for (size_t i = 0; i < ids.size(); ++i) {
            by_event[ids[i]].push_back(base + i);
        }
    }
    return by_event;
}

std::unordered_map<int64_t, std::vector<size_t>> build_event_index_from_rows(const RawInteractionReader& reader) {
    std::unordered_map<int64_t, std::vector<size_t>> by_event;
    static constexpr size_t kChunkRows = 65536;
    std::vector<int64_t> ids;
    for (size_t base = 0; base < reader.row_count; base += kChunkRows) {
        const size_t n = std::min(kChunkRows, reader.row_count - base);
        if (!reader.read_event_ids(base, n, ids)) {
            continue;
        }
        for (size_t i = 0; i < ids.size(); ++i) {
            by_event[ids[i]].push_back(base + i);
        }
    }
    return by_event;
}

int32_t select_trigger_id_stream(const StreamingContext& ctx, size_t event_index) {
    if (event_index >= ctx.event_to_exttrig_reg.size()) {
        return kInvalidTrigger;
    }

    const RefRegion trig_bounds = ctx.event_to_exttrig_reg[event_index];
    if (!is_valid_region(trig_bounds)) {
        return kInvalidTrigger;
    }

    const size_t start = static_cast<size_t>(trig_bounds.start);
    const size_t stop = static_cast<size_t>(trig_bounds.stop);
    if (start >= ctx.event_to_exttrig_ref.size()) {
        return kInvalidTrigger;
    }

    const size_t clamped_stop = std::min(stop, ctx.event_to_exttrig_ref.size());
    int32_t first = kInvalidTrigger;
    bool any = false;
    for (size_t idx = start; idx < clamped_stop; ++idx) {
        const uint32_t ext_idx = ctx.event_to_exttrig_ref[idx][1];
        if (ext_idx >= ctx.ext_trigs.size()) {
            continue;
        }
        const int32_t iogroup = ctx.ext_trigs[ext_idx].iogroup;
        if (first == kInvalidTrigger) {
            first = iogroup;
        }
        any = true;
        if (iogroup == 5) {
            return 5;
        }
    }

    return any ? first : kInvalidTrigger;
}

EventProducts collect_event_products_stream(
    StreamingContext& ctx,
    size_t event_index,
    const RawPacketFractionReader& frac_reader,
    const RawTrajectoryReader& traj_reader,
    const RawInteractionReader& int_reader) {
    EventProducts out;
    out.trigger_id = select_trigger_id_stream(ctx, event_index);

    if (event_index >= ctx.hit_event_bounds.size()) {
        return out;
    }

    const RefRegion event_bounds = ctx.hit_event_bounds[event_index];
    if (!is_valid_region(event_bounds)) {
        return out;
    }

    const size_t start = static_cast<size_t>(event_bounds.start);
    const size_t stop = static_cast<size_t>(event_bounds.stop);
    const size_t hit_count = stop - start;

    std::vector<PromptHit> event_hits;
    ctx.dset_hits.select({start}, {hit_count}).read(event_hits);

    if (event_hits.empty()) {
        return out;
    }

    uint32_t min_hit_id = UINT32_MAX;
    uint32_t max_hit_id = 0;
    for (const PromptHit& hit : event_hits) {
        if (hit.id < ctx.n_hit_to_pkt_reg) {
            min_hit_id = std::min(min_hit_id, hit.id);
            max_hit_id = std::max(max_hit_id, hit.id);
        }
    }

    std::vector<RefRegion> hit_pkt_regs;
    size_t hit_reg_base = 0;
    if (min_hit_id != UINT32_MAX) {
        hit_reg_base = static_cast<size_t>(min_hit_id);
        const size_t n = static_cast<size_t>(max_hit_id - min_hit_id + 1);
        read_rows_h5<RefRegion>(ctx.dset_hit_to_pkt_reg, hit_reg_base, n, ctx.n_hit_to_pkt_reg, hit_pkt_regs);
    }

    int32_t min_hit_ref = INT32_MAX;
    int32_t max_hit_ref_stop = INT32_MIN;
    for (const RefRegion& r : hit_pkt_regs) {
        if (!is_valid_region(r)) {
            continue;
        }
        min_hit_ref = std::min(min_hit_ref, r.start);
        max_hit_ref_stop = std::max(max_hit_ref_stop, r.stop);
    }

    std::vector<RefPair> hit_pkt_refs;
    size_t hit_ref_base = 0;
    if (min_hit_ref != INT32_MAX && max_hit_ref_stop > min_hit_ref) {
        hit_ref_base = static_cast<size_t>(min_hit_ref);
        const size_t n = static_cast<size_t>(max_hit_ref_stop - min_hit_ref);
        read_refpair_rows_h5(ctx.dset_hit_to_pkt_ref, hit_ref_base, n, ctx.n_hit_to_pkt_ref, hit_pkt_refs);
    }

    std::vector<uint32_t> pkt_ids(event_hits.size(), UINT32_MAX);
    uint32_t min_pkt_id = UINT32_MAX;
    uint32_t max_pkt_id = 0;
    for (size_t i = 0; i < event_hits.size(); ++i) {
        const uint32_t hit_id = event_hits[i].id;
        if (hit_id < hit_reg_base || static_cast<size_t>(hit_id - hit_reg_base) >= hit_pkt_regs.size()) {
            continue;
        }

        const RefRegion r = hit_pkt_regs[static_cast<size_t>(hit_id - hit_reg_base)];
        if (!is_valid_region(r) || r.start < static_cast<int32_t>(hit_ref_base)) {
            continue;
        }

        const size_t ref_idx = static_cast<size_t>(r.start - static_cast<int32_t>(hit_ref_base));
        if (ref_idx >= hit_pkt_refs.size()) {
            continue;
        }

        const uint32_t pkt_id = hit_pkt_refs[ref_idx][1];
        if (pkt_id >= ctx.n_pkt_to_seg_reg || pkt_id >= ctx.n_pkt_to_frac_reg) {
            continue;
        }

        pkt_ids[i] = pkt_id;
        min_pkt_id = std::min(min_pkt_id, pkt_id);
        max_pkt_id = std::max(max_pkt_id, pkt_id);
    }

    std::vector<RefRegion> pkt_seg_regs;
    std::vector<RefRegion> pkt_frac_regs;
    size_t pkt_reg_base = 0;
    if (min_pkt_id != UINT32_MAX) {
        pkt_reg_base = static_cast<size_t>(min_pkt_id);
        const size_t n = static_cast<size_t>(max_pkt_id - min_pkt_id + 1);
        read_rows_h5<RefRegion>(ctx.dset_pkt_to_seg_reg, pkt_reg_base, n, ctx.n_pkt_to_seg_reg, pkt_seg_regs);
        read_rows_h5<RefRegion>(ctx.dset_pkt_to_frac_reg, pkt_reg_base, n, ctx.n_pkt_to_frac_reg, pkt_frac_regs);
    }

    int32_t min_seg_ref = INT32_MAX;
    int32_t max_seg_ref_stop = INT32_MIN;
    int32_t min_frac_ref = INT32_MAX;
    int32_t max_frac_ref_stop = INT32_MIN;
    for (const RefRegion& r : pkt_seg_regs) {
        if (!is_valid_region(r)) continue;
        min_seg_ref = std::min(min_seg_ref, r.start);
        max_seg_ref_stop = std::max(max_seg_ref_stop, r.stop);
    }
    for (const RefRegion& r : pkt_frac_regs) {
        if (!is_valid_region(r)) continue;
        min_frac_ref = std::min(min_frac_ref, r.start);
        max_frac_ref_stop = std::max(max_frac_ref_stop, r.stop);
    }

    std::vector<RefPair> pkt_seg_refs;
    std::vector<RefPair> pkt_frac_refs;
    size_t seg_ref_base = 0;
    size_t frac_ref_base = 0;
    if (min_seg_ref != INT32_MAX && max_seg_ref_stop > min_seg_ref) {
        seg_ref_base = static_cast<size_t>(min_seg_ref);
        const size_t n = static_cast<size_t>(max_seg_ref_stop - min_seg_ref);
        read_refpair_rows_h5(ctx.dset_pkt_to_seg_ref, seg_ref_base, n, ctx.n_pkt_to_seg_ref, pkt_seg_refs);
    }
    if (min_frac_ref != INT32_MAX && max_frac_ref_stop > min_frac_ref) {
        frac_ref_base = static_cast<size_t>(min_frac_ref);
        const size_t n = static_cast<size_t>(max_frac_ref_stop - min_frac_ref);
        read_refpair_rows_h5(ctx.dset_pkt_to_frac_ref, frac_ref_base, n, ctx.n_pkt_to_frac_ref, pkt_frac_refs);
    }

    std::vector<uint32_t> seg_ids(event_hits.size(), UINT32_MAX);
    std::vector<uint32_t> frac_ids(event_hits.size(), UINT32_MAX);
    std::vector<uint16_t> match_counts(event_hits.size(), 0);
    std::vector<size_t> needed_seg_ids;
    std::vector<size_t> needed_frac_ids;
    needed_seg_ids.reserve(event_hits.size());
    needed_frac_ids.reserve(event_hits.size());

    for (size_t i = 0; i < event_hits.size(); ++i) {
        const uint32_t pkt_id = pkt_ids[i];
        if (pkt_id == UINT32_MAX || pkt_id < pkt_reg_base) {
            continue;
        }

        const size_t pkt_off = static_cast<size_t>(pkt_id - pkt_reg_base);
        if (pkt_off >= pkt_seg_regs.size() || pkt_off >= pkt_frac_regs.size()) {
            continue;
        }

        const RefRegion seg_r = pkt_seg_regs[pkt_off];
        if (is_valid_region(seg_r)) {
            match_counts[i] = static_cast<uint16_t>(region_size(seg_r));
            if (seg_r.start >= static_cast<int32_t>(seg_ref_base)) {
                const size_t ref_idx = static_cast<size_t>(seg_r.start - static_cast<int32_t>(seg_ref_base));
                if (ref_idx < pkt_seg_refs.size()) {
                    seg_ids[i] = pkt_seg_refs[ref_idx][1];
                    needed_seg_ids.push_back(static_cast<size_t>(seg_ids[i]));
                }
            }
        }

        const RefRegion frac_r = pkt_frac_regs[pkt_off];
        if (is_valid_region(frac_r) && frac_r.start >= static_cast<int32_t>(frac_ref_base)) {
            const size_t ref_idx = static_cast<size_t>(frac_r.start - static_cast<int32_t>(frac_ref_base));
            if (ref_idx < pkt_frac_refs.size()) {
                frac_ids[i] = pkt_frac_refs[ref_idx][1];
                needed_frac_ids.push_back(static_cast<size_t>(frac_ids[i]));
            }
        }
    }

    std::sort(needed_seg_ids.begin(), needed_seg_ids.end());
    needed_seg_ids.erase(std::unique(needed_seg_ids.begin(), needed_seg_ids.end()), needed_seg_ids.end());
    std::sort(needed_frac_ids.begin(), needed_frac_ids.end());
    needed_frac_ids.erase(std::unique(needed_frac_ids.begin(), needed_frac_ids.end()), needed_frac_ids.end());

    const auto seg_spans = contiguous_spans(needed_seg_ids);
    for (const auto& span : seg_spans) {
        if (span[0] + span[1] > ctx.n_segments) {
            continue;
        }

        bool all_cached = true;
        for (size_t id = span[0]; id < span[0] + span[1]; ++id) {
            if (ctx.segment_cache.find(id) == ctx.segment_cache.end()) {
                all_cached = false;
                break;
            }
        }
        if (all_cached) {
            continue;
        }

        std::vector<TrueSegment> rows;
        if (!read_rows_h5<TrueSegment>(ctx.dset_segments, span[0], span[1], ctx.n_segments, rows)) {
            continue;
        }
        for (size_t i = 0; i < rows.size(); ++i) {
            ctx.segment_cache.emplace(span[0] + i, rows[i]);
        }
    }

    const auto frac_spans = contiguous_spans(needed_frac_ids);
    for (const auto& span : frac_spans) {
        if (span[0] + span[1] > frac_reader.row_count) {
            continue;
        }

        bool all_cached = true;
        for (size_t id = span[0]; id < span[0] + span[1]; ++id) {
            if (ctx.fraction_cache.find(id) == ctx.fraction_cache.end()) {
                all_cached = false;
                break;
            }
        }
        if (all_cached) {
            continue;
        }

        std::vector<PacketFraction> rows;
        if (!frac_reader.read_rows(span[0], span[1], rows)) {
            continue;
        }
        for (size_t i = 0; i < rows.size(); ++i) {
            ctx.fraction_cache.emplace(span[0] + i, rows[i]);
        }
    }

    for (size_t i = 0; i < event_hits.size(); ++i) {
        const PromptHit& hit = event_hits[i];
        out.hit_x.push_back(hit.x);
        out.hit_y.push_back(hit.y);
        out.hit_z.push_back(hit.z);
        out.hit_charge.push_back(hit.Q);
        out.hit_E.push_back(hit.E);
        out.hit_ts.push_back(hit.ts_pps);

        uint16_t matches = 0;
        float packet_fraction = 0.0f;
        int32_t pdg = 0;
        int32_t segment_id = 0;
        int64_t file_traj_id = 0;
        int64_t traj_id = 0;
        int64_t vertex_id = 0;

        matches = match_counts[i];
        const uint32_t seg_id = seg_ids[i];
        if (seg_id != UINT32_MAX) {
            const auto seg_it = ctx.segment_cache.find(static_cast<size_t>(seg_id));
            if (seg_it != ctx.segment_cache.end()) {
                const TrueSegment& true_seg = seg_it->second;
                pdg = true_seg.pdg_id;
                segment_id = static_cast<int32_t>(true_seg.segment_id);
                file_traj_id = static_cast<int64_t>(true_seg.file_traj_id);
                traj_id = static_cast<int64_t>(true_seg.traj_id);
                vertex_id = static_cast<int64_t>(true_seg.vertex_id);

                if (out.spill_id < 0) {
                    out.spill_id = true_seg.event_id;
                }

                const uint32_t frac_id = frac_ids[i];
                if (frac_id != UINT32_MAX) {
                    const auto frac_it = ctx.fraction_cache.find(static_cast<size_t>(frac_id));
                    if (frac_it != ctx.fraction_cache.end()) {
                        packet_fraction = resolve_packet_fraction(frac_it->second, true_seg.segment_id);
                    }
                }
            }
        }

        out.hit_matches.push_back(matches);
        out.hit_packetFrac.push_back(packet_fraction);
        out.hit_pdg.push_back(pdg);
        out.hit_segmentID.push_back(segment_id);
        out.hit_particleID.push_back(file_traj_id);
        out.hit_particleIDLocal.push_back(traj_id);
        out.hit_vertexID.push_back(vertex_id);
    }

    if (out.spill_id >= 0) {
        const auto traj_it = ctx.traj_rows_by_event.find(out.spill_id);
        if (traj_it != ctx.traj_rows_by_event.end()) {
            const auto spans = contiguous_spans(traj_it->second);
            for (const auto& span : spans) {
                std::vector<Trajectory> rows;
                if (!traj_reader.read_rows(span[0], span[1], rows)) {
                    continue;
                }

                for (const Trajectory& t : rows) {
                    out.mcp_startx.push_back(t.xyz_start[0]);
                    out.mcp_starty.push_back(t.xyz_start[1]);
                    out.mcp_startz.push_back(t.xyz_start[2]);
                    out.mcp_endx.push_back(t.xyz_end[0]);
                    out.mcp_endy.push_back(t.xyz_end[1]);
                    out.mcp_endz.push_back(t.xyz_end[2]);
                    out.mcp_id.push_back(static_cast<int64_t>(t.file_traj_id));
                    out.mcp_idLocal.push_back(static_cast<int64_t>(t.traj_id));
                    out.mcp_pdg.push_back(t.pdg_id);
                    out.mcp_energy.push_back(t.E_start * kMeVToGeV);
                    out.mcp_px.push_back(t.pxyz_start[0] * kMeVToGeV);
                    out.mcp_py.push_back(t.pxyz_start[1] * kMeVToGeV);
                    out.mcp_pz.push_back(t.pxyz_start[2] * kMeVToGeV);
                    out.mcp_nuid.push_back(static_cast<int64_t>(t.vertex_id));
                    out.mcp_mother.push_back(t.parent_id);
                }
            }
        }

        const auto int_it = ctx.int_rows_by_event.find(out.spill_id);
        if (int_it != ctx.int_rows_by_event.end()) {
            const auto spans = contiguous_spans(int_it->second);
            for (const auto& span : spans) {
                std::vector<Interaction> rows;
                if (!int_reader.read_rows(span[0], span[1], rows)) {
                    continue;
                }

                for (const Interaction& in : rows) {
                    out.nuID.push_back(static_cast<int64_t>(in.vertex_id));
                    out.nue.push_back(in.Enu * kMeVToGeV);
                    out.nuPDG.push_back(in.nu_pdg);
                    out.nupx.push_back(in.nu_4mom[0] * kMeVToGeV);
                    out.nupy.push_back(in.nu_4mom[1] * kMeVToGeV);
                    out.nupz.push_back(in.nu_4mom[2] * kMeVToGeV);
                    out.nuvtxx.push_back(in.x_vert);
                    out.nuvtxy.push_back(in.y_vert);
                    out.nuvtxz.push_back(in.z_vert);
                    out.ccnc.push_back(in.isCC ? 0 : 1);

                    int32_t mode = 1000;
                    if (in.isQES) mode = 0;
                    if (in.isRES) mode = 1;
                    if (in.isDIS) mode = 2;
                    if (in.isCOH) mode = 3;
                    if (in.isCOH && in.isQES) mode = 4;
                    if (in.isMEC) mode = 10;
                    out.mode.push_back(mode);
                }
            }
        }
    }

    return out;
}

HepEVD::DetectorGeometry get_ndlar_geometry() {
    HepEVD::Volumes volumes;
    for (const auto& x : NDLAR_X_VALUES) {
        for (const auto& z : NDLAR_Z_VALUES) {
            HepEVD::BoxVolume larTpc({x, NDLAR_Y_CONST, z}, NDLAR_X_WIDTH, NDLAR_Y_WIDTH, NDLAR_Z_WIDTH);
            volumes.push_back(larTpc);
        }
    }
    return HepEVD::DetectorGeometry(volumes);
}

HighFive::CompoundType create_compound_PromptHit() {
    return {{"id", HighFive::AtomicType<uint32_t>{}}, {"x", HighFive::AtomicType<float>{}}, {"y", HighFive::AtomicType<float>{}},
            {"z", HighFive::AtomicType<float>{}}, {"Q", HighFive::AtomicType<float>{}}, {"E", HighFive::AtomicType<float>{}},
            {"ts_pps", HighFive::AtomicType<float>{}}};
}

HighFive::CompoundType create_compound_EventRow() {
    return {{"id", HighFive::AtomicType<int64_t>{}},          {"ts_start", HighFive::AtomicType<int64_t>{}},
            {"ts_end", HighFive::AtomicType<int64_t>{}},      {"unix_ts", HighFive::AtomicType<int64_t>{}},
            {"unix_ts_usec", HighFive::AtomicType<int64_t>{}}};
}

HighFive::CompoundType create_compound_ExtTrig() {
    return {{"iogroup", HighFive::AtomicType<int32_t>{}}};
}

HighFive::CompoundType create_compound_RefRegion() {
    return {{"start", HighFive::AtomicType<int32_t>{}}, {"stop", HighFive::AtomicType<int32_t>{}}};
}

HighFive::CompoundType create_compound_TrueSegment() {
    return {{"segment_id", HighFive::AtomicType<uint32_t>{}}, {"pdg_id", HighFive::AtomicType<int32_t>{}},
            {"file_traj_id", HighFive::AtomicType<uint32_t>{}}, {"traj_id", HighFive::AtomicType<uint32_t>{}},
            {"vertex_id", HighFive::AtomicType<uint64_t>{}}, {"event_id", HighFive::AtomicType<int64_t>{}}};
}

}  // namespace

HIGHFIVE_REGISTER_TYPE(PromptHit, create_compound_PromptHit)
HIGHFIVE_REGISTER_TYPE(EventRow, create_compound_EventRow)
HIGHFIVE_REGISTER_TYPE(ExtTrig, create_compound_ExtTrig)
HIGHFIVE_REGISTER_TYPE(RefRegion, create_compound_RefRegion)
HIGHFIVE_REGISTER_TYPE(TrueSegment, create_compound_TrueSegment)

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <h5flow_file.hdf5>\n";
        return 1;
    }

    const auto geometry = get_ndlar_geometry();
    auto server = std::make_unique<HepEVD::HepEVDServer>(geometry);

    try {
        const auto t0_total = SteadyClock::now();
        HighFive::File file(argv[1], HighFive::File::ReadOnly);

        const auto t0_meta = SteadyClock::now();
        StreamingContext ctx;
        ctx.dset_hits = file.getDataSet("charge/calib_prompt_hits/data");
        ctx.dset_hit_to_pkt_reg = file.getDataSet("charge/calib_prompt_hits/ref/charge/packets/ref_region");
        ctx.dset_hit_to_pkt_ref = file.getDataSet("charge/calib_prompt_hits/ref/charge/packets/ref");
        ctx.dset_pkt_to_seg_reg = file.getDataSet("charge/packets/ref/mc_truth/segments/ref_region");
        ctx.dset_pkt_to_seg_ref = file.getDataSet("charge/packets/ref/mc_truth/segments/ref");
        ctx.dset_pkt_to_frac_reg = file.getDataSet("charge/packets/ref/mc_truth/packet_fraction/ref_region");
        ctx.dset_pkt_to_frac_ref = file.getDataSet("charge/packets/ref/mc_truth/packet_fraction/ref");
        ctx.dset_segments = file.getDataSet("mc_truth/segments/data");

        file.getDataSet("charge/events/data").read(ctx.events);
        file.getDataSet("charge/ext_trigs/data").read(ctx.ext_trigs);
        file.getDataSet("charge/events/ref/charge/calib_prompt_hits/ref_region").read(ctx.hit_event_bounds);
        file.getDataSet("charge/events/ref/charge/ext_trigs/ref_region").read(ctx.event_to_exttrig_reg);
        file.getDataSet("charge/events/ref/charge/ext_trigs/ref").read(ctx.event_to_exttrig_ref);

        ctx.n_hit_to_pkt_reg = ctx.dset_hit_to_pkt_reg.getElementCount();
        ctx.n_hit_to_pkt_ref = ctx.dset_hit_to_pkt_ref.getElementCount();
        ctx.n_pkt_to_seg_reg = ctx.dset_pkt_to_seg_reg.getElementCount();
        ctx.n_pkt_to_seg_ref = ctx.dset_pkt_to_seg_ref.getElementCount();
        ctx.n_pkt_to_frac_reg = ctx.dset_pkt_to_frac_reg.getElementCount();
        ctx.n_pkt_to_frac_ref = ctx.dset_pkt_to_frac_ref.getElementCount();
        ctx.n_segments = ctx.dset_segments.getElementCount();
        const auto t1_meta = SteadyClock::now();

        const auto t0_index = SteadyClock::now();
        RawPacketFractionReader frac_reader(file.getId());
        RawTrajectoryReader traj_reader(file.getId());
        RawInteractionReader int_reader(file.getId());

        ctx.traj_rows_by_event = build_event_index_from_rows(traj_reader);
        ctx.int_rows_by_event = build_event_index_from_rows(int_reader);
        const auto t1_index = SteadyClock::now();

        const size_t num_events = ctx.events.size();
        std::cout << "File contains " << num_events << " events.\n";
        std::cout << "-------------------------------------------\n";
        std::cout << "Timing: metadata_load_ms=" << elapsed_ms(t0_meta, t1_meta)
                  << ", index_build_ms=" << elapsed_ms(t0_index, t1_index) << "\n";

        double total_collect_ms = 0.0;
        double total_add_hits_ms = 0.0;
        double total_add_mc_hits_ms = 0.0;
        double total_server_cycle_ms = 0.0;

        for (size_t i = 0; i < num_events; ++i) {
            const auto t0_collect = SteadyClock::now();
            EventProducts ev = collect_event_products_stream(ctx, i, frac_reader, traj_reader, int_reader);
            const auto t1_collect = SteadyClock::now();

            HepEVD::Hits evd_hits;
            HepEVD::MCHits evd_mc_hits;
            for (size_t j = 0; j < ev.hit_x.size(); ++j) {
                evd_hits.push_back(new HepEVD::Hit({ev.hit_x[j], ev.hit_y[j], ev.hit_z[j]}, ev.hit_E[j]));
                evd_mc_hits.push_back(new HepEVD::MCHit({ev.hit_x[j], ev.hit_y[j], ev.hit_z[j]}, ev.hit_pdg[j], ev.hit_E[j]));
            }

            int debug_printed = 0;
            for (size_t j = 0; j < ev.hit_pdg.size() && debug_printed < kDebugMatchPrintLimit; ++j) {
                if (ev.hit_segmentID[j] != 0) {
                    std::cout << "  Hit " << j << " -> PDG: " << ev.hit_pdg[j] << ", frac: " << ev.hit_packetFrac[j] << "\n";
                    ++debug_printed;
                }
            }

            std::cout << "Event " << i
                      << ": trigger=" << ev.trigger_id
                      << ", hits=" << ev.hit_x.size()
                      << ", matched_values=" << ev.hit_pdg.size()
                      << ", spill=" << ev.spill_id
                      << ", traj=" << ev.mcp_id.size()
                      << ", nu_vtx=" << ev.nuID.size() << "\n";

            const auto t0_add_hits = SteadyClock::now();
            server->addHits(evd_hits);
            const auto t1_add_hits = SteadyClock::now();

            const auto t0_add_mc_hits = SteadyClock::now();
            server->addMCHits(evd_mc_hits);
            const auto t1_add_mc_hits = SteadyClock::now();

            const auto t0_server_cycle = SteadyClock::now();
            server->startServer();
            server->resetServer();
            const auto t1_server_cycle = SteadyClock::now();

            const double collect_ms = elapsed_ms(t0_collect, t1_collect);
            const double add_hits_ms = elapsed_ms(t0_add_hits, t1_add_hits);
            const double add_mc_hits_ms = elapsed_ms(t0_add_mc_hits, t1_add_mc_hits);
            const double server_cycle_ms = elapsed_ms(t0_server_cycle, t1_server_cycle);

            total_collect_ms += collect_ms;
            total_add_hits_ms += add_hits_ms;
            total_add_mc_hits_ms += add_mc_hits_ms;
            total_server_cycle_ms += server_cycle_ms;

            std::cout << "  timing_ms: collect=" << collect_ms
                      << ", addHits=" << add_hits_ms
                      << ", addMCHits=" << add_mc_hits_ms
                      << ", serverCycle=" << server_cycle_ms << "\n";
        }

        const auto t1_total = SteadyClock::now();
        const double n = num_events > 0 ? static_cast<double>(num_events) : 1.0;
        std::cout << "-------------------------------------------\n";
        std::cout << "Timing summary: total_ms=" << elapsed_ms(t0_total, t1_total)
                  << ", avg_collect_ms=" << (total_collect_ms / n)
                  << ", avg_addHits_ms=" << (total_add_hits_ms / n)
                  << ", avg_addMCHits_ms=" << (total_add_mc_hits_ms / n)
                  << ", avg_serverCycle_ms=" << (total_server_cycle_ms / n) << "\n";
    } catch (const HighFive::Exception& err) {
        std::cerr << "HDF5 Error: " << err.what() << "\n";
        return 1;
    } catch (const std::exception& err) {
        std::cerr << "Runtime error: " << err.what() << "\n";
        return 1;
    }

    return 0;
}
