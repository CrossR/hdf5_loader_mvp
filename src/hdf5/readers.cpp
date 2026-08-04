#include "ndlar/hdf5/readers.hpp"

#include <algorithm>
#include <array>
#include <stdexcept>
#include <string>

#include "ndlar/hdf5/paths.hpp"

namespace ndlar::hdf5 {

namespace {

hid_t open_dataset_or_throw(hid_t file_id, const char* dataset_path, const char* error_context) {
    const hid_t dset = H5Dopen2(file_id, dataset_path, H5P_DEFAULT);
    if (dset < 0) {
        throw std::runtime_error(std::string("Failed to open ") + error_context + ": " + dataset_path);
    }
    return dset;
}

size_t dataset_row_count_or_throw(hid_t dset, const char* error_context) {
    hid_t space = H5Dget_space(dset);
    if (space < 0) {
        throw std::runtime_error(std::string("Failed to get dataspace for ") + error_context);
    }
    const hssize_t nrows = H5Sget_simple_extent_npoints(space);
    H5Sclose(space);
    if (nrows < 0) {
        throw std::runtime_error(std::string("Failed to get row count for ") + error_context);
    }
    return static_cast<size_t>(nrows);
}

template <typename Reader>
std::unordered_map<int64_t, std::vector<size_t>> build_event_index_from_rows_impl(const Reader& reader) {
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

}  // namespace

template <typename T>
bool read_rows_1d_hyperslab(
    hid_t dset,
    hid_t mem_type,
    size_t row_count,
    size_t first_idx,
    size_t count,
    std::vector<T>& out,
    const char* error_context) {
    if (count == 0 || first_idx >= row_count || first_idx + count > row_count) {
        out.clear();
        return false;
    }

    hid_t filespace = H5Dget_space(dset);
    if (filespace < 0) {
        throw std::runtime_error(std::string("Failed to get filespace for ") + error_context);
    }

    const hsize_t start[1] = {first_idx};
    const hsize_t hcount[1] = {count};
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, nullptr, hcount, nullptr);
    hid_t memspace = H5Screate_simple(1, hcount, nullptr);

    out.resize(count);
    const herr_t status = H5Dread(dset, mem_type, memspace, filespace, H5P_DEFAULT, out.data());
    H5Sclose(memspace);
    H5Sclose(filespace);

    if (status < 0) {
        out.clear();
        return false;
    }
    return true;
}

RawPacketFractionReader::RawPacketFractionReader(hid_t file_id) {
    dset = open_dataset_or_throw(file_id, paths::dataset::kPacketFraction, "packet_fraction dataset");
    row_count = dataset_row_count_or_throw(dset, "packet_fraction dataset");

    hsize_t dims[1] = {20};
    seg_array = H5Tarray_create2(H5T_NATIVE_LLONG, 1, dims);
    frac_array = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, dims);
    mem_type = H5Tcreate(H5T_COMPOUND, sizeof(PacketFraction));
    H5Tinsert(mem_type, "segment_ids", HOFFSET(PacketFraction, segment_ids), seg_array);
    H5Tinsert(mem_type, "fraction", HOFFSET(PacketFraction, fraction), frac_array);
}

RawPacketFractionReader::~RawPacketFractionReader() {
    if (mem_type >= 0) H5Tclose(mem_type);
    if (frac_array >= 0) H5Tclose(frac_array);
    if (seg_array >= 0) H5Tclose(seg_array);
    if (dset >= 0) H5Dclose(dset);
}

bool RawPacketFractionReader::read_rows(size_t first_idx, size_t count, std::vector<PacketFraction>& out) const {
    return read_rows_1d_hyperslab(dset, mem_type, row_count, first_idx, count, out, "packet_fraction");
}

RawRefRegionReader::RawRefRegionReader(hid_t file_id, const char* dataset_path) {
    dset = open_dataset_or_throw(file_id, dataset_path, "ref_region dataset");
    row_count = dataset_row_count_or_throw(dset, dataset_path);

    mem_type = H5Tcreate(H5T_COMPOUND, sizeof(RefRegion));
    H5Tinsert(mem_type, "start", HOFFSET(RefRegion, start), H5T_NATIVE_INT);
    H5Tinsert(mem_type, "stop", HOFFSET(RefRegion, stop), H5T_NATIVE_INT);
}

RawRefRegionReader::~RawRefRegionReader() {
    if (mem_type >= 0) H5Tclose(mem_type);
    if (dset >= 0) H5Dclose(dset);
}

bool RawRefRegionReader::read_rows(size_t first_idx, size_t count, std::vector<RefRegion>& out) const {
    return read_rows_1d_hyperslab(dset, mem_type, row_count, first_idx, count, out, "ref_region");
}

RawRefPairReader::RawRefPairReader(hid_t file_id, const char* dataset_path) {
    dset = open_dataset_or_throw(file_id, dataset_path, "ref_pair dataset");

    hid_t space = H5Dget_space(dset);
    if (space < 0) {
        throw std::runtime_error(std::string("Failed to get dataspace: ") + dataset_path);
    }

    const int ndims = H5Sget_simple_extent_ndims(space);
    if (ndims == 2) {
        hsize_t dims[2] = {0, 0};
        H5Sget_simple_extent_dims(space, dims, nullptr);
        is_2d = (dims[1] == 2);
        row_count = static_cast<size_t>(dims[0]);
    } else {
        const hssize_t nrows = H5Sget_simple_extent_npoints(space);
        if (nrows < 0) {
            H5Sclose(space);
            throw std::runtime_error(std::string("Failed to get row count: ") + dataset_path);
        }
        row_count = static_cast<size_t>(nrows);
    }
    H5Sclose(space);

    hsize_t dims[1] = {2};
    pair_array_type = H5Tarray_create2(H5T_NATIVE_UINT, 1, dims);
}

RawRefPairReader::~RawRefPairReader() {
    if (pair_array_type >= 0) H5Tclose(pair_array_type);
    if (dset >= 0) H5Dclose(dset);
}

bool RawRefPairReader::read_rows(size_t first_idx, size_t count, std::vector<RefPair>& out) const {
    if (count == 0 || first_idx >= row_count || first_idx + count > row_count) {
        out.clear();
        return false;
    }

    hid_t filespace = H5Dget_space(dset);
    if (filespace < 0) {
        throw std::runtime_error("Failed to get ref-pair filespace");
    }

    if (is_2d) {
        const hsize_t start[2] = {first_idx, 0};
        const hsize_t hcount[2] = {count, 2};
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, nullptr, hcount, nullptr);
        hid_t memspace = H5Screate_simple(2, hcount, nullptr);

        out.resize(count);
        const herr_t status = H5Dread(dset, H5T_NATIVE_UINT, memspace, filespace, H5P_DEFAULT, out.data());
        H5Sclose(memspace);
        H5Sclose(filespace);

        if (status < 0) {
            out.clear();
            return false;
        }
        return true;
    }

    const hsize_t start[1] = {first_idx};
    const hsize_t hcount[1] = {count};
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, nullptr, hcount, nullptr);
    hid_t memspace = H5Screate_simple(1, hcount, nullptr);

    out.resize(count);
    const herr_t status = H5Dread(dset, pair_array_type, memspace, filespace, H5P_DEFAULT, out.data());
    H5Sclose(memspace);
    H5Sclose(filespace);

    if (status < 0) {
        out.clear();
        return false;
    }
    return true;
}

RawTrueSegmentReader::RawTrueSegmentReader(hid_t file_id, const char* dataset_path) {
    dset = open_dataset_or_throw(file_id, dataset_path, "true_segment dataset");
    row_count = dataset_row_count_or_throw(dset, dataset_path);

    mem_type = H5Tcreate(H5T_COMPOUND, sizeof(TrueSegment));
    H5Tinsert(mem_type, "segment_id", HOFFSET(TrueSegment, segment_id), H5T_NATIVE_UINT);
    H5Tinsert(mem_type, "pdg_id", HOFFSET(TrueSegment, pdg_id), H5T_NATIVE_INT);
    H5Tinsert(mem_type, "file_traj_id", HOFFSET(TrueSegment, file_traj_id), H5T_NATIVE_UINT);
    H5Tinsert(mem_type, "traj_id", HOFFSET(TrueSegment, traj_id), H5T_NATIVE_UINT);
    H5Tinsert(mem_type, "vertex_id", HOFFSET(TrueSegment, vertex_id), H5T_NATIVE_ULLONG);
    H5Tinsert(mem_type, "event_id", HOFFSET(TrueSegment, event_id), H5T_NATIVE_LLONG);
}

RawTrueSegmentReader::~RawTrueSegmentReader() {
    if (mem_type >= 0) H5Tclose(mem_type);
    if (dset >= 0) H5Dclose(dset);
}

bool RawTrueSegmentReader::read_rows(size_t first_idx, size_t count, std::vector<TrueSegment>& out) const {
    return read_rows_1d_hyperslab(dset, mem_type, row_count, first_idx, count, out, "true_segment");
}

RawTrajectoryReader::RawTrajectoryReader(hid_t file_id) {
    dset = open_dataset_or_throw(file_id, paths::dataset::kTrajectories, "trajectories dataset");
    row_count = dataset_row_count_or_throw(dset, "trajectories dataset");

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

RawTrajectoryReader::~RawTrajectoryReader() {
    if (event_id_mem_type >= 0) H5Tclose(event_id_mem_type);
    if (mem_type >= 0) H5Tclose(mem_type);
    if (vec3_type >= 0) H5Tclose(vec3_type);
    if (dset >= 0) H5Dclose(dset);
}

bool RawTrajectoryReader::read_rows(size_t first_idx, size_t count, std::vector<Trajectory>& out) const {
    return read_rows_1d_hyperslab(dset, mem_type, row_count, first_idx, count, out, "trajectories");
}

bool RawTrajectoryReader::read_event_ids(size_t first_idx, size_t count, std::vector<int64_t>& out) const {
    return read_rows_1d_hyperslab(dset, event_id_mem_type, row_count, first_idx, count, out, "trajectory_ids");
}

RawInteractionReader::RawInteractionReader(hid_t file_id) {
    dset = open_dataset_or_throw(file_id, paths::dataset::kInteractions, "interactions dataset");
    row_count = dataset_row_count_or_throw(dset, "interactions dataset");

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

RawInteractionReader::~RawInteractionReader() {
    if (event_id_mem_type >= 0) H5Tclose(event_id_mem_type);
    if (mem_type >= 0) H5Tclose(mem_type);
    if (vec4_type >= 0) H5Tclose(vec4_type);
    if (dset >= 0) H5Dclose(dset);
}

bool RawInteractionReader::read_rows(size_t first_idx, size_t count, std::vector<Interaction>& out) const {
    return read_rows_1d_hyperslab(dset, mem_type, row_count, first_idx, count, out, "interactions");
}

bool RawInteractionReader::read_event_ids(size_t first_idx, size_t count, std::vector<int64_t>& out) const {
    return read_rows_1d_hyperslab(dset, event_id_mem_type, row_count, first_idx, count, out, "interaction_ids");
}

std::vector<std::array<size_t, 2>> contiguous_spans(const std::vector<size_t>& indices, size_t max_gap) {
    std::vector<std::array<size_t, 2>> spans;
    if (indices.empty()) {
        return spans;
    }

    size_t span_start = indices[0];
    size_t prev = indices[0];
    for (size_t i = 1; i < indices.size(); ++i) {
        const size_t cur = indices[i];
        if (cur <= prev + max_gap + 1) {
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

std::unordered_map<int64_t, std::vector<size_t>> build_event_index_from_rows(const RawTrajectoryReader& reader) {
    return build_event_index_from_rows_impl(reader);
}

std::unordered_map<int64_t, std::vector<size_t>> build_event_index_from_rows(const RawInteractionReader& reader) {
    return build_event_index_from_rows_impl(reader);
}

}  // namespace ndlar::hdf5
