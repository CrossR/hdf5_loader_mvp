#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <hdf5.h>

#include "ndlar/hdf5/unique_hid.hpp"

namespace ndlar::hdf5
{

// Forward declarations of raw dataset readers.
class RawTrajectoryReader;
class RawInteractionReader;

// Group sorted indices into [start, length] spans, allowing optional small gaps.
std::vector<std::array<size_t, 2>> contiguous_spans(const std::vector<size_t> &indices, size_t max_gap = 0);

// Build event_id -> row_indices map from trajectory rows.
std::unordered_map<int64_t, std::vector<size_t>> build_event_index_from_rows(const RawTrajectoryReader &reader);

// Build event_id -> row_indices map from interaction rows.
std::unordered_map<int64_t, std::vector<size_t>> build_event_index_from_rows(const RawInteractionReader &reader);

// Open a dataset and return its hid_t, or throw an exception if it fails.
hid_t open_dataset_or_throw(hid_t file_id, const char *dataset_path, const char *error_context);

// Open a dataset and return its filespace hid_t, or throw an exception if it fails.
hid_t get_filespace_or_throw(hid_t dset, const char *error_context);

// Get the number of rows in a dataset, or throw an exception if it fails.
size_t dataset_row_count_or_throw(hid_t dset, const char *error_context);

// Build an event index from rows of a dataset, using the provided Reader type.
template <typename Reader>
std::unordered_map<int64_t, std::vector<size_t>> build_event_index_from_rows_impl(const Reader &reader);

// Read a contiguous hyperslab of rows from a 1D dataset into a vector, or
// return false if the requested range is invalid.
template <typename T>
bool read_rows_1d_hyperslab(hid_t dset, hid_t mem_type, hid_t filespace, size_t row_count, size_t first_idx, size_t count,
    std::vector<T> &out, const char *error_context)
{
    if (count == 0 || first_idx >= row_count || first_idx + count > row_count)
    {
        out.clear();
        return false;
    }

    const hsize_t start[1] = {first_idx};
    const hsize_t hcount[1] = {count};
    const herr_t select_status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, nullptr, hcount, nullptr);
    if (select_status < 0)
    {
        throw std::runtime_error(std::string("Failed to select hyperslab for ") + error_context);
    }

    hid_t memspace = H5Screate_simple(1, hcount, nullptr);
    if (memspace < 0)
    {
        throw std::runtime_error(std::string("Failed to create memspace for ") + error_context);
    }

    out.resize(count);
    const herr_t status = H5Dread(dset, mem_type, memspace, filespace, H5P_DEFAULT, out.data());
    H5Sclose(memspace);

    if (status < 0)
    {
        out.clear();
        return false;
    }
    return true;
}

/*
 * Base class for raw HDF5 dataset readers.
 */
class RawReaderBase
{
public:
    RawReaderBase(hid_t file_id, const char *dataset_path, const char *error_context)
    {
        dset_.reset(open_dataset_or_throw(file_id, dataset_path, error_context), H5Dclose);
        filespace_.reset(get_filespace_or_throw(dset_.get(), error_context), H5Sclose);
        row_count = dataset_row_count_or_throw(dset_.get(), error_context);
    }

    virtual ~RawReaderBase() = default;
    // ^ No manual closes needed anymore! UniqueHID handles it.

    size_t row_count = 0;

protected:
    UniqueHID dset_;
    UniqueHID filespace_;
    UniqueHID mem_type_;
};

} // namespace ndlar::hdf5
