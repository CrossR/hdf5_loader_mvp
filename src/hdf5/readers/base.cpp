
#include <stdexcept>
#include <string>

#include "ndlar/hdf5/readers.hpp"
#include "ndlar/hdf5/readers/base.hpp"

namespace ndlar::hdf5
{

namespace
{

hid_t open_dataset_or_throw(hid_t file_id, const char *dataset_path, const char *error_context)
{
    const hid_t dset = H5Dopen2(file_id, dataset_path, H5P_DEFAULT);
    if (dset < 0)
    {
        throw std::runtime_error(std::string("Failed to open ") + error_context + ": " + dataset_path);
    }
    return dset;
}

hid_t get_filespace_or_throw(hid_t dset, const char *error_context)
{
    const hid_t filespace = H5Dget_space(dset);
    if (filespace < 0)
    {
        throw std::runtime_error(std::string("Failed to get filespace for ") + error_context);
    }
    return filespace;
}

size_t dataset_row_count_or_throw(hid_t dset, const char *error_context)
{
    hid_t space = H5Dget_space(dset);
    if (space < 0)
    {
        throw std::runtime_error(std::string("Failed to get dataspace for ") + error_context);
    }
    const hssize_t nrows = H5Sget_simple_extent_npoints(space);
    H5Sclose(space);
    if (nrows < 0)
    {
        throw std::runtime_error(std::string("Failed to get row count for ") + error_context);
    }
    return static_cast<size_t>(nrows);
}

} // anonymous namespace

RawReaderBase::RawReaderBase(hid_t file_id, const char *dataset_path, const char *error_context)
{
    dset_.reset(open_dataset_or_throw(file_id, dataset_path, error_context), H5Dclose);
    filespace_.reset(get_filespace_or_throw(dset_.get(), error_context), H5Sclose);
    row_count = dataset_row_count_or_throw(dset_.get(), error_context);
}

} // namespace ndlar::hdf5
