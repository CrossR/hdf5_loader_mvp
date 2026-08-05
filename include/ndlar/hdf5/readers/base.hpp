#pragma once

#include <hdf5.h>
#include <vector>

#include "ndlar/hdf5/unique_hid.hpp"

namespace ndlar::hdf5
{

class RawReaderBase
{
public:
    RawReaderBase(hid_t file_id, const char *dataset_path, const char *error_context);
    virtual ~RawReaderBase() = default;

    size_t row_count = 0;

protected:
    UniqueHID dset_;
    UniqueHID filespace_;
    UniqueHID mem_type_;
};

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
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, nullptr, hcount, nullptr);

    // Safely manage memspace with UniqueHID
    UniqueHID memspace(H5Screate_simple(1, hcount, nullptr), H5Sclose);
    out.resize(count);

    const herr_t status = H5Dread(dset, mem_type, memspace.get(), filespace, H5P_DEFAULT, out.data());

    if (status < 0)
    {
        out.clear();
        return false;
    }
    return true;
}

} // namespace ndlar::hdf5
