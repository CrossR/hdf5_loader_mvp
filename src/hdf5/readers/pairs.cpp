
#include "ndlar/hdf5/readers.hpp"
#include "ndlar/hdf5/readers/pairs.hpp"

namespace ndlar::hdf5
{

RawRefPairReader::RawRefPairReader(hid_t file_id, const char *dataset_path)
{
    dset = open_dataset_or_throw(file_id, dataset_path, "ref_pair dataset");
    filespace = get_filespace_or_throw(dset, dataset_path);

    hid_t space = H5Dget_space(dset);
    if (space < 0)
    {
        throw std::runtime_error(std::string("Failed to get dataspace: ") + dataset_path);
    }

    const int ndims = H5Sget_simple_extent_ndims(space);
    if (ndims == 2)
    {
        hsize_t dims[2] = {0, 0};
        H5Sget_simple_extent_dims(space, dims, nullptr);
        is_2d = (dims[1] == 2);
        row_count = static_cast<size_t>(dims[0]);
    }
    else
    {
        const hssize_t nrows = H5Sget_simple_extent_npoints(space);
        if (nrows < 0)
        {
            H5Sclose(space);
            throw std::runtime_error(std::string("Failed to get row count: ") + dataset_path);
        }
        row_count = static_cast<size_t>(nrows);
    }
    H5Sclose(space);

    hsize_t dims[1] = {2};
    pair_array_type = H5Tarray_create2(H5T_NATIVE_UINT, 1, dims);
}

RawRefPairReader::~RawRefPairReader()
{
    if (pair_array_type >= 0)
        H5Tclose(pair_array_type);
    if (dset >= 0)
        H5Dclose(dset);
    if (filespace >= 0)
        H5Sclose(filespace);
}

bool RawRefPairReader::read_rows(size_t first_idx, size_t count, std::vector<RefPair> &out) const
{
    if (count == 0 || first_idx >= row_count || first_idx + count > row_count)
    {
        out.clear();
        return false;
    }

    if (is_2d)
    {
        const hsize_t start[2] = {first_idx, 0};
        const hsize_t hcount[2] = {count, 2};
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, nullptr, hcount, nullptr);
        hid_t memspace = H5Screate_simple(2, hcount, nullptr);

        out.resize(count);
        const herr_t status = H5Dread(dset, H5T_NATIVE_UINT, memspace, filespace, H5P_DEFAULT, out.data());
        H5Sclose(memspace);

        if (status < 0)
        {
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

    if (status < 0)
    {
        out.clear();
        return false;
    }
    return true;
}

} // namespace ndlar::hdf5
