#include "ndlar/hdf5/readers/pairs.hpp"

namespace ndlar::hdf5
{

RawRefPairReader::RawRefPairReader(hid_t file_id, const char *dataset_path) :
    RawReaderBase(file_id, dataset_path, "ref_pair dataset")
{
    // Reuse the filespace_ from the base class instead of opening a new one
    const int ndims = H5Sget_simple_extent_ndims(filespace_.get());

    if (ndims == 2)
    {
        hsize_t dims[2] = {0, 0};
        H5Sget_simple_extent_dims(filespace_.get(), dims, nullptr);
        is_2d = (dims[1] == 2);

        // Override the base class row_count
        row_count = static_cast<size_t>(dims[0]);
    }
    else
        is_2d = false;

    hsize_t dims[1] = {2};
    pair_array_type_.reset(H5Tarray_create2(H5T_NATIVE_UINT, 1, dims), H5Tclose);
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
        H5Sselect_hyperslab(filespace_.get(), H5S_SELECT_SET, start, nullptr, hcount, nullptr);

        UniqueHID memspace(H5Screate_simple(2, hcount, nullptr), H5Sclose);
        out.resize(count);

        const herr_t status = H5Dread(dset_.get(), H5T_NATIVE_UINT, memspace.get(), filespace_.get(), H5P_DEFAULT, out.data());
        if (status < 0)
        {
            out.clear();
            return false;
        }
        return true;
    }

    const hsize_t start[1] = {first_idx};
    const hsize_t hcount[1] = {count};
    H5Sselect_hyperslab(filespace_.get(), H5S_SELECT_SET, start, nullptr, hcount, nullptr);

    UniqueHID memspace(H5Screate_simple(1, hcount, nullptr), H5Sclose);
    out.resize(count);

    const herr_t status = H5Dread(dset_.get(), pair_array_type_.get(), memspace.get(), filespace_.get(), H5P_DEFAULT, out.data());
    if (status < 0)
    {
        out.clear();
        return false;
    }
    return true;
}

} // namespace ndlar::hdf5
