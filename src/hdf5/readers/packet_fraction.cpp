#include "ndlar/hdf5/paths.hpp"

#include "ndlar/hdf5/readers.hpp"
#include "ndlar/hdf5/readers/packet_fraction.hpp"

namespace ndlar::hdf5
{

RawPacketFractionReader::RawPacketFractionReader(hid_t file_id)
{
    dset = open_dataset_or_throw(file_id, paths::dataset::kPacketFraction, "packet_fraction dataset");
    filespace = get_filespace_or_throw(dset, "packet_fraction dataset");
    row_count = dataset_row_count_or_throw(dset, "packet_fraction dataset");

    hsize_t dims[1] = {20};
    seg_array = H5Tarray_create2(H5T_NATIVE_LLONG, 1, dims);
    frac_array = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, dims);
    mem_type = H5Tcreate(H5T_COMPOUND, sizeof(PacketFraction));
    H5Tinsert(mem_type, "segment_ids", HOFFSET(PacketFraction, segment_ids), seg_array);
    H5Tinsert(mem_type, "fraction", HOFFSET(PacketFraction, fraction), frac_array);
}

RawPacketFractionReader::~RawPacketFractionReader()
{
    if (mem_type >= 0)
        H5Tclose(mem_type);
    if (frac_array >= 0)
        H5Tclose(frac_array);
    if (seg_array >= 0)
        H5Tclose(seg_array);
    if (dset >= 0)
        H5Dclose(dset);
    if (filespace >= 0)
        H5Sclose(filespace);
}

bool RawPacketFractionReader::read_rows(size_t first_idx, size_t count, std::vector<PacketFraction> &out) const
{
    return read_rows_1d_hyperslab(dset, mem_type, filespace, row_count, first_idx, count, out, "packet_fraction");
}

} // namespace ndlar::hdf5
