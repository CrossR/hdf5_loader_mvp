#include "ndlar/hdf5/paths.hpp"

#include "ndlar/hdf5/readers.hpp"
#include "ndlar/hdf5/readers/packet_fraction.hpp"
#include "ndlar/hdf5/types.hpp"

namespace ndlar::hdf5
{

RawPacketFractionReader::RawPacketFractionReader(hid_t file_id) :
    RawReaderBase(file_id, paths::dataset::kPacketFraction, "packet_fraction dataset")
{
    // Dynamically calculate the array length from the struct.
    constexpr size_t kSegArraySize = sizeof(PacketFraction::segment_ids) / sizeof(PacketFraction::segment_ids[0]);
    hsize_t dims[1] = {kSegArraySize};

    // Define the arrays as unique HIDs to ensure proper cleanup automatically.
    UniqueHID seg_array(H5Tarray_create2(H5T_NATIVE_LLONG, 1, dims), H5Tclose);
    UniqueHID frac_array(H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, dims), H5Tclose);

    mem_type_.reset(H5Tcreate(H5T_COMPOUND, sizeof(PacketFraction)), H5Tclose);
    H5Tinsert(mem_type_.get(), "segment_ids", HOFFSET(PacketFraction, segment_ids), seg_array.get());
    H5Tinsert(mem_type_.get(), "fraction", HOFFSET(PacketFraction, fraction), frac_array.get());
}

bool RawPacketFractionReader::read_rows(size_t first_idx, size_t count, std::vector<PacketFraction> &out) const
{
    return read_rows_1d_hyperslab(dset_.get(), mem_type_.get(), filespace_.get(), row_count, first_idx, count, out, "packet_fraction");
}

} // namespace ndlar::hdf5
