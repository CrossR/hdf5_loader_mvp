#include <iostream>

#include <highfive/H5File.hpp>

#include "ndlar/hdf5/paths.hpp"
#include "ndlar/hdf5/readers.hpp"
#include "ndlar/interface.hpp"

namespace ndlar::hdf5
{

HDF5EventProvider::HDF5EventProvider(const std::string &filepath)
{
    // Open the file in read-only mode
    file_ = std::make_unique<HighFive::File>(filepath, HighFive::File::ReadOnly);

    // Initialize the streaming context (readers and ref tables)
    initialize_streaming_context(*file_, ctx_);

    // Initialize the truth readers
    frac_reader_ = std::make_unique<RawPacketFractionReader>(*file_, paths::dataset::kHitBacktrack);
    traj_reader_ = std::make_unique<RawTrajectoryReader>(*file_, paths::dataset::kTrajectories);
    int_reader_ = std::make_unique<RawInteractionReader>(*file_, paths::dataset::kInteractions);

    // Build the trajectory / interaction / segment indices across the whole file
    ctx_.traj_rows_by_event = build_event_index_from_rows(*traj_reader_);
    ctx_.int_rows_by_event = build_event_index_from_rows(*int_reader_);
}

HDF5EventProvider::~HDF5EventProvider() = default;

size_t HDF5EventProvider::get_num_events() const
{
    return ctx_.events.size();
}

EventProducts HDF5EventProvider::get_event(size_t event_index)
{
    if (event_index >= ctx_.events.size())
    {
        throw std::out_of_range("Event index out of range");
    }

    return collect_event_products_stream(ctx_, event_index, *frac_reader_, *traj_reader_, *int_reader_);
}

void HDF5EventProvider::clear_caches()
{
    ctx_.fraction_blocks.clear();
    ctx_.last_fraction_block_ptr = nullptr;
    ctx_.last_fraction_block_base = SIZE_MAX;

    ctx_.needed_seg_ids.clear();
    ctx_.needed_frac_ids.clear();

    // TODO: Is this worth doing?
    //       Check RAM floor vs timing.
    // ctx_.segment_cache.clear();
    // ctx_.segment_cache_valid.clear();
}

} // namespace ndlar::hdf5
