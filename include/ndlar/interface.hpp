#pragma once

#include <memory>
#include <string>

#include "ndlar/hdf5/collector.hpp"
#include "ndlar/hdf5/types.hpp"
#include "ndlar/hdf5/paths.hpp"

namespace HighFive
{
class File;
}

namespace ndlar::hdf5
{

class HDF5EventProvider
{
public:
    // Opens the HDF5 file and builds the event indices
    explicit HDF5EventProvider(const std::string &filepath, paths::HitType hit_type = paths::HitType::Prompt);
    ~HDF5EventProvider();

    // Returns the total number of events in the file
    size_t get_num_events() const;

    // Collects and returns all event info (hits, truth trajectories, and interactions) for a single event
    EventProducts get_event(size_t event_index);

    // Clears the RAM caches
    void clear_caches();

private:
    std::unique_ptr<HighFive::File> file_;
    StreamingContext ctx_;

    std::unique_ptr<RawPacketFractionReader> frac_reader_;
    std::unique_ptr<RawTrajectoryReader> traj_reader_;
    std::unique_ptr<RawInteractionReader> int_reader_;
};

} // namespace ndlar::hdf5
