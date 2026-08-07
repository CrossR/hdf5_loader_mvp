#pragma once

#include <string>

namespace ndlar::hdf5::paths
{

enum class HitType
{
    Prompt,
    Merged,
    Final
};

inline std::string to_string(HitType type)
{
    switch (type)
    {
        case HitType::Merged:
            return "merged";
        case HitType::Final:
            return "final";
        default:
            return "prompt";
    }
}

// Dynamic paths that depend on the hit type
struct PathResolver
{
    std::string hit_key;

    explicit PathResolver(HitType type = HitType::Prompt) :
        hit_key(to_string(type))
    {
    }

    // Dynamic Datasets
    std::string hits() const
    {
        return "charge/calib_" + hit_key + "_hits/data";
    }
    std::string hit_backtrack() const
    {
        return "mc_truth/calib_" + hit_key + "_hit_backtrack/data";
    }

    // Dynamic Ref Regions
    std::string event_to_hits_reg() const
    {
        return "charge/events/ref/charge/calib_" + hit_key + "_hits/ref_region";
    }
    std::string hit_to_packet_reg() const
    {
        return "charge/calib_" + hit_key + "_hits/ref/charge/packets/ref_region";
    }
    std::string hit_to_backtrack_reg() const
    {
        return "charge/calib_" + hit_key + "_hits/ref/mc_truth/calib_" + hit_key + "_hit_backtrack/ref_region";
    }

    // Dynamic Ref Data
    std::string hit_to_packet_ref() const
    {
        return "charge/calib_" + hit_key + "_hits/ref/charge/packets/ref";
    }
    std::string hit_to_backtrack_ref() const
    {
        return "charge/calib_" + hit_key + "_hits/ref/mc_truth/calib_" + hit_key + "_hit_backtrack/ref";
    }
};

// Static dataset paths that do not depend on the hit type
namespace dataset
{
inline constexpr const char *kEvents = "charge/events/data";
inline constexpr const char *kExtTrigs = "charge/ext_trigs/data";
inline constexpr const char *kPacketFraction = "mc_truth/packet_fraction/data";
inline constexpr const char *kSegments = "mc_truth/segments/data";
inline constexpr const char *kTrajectories = "mc_truth/trajectories/data";
inline constexpr const char *kInteractions = "mc_truth/interactions/data";
} // namespace dataset

namespace ref_region
{
inline constexpr const char *kPacketToSegment = "charge/packets/ref/mc_truth/segments/ref_region";
inline constexpr const char *kEventToExtTrigs = "charge/events/ref/charge/ext_trigs/ref_region";
} // namespace ref_region

namespace ref_data
{
inline constexpr const char *kPacketToSegment = "charge/packets/ref/mc_truth/segments/ref";
inline constexpr const char *kEventToExtTrigs = "charge/events/ref/charge/ext_trigs/ref";
} // namespace ref_data

} // namespace ndlar::hdf5::paths
