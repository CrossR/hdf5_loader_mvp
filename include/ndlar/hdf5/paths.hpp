#pragma once

namespace ndlar::hdf5::paths
{

namespace dataset
{
inline constexpr const char *kPromptHits = "charge/calib_prompt_hits/data";
inline constexpr const char *kEvents = "charge/events/data";
inline constexpr const char *kExtTrigs = "charge/ext_trigs/data";
inline constexpr const char *kPacketFraction = "mc_truth/packet_fraction/data";
inline constexpr const char *kHitBacktrack = "mc_truth/calib_prompt_hit_backtrack/data";
inline constexpr const char *kSegments = "mc_truth/segments/data";
inline constexpr const char *kTrajectories = "mc_truth/trajectories/data";
inline constexpr const char *kInteractions = "mc_truth/interactions/data";
} // namespace dataset

namespace ref_region
{
inline constexpr const char *kHitToPacket = "charge/calib_prompt_hits/ref/charge/packets/ref_region";
inline constexpr const char *kPacketToSegment = "charge/packets/ref/mc_truth/segments/ref_region";
inline constexpr const char *kPacketToFraction = "charge/packets/ref/mc_truth/packet_fraction/ref_region";
inline constexpr const char *kEventToHits = "charge/events/ref/charge/calib_prompt_hits/ref_region";
inline constexpr const char *kEventToExtTrigs = "charge/events/ref/charge/ext_trigs/ref_region";
inline constexpr const char *kHitToBacktrack = "charge/calib_prompt_hits/ref/mc_truth/calib_prompt_hit_backtrack/ref_region";
} // namespace ref_region

namespace ref_data
{
inline constexpr const char *kHitToPacket = "charge/calib_prompt_hits/ref/charge/packets/ref";
inline constexpr const char *kPacketToSegment = "charge/packets/ref/mc_truth/segments/ref";
inline constexpr const char *kPacketToFraction = "charge/packets/ref/mc_truth/packet_fraction/ref";
inline constexpr const char *kEventToExtTrigs = "charge/events/ref/charge/ext_trigs/ref";
inline constexpr const char *kHitToBacktrack = "charge/calib_prompt_hits/ref/mc_truth/calib_prompt_hit_backtrack/ref";
} // namespace ref_data

} // namespace ndlar::hdf5::paths
