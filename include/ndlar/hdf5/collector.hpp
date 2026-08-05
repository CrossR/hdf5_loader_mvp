#pragma once

#include <memory>
#include <unordered_map>
#include <vector>

#include <highfive/H5File.hpp>

#include "ndlar/hdf5/readers/interaction.hpp"
#include "ndlar/hdf5/readers/packet_fraction.hpp"
#include "ndlar/hdf5/readers/pairs.hpp"
#include "ndlar/hdf5/readers/prompt_hit.hpp"
#include "ndlar/hdf5/readers/region.hpp"
#include "ndlar/hdf5/readers/segment.hpp"
#include "ndlar/hdf5/readers/trajectory.hpp"

namespace ndlar::hdf5
{

// Working buffers, readers, and caches reused while collecting one event at a time.
struct StreamingContext
{
    std::vector<EventRow> events;
    std::vector<ExtTrig> ext_trigs;
    std::vector<RefRegion> hit_event_bounds;
    std::vector<RefRegion> event_to_exttrig_reg;
    std::vector<RefPair> event_to_exttrig_ref;

    size_t last_fraction_block_base = SIZE_MAX;
    std::vector<PacketFraction> *last_fraction_block_ptr = nullptr;

    std::unique_ptr<RawPromptHitReader> prompt_hit_reader;
    std::unique_ptr<RawRefRegionReader> hit_to_pkt_reg_reader;
    std::unique_ptr<RawRefPairReader> hit_to_pkt_ref_reader;
    std::unique_ptr<RawRefRegionReader> pkt_to_seg_reg_reader;
    std::unique_ptr<RawRefPairReader> pkt_to_seg_ref_reader;
    std::unique_ptr<RawRefRegionReader> pkt_to_frac_reg_reader;
    std::unique_ptr<RawRefPairReader> pkt_to_frac_ref_reader;
    std::unique_ptr<RawTrueSegmentReader> segment_reader;

    std::unordered_map<int64_t, std::vector<size_t>> traj_rows_by_event;
    std::unordered_map<int64_t, std::vector<size_t>> int_rows_by_event;

    std::vector<TrueSegment> segment_cache;
    std::vector<uint8_t> segment_cache_valid;
    std::unordered_map<size_t, std::vector<PacketFraction>> fraction_blocks;

    std::vector<PromptHit> event_hits;
    std::vector<RefRegion> hit_pkt_regs;
    std::vector<RefPair> hit_pkt_refs;
    std::vector<uint32_t> pkt_ids;
    std::vector<RefRegion> pkt_seg_regs;
    std::vector<RefRegion> pkt_frac_regs;
    std::vector<RefPair> pkt_seg_refs;
    std::vector<RefPair> pkt_frac_refs;
    std::vector<std::vector<uint32_t>> seg_ids;
    std::vector<std::vector<uint32_t>> frac_ids;
    std::vector<uint16_t> match_counts;
    std::vector<size_t> needed_seg_ids;
    std::vector<size_t> needed_frac_ids;
    std::vector<Trajectory> trajectory_rows;
    std::vector<Interaction> interaction_rows;
    std::vector<PacketFraction> fraction_rows;
};

// Load file-level datasets and initialize readers/caches for streaming access.
void initialize_streaming_context(HighFive::File &file, StreamingContext &ctx);

// Resolve the trigger ID for an event; returns kInvalidTrigger if none is available.
int32_t select_trigger_id_stream(const StreamingContext &ctx, size_t event_index);

// Collect one event's hit, truth-match, trajectory, and interaction products.
EventProducts collect_event_products_stream(StreamingContext &ctx, size_t event_index, const RawPacketFractionReader &frac_reader,
    const RawTrajectoryReader &traj_reader, const RawInteractionReader &int_reader);

// Print a short debug summary of non-empty hit truth matches.
void print_debug_matches(const EventProducts &event_products);

} // namespace ndlar::hdf5
