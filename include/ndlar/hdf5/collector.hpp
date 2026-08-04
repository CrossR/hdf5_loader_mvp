#pragma once

#include <memory>
#include <unordered_map>
#include <vector>

#include <highfive/H5File.hpp>

#include "ndlar/hdf5/readers.hpp"

namespace ndlar::hdf5 {

struct StreamingContext {
    std::vector<EventRow> events;
    std::vector<ExtTrig> ext_trigs;
    std::vector<RefRegion> hit_event_bounds;
    std::vector<RefRegion> event_to_exttrig_reg;
    std::vector<RefPair> event_to_exttrig_ref;

    HighFive::DataSet dset_hits;
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
    std::vector<uint32_t> seg_ids;
    std::vector<uint32_t> frac_ids;
    std::vector<uint16_t> match_counts;
    std::vector<size_t> needed_seg_ids;
    std::vector<size_t> needed_frac_ids;
    std::vector<Trajectory> trajectory_rows;
    std::vector<Interaction> interaction_rows;
    std::vector<PacketFraction> fraction_rows;
};

void initialize_streaming_context(HighFive::File& file, StreamingContext& ctx);

int32_t select_trigger_id_stream(const StreamingContext& ctx, size_t event_index);

EventProducts collect_event_products_stream(
    StreamingContext& ctx,
    size_t event_index,
    const RawPacketFractionReader& frac_reader,
    const RawTrajectoryReader& traj_reader,
    const RawInteractionReader& int_reader);

void print_debug_matches(const EventProducts& event_products);
void print_collect_breakdown(const EventProducts& event_products);

}  // namespace ndlar::hdf5
