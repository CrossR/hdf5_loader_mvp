#pragma once

#include <memory>
#include <unordered_map>
#include <vector>

#include <highfive/H5File.hpp>

#include "ndlar/hdf5/paths.hpp"
#include "ndlar/hdf5/readers.hpp"

namespace ndlar::hdf5
{

// Working buffers, readers, and caches reused while collecting one event at a time.
struct StreamingContext
{
    bool has_mc = true;

    std::vector<EventRow> events;
    std::vector<ExtTrig> ext_trigs;
    std::vector<RefRegion> hit_event_bounds;
    std::vector<RefRegion> event_to_exttrig_reg;
    std::vector<RefPair> event_to_exttrig_ref;

    size_t last_fraction_block_base = SIZE_MAX;
    std::vector<PacketFraction> *last_fraction_block_ptr = nullptr;

    std::unique_ptr<RawCaloHitReader> calo_hit_reader;
    std::unique_ptr<RawTrueSegmentReader> segment_reader;
    std::unique_ptr<RawRefRegionReader> hit_to_pkt_reg_reader;
    std::unique_ptr<RawRefPairReader> hit_to_pkt_ref_reader;
    std::unique_ptr<RawRefRegionReader> pkt_to_seg_reg_reader;
    std::unique_ptr<RawRefPairReader> pkt_to_seg_ref_reader;
    std::unique_ptr<RawRefRegionReader> hit_to_btrk_reg_reader;
    std::unique_ptr<RawRefPairReader> hit_to_btrk_ref_reader;

    std::unordered_map<int64_t, std::vector<size_t>> seg_rows_by_event;
    std::unordered_map<int64_t, std::vector<size_t>> traj_rows_by_event;
    std::unordered_map<int64_t, std::vector<size_t>> int_rows_by_event;

    std::vector<TrueSegment> segment_cache;
    std::vector<uint8_t> segment_cache_valid;
    std::unordered_map<size_t, std::vector<PacketFraction>> fraction_blocks;

    std::vector<CaloHit> event_hits;
    std::vector<RefRegion> hit_pkt_regs;
    std::vector<RefPair> hit_pkt_refs;
    std::vector<RefRegion> pkt_seg_regs;
    std::vector<RefPair> pkt_seg_refs;
    std::vector<RefRegion> hit_btrk_regs;
    std::vector<RefPair> hit_btrk_refs;
    std::vector<uint32_t> hit_to_btrk_map;
    std::vector<uint32_t> pkt_ids;
    std::vector<RefRegion> pkt_frac_regs;
    std::vector<RefPair> pkt_frac_refs;
    std::vector<std::vector<uint32_t>> seg_ids;
    std::vector<std::vector<uint32_t>> frac_ids;
    std::vector<uint16_t> match_counts;
    std::vector<size_t> needed_seg_ids;
    std::vector<size_t> needed_frac_ids;
    std::vector<Trajectory> trajectory_rows;
    std::vector<Interaction> interaction_rows;
    std::vector<PacketFraction> fraction_rows;

    std::unique_ptr<RawPacketFractionReader> frac_reader;
    std::unique_ptr<RawTrajectoryReader> traj_reader;
    std::unique_ptr<RawInteractionReader> int_reader;

    bool is_setup() const
    {
        bool ok = true;

        // Must have the core charge data readers
        if (!calo_hit_reader)
        {
            std::cerr << "Missing or failed to load: calo_hit_reader (hits dataset)!\n";
            ok = false;
        }
        if (events.empty())
        {
            std::cerr << "Missing or empty: events dataset!\n";
            ok = false;
        }
        if (hit_event_bounds.empty())
        {
            std::cerr << "Missing or empty: hit_event_bounds (event -> hits references)!\n";
            ok = false;
        }

        // If we aren't using MC, stop here.
        if (!has_mc)
            return ok;

        // Otherwise, check the MC readers too, so the error messages are complete.
        if (!segment_reader)
        {
            std::cerr << "Missing: segment_reader! If this is a data file, use --data!\n";
            ok = false;
        }
        if (!hit_to_btrk_reg_reader || !hit_to_btrk_ref_reader)
        {
            std::cerr << "Missing: backtrack references! If this is a data file, use --data!\n";
            ok = false;
        }
        if (!hit_to_pkt_reg_reader || !hit_to_pkt_ref_reader)
        {
            std::cerr << "Missing: packet references! If this is a data file, use --data!\n";
            ok = false;
        }

        return ok;
    }
};

// Load file-level datasets and initialize readers/caches for streaming access.
void initialize_streaming_context(
    HighFive::File &file, StreamingContext &ctx, const paths::HitType hit_type = paths::HitType::Prompt, const bool has_mc = true);

// Clear the caches in the streaming context, to free memory between events or files.
void clear_caches(StreamingContext &ctx);

// Resolve the trigger ID for an event; returns kInvalidTrigger if none is available.
int32_t select_trigger_id_stream(const StreamingContext &ctx, size_t event_index);

// Collect one event's hit, truth-match, trajectory, and interaction products.
EventProducts collect_event_products_stream(StreamingContext &ctx, size_t event_index);

// Print a short debug summary of non-empty hit truth matches.
void print_debug_matches(const EventProducts &event_products);

} // namespace ndlar::hdf5
