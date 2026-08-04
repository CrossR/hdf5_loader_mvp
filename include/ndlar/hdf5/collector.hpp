#pragma once

#include <memory>
#include <unordered_map>
#include <vector>

#include <highfive/H5File.hpp>

#include "ndlar/hdf5/readers.hpp"

namespace ndlar::hdf5
{

// Context for streaming event products from an HDF5 file.
//
// Essentially, this struct holds all the necessary information to read event
// products from an HDF5 file in a streaming fashion. It contains references to
// datasets, readers for various data types, and caches for segments and
// fractions.
struct StreamingContext
{
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

/*
 * Initializes the streaming context by reading metadata from the HDF5 file.
 *
 * @param file The HDF5 file to read from.
 * @param ctx The streaming context to initialize.
*/
void initialize_streaming_context(HighFive::File &file, StreamingContext &ctx);

/*
 * Selects the trigger ID for a given event index from the streaming context.
 *
 * @param ctx The streaming context containing event information.
 * @param event_index The index of the event for which to select the trigger ID.
 * @return The selected trigger ID, or kInvalidTrigger if not found.
 */
int32_t select_trigger_id_stream(const StreamingContext &ctx, size_t event_index);

/*
 * Collects all event products for a given event index from the streaming context.
 *
 * @param ctx The streaming context containing event information.
 * @param event_index The index of the event for which to collect products.
 * @param frac_reader The reader for raw packet fractions.
 * @param traj_reader The reader for raw trajectories.
 * @param int_reader The reader for raw interactions.
 * @return An EventProducts struct containing all collected products for the event.
*/
EventProducts collect_event_products_stream(StreamingContext &ctx, size_t event_index, const RawPacketFractionReader &frac_reader,
    const RawTrajectoryReader &traj_reader, const RawInteractionReader &int_reader);

/*
 * Prints debug information about the matches found in the event products.
 *
 * @param event_products The event products containing match information.
 */
void print_debug_matches(const EventProducts &event_products);

/*
 * Prints a breakdown of the collection process for the given event products.
 *
 * @param event_products The event products for which to print the collection breakdown.
 */
void print_collect_breakdown(const EventProducts &event_products);

} // namespace ndlar::hdf5
