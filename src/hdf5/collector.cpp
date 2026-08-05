#include "ndlar/hdf5/collector.hpp"

#include <algorithm>
#include <iostream>

#include "ndlar/common.hpp"
#include "ndlar/hdf5/highfive_types.hpp"
#include "ndlar/hdf5/paths.hpp"
#include "ndlar/hdf5/readers.hpp"

namespace ndlar::hdf5
{

namespace
{

std::pair<int32_t, int32_t> min_start_max_stop(const std::vector<RefRegion> &regions)
{
    int32_t min_start = INT32_MAX;
    int32_t max_stop = INT32_MIN;
    for (const RefRegion &region : regions)
    {
        if (!is_valid_region(region))
        {
            continue;
        }
        min_start = std::min(min_start, region.start);
        max_stop = std::max(max_stop, region.stop);
    }
    return {min_start, max_stop};
}

void load_ref_pairs_window(const RawRefPairReader *reader, int32_t min_start, int32_t max_stop, size_t &out_base, std::vector<RefPair> &out_rows)
{
    out_rows.clear();
    out_base = 0;
    if (reader == nullptr || min_start == INT32_MAX || max_stop <= min_start)
    {
        return;
    }

    out_base = static_cast<size_t>(min_start);
    const size_t count = static_cast<size_t>(max_stop - min_start);
    reader->read_rows(out_base, count, out_rows);
}

bool read_event_hits(StreamingContext &ctx, size_t event_index)
{
    if (event_index >= ctx.hit_event_bounds.size())
    {
        return false;
    }
    const RefRegion event_bounds = ctx.hit_event_bounds[event_index];
    if (!is_valid_region(event_bounds))
    {
        return false;
    }
    const size_t start = static_cast<size_t>(event_bounds.start);
    const size_t stop = static_cast<size_t>(event_bounds.stop);
    const size_t hit_count = stop - start;

    ctx.event_hits.clear();
    ctx.prompt_hit_reader->read_rows(start, hit_count, ctx.event_hits);
    return !ctx.event_hits.empty();
}

void resolve_hit_references(StreamingContext &ctx)
{
    uint32_t min_hit_id = UINT32_MAX;
    uint32_t max_hit_id = 0;
    for (const PromptHit &hit : ctx.event_hits)
    {
        if (ctx.hit_to_pkt_reg_reader != nullptr && hit.id < ctx.hit_to_pkt_reg_reader->row_count)
        {
            min_hit_id = std::min(min_hit_id, hit.id);
            max_hit_id = std::max(max_hit_id, hit.id);
        }
    }

    ctx.hit_pkt_regs.clear();
    size_t hit_reg_base = 0;
    if (min_hit_id != UINT32_MAX)
    {
        hit_reg_base = static_cast<size_t>(min_hit_id);
        const size_t n = static_cast<size_t>(max_hit_id - min_hit_id + 1);
        ctx.hit_to_pkt_reg_reader->read_rows(hit_reg_base, n, ctx.hit_pkt_regs);
    }

    const auto [min_hit_ref, max_hit_ref_stop] = min_start_max_stop(ctx.hit_pkt_regs);
    size_t hit_ref_base = 0;
    load_ref_pairs_window(ctx.hit_to_pkt_ref_reader.get(), min_hit_ref, max_hit_ref_stop, hit_ref_base, ctx.hit_pkt_refs);

    ctx.pkt_ids.assign(ctx.event_hits.size(), UINT32_MAX);
    uint32_t min_pkt_id = UINT32_MAX;
    uint32_t max_pkt_id = 0;

    for (size_t i = 0; i < ctx.event_hits.size(); ++i)
    {
        const uint32_t hit_id = ctx.event_hits[i].id;
        if (hit_id < hit_reg_base || static_cast<size_t>(hit_id - hit_reg_base) >= ctx.hit_pkt_regs.size())
            continue;

        const RefRegion r = ctx.hit_pkt_regs[static_cast<size_t>(hit_id - hit_reg_base)];
        if (!is_valid_region(r) || r.start < static_cast<int32_t>(hit_ref_base))
            continue;

        const size_t ref_idx = static_cast<size_t>(r.start - static_cast<int32_t>(hit_ref_base));
        if (ref_idx >= ctx.hit_pkt_refs.size())
            continue;

        const uint32_t pkt_id = ctx.hit_pkt_refs[ref_idx][1];
        if (pkt_id >= ctx.pkt_to_seg_reg_reader->row_count || pkt_id >= ctx.pkt_to_frac_reg_reader->row_count)
            continue;

        ctx.pkt_ids[i] = pkt_id;
        min_pkt_id = std::min(min_pkt_id, pkt_id);
        max_pkt_id = std::max(max_pkt_id, pkt_id);
    }

    ctx.pkt_seg_regs.clear();
    ctx.pkt_frac_regs.clear();
    size_t pkt_reg_base = 0;

    if (min_pkt_id != UINT32_MAX)
    {
        pkt_reg_base = static_cast<size_t>(min_pkt_id);
        const size_t n = static_cast<size_t>(max_pkt_id - min_pkt_id + 1);
        ctx.pkt_to_seg_reg_reader->read_rows(pkt_reg_base, n, ctx.pkt_seg_regs);
        ctx.pkt_to_frac_reg_reader->read_rows(pkt_reg_base, n, ctx.pkt_frac_regs);
    }

    const auto [min_seg_ref, max_seg_ref_stop] = min_start_max_stop(ctx.pkt_seg_regs);
    const auto [min_frac_ref, max_frac_ref_stop] = min_start_max_stop(ctx.pkt_frac_regs);

    size_t seg_ref_base = 0;
    size_t frac_ref_base = 0;
    load_ref_pairs_window(ctx.pkt_to_seg_ref_reader.get(), min_seg_ref, max_seg_ref_stop, seg_ref_base, ctx.pkt_seg_refs);
    load_ref_pairs_window(ctx.pkt_to_frac_ref_reader.get(), min_frac_ref, max_frac_ref_stop, frac_ref_base, ctx.pkt_frac_refs);

    ctx.seg_ids.assign(ctx.event_hits.size(), {});
    ctx.frac_ids.assign(ctx.event_hits.size(), {});
    ctx.match_counts.assign(ctx.event_hits.size(), 0);
    ctx.needed_seg_ids.clear();
    ctx.needed_frac_ids.clear();
    ctx.needed_seg_ids.reserve(ctx.event_hits.size() * 2); // Estimate 2 matches per hit
    ctx.needed_frac_ids.reserve(ctx.event_hits.size() * 2);

    for (size_t i = 0; i < ctx.event_hits.size(); ++i)
    {
        const uint32_t pkt_id = ctx.pkt_ids[i];
        if (pkt_id == UINT32_MAX || pkt_id < pkt_reg_base)
            continue;

        const size_t pkt_off = static_cast<size_t>(pkt_id - pkt_reg_base);
        if (pkt_off >= ctx.pkt_seg_regs.size() || pkt_off >= ctx.pkt_frac_regs.size())
            continue;

        const RefRegion seg_r = ctx.pkt_seg_regs[pkt_off];
        const RefRegion frac_r = ctx.pkt_frac_regs[pkt_off];

        if (is_valid_region(seg_r))
        {
            const size_t num_matches = region_size(seg_r);
            ctx.match_counts[i] = static_cast<uint16_t>(num_matches);

            for (size_t m = 0; m < num_matches; ++m)
            {
                int32_t s_idx = seg_r.start + static_cast<int32_t>(m);
                if (s_idx >= static_cast<int32_t>(seg_ref_base))
                {
                    size_t ref_idx = static_cast<size_t>(s_idx - static_cast<int32_t>(seg_ref_base));
                    if (ref_idx < ctx.pkt_seg_refs.size())
                    {
                        uint32_t id = ctx.pkt_seg_refs[ref_idx][1];
                        ctx.seg_ids[i].push_back(id);
                        ctx.needed_seg_ids.push_back(static_cast<size_t>(id));
                    }
                }

                if (is_valid_region(frac_r) && m < static_cast<size_t>(region_size(frac_r)))
                {
                    int32_t f_idx = frac_r.start + static_cast<int32_t>(m);
                    if (f_idx >= static_cast<int32_t>(frac_ref_base))
                    {
                        size_t ref_idx = static_cast<size_t>(f_idx - static_cast<int32_t>(frac_ref_base));
                        if (ref_idx < ctx.pkt_frac_refs.size())
                        {
                            uint32_t id = ctx.pkt_frac_refs[ref_idx][1];
                            ctx.frac_ids[i].push_back(id);
                            ctx.needed_frac_ids.push_back(static_cast<size_t>(id));
                        }
                    }
                }
            }
        }
    }

    std::sort(ctx.needed_seg_ids.begin(), ctx.needed_seg_ids.end());
    ctx.needed_seg_ids.erase(std::unique(ctx.needed_seg_ids.begin(), ctx.needed_seg_ids.end()), ctx.needed_seg_ids.end());

    std::sort(ctx.needed_frac_ids.begin(), ctx.needed_frac_ids.end());
    ctx.needed_frac_ids.erase(std::unique(ctx.needed_frac_ids.begin(), ctx.needed_frac_ids.end()), ctx.needed_frac_ids.end());
}

void populate_truth_products(StreamingContext &ctx, EventProducts &out, const RawTrajectoryReader &traj_reader, const RawInteractionReader &int_reader)
{
    if (out.spill_id < 0)
        return;

    const auto traj_it = ctx.traj_rows_by_event.find(out.spill_id);
    if (traj_it != ctx.traj_rows_by_event.end())
    {
        out.reserve_trajectory_products(traj_it->second.size());
        const auto spans = contiguous_spans(traj_it->second, ndlar::kCacheReadGapTolerance);
        for (const auto &span : spans)
        {
            if (traj_reader.read_rows(span[0], span[1], ctx.trajectory_rows))
            {
                append_trajectory_products(ctx.trajectory_rows, out);
            }
        }
    }

    const auto int_it = ctx.int_rows_by_event.find(out.spill_id);
    if (int_it != ctx.int_rows_by_event.end())
    {
        out.reserve_interaction_products(int_it->second.size());
        const auto spans = contiguous_spans(int_it->second, ndlar::kCacheReadGapTolerance);
        for (const auto &span : spans)
        {
            if (int_reader.read_rows(span[0], span[1], ctx.interaction_rows))
            {
                append_interaction_products(ctx.interaction_rows, out);
            }
        }
    }
}

} // namespace

void print_debug_matches(const EventProducts &event_products)
{
    int debug_printed = 0;
    for (size_t hit_index = 0; hit_index < event_products.hit_pdg.size() && debug_printed < ndlar::kDebugMatchPrintLimit; ++hit_index)
    {
        // Skip hits with no matches or where the first match is segment 0
        if (event_products.hit_segmentID[hit_index].empty() || event_products.hit_segmentID[hit_index][0] == 0)
        {
            continue;
        }

        std::cout << "  Hit " << hit_index << " -> PDG: " << event_products.hit_pdg[hit_index][0]
                  << ", frac: " << event_products.hit_packetFrac[hit_index][0] << "\n";
        ++debug_printed;
    }
}

size_t fraction_block_base(size_t row_id)
{
    return (row_id / ndlar::kFractionBlockRows) * ndlar::kFractionBlockRows;
}

const PacketFraction *get_cached_fraction_row(StreamingContext &ctx, size_t row_id)
{
    const size_t block_base = fraction_block_base(row_id);

    // Check fast 1-element cache first
    if (ctx.last_fraction_block_base == block_base && ctx.last_fraction_block_ptr != nullptr)
    {
        const size_t offset = row_id - block_base;
        if (offset >= ctx.last_fraction_block_ptr->size())
            return nullptr;
        return &(*ctx.last_fraction_block_ptr)[offset];
    }

    // Fallback to the slow map lookup
    const auto block_it = ctx.fraction_blocks.find(block_base);
    if (block_it == ctx.fraction_blocks.end())
    {
        return nullptr;
    }

    // Update the fast cache
    ctx.last_fraction_block_base = block_base;
    ctx.last_fraction_block_ptr = &(block_it->second);

    const size_t offset = row_id - block_base;
    if (offset >= block_it->second.size())
    {
        return nullptr;
    }
    return &block_it->second[offset];
}

void populate_hit_products(StreamingContext &ctx, EventProducts &out)
{
    for (size_t i = 0; i < ctx.event_hits.size(); ++i)
    {
        const PromptHit &hit = ctx.event_hits[i];
        out.hit_x[i] = hit.x;
        out.hit_y[i] = hit.y;
        out.hit_z[i] = hit.z;
        out.hit_charge[i] = hit.Q;
        out.hit_E[i] = hit.E;
        out.hit_ts[i] = hit.ts_pps;

        size_t n_matches = ctx.seg_ids[i].size();
        out.hit_matches[i] = ctx.match_counts[i];

        // Resize the inner vectors for the current hit
        out.hit_pdg[i].resize(n_matches, 0);
        out.hit_segmentID[i].resize(n_matches, 0);
        out.hit_particleID[i].resize(n_matches, 0);
        out.hit_particleIDLocal[i].resize(n_matches, 0);
        out.hit_vertexID[i].resize(n_matches, 0);
        out.hit_packetFrac[i].resize(n_matches, 0.0f);

        for (size_t m = 0; m < n_matches; ++m)
        {
            const uint32_t seg_id = ctx.seg_ids[i][m];

            if (seg_id < ctx.segment_cache_valid.size() && ctx.segment_cache_valid[seg_id])
            {
                const TrueSegment &true_seg = ctx.segment_cache[seg_id];

                out.hit_pdg[i][m] = true_seg.pdg_id;
                out.hit_segmentID[i][m] = static_cast<int32_t>(true_seg.segment_id);
                out.hit_particleID[i][m] = static_cast<int64_t>(true_seg.file_traj_id);
                out.hit_particleIDLocal[i][m] = static_cast<int64_t>(true_seg.traj_id);
                out.hit_vertexID[i][m] = static_cast<int64_t>(true_seg.vertex_id);

                if (out.spill_id < 0)
                {
                    out.spill_id = true_seg.event_id;
                }

                // Safely grab fraction
                if (!ctx.frac_ids[i].empty()) {
                    const uint32_t frac_id = ctx.frac_ids[i][0];
                    const PacketFraction* frac_row = get_cached_fraction_row(ctx, static_cast<size_t>(frac_id));
                    if (frac_row != nullptr) {
                        out.hit_packetFrac[i][m] = resolve_packet_fraction(*frac_row, true_seg.segment_id);
                    }
                }
            }
        }
    }
}

void ensure_fraction_range_cached(StreamingContext &ctx, const RawPacketFractionReader &frac_reader, size_t first_id, size_t last_id_inclusive)
{
    const size_t first_block = fraction_block_base(first_id);
    const size_t last_block = fraction_block_base(last_id_inclusive);

    for (size_t block_base = first_block; block_base <= last_block; block_base += ndlar::kFractionBlockRows)
    {
        if (ctx.fraction_blocks.find(block_base) != ctx.fraction_blocks.end())
        {
            continue;
        }

        const size_t available = frac_reader.row_count - block_base;
        const size_t count = std::min(ndlar::kFractionBlockRows, available);
        auto &rows = ctx.fraction_rows;
        if (!frac_reader.read_rows(block_base, count, rows))
        {
            continue;
        }
        ctx.fraction_blocks.emplace(block_base, rows);
    }
}

void update_caches(StreamingContext &ctx, const RawPacketFractionReader &frac_reader)
{
    const auto seg_spans = contiguous_spans(ctx.needed_seg_ids, ndlar::kCacheReadGapTolerance);
    for (const auto &span : seg_spans)
    {
        if (span[0] + span[1] > ctx.segment_reader->row_count)
            continue;

        bool all_cached = true;
        for (size_t id = span[0]; id < span[0] + span[1]; ++id)
        {
            if (!ctx.segment_cache_valid[id])
            {
                all_cached = false;
                break;
            }
        }
        if (all_cached)
            continue;

        std::vector<TrueSegment> rows;
        if (!ctx.segment_reader->read_rows(span[0], span[1], rows))
            continue;

        for (size_t i = 0; i < rows.size(); ++i)
        {
            const size_t id = span[0] + i;
            ctx.segment_cache[id] = rows[i];
            ctx.segment_cache_valid[id] = 1;
        }
    }

    const auto frac_spans = contiguous_spans(ctx.needed_frac_ids, ndlar::kCacheReadGapTolerance);
    for (const auto &span : frac_spans)
    {
        if (span[0] + span[1] > frac_reader.row_count)
            continue;
        ensure_fraction_range_cached(ctx, frac_reader, span[0], span[0] + span[1] - 1);
    }
}

void initialize_streaming_context(HighFive::File &file, StreamingContext &ctx)
{
    ctx.prompt_hit_reader = std::make_unique<RawPromptHitReader>(file.getId());
    ctx.hit_to_pkt_reg_reader = std::make_unique<RawRefRegionReader>(file.getId(), paths::ref_region::kHitToPacket);
    ctx.hit_to_pkt_ref_reader = std::make_unique<RawRefPairReader>(file.getId(), paths::ref_data::kHitToPacket);
    ctx.pkt_to_seg_reg_reader = std::make_unique<RawRefRegionReader>(file.getId(), paths::ref_region::kPacketToSegment);
    ctx.pkt_to_seg_ref_reader = std::make_unique<RawRefPairReader>(file.getId(), paths::ref_data::kPacketToSegment);
    ctx.pkt_to_frac_reg_reader = std::make_unique<RawRefRegionReader>(file.getId(), paths::ref_region::kPacketToFraction);
    ctx.pkt_to_frac_ref_reader = std::make_unique<RawRefPairReader>(file.getId(), paths::ref_data::kPacketToFraction);
    ctx.segment_reader = std::make_unique<RawTrueSegmentReader>(file.getId(), paths::dataset::kSegments);

    file.getDataSet(paths::dataset::kEvents).read(ctx.events);
    file.getDataSet(paths::dataset::kExtTrigs).read(ctx.ext_trigs);
    file.getDataSet(paths::ref_region::kEventToHits).read(ctx.hit_event_bounds);
    file.getDataSet(paths::ref_region::kEventToExtTrigs).read(ctx.event_to_exttrig_reg);
    file.getDataSet(paths::ref_data::kEventToExtTrigs).read(ctx.event_to_exttrig_ref);

    ctx.segment_cache.resize(ctx.segment_reader->row_count);
    ctx.segment_cache_valid.assign(ctx.segment_reader->row_count, 0);
}

int32_t select_trigger_id_stream(const StreamingContext &ctx, size_t event_index)
{
    if (event_index >= ctx.event_to_exttrig_reg.size())
    {
        return ndlar::kInvalidTrigger;
    }

    const RefRegion trig_bounds = ctx.event_to_exttrig_reg[event_index];
    if (!is_valid_region(trig_bounds))
    {
        return ndlar::kInvalidTrigger;
    }

    const size_t start = static_cast<size_t>(trig_bounds.start);
    const size_t stop = static_cast<size_t>(trig_bounds.stop);
    if (start >= ctx.event_to_exttrig_ref.size())
    {
        return ndlar::kInvalidTrigger;
    }

    const size_t clamped_stop = std::min(stop, ctx.event_to_exttrig_ref.size());
    int32_t first = ndlar::kInvalidTrigger;
    bool any = false;
    for (size_t idx = start; idx < clamped_stop; ++idx)
    {
        const uint32_t ext_idx = ctx.event_to_exttrig_ref[idx][1];
        if (ext_idx >= ctx.ext_trigs.size())
        {
            continue;
        }
        const int32_t iogroup = ctx.ext_trigs[ext_idx].iogroup;
        if (first == ndlar::kInvalidTrigger)
        {
            first = iogroup;
        }
        any = true;
        if (iogroup == 5)
        {
            return 5;
        }
    }

    return any ? first : ndlar::kInvalidTrigger;
}

EventProducts collect_event_products_stream(StreamingContext &ctx, size_t event_index, const RawPacketFractionReader &frac_reader,
    const RawTrajectoryReader &traj_reader, const RawInteractionReader &int_reader)
{
    EventProducts out;
    const auto &evt = ctx.events[event_index];
    out.trigger_id = select_trigger_id_stream(ctx, event_index);
    out.event_start_t = evt.ts_start;
    out.event_end_t = evt.ts_end;
    out.unix_ts = evt.unix_ts;
    out.unix_ts_usec = evt.unix_ts_usec;

    if (!read_event_hits(ctx, event_index))
    {
        // Empty event
        return out;
    }

    out.reserve_hit_products(ctx.event_hits.size());

    resolve_hit_references(ctx);
    update_caches(ctx, frac_reader);

    populate_hit_products(ctx, out);
    populate_truth_products(ctx, out, traj_reader, int_reader);

    return out;
}

} // namespace ndlar::hdf5
