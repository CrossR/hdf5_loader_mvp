#include "ndlar/hdf5/collector.hpp"

#include <algorithm>
#include <iostream>

#include "ndlar/common.hpp"
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
        return false;

    const RefRegion event_bounds = ctx.hit_event_bounds[event_index];

    if (!is_valid_region(event_bounds))
        return false;

    const size_t start = static_cast<size_t>(event_bounds.start);
    const size_t stop = static_cast<size_t>(event_bounds.stop);
    const size_t hit_count = stop - start;

    ctx.event_hits.clear();
    ctx.calo_hit_reader->read_rows(start, hit_count, ctx.event_hits);

    return !ctx.event_hits.empty();
}

void resolve_hit_references(StreamingContext &ctx, EventProducts &out)
{
    uint32_t min_hit = UINT32_MAX, max_hit = 0;
    for (const auto &hit : ctx.event_hits)
    {
        min_hit = std::min(min_hit, hit.id);
        max_hit = std::max(max_hit, hit.id);
    }
    size_t hit_reg_base = min_hit;
    size_t hit_count = (min_hit == UINT32_MAX) ? 0 : (max_hit - min_hit + 1);

    if (hit_count > 0)
    {
        ctx.hit_to_pkt_reg_reader->read_rows(hit_reg_base, hit_count, ctx.hit_pkt_regs);
        ctx.hit_to_btrk_reg_reader->read_rows(hit_reg_base, hit_count, ctx.hit_btrk_regs);
    }
    else
    {
        ctx.hit_pkt_regs.clear();
        ctx.hit_btrk_regs.clear();
    }

    auto [min_hp_ref, max_hp_ref] = min_start_max_stop(ctx.hit_pkt_regs);
    size_t hit_pkt_ref_base = 0;
    load_ref_pairs_window(ctx.hit_to_pkt_ref_reader.get(), min_hp_ref, max_hp_ref, hit_pkt_ref_base, ctx.hit_pkt_refs);

    auto [min_hb_ref, max_hb_ref] = min_start_max_stop(ctx.hit_btrk_regs);
    size_t hit_btrk_ref_base = 0;
    load_ref_pairs_window(ctx.hit_to_btrk_ref_reader.get(), min_hb_ref, max_hb_ref, hit_btrk_ref_base, ctx.hit_btrk_refs);

    uint32_t min_pkt = UINT32_MAX, max_pkt = 0;
    for (const auto &r : ctx.hit_pkt_refs)
    {
        min_pkt = std::min(min_pkt, r[1]);
        max_pkt = std::max(max_pkt, r[1]);
    }
    size_t pkt_reg_base = min_pkt;
    size_t pkt_count = (min_pkt == UINT32_MAX) ? 0 : (max_pkt - min_pkt + 1);
    if (pkt_count > 0)
        ctx.pkt_to_seg_reg_reader->read_rows(pkt_reg_base, pkt_count, ctx.pkt_seg_regs);
    else
        ctx.pkt_seg_regs.clear();

    auto [min_ps_ref, max_ps_ref] = min_start_max_stop(ctx.pkt_seg_regs);
    size_t pkt_seg_ref_base = 0;
    load_ref_pairs_window(ctx.pkt_to_seg_ref_reader.get(), min_ps_ref, max_ps_ref, pkt_seg_ref_base, ctx.pkt_seg_refs);

    // Find first valid event_id by traversing hits -> packets -> segments
    out.spill_id = -1;
    for (const auto &hit : ctx.event_hits)
    {
        if (out.spill_id >= 0)
            break;
        if (hit.id >= hit_reg_base && hit.id - hit_reg_base < ctx.hit_pkt_regs.size())
        {
            const auto &hr = ctx.hit_pkt_regs[hit.id - hit_reg_base];
            if (is_valid_region(hr) && hr.start >= hit_pkt_ref_base)
            {
                for (int hp = hr.start; hp < hr.stop && out.spill_id < 0; ++hp)
                {
                    uint32_t pkt = ctx.hit_pkt_refs[hp - hit_pkt_ref_base][1];
                    if (pkt >= pkt_reg_base && pkt - pkt_reg_base < ctx.pkt_seg_regs.size())
                    {
                        const auto &pr = ctx.pkt_seg_regs[pkt - pkt_reg_base];
                        if (is_valid_region(pr) && pr.start >= pkt_seg_ref_base)
                        {
                            for (int ps = pr.start; ps < pr.stop && out.spill_id < 0; ++ps)
                            {
                                uint32_t seg = ctx.pkt_seg_refs[ps - pkt_seg_ref_base][1];
                                std::vector<TrueSegment> s;
                                if (ctx.segment_reader->read_rows(seg, 1, s) && !s.empty())
                                {
                                    out.spill_id = s[0].event_id;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Queue up the Backtrack rows via H5Flow reference links
    ctx.hit_to_btrk_map.assign(ctx.event_hits.size(), UINT32_MAX);
    ctx.needed_frac_ids.clear();
    for (size_t i = 0; i < ctx.event_hits.size(); ++i)
    {
        uint32_t hit_id = ctx.event_hits[i].id;
        if (hit_id >= hit_reg_base && hit_id - hit_reg_base < ctx.hit_btrk_regs.size())
        {
            const auto &r = ctx.hit_btrk_regs[hit_id - hit_reg_base];
            if (is_valid_region(r) && r.start >= hit_btrk_ref_base)
            {
                uint32_t btrk_id = ctx.hit_btrk_refs[r.start - hit_btrk_ref_base][1]; // Python takes [:,0]
                ctx.hit_to_btrk_map[i] = btrk_id;
                ctx.needed_frac_ids.push_back(btrk_id);
            }
        }
    }
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

    if (debug_printed == 0)
        std::cout << "  No matches found for this event.\n";
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
    // Build segment lookup directly for current spill_id
    std::unordered_map<uint32_t, TrueSegment> segment_lookup;
    const auto seg_it = ctx.seg_rows_by_event.find(out.spill_id);
    if (seg_it != ctx.seg_rows_by_event.end())
    {
        const auto seg_spans = contiguous_spans(seg_it->second, ndlar::kCacheReadGapTolerance);
        for (const auto &span : seg_spans)
        {
            std::vector<TrueSegment> rows;
            if (ctx.segment_reader->read_rows(span[0], span[1], rows))
            {
                for (const auto &seg : rows)
                {
                    if (seg.event_id == out.spill_id)
                    {
                        segment_lookup.emplace(seg.segment_id, seg);
                    }
                }
            }
        }
    }

    // Process hits and backtrack rows
    for (size_t i = 0; i < ctx.event_hits.size(); ++i)
    {
        const CaloHit &hit = ctx.event_hits[i];
        out.hit_x[i] = hit.x;
        out.hit_y[i] = hit.y;
        out.hit_z[i] = hit.z;
        out.hit_charge[i] = hit.Q;
        out.hit_E[i] = hit.E;
        out.hit_ts[i] = hit.ts_pps;
        out.hit_io_group[i] = hit.io_group;
        out.hit_io_channel[i] = hit.io_channel;
        out.hit_chip_id[i] = hit.chip_id;
        out.hit_channel_id[i] = hit.channel_id;

        uint32_t btrk_id = ctx.hit_to_btrk_map[i];
        const PacketFraction *b_row = (btrk_id != UINT32_MAX) ? get_cached_fraction_row(ctx, btrk_id) : nullptr;
        if (!b_row)
        {
            out.hit_matches[i] = 0;
            continue;
        }

        size_t n_matches = 0;
        for (size_t k = 0; k < 20; ++k)
        {
            if (b_row->fraction[k] != 0.0)
            {
                n_matches++;
            }
        }

        out.hit_matches[i] = static_cast<uint16_t>(n_matches);
        out.hit_pdg[i].resize(n_matches, 0);
        out.hit_segmentID[i].resize(n_matches, 0);
        out.hit_particleID[i].resize(n_matches, 0);
        out.hit_particleIDLocal[i].resize(n_matches, 0);
        out.hit_vertexID[i].resize(n_matches, 0);
        out.hit_packetFrac[i].resize(n_matches, 0.0f);

        size_t m = 0;
        for (size_t k = 0; k < 20; ++k)
        {
            if (b_row->fraction[k] != 0.0)
            {
                out.hit_packetFrac[i][m] = static_cast<float>(b_row->fraction[k]);

                uint32_t seg_val = static_cast<uint32_t>(b_row->segment_ids[k]);
                auto it = segment_lookup.find(seg_val);
                if (it != segment_lookup.end())
                {
                    const TrueSegment &true_seg = it->second;
                    out.hit_pdg[i][m] = true_seg.pdg_id;
                    out.hit_segmentID[i][m] = static_cast<int32_t>(true_seg.segment_id);
                    out.hit_particleID[i][m] = static_cast<int64_t>(true_seg.file_traj_id);
                    out.hit_particleIDLocal[i][m] = static_cast<int64_t>(true_seg.traj_id);
                    out.hit_vertexID[i][m] = static_cast<int64_t>(true_seg.vertex_id);
                }
                m++;
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

        // Invalidate fast cache pointer before map modification
        ctx.last_fraction_block_ptr = nullptr;
        ctx.last_fraction_block_base = SIZE_MAX;

        ctx.fraction_blocks.emplace(block_base, rows);
    }
}

void update_caches(StreamingContext &ctx, const RawPacketFractionReader &frac_reader)
{
    const auto frac_spans = contiguous_spans(ctx.needed_frac_ids, ndlar::kCacheReadGapTolerance);
    for (const auto &span : frac_spans)
    {
        if (span[0] + span[1] > frac_reader.row_count)
            continue;
        ensure_fraction_range_cached(ctx, frac_reader, span[0], span[0] + span[1] - 1);
    }
}

void initialize_streaming_context(HighFive::File &file, StreamingContext &ctx, paths::HitType hit_type)
{
    // TODO: We should add an "is_mc" flag somewhere, and use it here and elsewhere to skip MC for data.

    const paths::PathResolver resolver(hit_type);

    ctx.calo_hit_reader = std::make_unique<RawCaloHitReader>(file, resolver.hits());
    ctx.segment_reader = std::make_unique<RawTrueSegmentReader>(file, paths::dataset::kSegments);
    ctx.seg_rows_by_event = build_event_index_from_rows(*ctx.segment_reader);

    ctx.hit_to_pkt_reg_reader = std::make_unique<RawRefRegionReader>(file, resolver.hit_to_packet_reg());
    ctx.hit_to_pkt_ref_reader = std::make_unique<RawRefPairReader>(file, resolver.hit_to_packet_ref());
    ctx.pkt_to_seg_reg_reader = std::make_unique<RawRefRegionReader>(file, paths::ref_region::kPacketToSegment);
    ctx.pkt_to_seg_ref_reader = std::make_unique<RawRefPairReader>(file, paths::ref_data::kPacketToSegment);

    ctx.hit_to_btrk_reg_reader = std::make_unique<RawRefRegionReader>(file, resolver.hit_to_backtrack_reg());
    ctx.hit_to_btrk_ref_reader = std::make_unique<RawRefPairReader>(file, resolver.hit_to_backtrack_ref());

    file.getDataSet(paths::dataset::kEvents).read(ctx.events);
    file.getDataSet(paths::dataset::kExtTrigs).read(ctx.ext_trigs);
    file.getDataSet(resolver.event_to_hits_reg()).read(ctx.hit_event_bounds);
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
        return out;

    out.reserve_hit_products(ctx.event_hits.size());

    resolve_hit_references(ctx, out);
    update_caches(ctx, frac_reader);

    populate_hit_products(ctx, out);
    populate_truth_products(ctx, out, traj_reader, int_reader);

    return out;
}

} // namespace ndlar::hdf5
