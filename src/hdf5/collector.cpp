#include "ndlar/hdf5/collector.hpp"

#include <algorithm>
#include <climits>
#include <iostream>

#include "ndlar/common.hpp"
#include "ndlar/hdf5/highfive_types.hpp"

namespace ndlar::hdf5 {

void print_debug_matches(const EventProducts& event_products) {
    int debug_printed = 0;
    for (size_t hit_index = 0; hit_index < event_products.hit_pdg.size() && debug_printed < ndlar::kDebugMatchPrintLimit; ++hit_index) {
        if (event_products.hit_segmentID[hit_index] == 0) {
            continue;
        }

        std::cout << "  Hit " << hit_index
                  << " -> PDG: " << event_products.hit_pdg[hit_index]
                  << ", frac: " << event_products.hit_packetFrac[hit_index] << "\n";
        ++debug_printed;
    }
}

size_t fraction_block_base(size_t row_id) {
    return (row_id / ndlar::kFractionBlockRows) * ndlar::kFractionBlockRows;
}

const PacketFraction* get_cached_fraction_row(const StreamingContext& ctx, size_t row_id) {
    const size_t block_base = fraction_block_base(row_id);
    const auto block_it = ctx.fraction_blocks.find(block_base);
    if (block_it == ctx.fraction_blocks.end()) {
        return nullptr;
    }

    const size_t offset = row_id - block_base;
    if (offset >= block_it->second.size()) {
        return nullptr;
    }
    return &block_it->second[offset];
}

void ensure_fraction_range_cached(
    StreamingContext& ctx,
    const RawPacketFractionReader& frac_reader,
    size_t first_id,
    size_t last_id_inclusive) {
    const size_t first_block = fraction_block_base(first_id);
    const size_t last_block = fraction_block_base(last_id_inclusive);

    for (size_t block_base = first_block; block_base <= last_block; block_base += ndlar::kFractionBlockRows) {
        if (ctx.fraction_blocks.find(block_base) != ctx.fraction_blocks.end()) {
            continue;
        }

        const size_t available = frac_reader.row_count - block_base;
        const size_t count = std::min(ndlar::kFractionBlockRows, available);
        auto& rows = ctx.fraction_rows;
        if (!frac_reader.read_rows(block_base, count, rows)) {
            continue;
        }
        ctx.fraction_blocks.emplace(block_base, rows);
    }
}

void initialize_streaming_context(HighFive::File& file, StreamingContext& ctx) {
    ctx.dset_hits = file.getDataSet("charge/calib_prompt_hits/data");
    ctx.hit_to_pkt_reg_reader = std::make_unique<RawRefRegionReader>(
        file.getId(), "charge/calib_prompt_hits/ref/charge/packets/ref_region");
    ctx.hit_to_pkt_ref_reader = std::make_unique<RawRefPairReader>(
        file.getId(), "charge/calib_prompt_hits/ref/charge/packets/ref");
    ctx.pkt_to_seg_reg_reader = std::make_unique<RawRefRegionReader>(
        file.getId(), "charge/packets/ref/mc_truth/segments/ref_region");
    ctx.pkt_to_seg_ref_reader = std::make_unique<RawRefPairReader>(
        file.getId(), "charge/packets/ref/mc_truth/segments/ref");
    ctx.pkt_to_frac_reg_reader = std::make_unique<RawRefRegionReader>(
        file.getId(), "charge/packets/ref/mc_truth/packet_fraction/ref_region");
    ctx.pkt_to_frac_ref_reader = std::make_unique<RawRefPairReader>(
        file.getId(), "charge/packets/ref/mc_truth/packet_fraction/ref");
    ctx.segment_reader = std::make_unique<RawTrueSegmentReader>(file.getId(), "mc_truth/segments/data");

    file.getDataSet("charge/events/data").read(ctx.events);
    file.getDataSet("charge/ext_trigs/data").read(ctx.ext_trigs);
    file.getDataSet("charge/events/ref/charge/calib_prompt_hits/ref_region").read(ctx.hit_event_bounds);
    file.getDataSet("charge/events/ref/charge/ext_trigs/ref_region").read(ctx.event_to_exttrig_reg);
    file.getDataSet("charge/events/ref/charge/ext_trigs/ref").read(ctx.event_to_exttrig_ref);

    ctx.segment_cache.resize(ctx.segment_reader->row_count);
    ctx.segment_cache_valid.assign(ctx.segment_reader->row_count, 0);
}

int32_t select_trigger_id_stream(const StreamingContext& ctx, size_t event_index) {
    if (event_index >= ctx.event_to_exttrig_reg.size()) {
        return ndlar::kInvalidTrigger;
    }

    const RefRegion trig_bounds = ctx.event_to_exttrig_reg[event_index];
    if (!is_valid_region(trig_bounds)) {
        return ndlar::kInvalidTrigger;
    }

    const size_t start = static_cast<size_t>(trig_bounds.start);
    const size_t stop = static_cast<size_t>(trig_bounds.stop);
    if (start >= ctx.event_to_exttrig_ref.size()) {
        return ndlar::kInvalidTrigger;
    }

    const size_t clamped_stop = std::min(stop, ctx.event_to_exttrig_ref.size());
    int32_t first = ndlar::kInvalidTrigger;
    bool any = false;
    for (size_t idx = start; idx < clamped_stop; ++idx) {
        const uint32_t ext_idx = ctx.event_to_exttrig_ref[idx][1];
        if (ext_idx >= ctx.ext_trigs.size()) {
            continue;
        }
        const int32_t iogroup = ctx.ext_trigs[ext_idx].iogroup;
        if (first == ndlar::kInvalidTrigger) {
            first = iogroup;
        }
        any = true;
        if (iogroup == 5) {
            return 5;
        }
    }

    return any ? first : ndlar::kInvalidTrigger;
}

EventProducts collect_event_products_stream(
    StreamingContext& ctx,
    size_t event_index,
    const RawPacketFractionReader& frac_reader,
    const RawTrajectoryReader& traj_reader,
    const RawInteractionReader& int_reader) {
    EventProducts out;
    out.trigger_id = select_trigger_id_stream(ctx, event_index);

    if (event_index >= ctx.hit_event_bounds.size()) {
        return out;
    }

    const RefRegion event_bounds = ctx.hit_event_bounds[event_index];
    if (!is_valid_region(event_bounds)) {
        return out;
    }

    const size_t start = static_cast<size_t>(event_bounds.start);
    const size_t stop = static_cast<size_t>(event_bounds.stop);
    const size_t hit_count = stop - start;

    out.reserve_hit_products(hit_count);

    auto& event_hits = ctx.event_hits;
    event_hits.clear();
    ctx.dset_hits.select({start}, {hit_count}).read(event_hits);

    if (event_hits.empty()) {
        return out;
    }

    uint32_t min_hit_id = UINT32_MAX;
    uint32_t max_hit_id = 0;
    for (const PromptHit& hit : event_hits) {
        if (ctx.hit_to_pkt_reg_reader != nullptr && hit.id < ctx.hit_to_pkt_reg_reader->row_count) {
            min_hit_id = std::min(min_hit_id, hit.id);
            max_hit_id = std::max(max_hit_id, hit.id);
        }
    }

    auto& hit_pkt_regs = ctx.hit_pkt_regs;
    hit_pkt_regs.clear();
    size_t hit_reg_base = 0;
    if (min_hit_id != UINT32_MAX) {
        hit_reg_base = static_cast<size_t>(min_hit_id);
        const size_t n = static_cast<size_t>(max_hit_id - min_hit_id + 1);
        ctx.hit_to_pkt_reg_reader->read_rows(hit_reg_base, n, hit_pkt_regs);
    }

    int32_t min_hit_ref = INT32_MAX;
    int32_t max_hit_ref_stop = INT32_MIN;
    for (const RefRegion& r : hit_pkt_regs) {
        if (!is_valid_region(r)) {
            continue;
        }
        min_hit_ref = std::min(min_hit_ref, r.start);
        max_hit_ref_stop = std::max(max_hit_ref_stop, r.stop);
    }

    auto& hit_pkt_refs = ctx.hit_pkt_refs;
    hit_pkt_refs.clear();
    size_t hit_ref_base = 0;
    if (min_hit_ref != INT32_MAX && max_hit_ref_stop > min_hit_ref) {
        hit_ref_base = static_cast<size_t>(min_hit_ref);
        const size_t n = static_cast<size_t>(max_hit_ref_stop - min_hit_ref);
        ctx.hit_to_pkt_ref_reader->read_rows(hit_ref_base, n, hit_pkt_refs);
    }

    auto& pkt_ids = ctx.pkt_ids;
    pkt_ids.assign(event_hits.size(), UINT32_MAX);
    uint32_t min_pkt_id = UINT32_MAX;
    uint32_t max_pkt_id = 0;
    for (size_t i = 0; i < event_hits.size(); ++i) {
        const uint32_t hit_id = event_hits[i].id;
        if (hit_id < hit_reg_base || static_cast<size_t>(hit_id - hit_reg_base) >= hit_pkt_regs.size()) {
            continue;
        }

        const RefRegion r = hit_pkt_regs[static_cast<size_t>(hit_id - hit_reg_base)];
        if (!is_valid_region(r) || r.start < static_cast<int32_t>(hit_ref_base)) {
            continue;
        }

        const size_t ref_idx = static_cast<size_t>(r.start - static_cast<int32_t>(hit_ref_base));
        if (ref_idx >= hit_pkt_refs.size()) {
            continue;
        }

        const uint32_t pkt_id = hit_pkt_refs[ref_idx][1];
        if (pkt_id >= ctx.pkt_to_seg_reg_reader->row_count || pkt_id >= ctx.pkt_to_frac_reg_reader->row_count) {
            continue;
        }

        pkt_ids[i] = pkt_id;
        min_pkt_id = std::min(min_pkt_id, pkt_id);
        max_pkt_id = std::max(max_pkt_id, pkt_id);
    }

    auto& pkt_seg_regs = ctx.pkt_seg_regs;
    auto& pkt_frac_regs = ctx.pkt_frac_regs;
    pkt_seg_regs.clear();
    pkt_frac_regs.clear();
    size_t pkt_reg_base = 0;
    if (min_pkt_id != UINT32_MAX) {
        pkt_reg_base = static_cast<size_t>(min_pkt_id);
        const size_t n = static_cast<size_t>(max_pkt_id - min_pkt_id + 1);
        ctx.pkt_to_seg_reg_reader->read_rows(pkt_reg_base, n, pkt_seg_regs);
        ctx.pkt_to_frac_reg_reader->read_rows(pkt_reg_base, n, pkt_frac_regs);
    }

    int32_t min_seg_ref = INT32_MAX;
    int32_t max_seg_ref_stop = INT32_MIN;
    int32_t min_frac_ref = INT32_MAX;
    int32_t max_frac_ref_stop = INT32_MIN;
    for (const RefRegion& r : pkt_seg_regs) {
        if (!is_valid_region(r)) continue;
        min_seg_ref = std::min(min_seg_ref, r.start);
        max_seg_ref_stop = std::max(max_seg_ref_stop, r.stop);
    }
    for (const RefRegion& r : pkt_frac_regs) {
        if (!is_valid_region(r)) continue;
        min_frac_ref = std::min(min_frac_ref, r.start);
        max_frac_ref_stop = std::max(max_frac_ref_stop, r.stop);
    }

    auto& pkt_seg_refs = ctx.pkt_seg_refs;
    auto& pkt_frac_refs = ctx.pkt_frac_refs;
    pkt_seg_refs.clear();
    pkt_frac_refs.clear();
    size_t seg_ref_base = 0;
    size_t frac_ref_base = 0;
    if (min_seg_ref != INT32_MAX && max_seg_ref_stop > min_seg_ref) {
        seg_ref_base = static_cast<size_t>(min_seg_ref);
        const size_t n = static_cast<size_t>(max_seg_ref_stop - min_seg_ref);
        ctx.pkt_to_seg_ref_reader->read_rows(seg_ref_base, n, pkt_seg_refs);
    }
    if (min_frac_ref != INT32_MAX && max_frac_ref_stop > min_frac_ref) {
        frac_ref_base = static_cast<size_t>(min_frac_ref);
        const size_t n = static_cast<size_t>(max_frac_ref_stop - min_frac_ref);
        ctx.pkt_to_frac_ref_reader->read_rows(frac_ref_base, n, pkt_frac_refs);
    }

    auto& seg_ids = ctx.seg_ids;
    auto& frac_ids = ctx.frac_ids;
    auto& match_counts = ctx.match_counts;
    auto& needed_seg_ids = ctx.needed_seg_ids;
    auto& needed_frac_ids = ctx.needed_frac_ids;
    seg_ids.assign(event_hits.size(), UINT32_MAX);
    frac_ids.assign(event_hits.size(), UINT32_MAX);
    match_counts.assign(event_hits.size(), 0);
    needed_seg_ids.clear();
    needed_frac_ids.clear();
    needed_seg_ids.reserve(event_hits.size());
    needed_frac_ids.reserve(event_hits.size());

    for (size_t i = 0; i < event_hits.size(); ++i) {
        const uint32_t pkt_id = pkt_ids[i];
        if (pkt_id == UINT32_MAX || pkt_id < pkt_reg_base) {
            continue;
        }

        const size_t pkt_off = static_cast<size_t>(pkt_id - pkt_reg_base);
        if (pkt_off >= pkt_seg_regs.size() || pkt_off >= pkt_frac_regs.size()) {
            continue;
        }

        const RefRegion seg_r = pkt_seg_regs[pkt_off];
        if (is_valid_region(seg_r)) {
            match_counts[i] = static_cast<uint16_t>(region_size(seg_r));
            if (seg_r.start >= static_cast<int32_t>(seg_ref_base)) {
                const size_t ref_idx = static_cast<size_t>(seg_r.start - static_cast<int32_t>(seg_ref_base));
                if (ref_idx < pkt_seg_refs.size()) {
                    seg_ids[i] = pkt_seg_refs[ref_idx][1];
                    needed_seg_ids.push_back(static_cast<size_t>(seg_ids[i]));
                }
            }
        }

        const RefRegion frac_r = pkt_frac_regs[pkt_off];
        if (is_valid_region(frac_r) && frac_r.start >= static_cast<int32_t>(frac_ref_base)) {
            const size_t ref_idx = static_cast<size_t>(frac_r.start - static_cast<int32_t>(frac_ref_base));
            if (ref_idx < pkt_frac_refs.size()) {
                frac_ids[i] = pkt_frac_refs[ref_idx][1];
                needed_frac_ids.push_back(static_cast<size_t>(frac_ids[i]));
            }
        }
    }

    std::sort(needed_seg_ids.begin(), needed_seg_ids.end());
    needed_seg_ids.erase(std::unique(needed_seg_ids.begin(), needed_seg_ids.end()), needed_seg_ids.end());
    std::sort(needed_frac_ids.begin(), needed_frac_ids.end());
    needed_frac_ids.erase(std::unique(needed_frac_ids.begin(), needed_frac_ids.end()), needed_frac_ids.end());

    const auto seg_spans = contiguous_spans(needed_seg_ids, ndlar::kCacheReadGapTolerance);
    for (const auto& span : seg_spans) {
        if (span[0] + span[1] > ctx.segment_reader->row_count) {
            continue;
        }

        bool all_cached = true;
        for (size_t id = span[0]; id < span[0] + span[1]; ++id) {
            if (!ctx.segment_cache_valid[id]) {
                all_cached = false;
                break;
            }
        }
        if (all_cached) {
            continue;
        }

        std::vector<TrueSegment> rows;
        if (!ctx.segment_reader->read_rows(span[0], span[1], rows)) {
            continue;
        }
        for (size_t i = 0; i < rows.size(); ++i) {
            const size_t id = span[0] + i;
            ctx.segment_cache[id] = rows[i];
            ctx.segment_cache_valid[id] = 1;
        }
    }

    const auto frac_spans = contiguous_spans(needed_frac_ids, ndlar::kCacheReadGapTolerance);
    for (const auto& span : frac_spans) {
        if (span[0] + span[1] > frac_reader.row_count) {
            continue;
        }

        ensure_fraction_range_cached(ctx, frac_reader, span[0], span[0] + span[1] - 1);
    }

    for (size_t i = 0; i < event_hits.size(); ++i) {
        const PromptHit& hit = event_hits[i];
        out.hit_x.push_back(hit.x);
        out.hit_y.push_back(hit.y);
        out.hit_z.push_back(hit.z);
        out.hit_charge.push_back(hit.Q);
        out.hit_E.push_back(hit.E);
        out.hit_ts.push_back(hit.ts_pps);

        uint16_t matches = 0;
        float packet_fraction = 0.0f;
        int32_t pdg = 0;
        int32_t segment_id = 0;
        int64_t file_traj_id = 0;
        int64_t traj_id = 0;
        int64_t vertex_id = 0;

        matches = match_counts[i];
        const uint32_t seg_id = seg_ids[i];
        if (seg_id != UINT32_MAX) {
            const size_t seg_index = static_cast<size_t>(seg_id);
            if (seg_index < ctx.segment_cache_valid.size() && ctx.segment_cache_valid[seg_index]) {
                const TrueSegment& true_seg = ctx.segment_cache[seg_index];
                pdg = true_seg.pdg_id;
                segment_id = static_cast<int32_t>(true_seg.segment_id);
                file_traj_id = static_cast<int64_t>(true_seg.file_traj_id);
                traj_id = static_cast<int64_t>(true_seg.traj_id);
                vertex_id = static_cast<int64_t>(true_seg.vertex_id);

                if (out.spill_id < 0) {
                    out.spill_id = true_seg.event_id;
                }

                const uint32_t frac_id = frac_ids[i];
                if (frac_id != UINT32_MAX) {
                    const PacketFraction* frac_row = get_cached_fraction_row(ctx, static_cast<size_t>(frac_id));
                    if (frac_row != nullptr) {
                        packet_fraction = resolve_packet_fraction(*frac_row, true_seg.segment_id);
                    }
                }
            }
        }

        out.hit_matches.push_back(matches);
        out.hit_packetFrac.push_back(packet_fraction);
        out.hit_pdg.push_back(pdg);
        out.hit_segmentID.push_back(segment_id);
        out.hit_particleID.push_back(file_traj_id);
        out.hit_particleIDLocal.push_back(traj_id);
        out.hit_vertexID.push_back(vertex_id);
    }

    if (out.spill_id >= 0) {
        const auto traj_it = ctx.traj_rows_by_event.find(out.spill_id);
        if (traj_it != ctx.traj_rows_by_event.end()) {
            out.reserve_trajectory_products(traj_it->second.size());
            const auto spans = contiguous_spans(traj_it->second, ndlar::kCacheReadGapTolerance);
            for (const auto& span : spans) {
                auto& rows = ctx.trajectory_rows;
                if (!traj_reader.read_rows(span[0], span[1], rows)) {
                    continue;
                }

                append_trajectory_products(rows, out);
            }
        }

        const auto int_it = ctx.int_rows_by_event.find(out.spill_id);
        if (int_it != ctx.int_rows_by_event.end()) {
            out.reserve_interaction_products(int_it->second.size());
            const auto spans = contiguous_spans(int_it->second, ndlar::kCacheReadGapTolerance);
            for (const auto& span : spans) {
                auto& rows = ctx.interaction_rows;
                if (!int_reader.read_rows(span[0], span[1], rows)) {
                    continue;
                }

                append_interaction_products(rows, out);
            }
        }
    }

    return out;
}

}  // namespace ndlar::hdf5
