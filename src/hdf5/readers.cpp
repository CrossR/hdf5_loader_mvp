#include <algorithm>
#include <array>
#include <stdexcept>
#include <string>

#include "ndlar/hdf5/readers.hpp"

namespace ndlar::hdf5
{

template <typename Reader>
std::unordered_map<int64_t, std::vector<size_t>> build_event_index_from_rows_impl(const Reader &reader)
{
    std::unordered_map<int64_t, std::vector<size_t>> by_event;
    static constexpr size_t kChunkRows = 65536;
    std::vector<int64_t> ids;
    for (size_t base = 0; base < reader.row_count; base += kChunkRows)
    {
        const size_t n = std::min(kChunkRows, reader.row_count - base);
        if (!reader.read_event_ids(base, n, ids))
        {
            continue;
        }
        for (size_t i = 0; i < ids.size(); ++i)
        {
            by_event[ids[i]].push_back(base + i);
        }
    }
    return by_event;
}

std::vector<std::array<size_t, 2>> contiguous_spans(const std::vector<size_t> &indices, size_t max_gap)
{
    std::vector<std::array<size_t, 2>> spans;
    if (indices.empty())
    {
        return spans;
    }

    size_t span_start = indices[0];
    size_t prev = indices[0];
    for (size_t i = 1; i < indices.size(); ++i)
    {
        const size_t cur = indices[i];
        if (cur <= prev + max_gap + 1)
        {
            prev = cur;
            continue;
        }

        spans.push_back({span_start, prev - span_start + 1});
        span_start = cur;
        prev = cur;
    }
    spans.push_back({span_start, prev - span_start + 1});

    return spans;
}

std::unordered_map<int64_t, std::vector<size_t>> build_event_index_from_rows(const RawTrajectoryReader &reader)
{
    return build_event_index_from_rows_impl(reader);
}

std::unordered_map<int64_t, std::vector<size_t>> build_event_index_from_rows(const RawInteractionReader &reader)
{
    return build_event_index_from_rows_impl(reader);
}

} // namespace ndlar::hdf5
