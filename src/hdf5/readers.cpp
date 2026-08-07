
#include "ndlar/hdf5/readers.hpp"

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
