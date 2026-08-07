#pragma once

#include <cstddef>
#include <vector>

#include <highfive/H5DataSet.hpp>
#include <highfive/H5File.hpp>

#include "ndlar/hdf5/highfive_types.hpp"
#include "ndlar/hdf5/types.hpp"

namespace ndlar::hdf5
{

// Basic Table Reader
//
// Reads the given datatype into a vector, using the HighFive library. This is a
// generic reader that can be used for any dataset, but does not have any
// special optimizations for event_id filtering.
template <typename T>
class TableReader
{
public:
    TableReader(HighFive::File &file, const std::string &dataset_path) :
        dset_(file.getDataSet(dataset_path)),
        row_count(dset_.getDimensions()[0])
    {
    }

    bool read_rows(size_t first_idx, size_t count, std::vector<T> &out) const
    {
        if (count == 0 || first_idx >= row_count || first_idx + count > row_count)
        {
            out.clear();
            return false;
        }

        dset_.select({first_idx}, {count}).read(out);
        return true;
    }

protected:
    HighFive::DataSet dset_;

public:
    size_t row_count = 0;
};

// Event Table Reader - Specialised version of the TableReader, that utilises
// the EventID to ensure we only read the relevant event, rather than the full
// dataset. This is a performance optimisation for large datasets with many
// events.
template <typename T>
class EventTableReader : public TableReader<T>
{
public:
    using TableReader<T>::TableReader;

    bool read_event_ids(size_t first_idx, size_t count, std::vector<int64_t> &out) const
    {
        if (count == 0 || first_idx >= this->row_count || first_idx + count > this->row_count)
        {
            out.clear();
            return false;
        }

        // Read only the event_id column using the dummy struct
        std::vector<EventIdOnly> temp;
        this->dset_.select({first_idx}, {count}).read(temp);

        // Unpack into the flat integer vector, to be used by the caller
        out.resize(temp.size());
        for (size_t i = 0; i < temp.size(); ++i)
        {
            out[i] = temp[i].event_id;
        }

        return true;
    }
};

// Specialized Pair Reader for RefPair datasets.
//
// Annoyingly, this can't be simply done with TableReader<> because the Pairs
// can be [N], or [N,2] depending on the dataset. So we need to handle both
// cases, which means a bit more wrapping.
class RawRefPairReader
{
public:
    RawRefPairReader(HighFive::File &file, const std::string &dataset_path) :
        dset_(file.getDataSet(dataset_path))
    {
        auto dims = dset_.getDimensions();
        row_count = dims[0];
        is_2d = (dims.size() == 2 && dims[1] == 2);
    }

    bool read_rows(size_t first_idx, size_t count, std::vector<RefPair> &out) const
    {
        if (count == 0 || first_idx >= row_count || first_idx + count > row_count)
        {
            out.clear();
            return false;
        }

        if (is_2d)
            dset_.select({first_idx, 0}, {count, 2}).read(out);
        else
            dset_.select({first_idx}, {count}).read(out);
        return true;
    }

    size_t row_count = 0;

private:
    HighFive::DataSet dset_;
    bool is_2d = false;
};

// Various utility functions for working with the HDF5 datasets and references
template <typename Reader>
std::unordered_map<int64_t, std::vector<size_t>> build_event_index_from_rows(const Reader &reader)
{
    std::unordered_map<int64_t, std::vector<size_t>> index;
    std::vector<int64_t> event_ids;

    if (reader.read_event_ids(0, reader.row_count, event_ids))
    {
        for (size_t i = 0; i < event_ids.size(); ++i)
        {
            index[event_ids[i]].push_back(i);
        }
    }
    return index;
}

/*
 * Returns a list of contiguous spans from a list of indices.
 * Each span is represented as a pair of (start_index, length).
 *
 * @param indices A sorted vector of indices.
 * @param max_gap The maximum allowed gap between indices to consider them contiguous.
 * @return A vector of pairs representing the contiguous spans.
 */
inline std::vector<std::array<size_t, 2>> contiguous_spans(const std::vector<size_t> &indices, size_t max_gap)
{
    std::vector<std::array<size_t, 2>> spans;
    if (indices.empty())
        return spans;

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

// Finally, define all the specific readers for each dataset type, using the appropriate struct
//
// First, basic readers for the datasets that don't require event_id filtering
using RawPromptHitReader = TableReader<PromptHit>;
using RawPacketFractionReader = TableReader<PacketFraction>;
using RawRefRegionReader = TableReader<RefRegion>;

// Smart Readers
using RawTrueSegmentReader = EventTableReader<TrueSegment>;
using RawInteractionReader = EventTableReader<Interaction>;
using RawTrajectoryReader = EventTableReader<Trajectory>;

} // namespace ndlar::hdf5
