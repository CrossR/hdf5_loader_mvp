#pragma once

#include <chrono>
#include <climits>
#include <cstdint>

namespace ndlar
{

using SteadyClock = std::chrono::steady_clock;

// Number of hits to print out debug information for when printing out event products.
static constexpr int kDebugMatchPrintLimit = 3;
// Invalid trigger ID value to use when an event has no trigger ID.
static constexpr int32_t kInvalidTrigger = INT32_MAX;
// Conversion factor from MeV to GeV.
static constexpr float kMeVToGeV = 1.0e-3f;
// Tolerance for cache read gaps.
static constexpr size_t kCacheReadGapTolerance = 256;
// Number of rows in each fraction block.
static constexpr size_t kFractionBlockRows = 256;

/**
 * Returns the elapsed time in milliseconds between two time points.
 *
 * @param start The starting time point.
 * @param end The ending time point.
 * @return The elapsed time in milliseconds as a double.
 */
inline double elapsed_ms(const SteadyClock::time_point &start, const SteadyClock::time_point &end)
{
    return std::chrono::duration<double, std::milli>(end - start).count();
}

} // namespace ndlar
