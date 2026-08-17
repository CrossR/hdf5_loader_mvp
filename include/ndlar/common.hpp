#pragma once

#include <climits>
#include <cstdint>

namespace ndlar
{
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

} // namespace ndlar
