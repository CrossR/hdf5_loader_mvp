#pragma once

#include <chrono>
#include <climits>
#include <cstdint>

namespace ndlar {

using SteadyClock = std::chrono::steady_clock;

static constexpr int kDebugMatchPrintLimit = 3;
static constexpr int32_t kInvalidTrigger = INT32_MAX;
static constexpr float kMeVToGeV = 1.0e-3f;
static constexpr size_t kCacheReadGapTolerance = 256;
static constexpr size_t kFractionBlockRows = 256;

inline double elapsed_ms(const SteadyClock::time_point& start, const SteadyClock::time_point& end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

}  // namespace ndlar
