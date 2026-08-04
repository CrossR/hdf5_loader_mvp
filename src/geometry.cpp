#include "ndlar/geometry.hpp"

#include <vector>

namespace {

static const std::vector<float> NDLAR_X_VALUES = {
    -323.7130126953125,
    -276.2869873046875,
    -223.71299743652344,
    -176.28700256347656,
    -123.71299743652344,
    -76.28700256347656,
    -23.71299934387207,
    23.71299934387207,
    76.28700256347656,
    123.71299743652344,
    176.28700256347656,
    223.71299743652344,
    276.2869873046875,
    323.7130126953125,
};

static const std::vector<float> NDLAR_Z_VALUES = {
    465.7558898925781,
    565.7559204101563,
    665.7559204101563,
    765.7559204101563,
    865.7559204101563,
};

static const float NDLAR_Y_CONST = -66.87129974365234;
static const float NDLAR_X_WIDTH = 46.79100036621094;
static const float NDLAR_Y_WIDTH = 299.5989990234375;
static const float NDLAR_Z_WIDTH = 95.6635971069336;

}  // namespace

namespace ndlar {

HepEVD::DetectorGeometry get_ndlar_geometry() {
    HepEVD::Volumes volumes;
    volumes.reserve(NDLAR_X_VALUES.size() * NDLAR_Z_VALUES.size());
    for (const auto& x : NDLAR_X_VALUES) {
        for (const auto& z : NDLAR_Z_VALUES) {
            HepEVD::BoxVolume lar_tpc({x, NDLAR_Y_CONST, z}, NDLAR_X_WIDTH, NDLAR_Y_WIDTH, NDLAR_Z_WIDTH);
            volumes.push_back(lar_tpc);
        }
    }
    return HepEVD::DetectorGeometry(volumes);
}

}  // namespace ndlar
