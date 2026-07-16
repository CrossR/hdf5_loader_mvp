#include <iostream>
#include <vector>
#include <highfive/H5File.hpp>
#include "hep_evd.h"

// Hard code the NDLAr Geometry here, for use in HepEVD.
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

struct PromptHit {
    float x, y, z, E;
};

// Link the PromptHit struct to a HighFive compound type so we can read it from
// HDF5
HighFive::CompoundType create_compound_PromptHit() {
    return {
        {"x", HighFive::AtomicType<float>{}},
        {"y", HighFive::AtomicType<float>{}},
        {"z", HighFive::AtomicType<float>{}},
        {"E", HighFive::AtomicType<float>{}}
    };
}
HIGHFIVE_REGISTER_TYPE(PromptHit, create_compound_PromptHit)

// Link the RefRegion struct to a HighFive compound type so we can read it from
// HDF5
struct RefRegion {
    uint64_t start;
    uint64_t stop;
};

HighFive::CompoundType create_compound_RefRegion() {
    return {
        {"start", HighFive::AtomicType<uint64_t>{}},
        {"stop", HighFive::AtomicType<uint64_t>{}}
    };
}
HIGHFIVE_REGISTER_TYPE(RefRegion, create_compound_RefRegion)

// HepEVD Geometry Helper Functions
HepEVD::DetectorGeometry get_ndlar_geometry() {
    HepEVD::Volumes volumes;
    for (const auto& x : NDLAR_X_VALUES) {
        for (const auto& z : NDLAR_Z_VALUES) {
            HepEVD::BoxVolume larTpc({x, NDLAR_Y_CONST, z}, NDLAR_X_WIDTH, NDLAR_Y_WIDTH, NDLAR_Z_WIDTH);
            volumes.push_back(larTpc);
        }
    }

    return HepEVD::DetectorGeometry(volumes);
}

// Main Event Processing Loop
int main(int argc, char** argv) {

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <h5flow_file.hdf5>" << std::endl;
        return 1;
    }

    const auto geometry = get_ndlar_geometry();
    const auto server = new HepEVD::HepEVDServer(geometry);

    try {
        HighFive::File file(argv[1], HighFive::File::ReadOnly);

        // Open the hits dataset
        auto dset_hits = file.getDataSet("charge/calib_prompt_hits/data");

        // Open the h5flow relational bounds dataset
        auto dset_regions = file.getDataSet("charge/events/ref/charge/calib_prompt_hits/ref_region");

        // Read ALL the event bounding boxes into memory at once
        std::vector<RefRegion> regions;
        dset_regions.read(regions);

        size_t num_events = regions.size();
        std::cout << "File contains " << num_events << " events." << std::endl;
        std::cout << "-------------------------------------------" << std::endl;

        if (num_events == 0) {
            std::cerr << "No events found in the file." << std::endl;
            return 1;
        }

        // Per event loop
        for (size_t i = 0; i < num_events; ++i) {

            uint64_t start = regions[i].start;
            uint64_t stop = regions[i].stop;
            uint64_t count = stop - start;

            std::cout << "Event " << i << " contains " << count << " hits." << std::endl;

            // Skip empty events.
            if (count == 0) {
                continue;
            }

            // Slice the hits dataset dynamically for this specific event
            std::vector<PromptHit> event_hits;
            dset_hits.select({start}, {count}).read(event_hits);

            // Propagate over to HepEVD.
            HepEVD::Hits hits;

            for (const auto& hit : event_hits) {
                HepEVD::Hit* evd_hit = new HepEVD::Hit({hit.x, hit.y, hit.z}, hit.E);
                hits.push_back(evd_hit);
            }

            server->addHits(hits);
            server->startServer();
            server->resetServer();
        }

    } catch (const HighFive::Exception& err) {
        std::cerr << "HDF5 Error: " << err.what() << "\n";
        return 1;
    }

    return 0;
}
