
#include <iostream>
#include <memory>
#include <stdexcept>

#include "hep_evd.h"
#include <highfive/H5File.hpp>

#include "ndlar/common.hpp"
#include "ndlar/geometry.hpp"
#include "ndlar/hdf5/collector.hpp"
#include "ndlar/hdf5/paths.hpp"
#include "ndlar/hdf5/readers.hpp"
#include "ndlar/reader_app.hpp"

#if defined(ROOT_FOUND)
#include "ndlar/to_root.hpp"
#endif

// Include the LZF registration header from the submodule
extern "C"
{
    int register_lzf(void);
}

struct CLIArgs
{
    std::string h5flow_file = "";
    ndlar::hdf5::paths::HitType hit_type = ndlar::hdf5::paths::HitType::Prompt;
    bool has_mc = true;
};

void print_usage(const char *prog_name)
{
    std::cerr << "Usage: " << prog_name << " <h5flow_file.hdf5> [hit_type]\n";
    std::cerr << "  <h5flow_file.hdf5> : Path to the HDF5 file containing the event data.\n";
    std::cerr << "  --hit-type         : Optional. Specify the hit type to read. Valid options are:\n";
    std::cerr << "                       prompt (default), merged, final\n";
    std::cerr << "  --data             : Optional. Specify if the data is real (not MC). Default is MC.\n";
}

CLIArgs parse_cli_args(int argc, char **argv)
{
    CLIArgs args;

    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help")
        {
            print_usage(argv[0]);
            break;
        }
        else if (arg[0] != '-' && args.h5flow_file.empty())
        {
            args.h5flow_file = arg;
        }
        else if (arg == "--hit-type" && i + 1 < argc)
        {
            std::string hit_type_arg = argv[++i];

            if (hit_type_arg == "prompt")
                args.hit_type = ndlar::hdf5::paths::HitType::Prompt;
            else if (hit_type_arg == "merged")
                args.hit_type = ndlar::hdf5::paths::HitType::Merged;
            else if (hit_type_arg == "final")
                args.hit_type = ndlar::hdf5::paths::HitType::Final;
            else
                throw std::invalid_argument("Invalid hit type argument: " + hit_type_arg);

        }
        else if (arg == "--data")
        {
            args.has_mc = false;
        }
        else
        {
            throw std::invalid_argument("Unexpected argument: " + arg);
        }
    }

    if (args.h5flow_file.empty())
        throw std::runtime_error("Missing required argument: <h5flow_file.hdf5>");

    return args;
}

int ndlar::run_reader_app(int argc, char **argv)
{

    // Register the LZF filter first, before anything else.
    // Needed for 2x2 data and likely other bits of real data.
    if (register_lzf() < 0)
    {
        std::cerr << "Failed to register HDF5 LZF filter!" << std::endl;
        return 1;
    }

    CLIArgs args{parse_cli_args(argc, argv)};

    const auto file_path = args.h5flow_file;
    const auto hit_type = args.hit_type;
    const auto has_mc = args.has_mc;

    const auto geometry = ndlar::get_ndlar_geometry();
    auto server = std::make_unique<HepEVD::HepEVDServer>(geometry);

#if defined(ROOT_FOUND)
    std::cout << "ROOT found. Writing event products to ROOT file 'event_products.root'.\n";
    auto root_writer = std::make_unique<ndlar::RootWriter>("event_products.root");
#endif

    try
    {
        const auto t0_total = ndlar::SteadyClock::now();
        HighFive::File file(file_path, HighFive::File::ReadOnly);

        const auto t0_meta = ndlar::SteadyClock::now();
        ndlar::hdf5::StreamingContext ctx;
        ndlar::hdf5::initialize_streaming_context(file, ctx, hit_type, has_mc);

        if (!ctx.is_setup())
        {
            std::cerr << "Warning: File is missing required datasets for the requested configuration. Skipping..." << std::endl;
            return 1;
        }

        const auto t1_meta = ndlar::SteadyClock::now();

        const auto t0_index = ndlar::SteadyClock::now();
        const ndlar::hdf5::paths::PathResolver resolver(hit_type);
        const auto t1_index = ndlar::SteadyClock::now();

        const size_t num_events = ctx.events.size();
        std::cout << "File contains " << num_events << " events.\n";
        std::cout << "-------------------------------------------\n";
        std::cout << "Timing: metadata_load_ms=" << ndlar::elapsed_ms(t0_meta, t1_meta)
                  << ", index_build_ms=" << ndlar::elapsed_ms(t0_index, t1_index) << "\n";

        double total_collect_ms = 0.0;

        for (size_t i = 0; i < num_events; ++i)
        {
            // Grab all the products for the current event.
            const auto t0_collect = ndlar::SteadyClock::now();
            auto ev = ndlar::hdf5::collect_event_products_stream(ctx, i);
            const auto t1_collect = ndlar::SteadyClock::now();

#if defined(ROOT_FOUND)
            // Fill the ROOT tree with the collected event products.
            root_writer->Fill(ev, 0, 0, i);
#endif

            // Print how long it took to gather all the info for the current event.
            const double collect_ms = ndlar::elapsed_ms(t0_collect, t1_collect);
            total_collect_ms += collect_ms;
            std::cout << "  timing_ms: collect=" << collect_ms << "\n";

            // Print some debug info about the event products we just collected.
            ndlar::hdf5::print_debug_matches(ev);

            std::cout << "Event " << i << ": trigger=" << ev.trigger_id << ", hits=" << ev.hit_x.size() << ", matched_values=" << ev.hit_pdg.size()
                      << ", spill=" << ev.spill_id << ", traj=" << ev.mcp_id.size() << ", nu_vtx=" << ev.nuID.size() << "\n";

            // INFO: Since this is just an MVP...just run HepEVD, not full Pandora reco.
            HepEVD::Hits evd_hits;
            HepEVD::MCHits evd_mc_hits;
            for (size_t j = 0; j < ev.hit_x.size(); ++j)
            {
                evd_hits.push_back(new HepEVD::Hit({ev.hit_x[j], ev.hit_y[j], ev.hit_z[j]}, ev.hit_E[j]));

                int best_pdg = 0;
                if (ev.hit_pdg[j].empty())
                    continue; // No matches for this hit, skip it.

                // Find the match with the highest packet fraction
                size_t best_idx = 0;
                float max_frac = -1.0f;
                for (size_t m = 0; m < ev.hit_packetFrac[j].size(); ++m)
                {
                    if (ev.hit_packetFrac[j][m] > max_frac)
                    {
                        max_frac = ev.hit_packetFrac[j][m];
                        best_idx = m;
                    }
                }
                best_pdg = ev.hit_pdg[j][best_idx];

                evd_mc_hits.push_back(new HepEVD::MCHit({ev.hit_x[j], ev.hit_y[j], ev.hit_z[j]}, best_pdg, ev.hit_E[j]));
            }

            server->addHits(evd_hits);
            server->addMCHits(evd_mc_hits);

            server->startServer();
            server->resetServer();
        }

        const auto t1_total = ndlar::SteadyClock::now();
        const double n = num_events > 0 ? static_cast<double>(num_events) : 1.0;
        std::cout << "-------------------------------------------\n";
        std::cout << "Timing summary: total_ms=" << ndlar::elapsed_ms(t0_total, t1_total) << ", avg_collect_ms=" << (total_collect_ms / n) << "\n";

#if defined(ROOT_FOUND)
        std::cout << "Writing event products to ROOT file 'event_products.root'.\n";
        root_writer->Write();
        std::cout << "ROOT file write complete.\n";
#endif
    }
    catch (const HighFive::Exception &err)
    {
        std::cerr << "HDF5 Error: " << err.what() << "\n";
        return 1;
    }
    catch (const std::exception &err)
    {
        std::cerr << "Runtime error: " << err.what() << "\n";
        return 1;
    }

    return 0;
}
