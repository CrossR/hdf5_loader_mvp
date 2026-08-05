#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include "ndlar/common.hpp"

namespace ndlar::hdf5
{

using RefPair = std::array<uint32_t, 2>;

// Row schema for charge/calib_prompt_hits/data.
struct PromptHit
{
    uint32_t id;
    float x, y, z, Q, E, ts_pps;
};

// Row schema for charge/events/data.
struct EventRow
{
    int64_t id;
    int64_t ts_start;
    int64_t ts_end;
    int64_t unix_ts;
    int64_t unix_ts_usec;
};

// Row schema for charge/ext_trigs/data.
struct ExtTrig
{
    int32_t iogroup;
};

// Row schema for mc_truth/trajectories/data.
struct Trajectory
{
    int64_t event_id;
    std::array<float, 3> xyz_start;
    std::array<float, 3> xyz_end;
    uint32_t file_traj_id;
    uint32_t traj_id;
    int32_t pdg_id;
    float E_start;
    std::array<float, 3> pxyz_start;
    uint64_t vertex_id;
    int64_t parent_id;
};

// Row schema for mc_truth/interactions/data.
struct Interaction
{
    int64_t event_id;
    uint64_t vertex_id;
    float x_vert;
    float y_vert;
    float z_vert;
    float Enu;
    int32_t nu_pdg;
    std::array<float, 4> nu_4mom;
    uint8_t isCC;
    uint8_t isQES;
    uint8_t isRES;
    uint8_t isDIS;
    uint8_t isMEC;
    uint8_t isCOH;
};

// Row schema for mc_truth/packet_fraction/data.
struct PacketFraction
{
    std::array<int64_t, 20> segment_ids;
    std::array<double, 20> fraction;
};

// [start, stop) row interval stored in ref_region datasets.
struct RefRegion
{
    int32_t start;
    int32_t stop;
};

// Row schema for mc_truth/segments/data.
struct TrueSegment
{
    uint32_t segment_id;
    int32_t pdg_id;
    uint32_t file_traj_id;
    uint32_t traj_id;
    uint64_t vertex_id;
    int64_t event_id;
};

// Event payload consumed by HepEVD.
struct EventProducts
{
    int32_t trigger_id = ndlar::kInvalidTrigger;
    int64_t spill_id = -1;
    int32_t event_start_t = -5;
    int32_t event_end_t = -5;
    int32_t unix_ts = -5;
    int32_t unix_ts_usec = -5;

    std::vector<float> hit_x;
    std::vector<float> hit_y;
    std::vector<float> hit_z;
    std::vector<float> hit_charge;
    std::vector<float> hit_E;
    std::vector<float> hit_ts;
    std::vector<uint16_t> hit_matches;

    std::vector<std::vector<float>> hit_packetFrac;
    std::vector<std::vector<int32_t>> hit_pdg;
    std::vector<std::vector<int32_t>> hit_segmentID;
    std::vector<std::vector<int64_t>> hit_particleID;
    std::vector<std::vector<int64_t>> hit_particleIDLocal;
    std::vector<std::vector<int64_t>> hit_vertexID;

    std::vector<float> mcp_startx;
    std::vector<float> mcp_starty;
    std::vector<float> mcp_startz;
    std::vector<float> mcp_endx;
    std::vector<float> mcp_endy;
    std::vector<float> mcp_endz;
    std::vector<int64_t> mcp_id;
    std::vector<int64_t> mcp_idLocal;
    std::vector<int32_t> mcp_pdg;
    std::vector<float> mcp_energy;
    std::vector<float> mcp_px;
    std::vector<float> mcp_py;
    std::vector<float> mcp_pz;
    std::vector<int64_t> mcp_nuid;
    std::vector<int64_t> mcp_mother;

    std::vector<int64_t> nuID;
    std::vector<float> nue;
    std::vector<int32_t> nuPDG;
    std::vector<float> nupx;
    std::vector<float> nupy;
    std::vector<float> nupz;
    std::vector<float> nuvtxx;
    std::vector<float> nuvtxy;
    std::vector<float> nuvtxz;
    std::vector<int32_t> ccnc;
    std::vector<int32_t> mode;

    void reserve_hit_products(size_t hit_count);
    void reserve_trajectory_products(size_t trajectory_count);
    void reserve_interaction_products(size_t interaction_count);
};

// True when region encodes a non-empty half-open interval [start, stop).
bool is_valid_region(const RefRegion &region);

// Number of rows in a valid region, else zero.
int region_size(const RefRegion &region);

// Return the fraction corresponding to segment_id in one packet-fraction row.
float resolve_packet_fraction(const PacketFraction &row, uint32_t segment_id);

// Convert interaction flags into the mode encoding expected by output consumers.
int32_t interaction_mode(const Interaction &interaction);

// Append trajectory rows to event output arrays.
void append_trajectory_products(const std::vector<Trajectory> &rows, EventProducts &out);

// Append interaction rows to event output arrays.
void append_interaction_products(const std::vector<Interaction> &rows, EventProducts &out);

} // namespace ndlar::hdf5
