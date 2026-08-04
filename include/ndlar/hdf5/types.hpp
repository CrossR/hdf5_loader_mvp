#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include "ndlar/common.hpp"

namespace ndlar::hdf5
{

using RefPair = std::array<uint32_t, 2>;

/** A struct representing a prompt hit in the detector. */
struct PromptHit
{
    uint32_t id;
    float x, y, z, Q, E, ts_pps;
};

/** A struct representing an event row in the detector. */
struct EventRow
{
    int64_t id;
    int64_t ts_start;
    int64_t ts_end;
    int64_t unix_ts;
    int64_t unix_ts_usec;
};

/** A struct representing an external trigger in the detector. */
struct ExtTrig
{
    int32_t iogroup;
};

/** A struct representing a trajectory in the detector. */
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

/** A struct representing an interaction in the detector. */
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

/** A struct representing a packet fraction in the detector. */
struct PacketFraction
{
    std::array<int64_t, 20> segment_ids;
    std::array<double, 20> fraction;
};

/** A struct representing a reference region in the detector. */
struct RefRegion
{
    int32_t start;
    int32_t stop;
};

/** A struct representing a true segment in the detector. */
struct TrueSegment
{
    uint32_t segment_id;
    int32_t pdg_id;
    uint32_t file_traj_id;
    uint32_t traj_id;
    uint64_t vertex_id;
    int64_t event_id;
};

/** A struct representing the complete event products. */
struct EventProducts
{
    int32_t trigger_id = ndlar::kInvalidTrigger;
    int64_t spill_id = -1;

    std::vector<float> hit_x;
    std::vector<float> hit_y;
    std::vector<float> hit_z;
    std::vector<float> hit_charge;
    std::vector<float> hit_E;
    std::vector<float> hit_ts;
    std::vector<uint16_t> hit_matches;

    std::vector<float> hit_packetFrac;
    std::vector<int32_t> hit_pdg;
    std::vector<int32_t> hit_segmentID;
    std::vector<int64_t> hit_particleID;
    std::vector<int64_t> hit_particleIDLocal;
    std::vector<int64_t> hit_vertexID;

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

/**
 * A struct representing the context for streaming event products from an HDF5 file.
 *
 * @param region The reference region for the current event.
 * @return bool indicating whether the region is valid (true) or not (false).
 */
bool is_valid_region(const RefRegion &region);

/**
 * Returns the size of a reference region.
 *
 * @param region The reference region.
 * @return The size of the region.
 */
int region_size(const RefRegion &region);

/**
 * Resolves the packet fraction for a given segment ID from a PacketFraction row.
 *
 * @param row The PacketFraction row containing segment IDs and fractions.
 * @param segment_id The segment ID for which to resolve the fraction.
 * @return The resolved fraction as a float. Returns 0.0 if the segment ID is not found.
 */
float resolve_packet_fraction(const PacketFraction &row, uint32_t segment_id);

/**
 * Determines the interaction mode based on the provided Interaction struct.
 *
 * @param interaction The Interaction struct containing interaction details.
 * @return An int32_t representing the interaction mode.
 */
int32_t interaction_mode(const Interaction &interaction);

/**
 * Appends trajectory products to the provided EventProducts struct.
 *
 * @param rows The vector of Trajectory rows to append.
 * @param out The EventProducts struct to which the trajectory products will be appended.
 */
void append_trajectory_products(const std::vector<Trajectory> &rows, EventProducts &out);

/**
 * Appends interaction products to the provided EventProducts struct.
 *
 * @param rows The vector of Interaction rows to append.
 * @param out The EventProducts struct to which the interaction products will be appended.
 */
void append_interaction_products(const std::vector<Interaction> &rows, EventProducts &out);

} // namespace ndlar::hdf5
