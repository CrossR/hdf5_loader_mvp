#include "ndlar/hdf5/types.hpp"

namespace ndlar::hdf5 {

void EventProducts::reserve_hit_products(size_t hit_count) {
    hit_x.reserve(hit_count);
    hit_y.reserve(hit_count);
    hit_z.reserve(hit_count);
    hit_charge.reserve(hit_count);
    hit_E.reserve(hit_count);
    hit_ts.reserve(hit_count);
    hit_matches.reserve(hit_count);
    hit_packetFrac.reserve(hit_count);
    hit_pdg.reserve(hit_count);
    hit_segmentID.reserve(hit_count);
    hit_particleID.reserve(hit_count);
    hit_particleIDLocal.reserve(hit_count);
    hit_vertexID.reserve(hit_count);
}

void EventProducts::reserve_trajectory_products(size_t trajectory_count) {
    mcp_startx.reserve(trajectory_count);
    mcp_starty.reserve(trajectory_count);
    mcp_startz.reserve(trajectory_count);
    mcp_endx.reserve(trajectory_count);
    mcp_endy.reserve(trajectory_count);
    mcp_endz.reserve(trajectory_count);
    mcp_id.reserve(trajectory_count);
    mcp_idLocal.reserve(trajectory_count);
    mcp_pdg.reserve(trajectory_count);
    mcp_energy.reserve(trajectory_count);
    mcp_px.reserve(trajectory_count);
    mcp_py.reserve(trajectory_count);
    mcp_pz.reserve(trajectory_count);
    mcp_nuid.reserve(trajectory_count);
    mcp_mother.reserve(trajectory_count);
}

void EventProducts::reserve_interaction_products(size_t interaction_count) {
    nuID.reserve(interaction_count);
    nue.reserve(interaction_count);
    nuPDG.reserve(interaction_count);
    nupx.reserve(interaction_count);
    nupy.reserve(interaction_count);
    nupz.reserve(interaction_count);
    nuvtxx.reserve(interaction_count);
    nuvtxy.reserve(interaction_count);
    nuvtxz.reserve(interaction_count);
    ccnc.reserve(interaction_count);
    mode.reserve(interaction_count);
}

bool is_valid_region(const RefRegion& region) {
    return region.start >= 0 && region.stop > region.start;
}

int region_size(const RefRegion& region) {
    if (!is_valid_region(region)) {
        return 0;
    }
    return region.stop - region.start;
}

float resolve_packet_fraction(const PacketFraction& row, uint32_t segment_id) {
    for (size_t i = 0; i < row.segment_ids.size(); ++i) {
        if (row.segment_ids[i] == static_cast<int64_t>(segment_id)) {
            return static_cast<float>(row.fraction[i]);
        }
    }
    return 0.0f;
}

int32_t interaction_mode(const Interaction& interaction) {
    int32_t mode = 1000;
    if (interaction.isQES) mode = 0;
    if (interaction.isRES) mode = 1;
    if (interaction.isDIS) mode = 2;
    if (interaction.isCOH) mode = 3;
    if (interaction.isCOH && interaction.isQES) mode = 4;
    if (interaction.isMEC) mode = 10;
    return mode;
}

void append_trajectory_products(const std::vector<Trajectory>& rows, EventProducts& out) {
    for (const Trajectory& row : rows) {
        out.mcp_startx.push_back(row.xyz_start[0]);
        out.mcp_starty.push_back(row.xyz_start[1]);
        out.mcp_startz.push_back(row.xyz_start[2]);
        out.mcp_endx.push_back(row.xyz_end[0]);
        out.mcp_endy.push_back(row.xyz_end[1]);
        out.mcp_endz.push_back(row.xyz_end[2]);
        out.mcp_id.push_back(static_cast<int64_t>(row.file_traj_id));
        out.mcp_idLocal.push_back(static_cast<int64_t>(row.traj_id));
        out.mcp_pdg.push_back(row.pdg_id);
        out.mcp_energy.push_back(row.E_start * ndlar::kMeVToGeV);
        out.mcp_px.push_back(row.pxyz_start[0] * ndlar::kMeVToGeV);
        out.mcp_py.push_back(row.pxyz_start[1] * ndlar::kMeVToGeV);
        out.mcp_pz.push_back(row.pxyz_start[2] * ndlar::kMeVToGeV);
        out.mcp_nuid.push_back(static_cast<int64_t>(row.vertex_id));
        out.mcp_mother.push_back(row.parent_id);
    }
}

void append_interaction_products(const std::vector<Interaction>& rows, EventProducts& out) {
    for (const Interaction& row : rows) {
        out.nuID.push_back(static_cast<int64_t>(row.vertex_id));
        out.nue.push_back(row.Enu * ndlar::kMeVToGeV);
        out.nuPDG.push_back(row.nu_pdg);
        out.nupx.push_back(row.nu_4mom[0] * ndlar::kMeVToGeV);
        out.nupy.push_back(row.nu_4mom[1] * ndlar::kMeVToGeV);
        out.nupz.push_back(row.nu_4mom[2] * ndlar::kMeVToGeV);
        out.nuvtxx.push_back(row.x_vert);
        out.nuvtxy.push_back(row.y_vert);
        out.nuvtxz.push_back(row.z_vert);
        out.ccnc.push_back(row.isCC ? 0 : 1);
        out.mode.push_back(interaction_mode(row));
    }
}

}  // namespace ndlar::hdf5
