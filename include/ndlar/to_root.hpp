#pragma once

#include "ndlar/hdf5/types.hpp"
#include <string>
#include <vector>

class TFile;
class TTree;

namespace ndlar
{

class RootWriter
{
public:
    RootWriter(const std::string &output_filename, bool has_mc = true);
    ~RootWriter();

    /**
     * @brief Fill the ROOT tree with the given event data.
     *
     * @param ev The event data to fill into the tree.
     * @param run The run number of the event.
     * @param subrun The subrun number of the event.
     * @param event_id The event number of the event.
     */
    void Fill(const hdf5::EventProducts &ev, int run, int subrun, int event_id);

    /**
     * @brief Write the ROOT file to disk.
     */
    void Write();

private:
    /**
     * @brief Create and register a vector for a branch in the ROOT tree.
     *
     * @tparam T The type of the vector elements.
     * @param name The name of the branch.
     * @return A pointer to the created vector.
     */
    template <typename T>
    void CreateBranchVector(const std::string &name, std::vector<T> *&member_ptr);

    TFile *file_ = nullptr;
    TTree *tree_ = nullptr;
    bool has_mc_;

    // Branch Variables
    int run_ = 0, subrun_ = 0, event_ = 0;
    int event_start_t_ = -5, event_end_t_ = -5;
    int triggers_ = 0, unix_ts_ = -5, unix_ts_usec_ = -5, nhits_ = 0;

    // Hits
    std::vector<float> *x_ = nullptr, *y_ = nullptr, *z_ = nullptr;
    std::vector<float> *ts_ = nullptr, *E_ = nullptr, *charge_ = nullptr;
    std::vector<unsigned char> *io_group_ = nullptr, *io_channel_ = nullptr, *chip_id_ = nullptr, *channel_id_ = nullptr;

    // Hit-level MC Truth
    std::vector<std::vector<int>> *hit_pdg_ = nullptr;
    std::vector<std::vector<long>> *hit_segmentID_ = nullptr;
    std::vector<std::vector<long>> *hit_particleID_ = nullptr;
    std::vector<std::vector<long>> *hit_particleIDLocal_ = nullptr;
    std::vector<std::vector<long>> *hit_vertexID_ = nullptr;
    std::vector<std::vector<float>> *hit_packetFrac_ = nullptr;

    // MCParticles
    std::vector<float> *mcp_px_ = nullptr, *mcp_py_ = nullptr, *mcp_pz_ = nullptr, *mcp_energy_ = nullptr;
    std::vector<float> *mcp_startx_ = nullptr, *mcp_starty_ = nullptr, *mcp_startz_ = nullptr;
    std::vector<float> *mcp_endx_ = nullptr, *mcp_endy_ = nullptr, *mcp_endz_ = nullptr;
    std::vector<long> *mcp_id_ = nullptr, *mcp_idLocal_ = nullptr, *mcp_nuid_ = nullptr, *mcp_vertex_id_ = nullptr, *mcp_mother_ = nullptr;
    std::vector<int> *mcp_pdg_ = nullptr;

    // Neutrinos
    std::vector<float> *nuvtxx_ = nullptr, *nuvtxy_ = nullptr, *nuvtxz_ = nullptr;
    std::vector<float> *nupx_ = nullptr, *nupy_ = nullptr, *nupz_ = nullptr, *nue_ = nullptr;
    std::vector<long> *nuID_ = nullptr, *vertex_id_ = nullptr;
    std::vector<int> *nuPDG_ = nullptr, *mode_ = nullptr, *ccnc_ = nullptr;
};

} // namespace ndlar
