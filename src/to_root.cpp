#include <stdexcept>

#include <TFile.h>
#include <TInterpreter.h>
#include <TTree.h>

#include "ndlar/to_root.hpp"

namespace ndlar
{

RootWriter::RootWriter(const std::string &output_filename, bool has_mc) :
    has_mc_(has_mc)
{
    // Ensure that ROOT knows about the vector types we are using for nested vectors.
    gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
    gInterpreter->GenerateDictionary("vector<vector<long> >", "vector");
    gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");

    file_ = new TFile(output_filename.c_str(), "RECREATE");

    if (!file_ || file_->IsZombie())
        throw std::runtime_error("Failed to create ROOT file: " + output_filename);

    tree_ = new TTree("events", "events");

    // Initialize pointers
    x_ = new std::vector<float>();
    y_ = new std::vector<float>();
    z_ = new std::vector<float>();
    ts_ = new std::vector<float>();
    E_ = new std::vector<float>();
    charge_ = new std::vector<float>();
    io_group_ = new std::vector<unsigned char>();
    io_channel_ = new std::vector<unsigned char>();
    chip_id_ = new std::vector<unsigned char>();
    channel_id_ = new std::vector<unsigned char>();

    // Basic Event Info
    tree_->Branch("run", &run_);
    tree_->Branch("subrun", &subrun_);
    tree_->Branch("event", &event_);
    tree_->Branch("event_start_t", &event_start_t_);
    tree_->Branch("event_end_t", &event_end_t_);
    tree_->Branch("triggers", &triggers_);
    tree_->Branch("unix_ts", &unix_ts_);
    tree_->Branch("unix_ts_usec", &unix_ts_usec_);
    tree_->Branch("nhits", &nhits_);

    // Hit Arrays
    tree_->Branch("x", &x_);
    tree_->Branch("y", &y_);
    tree_->Branch("z", &z_);
    tree_->Branch("ts", &ts_);
    tree_->Branch("E", &E_);
    tree_->Branch("charge", &charge_);
    tree_->Branch("io_group", &io_group_);
    tree_->Branch("io_channel", &io_channel_);
    tree_->Branch("chip_id", &chip_id_);
    tree_->Branch("channel_id", &channel_id_);

    // Early return if this is not MC, since the rest of the branches are MC-only.
    if (!has_mc_)
        return;

    // Initialize MC Pointers
    hit_pdg_ = new std::vector<std::vector<int>>();
    hit_segmentID_ = new std::vector<std::vector<long>>();
    hit_particleID_ = new std::vector<std::vector<long>>();
    hit_particleIDLocal_ = new std::vector<std::vector<long>>();
    hit_vertexID_ = new std::vector<std::vector<long>>();
    hit_packetFrac_ = new std::vector<std::vector<float>>();

    mcp_px_ = new std::vector<float>();
    mcp_py_ = new std::vector<float>();
    mcp_pz_ = new std::vector<float>();
    mcp_energy_ = new std::vector<float>();
    mcp_startx_ = new std::vector<float>();
    mcp_starty_ = new std::vector<float>();
    mcp_startz_ = new std::vector<float>();
    mcp_endx_ = new std::vector<float>();
    mcp_endy_ = new std::vector<float>();
    mcp_endz_ = new std::vector<float>();
    mcp_id_ = new std::vector<long>();
    mcp_idLocal_ = new std::vector<long>();
    mcp_nuid_ = new std::vector<long>();
    mcp_vertex_id_ = new std::vector<long>();
    mcp_mother_ = new std::vector<long>();
    mcp_pdg_ = new std::vector<int>();

    nuvtxx_ = new std::vector<float>();
    nuvtxy_ = new std::vector<float>();
    nuvtxz_ = new std::vector<float>();
    nupx_ = new std::vector<float>();
    nupy_ = new std::vector<float>();
    nupz_ = new std::vector<float>();
    nue_ = new std::vector<float>();
    nuID_ = new std::vector<long>();
    vertex_id_ = new std::vector<long>();
    nuPDG_ = new std::vector<int>();
    mode_ = new std::vector<int>();
    ccnc_ = new std::vector<int>();

    // MC Branches
    tree_->Branch("hit_pdg", &hit_pdg_);
    tree_->Branch("hit_segmentID", &hit_segmentID_);
    tree_->Branch("hit_particleID", &hit_particleID_);
    tree_->Branch("hit_particleIDLocal", &hit_particleIDLocal_);
    tree_->Branch("hit_vertexID", &hit_vertexID_);
    tree_->Branch("hit_packetFrac", &hit_packetFrac_);

    tree_->Branch("mcp_px", &mcp_px_);
    tree_->Branch("mcp_py", &mcp_py_);
    tree_->Branch("mcp_pz", &mcp_pz_);
    tree_->Branch("mcp_energy", &mcp_energy_);
    tree_->Branch("mcp_startx", &mcp_startx_);
    tree_->Branch("mcp_starty", &mcp_starty_);
    tree_->Branch("mcp_startz", &mcp_startz_);
    tree_->Branch("mcp_endx", &mcp_endx_);
    tree_->Branch("mcp_endy", &mcp_endy_);
    tree_->Branch("mcp_endz", &mcp_endz_);
    tree_->Branch("mcp_id", &mcp_id_);
    tree_->Branch("mcp_idLocal", &mcp_idLocal_);
    tree_->Branch("mcp_nuid", &mcp_nuid_);
    tree_->Branch("mcp_vertex_id", &mcp_vertex_id_);
    tree_->Branch("mcp_mother", &mcp_mother_);
    tree_->Branch("mcp_pdg", &mcp_pdg_);

    tree_->Branch("nuvtxx", &nuvtxx_);
    tree_->Branch("nuvtxy", &nuvtxy_);
    tree_->Branch("nuvtxz", &nuvtxz_);
    tree_->Branch("nupx", &nupx_);
    tree_->Branch("nupy", &nupy_);
    tree_->Branch("nupz", &nupz_);
    tree_->Branch("nue", &nue_);
    tree_->Branch("nuID", &nuID_);
    tree_->Branch("vertex_id", &vertex_id_);
    tree_->Branch("nuPDG", &nuPDG_);
    tree_->Branch("mode", &mode_);
    tree_->Branch("ccnc", &ccnc_);
}

RootWriter::~RootWriter()
{
    // Write automatically cleans up the file, but we should clean up vectors
    delete x_;
    delete y_;
    delete z_;
    delete ts_;
    delete E_;
    delete charge_;
    delete io_group_;
    delete io_channel_;
    delete chip_id_;
    delete channel_id_;

    if (has_mc_)
    {
        delete hit_pdg_;
        delete hit_segmentID_;
        delete hit_particleID_;
        delete hit_particleIDLocal_;
        delete hit_vertexID_;
        delete hit_packetFrac_;
        delete mcp_px_;
        delete mcp_py_;
        delete mcp_pz_;
        delete mcp_energy_;
        delete mcp_startx_;
        delete mcp_starty_;
        delete mcp_startz_;
        delete mcp_endx_;
        delete mcp_endy_;
        delete mcp_endz_;
        delete mcp_id_;
        delete mcp_idLocal_;
        delete mcp_nuid_;
        delete mcp_vertex_id_;
        delete mcp_mother_;
        delete mcp_pdg_;
        delete nuvtxx_;
        delete nuvtxy_;
        delete nuvtxz_;
        delete nupx_;
        delete nupy_;
        delete nupz_;
        delete nue_;
        delete nuID_;
        delete vertex_id_;
        delete nuPDG_;
        delete mode_;
        delete ccnc_;
    }
    if (file_)
    {
        file_->Close();
        delete file_;
    }
}

void RootWriter::Fill(const hdf5::EventProducts &ev, int run, int subrun, int event_id)
{
    // Set basic event info
    run_ = run;
    subrun_ = subrun;
    event_ = event_id;
    triggers_ = ev.trigger_id;
    nhits_ = static_cast<int>(ev.hit_x.size());
    event_start_t_ = ev.event_start_t;
    event_end_t_ = ev.event_end_t;
    unix_ts_ = ev.unix_ts;
    unix_ts_usec_ = ev.unix_ts_usec;

    // Copy hit data
    *x_ = ev.hit_x;
    *y_ = ev.hit_y;
    *z_ = ev.hit_z;
    *ts_ = ev.hit_ts;
    *E_ = ev.hit_E;
    *charge_ = ev.hit_charge;
    *io_group_ = ev.hit_io_group;
    *io_channel_ = ev.hit_io_channel;
    *chip_id_ = ev.hit_chip_id;
    *channel_id_ = ev.hit_channel_id;

    if (!has_mc_)
    {
        tree_->Fill();
        return;
    }

    // Group hit truth into nested vectors
    hit_pdg_->clear();
    hit_segmentID_->clear();
    hit_particleID_->clear();
    hit_particleIDLocal_->clear();
    hit_vertexID_->clear();
    hit_packetFrac_->clear();

    hit_pdg_->reserve(nhits_);
    hit_segmentID_->reserve(nhits_);
    hit_particleID_->reserve(nhits_);
    hit_particleIDLocal_->reserve(nhits_);
    hit_vertexID_->reserve(nhits_);
    hit_packetFrac_->reserve(nhits_);

    for (size_t i = 0; i < ev.hit_x.size(); ++i)
    {
        hit_pdg_->push_back(ev.hit_pdg[i]);
        hit_packetFrac_->push_back(ev.hit_packetFrac[i]);

        // Cast 64-bit ints to ROOT 'long' for safety across platforms
        std::vector<long> segs(ev.hit_segmentID[i].begin(), ev.hit_segmentID[i].end());
        std::vector<long> pids(ev.hit_particleID[i].begin(), ev.hit_particleID[i].end());
        std::vector<long> locals(ev.hit_particleIDLocal[i].begin(), ev.hit_particleIDLocal[i].end());
        std::vector<long> vtxs(ev.hit_vertexID[i].begin(), ev.hit_vertexID[i].end());

        hit_segmentID_->push_back(std::move(segs));
        hit_particleID_->push_back(std::move(pids));
        hit_particleIDLocal_->push_back(std::move(locals));
        hit_vertexID_->push_back(std::move(vtxs));
    }

    // Copy MCParticle arrays
    *mcp_px_ = ev.mcp_px;
    *mcp_py_ = ev.mcp_py;
    *mcp_pz_ = ev.mcp_pz;
    *mcp_energy_ = ev.mcp_energy;
    *mcp_startx_ = ev.mcp_startx;
    *mcp_starty_ = ev.mcp_starty;
    *mcp_startz_ = ev.mcp_startz;
    *mcp_endx_ = ev.mcp_endx;
    *mcp_endy_ = ev.mcp_endy;
    *mcp_endz_ = ev.mcp_endz;
    *mcp_pdg_ = ev.mcp_pdg;

    // Assign the various ID arrays
    mcp_id_->assign(ev.mcp_id.begin(), ev.mcp_id.end());
    mcp_idLocal_->assign(ev.mcp_idLocal.begin(), ev.mcp_idLocal.end());
    mcp_nuid_->assign(ev.mcp_nuid.begin(), ev.mcp_nuid.end());
    mcp_vertex_id_->assign(ev.mcp_nuid.begin(), ev.mcp_nuid.end());
    mcp_mother_->assign(ev.mcp_mother.begin(), ev.mcp_mother.end());

    // Copy Neutrino arrays
    *nuvtxx_ = ev.nuvtxx;
    *nuvtxy_ = ev.nuvtxy;
    *nuvtxz_ = ev.nuvtxz;
    *nupx_ = ev.nupx;
    *nupy_ = ev.nupy;
    *nupz_ = ev.nupz;
    *nue_ = ev.nue;
    *nuPDG_ = ev.nuPDG;
    *mode_ = ev.mode;
    *ccnc_ = ev.ccnc;

    nuID_->assign(ev.nuID.begin(), ev.nuID.end());
    vertex_id_->assign(ev.nuID.begin(), ev.nuID.end());

    tree_->Fill();
}

void RootWriter::Write()
{
    if (file_)
    {
        file_->cd();
        tree_->Write();
    }
}

} // namespace ndlar
