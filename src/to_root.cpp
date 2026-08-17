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
    tree_->SetAutoFlush(100);

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

    // Hit-level Variables
    this->CreateBranchVector("x", x_);
    this->CreateBranchVector("y", y_);
    this->CreateBranchVector("z", z_);
    this->CreateBranchVector("ts", ts_);
    this->CreateBranchVector("E", E_);
    this->CreateBranchVector("charge", charge_);
    this->CreateBranchVector("io_group", io_group_);
    this->CreateBranchVector("io_channel", io_channel_);
    this->CreateBranchVector("chip_id", chip_id_);
    this->CreateBranchVector("channel_id", channel_id_);

    // Early return if this is not MC, since the rest of the branches are MC-only.
    if (!has_mc_)
        return;

    // Initialize MC Pointers
    this->CreateBranchVector("hit_pdg", hit_pdg_);
    this->CreateBranchVector("hit_segmentID", hit_segmentID_);
    this->CreateBranchVector("hit_particleID", hit_particleID_);
    this->CreateBranchVector("hit_particleIDLocal", hit_particleIDLocal_);
    this->CreateBranchVector("hit_vertexID", hit_vertexID_);
    this->CreateBranchVector("hit_packetFrac", hit_packetFrac_);

    this->CreateBranchVector("mcp_px", mcp_px_);
    this->CreateBranchVector("mcp_py", mcp_py_);
    this->CreateBranchVector("mcp_pz", mcp_pz_);
    this->CreateBranchVector("mcp_energy", mcp_energy_);
    this->CreateBranchVector("mcp_startx", mcp_startx_);
    this->CreateBranchVector("mcp_starty", mcp_starty_);
    this->CreateBranchVector("mcp_startz", mcp_startz_);
    this->CreateBranchVector("mcp_endx", mcp_endx_);
    this->CreateBranchVector("mcp_endy", mcp_endy_);
    this->CreateBranchVector("mcp_endz", mcp_endz_);
    this->CreateBranchVector("mcp_id", mcp_id_);
    this->CreateBranchVector("mcp_idLocal", mcp_idLocal_);
    this->CreateBranchVector("mcp_nuid", mcp_nuid_);
    this->CreateBranchVector("mcp_vertex_id", mcp_vertex_id_);
    this->CreateBranchVector("mcp_mother", mcp_mother_);
    this->CreateBranchVector("mcp_pdg", mcp_pdg_);

    this->CreateBranchVector("nuvtxx", nuvtxx_);
    this->CreateBranchVector("nuvtxy", nuvtxy_);
    this->CreateBranchVector("nuvtxz", nuvtxz_);
    this->CreateBranchVector("nupx", nupx_);
    this->CreateBranchVector("nupy", nupy_);
    this->CreateBranchVector("nupz", nupz_);
    this->CreateBranchVector("nue", nue_);
    this->CreateBranchVector("nuID", nuID_);
    this->CreateBranchVector("vertex_id", vertex_id_);
    this->CreateBranchVector("nuPDG", nuPDG_);
    this->CreateBranchVector("mode", mode_);
    this->CreateBranchVector("ccnc", ccnc_);
}

template <typename T>
void RootWriter::CreateBranchVector(const std::string &name, std::vector<T>*& member_ptr)
{
    member_ptr = new std::vector<T>();
    tree_->Branch(name.c_str(), &member_ptr);
}

RootWriter::~RootWriter()
{
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
    hit_pdg_->resize(nhits_);
    hit_segmentID_->resize(nhits_);
    hit_particleID_->resize(nhits_);
    hit_particleIDLocal_->resize(nhits_);
    hit_vertexID_->resize(nhits_);
    hit_packetFrac_->resize(nhits_);

    for (size_t i = 0; i < ev.hit_x.size(); ++i)
    {
        // Standard assignment reuses capacity automatically
        (*hit_pdg_)[i] = ev.hit_pdg[i];
        (*hit_packetFrac_)[i] = ev.hit_packetFrac[i];

        // Use assign() to safely cast the 64-bit ints to long while reusing capacity
        (*hit_segmentID_)[i].assign(ev.hit_segmentID[i].begin(), ev.hit_segmentID[i].end());
        (*hit_particleID_)[i].assign(ev.hit_particleID[i].begin(), ev.hit_particleID[i].end());
        (*hit_particleIDLocal_)[i].assign(ev.hit_particleIDLocal[i].begin(), ev.hit_particleIDLocal[i].end());
        (*hit_vertexID_)[i].assign(ev.hit_vertexID[i].begin(), ev.hit_vertexID[i].end());
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
