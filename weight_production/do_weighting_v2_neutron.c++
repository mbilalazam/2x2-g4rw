#include "TFile.h"
#include "TTree.h"
#include "EDepSim/TG4Event.h"
#include "EDepSim/TG4Trajectory.h"

#include "boost/program_options/options_description.hpp"
#include "boost/program_options/variables_map.hpp"
#include "boost/program_options/parsers.hpp"

#include "cetlib/filepath_maker.h"
#include "cetlib/filesystem.h"
#include "cetlib/search_path.h"
#include "fhiclcpp/intermediate_table.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "fhiclcpp/ParameterSet.h"

#include "geant4reweight/src/ReweightBase/G4ReweightManager.hh"
#include "geant4reweight/src/ReweightBase/G4MultiReweighter.hh"
#include "geant4reweight/src/ReweightBase/G4ReweightTraj.hh"
#include "geant4reweight/src/ReweightBase/G4ReweightStep.hh"

#include <string>
#include <iostream>
#include <set>
#include <stdexcept>

// ── helpers ───────────────────────────────────────────────────────────────

fhicl::ParameterSet MakeFCLPars(std::string fcl_file) {
    char const* fhicl_env = getenv("FHICL_FILE_PATH");
    std::string search_path;
    if (fhicl_env == nullptr) {
        std::cerr << "Expected environment variable FHICL_FILE_PATH is missing or empty: using \".\"\n";
        search_path = ".";
    } else {
        search_path = std::string{fhicl_env};
    }
    cet::filepath_first_absolute_or_lookup_with_dot lookupPolicy{search_path};
    return fhicl::ParameterSet::make(fcl_file, lookupPolicy);
}

TFile* OpenFile(const std::string filename) {
    TFile* theFile = nullptr;
    std::cout << "Searching for " << filename << std::endl;
    if (cet::file_exists(filename)) {
        theFile = new TFile(filename.c_str());
    } else {
        cet::search_path sp{"FW_SEARCH_PATH"};
        std::string found_filename;
        if (!sp.find_file(filename, found_filename))
            throw cet::exception("do_weighting") << "Could not find " << filename;
        theFile = new TFile(found_filename.c_str());
    }
    if (!theFile || theFile->IsZombie() || !theFile->IsOpen())
        throw cet::exception("do_weighting") << "Could not open " << filename;
    return theFile;
}

// PDG codes we want to reweight
std::vector<int> ValidPDGs { 2112 }; // for neutron

bool IsTarget(int pdg) {
    return std::find(ValidPDGs.begin(), ValidPDGs.end(), pdg) != ValidPDGs.end();
}

// Child PDGs needed for GetInteractionSubtype classification (proton, neutron)
bool IsChild(int pdg) {
    return (pdg == 2212 || pdg == 2112);
}

// ── Parse file index from filename ────────────────────────────────────────
// Expects a 7-digit number somewhere in the basename, e.g.:
//   MiniRun5_1E19_RHC.spill.0000003.EDEPSIM_SPILLS.root  ->  3
int ParseFileIdx(const std::string& filepath) {
    std::string base = filepath.substr(filepath.rfind('/') + 1);
    for (size_t i = 0; i + 7 <= base.size(); ++i) {
        bool all_digits = true;
        for (size_t j = i; j < i + 7; ++j) {
            if (!std::isdigit(base[j])) { all_digits = false; break; }
        }
        if (all_digits) {
            return std::stoi(base.substr(i, 7));
        }
    }
    std::cerr << "WARNING: could not parse 7-digit file index from: " << base
              << " -- defaulting to -1\n";
    return -1;
}

// ── main ──────────────────────────────────────────────────────────────────

int main(int argc, char** argv) {
    namespace po = boost::program_options;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help",                           "produce help message")
        ("i", po::value<std::string>(),    "Input EDepSim file")
        ("c", po::value<std::string>(),    "FCL file")
        ("o", po::value<std::string>(),    "Output file");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) { std::cout << desc << "\n"; return 1; }

    auto get_arg = [&](const char* key) -> std::string {
        if (!vm.count(key)) {
            std::cerr << "Need --" << key << std::endl; std::exit(1);
        }
        return vm[key].as<std::string>();
    };

    std::string input_file  = get_arg("i");
    std::string output_file = get_arg("o");
    std::string fcl_file    = get_arg("c");

    // Parse file index from input filename
    int file_idx = ParseFileIdx(input_file);
    std::cout << "File index parsed from filename: " << file_idx << std::endl;

    // ── FCL / fracs setup ─────────────────────────────────────────────────
    auto fcl_pars        = MakeFCLPars(fcl_file);
    auto material        = fcl_pars.get<fhicl::ParameterSet>("Material");
    auto fracs_file_name = fcl_pars.get<std::string>("FracsFile");
    auto g4rw_pars       = fcl_pars.get<std::vector<fhicl::ParameterSet>>("ParameterSet");
    // auto set_weight      = fcl_pars.get<double>("Weight", 2.);

    TFile* fracs_file = OpenFile(fracs_file_name);

    G4ReweightManager* RWManager = new G4ReweightManager({material});
    G4MultiReweighter* reweighter = new G4MultiReweighter(
        2112, *fracs_file, g4rw_pars, material, RWManager);
    
    // reweighter->SetAllParameterValues({set_weight});

    // ── Input tree ────────────────────────────────────────────────────────
    TFile f(input_file.c_str(), "open");
    TTree* tree    = (TTree*)f.Get("EDepSimEvents");
    TG4Event* event = nullptr;
    tree->SetBranchAddress("Event", &event);

    // ── Output tree ───────────────────────────────────────────────────────
    // One entry per (file, interaction, particle).
    // Lookup key in analysis:
    //   (file_idx, event_id, track_id)
    //   where event_id = InteractionNumber = (nu.id % 10000000) - 1
    //   and   track_id = prim.G4ID

    TFile fOut(output_file.c_str(), "recreate");
    TTree outtree("g4rwEvt", "G4RW weights per particle per event");

    int    out_file_idx  = file_idx;  // constant for this run
    int    out_spill_id  = -1;        // raw EDepSim spill index
    int    out_event_id  = -1;        // Primaries.InteractionNumber for this interaction
    int    out_run_id    = -1;        // RunId: <1e9=argon, 1e9=rock
    int    out_track_id  = -1;        // G4 TrackId (= CAF prim.G4ID)
    int    out_pdg       =  0;        // PDG code
    double out_weight_1p33 =  1.;     // +1sigma
    double out_weight_0p67 =  1.;     // -1sigma
    double out_weight_1p66 =  1.;     // +2sigma
    double out_weight_0p34 =  1.;     // -2sigma
    double out_weight_2p00 =  1.;     // +3sigma
    double out_len         =  0.;     // track length in cm

    outtree.Branch("file_idx",  &out_file_idx,  "file_idx/I");
    outtree.Branch("spill_id",  &out_spill_id,  "spill_id/I");
    outtree.Branch("event_id",  &out_event_id,  "event_id/I");
    outtree.Branch("run_id",    &out_run_id,    "run_id/I");
    outtree.Branch("track_id",  &out_track_id,  "track_id/I");
    outtree.Branch("pdg",       &out_pdg,       "pdg/I");
    outtree.Branch("weight_1p33", &out_weight_1p33, "weight_1p33/D");
    outtree.Branch("weight_0p67", &out_weight_0p67, "weight_0p67/D");
    outtree.Branch("weight_1p66", &out_weight_1p66, "weight_1p66/D");
    outtree.Branch("weight_0p34", &out_weight_0p34, "weight_0p34/D");
    outtree.Branch("weight_2p00", &out_weight_2p00, "weight_2p00/D");
    outtree.Branch("len",       &out_len,       "len/D");

    // ── Event loop ────────────────────────────────────────────────────────
    int n_written = 0;

    for (int i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);

        out_run_id  = event->RunId;
        out_spill_id = i;

        if (out_run_id >= 1000000000) continue;  // skip rock spills

        // ── Loop over primaries (neutrino interactions) in this spill ──────
        // A spill can have multiple neutrino interactions. Each primary
        // has its own InteractionNumber which matches CAF nu.id % 10000000 - 1.
        for (const auto& primary : event->Primaries) {

            // out_event_id = primary.InteractionNumber;
            // out_event_id = event->EventId;
            out_event_id = event->RunId * 1000000 + event->EventId;

            // Build set of TrackIds belonging to this primary's particles
            // These are the direct GENIE final-state particles for this interaction
            std::set<int> primaryTrackIds;
            for (const auto& part : primary.Particles) {
                primaryTrackIds.insert(part.TrackId);
            }

            // ── Loop over trajectories belonging to THIS primary ───────────
            for (const auto& traj : event->Trajectories) {

                // Only process direct primaries from neutrino vertex (ParentId==-1)
                if (traj.GetParentId() != -1) continue;

                // Only process target PDGs
                int pdg = traj.GetPDGCode();
                if (!IsTarget(pdg)) continue;

                if (traj.Points.empty()) continue;

                int track_id = traj.GetTrackId();

                // Only process trajectories belonging to THIS primary interaction
                if (primaryTrackIds.find(track_id) == primaryTrackIds.end()) continue;

                // ── Build G4RW trajectory ──────────────────────────────────
                G4ReweightTraj rw_traj(track_id, pdg, -1, i, {0,0});

                // Add child nucleons for GetInteractionSubtype classification (proton, neutron)
                for (const auto& child_traj : event->Trajectories) {
                    if (child_traj.GetParentId() != track_id) continue;
                    int cpid = child_traj.GetPDGCode();
                    if (!IsChild(cpid)) continue;
                    rw_traj.AddChild(new G4ReweightTraj(
                        child_traj.GetTrackId(), cpid, track_id, i, {0,0}));
                }

                // ── Build steps ───────────────────────────────────────────
                out_len = 0.;
                bool last_is_inel = false;
                if (traj.Points.size() > 1) {
                    const auto& last_pt = traj.Points.back();
                    last_is_inel = (last_pt.GetProcess() == 4 &&
                                    last_pt.GetSubprocess() == 121);
                }

                for (size_t j = 1; j < traj.Points.size(); ++j) {
                    std::string proc = "default";
                    if (j == traj.Points.size()-1 && last_is_inel)
                        proc = "neutronInelastic";

                    const auto& pt      = traj.Points[j];
                    const auto& prev_pt = traj.Points[j-1];
                    auto dist = pt.GetPosition() - prev_pt.GetPosition();
                    double len = sqrt(dist.X()*dist.X() +
                                      dist.Y()*dist.Y() +
                                      dist.Z()*dist.Z()) / 10.;
                    out_len += len;

                    double preP[3]  = { prev_pt.GetMomentum()[0],
                                        prev_pt.GetMomentum()[1],
                                        prev_pt.GetMomentum()[2] };
                    double postP[3] = { pt.GetMomentum()[0],
                                        pt.GetMomentum()[1],
                                        pt.GetMomentum()[2] };

                    G4ReweightStep* step = new G4ReweightStep(
                        track_id, pdg, 0, i, preP, postP, len, proc);
                    rw_traj.AddStep(step);
                }
    
                reweighter->SetAllParameterValues({1.33});
                out_weight_1p33 = reweighter->GetWeightFromSetParameters(rw_traj);
                reweighter->SetAllParameterValues({0.67});
                out_weight_0p67 = reweighter->GetWeightFromSetParameters(rw_traj);
                reweighter->SetAllParameterValues({1.66});
                out_weight_1p66 = reweighter->GetWeightFromSetParameters(rw_traj);
                reweighter->SetAllParameterValues({0.34});
                out_weight_0p34 = reweighter->GetWeightFromSetParameters(rw_traj);
                reweighter->SetAllParameterValues({2.00});
                out_weight_2p00 = reweighter->GetWeightFromSetParameters(rw_traj);
                
                out_track_id = track_id;
                out_pdg      = pdg;

                outtree.Fill();
                n_written++;
            }
        } // end primary loop
    } // end spill loop

    outtree.Write();
    fOut.Close();
    f.Close();

    std::cout << "Done. Written " << n_written << " entries to " << output_file << std::endl;
    return 0;
}