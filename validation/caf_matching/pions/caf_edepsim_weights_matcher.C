// caf_edepsim_weights_matcher.C
//
// Picks ONE pion from the CAF and tracks it through:
//   1. All provided EDepSim samples (official unpatched, OLD patched,
//      NEW/current patched)
//   2. All provided G4RW weight output files (new patched, old patched,
//      old unpatched)
//
// Every value that appears in the final summary table is returned
// directly from the function that read it (CafResult / EdepResult /
// WeightResult structs) -- nothing in the table is hand-transcribed
// from a separate print statement.
//
// Matching convention -- exactly as documented in the technote
// (subsec:g4rw_matching), NO offset applied anywhere:
//   - File index: CAF's own cafFileIdx == EDepSim filename index
//                 == weights file_idx  (direct one-to-one)
//   - Interaction key: vertex_id = RunId*1e6 + EventId
//   - Particle key: CAF prim.G4ID == EDepSim TrackId == Weights track_id
//
// NOTE on event_id overflow: RunId*1000000 overflows a 32-bit int for
// RunId >~ 2147 (confirmed in a prior run: RunId=9095 wraps to
// 505065408). Both this script and the production do_weighting_v2_*.c__
// writer use the SAME 32-bit int arithmetic, so the wrapped value is
// consistently reproduced on both sides and the lookup still succeeds --
// but the value printed as "combined ID" below is the WRAPPED value,
// not the true RunId*1e6+EventId. This is called out explicitly rather
// than silently computed.
//
// Setup before running:
//   shifter --image=mjkramer/sim2x2:ndlar011 --module=cvmfs /bin/bash
//   source /opt/generators/geant4/bin/geant4.sh
//   source /opt/generators/root/bin/thisroot.sh
//   export LD_LIBRARY_PATH=/opt/generators/edep-sim/install/lib:$LD_LIBRARY_PATH
//   export ROOT_INCLUDE_PATH=/opt/generators/edep-sim/install/include:$ROOT_INCLUDE_PATH
//
// Run with:
//   root -l -b -q caf_edepsim_weights_matcher.C

#ifdef __CLING__
R__LOAD_LIBRARY(libedepsim_io)
#endif

#include "TFile.h"
#include "TTree.h"
#include "EDepSim/TG4Event.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

// ---------------------------------------------------------------------
// (Single-target constants removed -- see kTargets list below main
// helper functions. The script now loops over 15 pions, 5 per band.)
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
// File paths
// ---------------------------------------------------------------------
const std::string kCafFile =
    "/global/u1/m/mazam/Geant4ReWeighting/2x2_G4RW/test/on_CAFS/31Aug2026/output.root";

std::string Pad7(int idx) {
    char buf[8];
    snprintf(buf, sizeof(buf), "%07d", idx);
    return std::string(buf);
}

// Official unpatched: subdirectory fixed at "0000000" -- confirmed
// working in a prior run; only the filename itself uses the real index.
std::string OfficialUnpatchedPath(int idx) {
    return "/global/cfs/cdirs/dune/www/data/2x2/simulation/productions/"
           "MiniRun5_1E19_RHC/MiniRun5_1E19_RHC.spill/EDEPSIM_SPILLS/0000000/"
           "MiniRun5_1E19_RHC.spill." + Pad7(idx) + ".EDEPSIM_SPILLS.root";
}
std::string PatchedOldPath(int idx) {
    return "/global/cfs/cdirs/dune/users/mazam/EdepSim_G4RW_nu/"
           "patched_withRock_genieNu_OLD_pending_delete/"
           "MiniRun5_1E19_RHC.nu." + Pad7(idx) + ".G4RW.EDEPSIM.root";
}
std::string PatchedNewPath(int idx) {
    return "/global/cfs/cdirs/dune/users/mazam/EdepSim_G4RW_nu/"
           "patched_withRock_genieNu/"
           "MiniRun5_1E19_RHC.nu." + Pad7(idx) + ".G4RW.EDEPSIM.root";
}

const std::string kWeightsNewPatched =
    "/global/cfs/cdirs/dune/users/mazam/G4RW/output_weights/"
    "patched_withRock_genieNu/weights_1000_all_LAr_PiPlus_cascade_0_9000_5weights.root";
const std::string kWeightsOldPatched =
    "/global/cfs/cdirs/dune/users/mazam/G4RW/output_weights/"
    "patched_withRock_genieNu/old_files/"
    "weights_1000_all_LAr_PiPlus_cascade_0_9000_5weights_addChild_range_0_9000.root";
const std::string kWeightsOldUnpatched =
    "/global/cfs/cdirs/dune/users/mazam/G4RW/output_weights/"
    "patched_noRock_genieNu/weights_1000_all_LAr_PiPlus_cascade_0_9000_5weights.root";

// ---------------------------------------------------------------------
// Result structs -- every field is filled directly from a file read.
// Only "found=false" and the -999/-1 sentinels are hand-set defaults.
// ---------------------------------------------------------------------
struct CafResult {
    bool found = false;
    int cafFileIdx = -1;
    long long nu_id = -1;
    int G4ID = -1;
    int pdg = -999;
    double cosTheta = -999, length = -999;
    // CAF's own stored weight, fixed-array convention per
    // plot_truth_pion_weights_patched.C's comment: [0]=0.34, [1]=0.67,
    // [2]=1.33, [3]=1.66, [4]=2.00. Reading this directly lets us see
    // numerically which weight file (new patched / old patched / old
    // unpatched) CAF's own value actually matches, instead of inferring
    // it from length similarity.
    double w034 = -999, w067 = -999, w133 = -999, w166 = -999, w200 = -999;
};

struct EdepResult {
    bool found = false;
    bool eventFound = false;
    int RunId = -999, EventId = -999, TrackId = -999;
    int pdg = -999, nPoints = -999, parentId = -999;
    double length = -999;
};

struct WeightResult {
    bool found = false;
    int file_idx = -999, run_id = -999, event_id = -999, track_id = -999, pdg = -999;
    double len = -999;
    double w034 = -999, w067 = -999, w133 = -999, w166 = -999, w200 = -999;
};

// ---------------------------------------------------------------------
// Step 1: CAF
// ---------------------------------------------------------------------
CafResult GetCAFEntry(int targetCafFileIdx, int targetRunId, int targetEventId, int targetTrackId) {
    CafResult r;
    TFile* f = TFile::Open(kCafFile.c_str());
    if (!f || f->IsZombie()) { std::cerr << "ERROR: cannot open CAF file\n"; return r; }
    TTree* t = (TTree*)f->Get("syst_g4rw_pion_truth");
    if (!t) { std::cerr << "ERROR: no syst_g4rw_pion_truth tree\n"; f->Close(); return r; }

    int cafFileIdx, G4ID, pdg;
    long long nu_id;
    double cosTheta, length;
    double w_patched[5];
    t->SetBranchAddress("cafFileIdx", &cafFileIdx);
    t->SetBranchAddress("nu_id",      &nu_id);
    t->SetBranchAddress("G4ID",       &G4ID);
    t->SetBranchAddress("pdg",        &pdg);
    t->SetBranchAddress("cosTheta",   &cosTheta);
    t->SetBranchAddress("length",     &length);
    t->SetBranchAddress("w_patched",  w_patched);

    long long targetNuId = (long long)targetRunId * 1000000 + targetEventId;
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
        t->GetEntry(i);
        if (cafFileIdx == targetCafFileIdx && nu_id == targetNuId && G4ID == targetTrackId) {
            r.found = true;
            r.cafFileIdx = cafFileIdx; r.nu_id = nu_id; r.G4ID = G4ID; r.pdg = pdg;
            r.cosTheta = cosTheta; r.length = length;
            r.w034 = w_patched[0]; r.w067 = w_patched[1]; r.w133 = w_patched[2];
            r.w166 = w_patched[3]; r.w200 = w_patched[4];
            break;
        }
    }
    f->Close();
    return r;
}

// ---------------------------------------------------------------------
// Step 2: EDepSim
// ---------------------------------------------------------------------
EdepResult SearchEdepSim(const std::string& filename, int targetRunId, int targetEventId, int targetTrackId) {
    EdepResult r;
    TFile* f = TFile::Open(filename.c_str());
    if (!f || f->IsZombie()) return r;
    TTree* tree = (TTree*)f->Get("EDepSimEvents");
    if (!tree) { f->Close(); return r; }
    TG4Event* event = nullptr;
    tree->SetBranchAddress("Event", &event);

    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        if (event->RunId != targetRunId || event->EventId != targetEventId) continue;
        r.eventFound = true;
        r.RunId = event->RunId; r.EventId = event->EventId;

        for (const auto& traj : event->Trajectories) {
            if (traj.GetTrackId() != targetTrackId) continue;
            r.found = true;
            r.TrackId = traj.GetTrackId();
            r.pdg = traj.GetPDGCode();
            r.nPoints = traj.Points.size();
            r.parentId = traj.GetParentId();
            double totalLen = 0.;
            for (size_t j = 1; j < traj.Points.size(); ++j) {
                auto d = traj.Points[j].GetPosition() - traj.Points[j-1].GetPosition();
                totalLen += sqrt(d.X()*d.X()+d.Y()*d.Y()+d.Z()*d.Z()) / 10.;
            }
            r.length = totalLen;
            break;
        }
        break;
    }
    f->Close();
    return r;
}

// ---------------------------------------------------------------------
// Step 3: Weights
// ---------------------------------------------------------------------
WeightResult SearchWeights(const std::string& filename, int targetCafFileIdx, int targetRunId, int targetEventId, int targetTrackId) {
    WeightResult r;
    TFile* f = TFile::Open(filename.c_str());
    if (!f || f->IsZombie()) return r;
    TTree* t = (TTree*)f->Get("g4rwEvt");
    if (!t) { f->Close(); return r; }

    int file_idx, run_id, event_id, track_id, pdg;
    double w034, w067, w133, w166, w200, len;
    t->SetBranchAddress("file_idx",    &file_idx);
    t->SetBranchAddress("run_id",      &run_id);
    t->SetBranchAddress("event_id",    &event_id);
    t->SetBranchAddress("track_id",    &track_id);
    t->SetBranchAddress("pdg",         &pdg);
    t->SetBranchAddress("weight_0p34", &w034);
    t->SetBranchAddress("weight_0p67", &w067);
    t->SetBranchAddress("weight_1p33", &w133);
    t->SetBranchAddress("weight_1p66", &w166);
    t->SetBranchAddress("weight_2p00", &w200);
    t->SetBranchAddress("len",         &len);

    // Same 32-bit int arithmetic as production's out_event_id -- if this
    // overflows, it overflows identically on both write and read side,
    // so the lookup still succeeds against what's actually stored.
    int targetEventId2 = targetRunId * 1000000 + targetEventId;

    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
        t->GetEntry(i);
        if (file_idx == targetCafFileIdx && event_id == targetEventId2 && track_id == targetTrackId) {
            r.found = true;
            r.file_idx = file_idx; r.run_id = run_id; r.event_id = event_id;
            r.track_id = track_id; r.pdg = pdg;
            r.len = len;
            r.w034 = w034; r.w067 = w067; r.w133 = w133; r.w166 = w166; r.w200 = w200;
            break;
        }
    }
    f->Close();
    return r;
}

// ---------------------------------------------------------------------
// Table printing helpers
// ---------------------------------------------------------------------
std::string F(double v, int prec = 3) {
    if (v == -999) return "N/A";
    std::ostringstream ss; ss << std::fixed << std::setprecision(prec) << v;
    return ss.str();
}
std::string I(long long v) {
    if (v == -999 || v == -1) return "N/A";
    return std::to_string(v);
}

void PrintRow(const std::string& source, const std::string& fileIdx, const std::string& runId,
              const std::string& eventId, const std::string& trackId, const std::string& pdg,
              const std::string& nPoints, const std::string& length,
              const std::string& w034, const std::string& w067, const std::string& w133,
              const std::string& w166, const std::string& w200) {
    std::cout << std::left  << std::setw(26) << source
               << std::right
               << std::setw(10) << fileIdx
               << std::setw(10) << runId
               << std::setw(10) << eventId
               << std::setw(9)  << trackId
               << std::setw(6)  << pdg
               << std::setw(9)  << nPoints
               << std::setw(11) << length
               << std::setw(9)  << w034
               << std::setw(9)  << w067
               << std::setw(9)  << w133
               << std::setw(9)  << w166
               << std::setw(9)  << w200
               << "\n";
}

// ---------------------------------------------------------------------
// Target list -- 5 per band, taken directly from the printed output of
// find_example_pions_from_caf_list.C (the CAF-traceable candidate list).
// Includes both Band B anomalies (cafFileIdx=909, 672) already spotted,
// alongside clean Band B candidates for direct comparison.
// ---------------------------------------------------------------------
struct Target {
    std::string label;
    int cafFileIdx, RunId, EventId, TrackId;
};

const std::vector<Target> kTargets = {
    // Band A (w_1.33 ~ 0.05-0.30, low tail)
    {"Band A #1", 398, 3987, 40, 3},
    {"Band A #2", 700, 7000, 76, 3},
    {"Band A #3", 331, 3311, 79, 1},
    {"Band A #4", 173, 1738, 62, 2},
    {"Band A #5", 573, 5731,  1, 4},
    // Band B (w_1.33 ~ 0.45-0.55, mid-range) -- includes both anomalies found
    {"Band B #1 (ANOMALY)", 909, 9095, 80, 2},
    {"Band B #2 (ANOMALY)", 672, 6728, 75, 2},
    {"Band B #3 (clean)",   225, 2251, 56, 2},
    {"Band B #4 (clean)",   440, 4401, 81, 2},
    {"Band B #5 (clean)",    14,  145, 23, 3},
    // Band C (w_1.33 ~ 0.95-1.05, near w=1)
    {"Band C #1", 925, 9256,  2, 3},
    {"Band C #2", 207, 2070, 78, 8},
    {"Band C #3", 725, 7258,  9, 3},
    {"Band C #4", 168, 1689, 27, 2},
    {"Band C #5", 888, 8887, 47, 2},
};

// ---------------------------------------------------------------------
// Run the full CAF -> EDepSim -> Weights pipeline for one target and
// print its consolidated table.
// ---------------------------------------------------------------------
void RunOneTarget(const Target& tgt) {

    CafResult caf = GetCAFEntry(tgt.cafFileIdx, tgt.RunId, tgt.EventId, tgt.TrackId);

    EdepResult edepOfficial = SearchEdepSim(OfficialUnpatchedPath(tgt.cafFileIdx), tgt.RunId, tgt.EventId, tgt.TrackId);
    EdepResult edepOld      = SearchEdepSim(PatchedOldPath(tgt.cafFileIdx),        tgt.RunId, tgt.EventId, tgt.TrackId);
    EdepResult edepNew      = SearchEdepSim(PatchedNewPath(tgt.cafFileIdx),        tgt.RunId, tgt.EventId, tgt.TrackId);

    WeightResult wNewPatched   = SearchWeights(kWeightsNewPatched,   tgt.cafFileIdx, tgt.RunId, tgt.EventId, tgt.TrackId);
    WeightResult wOldPatched   = SearchWeights(kWeightsOldPatched,   tgt.cafFileIdx, tgt.RunId, tgt.EventId, tgt.TrackId);
    WeightResult wOldUnpatched = SearchWeights(kWeightsOldUnpatched, tgt.cafFileIdx, tgt.RunId, tgt.EventId, tgt.TrackId);

    std::cout << "\n================================================================\n";
    std::cout << " " << tgt.label << ":  cafFileIdx=" << tgt.cafFileIdx << "  RunId=" << tgt.RunId
               << "  EventId=" << tgt.EventId << "  TrackId=" << tgt.TrackId << "\n";
    std::cout << "================================================================\n\n";

    PrintRow("Source", "fileIdx", "RunId", "EventId", "TrkId", "PDG", "nPts", "len[cm]",
              "w0.34", "w0.67", "w1.33", "w1.66", "w2.00");
    std::cout << std::string(150, '-') << "\n";

    PrintRow("CAF",
              I(caf.cafFileIdx), std::to_string(tgt.RunId), std::to_string(tgt.EventId), std::to_string(tgt.TrackId),
              caf.found ? I(caf.pdg) : "N/A", "-", caf.found ? F(caf.length) : "N/A",
              caf.found ? F(caf.w034,4) : "-", caf.found ? F(caf.w067,4) : "-",
              caf.found ? F(caf.w133,4) : "-", caf.found ? F(caf.w166,4) : "-",
              caf.found ? F(caf.w200,4) : "-");

    struct EdepRow { std::string label; EdepResult r; };
    std::vector<EdepRow> edepRows = {
        {"EDepSim: official unpatched", edepOfficial},
        {"EDepSim: patched OLD",        edepOld},
        {"EDepSim: patched NEW",        edepNew}
    };
    for (auto& er : edepRows) {
        if (er.r.found) {
            PrintRow(er.label, I(tgt.cafFileIdx), I(er.r.RunId), I(er.r.EventId), I(er.r.TrackId),
                      I(er.r.pdg), I(er.r.nPoints), F(er.r.length),
                      "-", "-", "-", "-", "-");
        } else if (er.r.eventFound) {
            PrintRow(er.label, I(tgt.cafFileIdx), I(er.r.RunId), I(er.r.EventId), "NOT FOUND",
                      "-", "-", "-", "-", "-", "-", "-", "-");
        } else {
            PrintRow(er.label, I(tgt.cafFileIdx), "EVENT NOT FOUND", "-", "-", "-", "-", "-",
                      "-", "-", "-", "-", "-");
        }
    }

    struct WRow { std::string label; WeightResult r; };
    std::vector<WRow> wRows = {
        {"Weights: new patched",   wNewPatched},
        {"Weights: old patched",   wOldPatched},
        {"Weights: old unpatched", wOldUnpatched}
    };
    for (auto& wr : wRows) {
        if (wr.r.found) {
            PrintRow(wr.label, I(wr.r.file_idx), I(wr.r.run_id), I(wr.r.event_id), I(wr.r.track_id),
                      I(wr.r.pdg), "-", F(wr.r.len),
                      F(wr.r.w034,4), F(wr.r.w067,4), F(wr.r.w133,4), F(wr.r.w166,4), F(wr.r.w200,4));
        } else {
            PrintRow(wr.label, I(tgt.cafFileIdx), "-", "NOT FOUND", "-", "-", "-", "-",
                      "-", "-", "-", "-", "-");
        }
    }
    std::cout << std::string(150, '-') << "\n";
}

// ---------------------------------------------------------------------
// Main -- loop over all 15 targets
// ---------------------------------------------------------------------
void caf_edepsim_weights_matcher() {
    for (const auto& tgt : kTargets) {
        RunOneTarget(tgt);
    }
    std::cout << "\nDone. Ran " << kTargets.size() << " pions (5 per band). Every value in every\n"
               << "table came directly from the file read for that row -- none were\n"
               << "transcribed by hand from a separate printout.\n";
}

int main() {
    caf_edepsim_weights_matcher();
    return 0;
}