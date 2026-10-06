// caf_edepsim_weights_matcher_proton.C
//
// Proton version of caf_edepsim_weights_matcher.C. Picks ONE proton from
// the CAF and tracks it through EDepSim (official unpatched, patched) and
// the G4RW weight output files, printing one consolidated table -- every
// value returned directly from the read for that row, nothing hand-
// transcribed.
//
// Currently targeting the track that failed the low-cafFileIdx EDepSim
// lookup: cafFileIdx=1, RunId=10, EventId=16, TrackId=2 -- to see, in one
// place, whether CAF/official-EDepSim/weights agree with each other even
// though patched EDepSim doesn't have it.
//
// Matching convention -- same as pion, direct one-to-one, no offset:
//   File index: CAF's own cafFileIdx == EDepSim filename index
//               == weights file_idx
//   Interaction key: vertex_id = RunId*1e6 + EventId
//   Particle key: CAF G4ID == EDepSim TrackId == Weights track_id
//
// Setup before running:
//   shifter --image=mjkramer/sim2x2:ndlar011 --module=cvmfs /bin/bash
//   source /opt/generators/geant4/bin/geant4.sh
//   source /opt/generators/root/bin/thisroot.sh
//   export LD_LIBRARY_PATH=/opt/generators/edep-sim/install/lib:$LD_LIBRARY_PATH
//   export ROOT_INCLUDE_PATH=/opt/generators/edep-sim/install/include:$ROOT_INCLUDE_PATH
//
// Run with:
//   root -l -b -q caf_edepsim_weights_matcher_proton.C

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
// TARGET
// ---------------------------------------------------------------------
const int kCafFileIdx   = 1;
const int kTargetRunId  = 10;
const int kTargetEventId= 16;
const int kTargetTrackId= 2;

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
// working for pion in a prior run; only the filename itself uses the
// real index.
std::string OfficialUnpatchedPath(int idx) {
    return "/global/cfs/cdirs/dune/www/data/2x2/simulation/productions/"
           "MiniRun5_1E19_RHC/MiniRun5_1E19_RHC.spill/EDEPSIM_SPILLS/0000000/"
           "MiniRun5_1E19_RHC.spill." + Pad7(idx) + ".EDEPSIM_SPILLS.root";
}
std::string PatchedNewPath(int idx) {
    return "/global/cfs/cdirs/dune/users/mazam/EdepSim_G4RW_nu/"
           "patched_withRock_genieNu/"
           "MiniRun5_1E19_RHC.nu." + Pad7(idx) + ".G4RW.EDEPSIM.root";
}
std::string PatchedOldPath(int idx) {
    return "/global/cfs/cdirs/dune/users/mazam/EdepSim_G4RW_nu/"
           "patched_withRock_genieNu_OLD_pending_delete/"
           "MiniRun5_1E19_RHC.nu." + Pad7(idx) + ".G4RW.EDEPSIM.root";
}

const std::string kWeightsPatched =
    "/global/cfs/cdirs/dune/users/mazam/G4RW/output_weights/"
    "patched_withRock_genieNu/weights_1000_all_LAr_P_cascade_0_4000_5weights.root";
const std::string kWeightsUnpatched =
    "/global/cfs/cdirs/dune/users/mazam/G4RW/output_weights/"
    "unpatched_withRock_genieNu/weights_1000_all_LAr_P_cascade_0_4000_5weights.root";

// ---------------------------------------------------------------------
// Result structs
// ---------------------------------------------------------------------
struct CafResult {
    bool found = false;
    int cafFileIdx = -1;
    long long nu_id = -1;
    int G4ID = -1;
    int pdg = -999;
    double cosTheta = -999, length = -999;
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
    TTree* t = (TTree*)f->Get("syst_g4rw_proton_truth");
    if (!t) { std::cerr << "ERROR: no syst_g4rw_proton_truth tree\n"; f->Close(); return r; }

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

    int targetEventId2 = targetRunId * 1000000 + targetEventId;   // same overflow convention as pion

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
// Main
// ---------------------------------------------------------------------
void caf_edepsim_weights_matcher_proton() {

    CafResult caf = GetCAFEntry(kCafFileIdx, kTargetRunId, kTargetEventId, kTargetTrackId);

    EdepResult edepOfficial = SearchEdepSim(OfficialUnpatchedPath(kCafFileIdx), kTargetRunId, kTargetEventId, kTargetTrackId);
    EdepResult edepPatchedOld = SearchEdepSim(PatchedOldPath(kCafFileIdx),      kTargetRunId, kTargetEventId, kTargetTrackId);
    EdepResult edepPatchedNew = SearchEdepSim(PatchedNewPath(kCafFileIdx),      kTargetRunId, kTargetEventId, kTargetTrackId);

    WeightResult wPatched   = SearchWeights(kWeightsPatched,   kCafFileIdx, kTargetRunId, kTargetEventId, kTargetTrackId);
    WeightResult wUnpatched = SearchWeights(kWeightsUnpatched, kCafFileIdx, kTargetRunId, kTargetEventId, kTargetTrackId);

    std::cout << "\n================================================================\n";
    std::cout << " TARGET: cafFileIdx=" << kCafFileIdx << "  RunId=" << kTargetRunId
               << "  EventId=" << kTargetEventId << "  TrackId=" << kTargetTrackId << "\n";
    std::cout << "================================================================\n\n";

    PrintRow("Source", "fileIdx", "RunId", "EventId", "TrkId", "PDG", "nPts", "len[cm]",
              "w0.34", "w0.67", "w1.33", "w1.66", "w2.00");
    std::cout << std::string(150, '-') << "\n";

    PrintRow("CAF",
              I(caf.cafFileIdx), std::to_string(kTargetRunId), std::to_string(kTargetEventId), std::to_string(kTargetTrackId),
              caf.found ? I(caf.pdg) : "N/A", "-", caf.found ? F(caf.length) : "N/A",
              caf.found ? F(caf.w034,4) : "-", caf.found ? F(caf.w067,4) : "-",
              caf.found ? F(caf.w133,4) : "-", caf.found ? F(caf.w166,4) : "-",
              caf.found ? F(caf.w200,4) : "-");

    struct EdepRow { std::string label; EdepResult r; };
    std::vector<EdepRow> edepRows = {
        {"EDepSim: official unpatched", edepOfficial},
        {"EDepSim: patched OLD",        edepPatchedOld},
        {"EDepSim: patched NEW",        edepPatchedNew}
    };
    for (auto& er : edepRows) {
        if (er.r.found) {
            PrintRow(er.label, I(kCafFileIdx), I(er.r.RunId), I(er.r.EventId), I(er.r.TrackId),
                      I(er.r.pdg), I(er.r.nPoints), F(er.r.length),
                      "-", "-", "-", "-", "-");
        } else if (er.r.eventFound) {
            PrintRow(er.label, I(kCafFileIdx), I(er.r.RunId), I(er.r.EventId), "NOT FOUND",
                      "-", "-", "-", "-", "-", "-", "-", "-");
        } else {
            PrintRow(er.label, I(kCafFileIdx), "EVENT NOT FOUND", "-", "-", "-", "-", "-",
                      "-", "-", "-", "-", "-");
        }
    }

    struct WRow { std::string label; WeightResult r; };
    std::vector<WRow> wRows = {
        {"Weights: patched",   wPatched},
        {"Weights: unpatched", wUnpatched}
    };
    for (auto& wr : wRows) {
        if (wr.r.found) {
            PrintRow(wr.label, I(wr.r.file_idx), I(wr.r.run_id), I(wr.r.event_id), I(wr.r.track_id),
                      I(wr.r.pdg), "-", F(wr.r.len),
                      F(wr.r.w034,4), F(wr.r.w067,4), F(wr.r.w133,4), F(wr.r.w166,4), F(wr.r.w200,4));
        } else {
            PrintRow(wr.label, I(kCafFileIdx), "-", "NOT FOUND", "-", "-", "-", "-",
                      "-", "-", "-", "-", "-");
        }
    }

    std::cout << std::string(150, '-') << "\n";
    std::cout << "\nDone. Every value above came directly from the file read for that row --\n"
               << "none were transcribed by hand from a separate printout.\n";
}

int main() {
    caf_edepsim_weights_matcher_proton();
    return 0;
}