// validate_proton_example_5.c++
//
// Proton analog of validate_pion_example_N.c++. Mirrors its structure
// exactly: opens the patched EDepSim trajectory directly, classifies via
// the reweighter's OWN GetInteractionSubtype() (not a hand-rolled
// function), reads real Sigma_I(p) at every step, computes the reduced-
// formula prediction, and compares against the saved production weight.
//
// *** THINGS INFERRED BY PARALLEL STRUCTURE, NOT INDEPENDENTLY CONFIRMED ***
// (flagged inline below with "UNCONFIRMED" -- check these first if this
// fails to compile or gives nonsensical output):
//   1. G4ProtonReweighter's constructor signature is assumed identical
//      to G4PiPlusReweighter's (FSInput, FSScales, material_pars,
//      rw_manager, inputElasticBiasHist, fix). Confirmed to EXIST as a
//      real custom-patched class (README_geant4reweight_custom_install.md),
//      but its exact constructor args were not directly inspected here.
//   2. The five exclusive channel KEYS are assumed to be
//      {"0n0p","0n1p","1n0p","1n1p","Other"} (matching reweight_proton.fcl's
//      older sibling reweight_P.fcl Cut names) -- "0n0p" is directly
//      confirmed present in G4ProtonReweighter.cc per the install
//      README's own grep check; the other four are inferred by pattern.
//   3. GetInteractionSubtype() is assumed to exist on G4ProtonReweighter
//      with the same signature as G4PiPlusReweighter's -- supported by
//      the README's note "Extended proton channel classification to
//      produce exclusive channels" but not directly viewed.
//   4. GetInelasticXSec()/GetElasticXSec() are base-class (G4Reweighter)
//      methods, so these SHOULD be available unchanged regardless of
//      species -- this one is on firmer ground than 1-3 above.
//
// fInelastic tag string ("protonInelastic") IS directly confirmed from
// do_weighting_v2_protons.c__ source (last_is_inel check, same 4/121
// numeric criterion as pion, just a different tag name).
//
// Target track and file identifiers are set by the block below -- change
// ONLY that block to adapt this template to a different track.
//
// Compile:
//   g++ -o validate_proton_example_5 validate_proton_example_5.c++ \
//     -I ${GEANT4REWEIGHT_INC} \
//     -I ${CETLIB_INC} -I ${CETLIB_EXCEPT_INC} -I ${FHICLCPP_INC} \
//     -I ${BOOST_INC} \
//     -I ${EDEPSIM_INC} -L ${EDEPSIM_LIB} -ledepsim_io -I ${ROOT_INC} \
//     -L ${GEANT4REWEIGHT_LIB} -lReweightBaseLib \
//     -L ${FHICLCPP_LIB} -lfhiclcpp \
//     -L ${CETLIB_EXCEPT_LIB} -lcetlib_except -L ${CETLIB_LIB} -lcetlib \
//     -L ${BOOST_LIB} -lboost_program_options \
//     -L ${ROOTSYS}/lib \
//     -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree \
//     -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread \
//     -lMultiProc -lROOTDataFrame -pthread -lm -ldl -rdynamic \
//     --std=c++17
//
// Run:
//   ./validate_proton_example_5 reweight_proton.fcl > validate_proton_example_5_output.txt 2>&1

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1D.h"
#include "EDepSim/TG4Event.h"
#include "EDepSim/TG4Trajectory.h"

#include "cetlib/filepath_maker.h"
#include "cetlib/filesystem.h"
#include "cetlib/search_path.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/intermediate_table.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "fhiclcpp/ParameterSet.h"

#include "geant4reweight/src/ReweightBase/G4ReweightManager.hh"
#include "geant4reweight/src/ReweightBase/G4ProtonReweighter.hh"   // UNCONFIRMED header path/name -- adjust if different
#include "geant4reweight/src/ReweightBase/G4ReweightTraj.hh"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <map>

// ---------------------------------------------------------------------
// TARGET TRACK -- the only block that changes between example scripts.
// ---------------------------------------------------------------------
const int    kTargetRunId   = 7783;
const int    kTargetEventId = 95;
const int    kTargetTrackId = 3;
const std::string kPatchedFile =
    "/global/cfs/cdirs/dune/users/mazam/EdepSim_G4RW_nu/patched_withRock_genieNu/MiniRun5_1E19_RHC.nu.0000778.G4RW.EDEPSIM.root";
const std::string kWeightFile =
    "/global/cfs/cdirs/dune/users/mazam/G4RW/output_weights/patched_withRock_genieNu/weights_1000_all_LAr_P_cascade_0_4000_5weights.root";

fhicl::ParameterSet MakeFCLPars(std::string fcl_file) {
    char const* fhicl_env = getenv("FHICL_FILE_PATH");
    std::string search_path = fhicl_env ? std::string{fhicl_env} : ".";
    cet::filepath_first_absolute_or_lookup_with_dot lookupPolicy{search_path};
    return fhicl::ParameterSet::make(fcl_file, lookupPolicy);
}

TFile* OpenFile(const std::string filename) {
    TFile* theFile = nullptr;
    if (cet::file_exists(filename)) {
        theFile = new TFile(filename.c_str());
    } else {
        cet::search_path sp{"FW_SEARCH_PATH"};
        std::string found_filename;
        if (!sp.find_file(filename, found_filename))
            throw cet::exception("validate_proton_example") << "Could not find " << filename;
        theFile = new TFile(found_filename.c_str());
    }
    if (!theFile || theFile->IsZombie() || !theFile->IsOpen())
        throw cet::exception("validate_proton_example") << "Could not open " << filename;
    return theFile;
}

// ---------------------------------------------------------------------
// Print one trajectory's step table.
// ---------------------------------------------------------------------
void printTrajectoryTable(const TG4Trajectory &traj, double &totalLen,
                           std::vector<double> &stepMomenta) {
    totalLen = 0.;
    stepMomenta.clear();

    std::cout << std::setw(5)  << "step"
              << std::setw(12) << "x[mm]" << std::setw(12) << "y[mm]" << std::setw(12) << "z[mm]"
              << std::setw(10) << "stepL[cm]" << std::setw(10) << "|p|[MeV]"
              << std::setw(8)  << "proc" << std::setw(10) << "subproc" << "\n";
    std::cout << std::string(80, '-') << "\n";

    for (size_t j = 0; j < traj.Points.size(); ++j) {
        const auto &pt = traj.Points[j];
        double stepL = 0.;
        if (j > 0) {
            const auto &prev = traj.Points[j-1];
            auto d = pt.GetPosition() - prev.GetPosition();
            stepL = d.Vect().Mag() / 10.0;
            totalLen += stepL;
        }
        double pmag = pt.GetMomentum().Mag();
        std::cout << std::setw(5)  << j
                  << std::setw(12) << std::fixed << std::setprecision(1) << pt.GetPosition().X()
                  << std::setw(12) << pt.GetPosition().Y()
                  << std::setw(12) << pt.GetPosition().Z()
                  << std::setw(10) << std::setprecision(2) << stepL
                  << std::setw(10) << std::setprecision(1) << pmag
                  << std::setw(8)  << pt.GetProcess()
                  << std::setw(10) << pt.GetSubprocess() << "\n";
        stepMomenta.push_back(pmag);
    }
    std::cout << "Total length (sum of steps): " << std::setprecision(2) << totalLen << " cm\n";
}

// ---------------------------------------------------------------------
// Find a specific (RunId, EventId, TrackId) primary trajectory directly.
// ---------------------------------------------------------------------
bool findTrack(const std::string& filename, int targetRunId, int targetEventId, int targetTrackId,
               TG4Event** outEvent, const TG4Trajectory** outTraj, TFile** outFile) {
    TFile* f = TFile::Open(filename.c_str());
    if (!f || f->IsZombie()) return false;
    TTree* tree = (TTree*)f->Get("EDepSimEvents");
    TG4Event* eventPtr = nullptr;
    tree->SetBranchAddress("Event", &eventPtr);

    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        if (eventPtr->RunId != targetRunId) continue;
        if (eventPtr->EventId != targetEventId) continue;
        for (const auto& traj : eventPtr->Trajectories) {
            if (traj.GetTrackId() != targetTrackId) continue;
            if (traj.GetPDGCode() != 2212) continue;
            *outEvent = eventPtr;
            *outTraj = &traj;
            *outFile = f;
            return true;
        }
    }
    return false;
}

int main(int argc, char** argv) {

    std::string fclFile = (argc > 1) ? argv[1] : "reweight_proton.fcl";

    std::cout << "================================================================\n";
    std::cout << " VALIDATE PROTON EXAMPLE\n";
    std::cout << " Target: File 0000778, RunId=" << kTargetRunId
               << ", EventId=" << kTargetEventId
               << ", TrackId=" << kTargetTrackId << ", PDG=2212\n";
    std::cout << "================================================================\n";

    auto fclPars   = MakeFCLPars(fclFile);
    auto material  = fclPars.get<fhicl::ParameterSet>("Material");
    auto fracsName = fclPars.get<std::string>("FracsFile");
    TFile* fracsFile = OpenFile(fracsName);

    G4ReweightManager* rwManager = new G4ReweightManager({material});

    // UNCONFIRMED: exclusive channel key set for proton. "0n0p" is
    // directly confirmed to exist in G4ProtonReweighter.cc; the other
    // four are inferred from reweight_P.fcl's Cut names by pattern.
    std::map<std::string, TH1D*> dummyScales;
    for (auto name : {"total"}) {
        TH1D* h = new TH1D(name, "", 1, 0., 5000.);
        h->SetBinContent(1, 1.0);
        dummyScales[name] = h;
    }
    // UNCONFIRMED: constructor signature assumed identical to G4PiPlusReweighter's.
    G4ProtonReweighter* reweighter = new G4ProtonReweighter(
        fracsFile, dummyScales, material, rwManager, nullptr, false);

    // ---- PATCHED trajectory ----
    std::cout << "\n================================================================\n";
    std::cout << " PATCHED trajectory\n";
    std::cout << "================================================================\n";

    TG4Event* evP = nullptr; const TG4Trajectory* trajP = nullptr; TFile* fP = nullptr;
    if (!findTrack(kPatchedFile, kTargetRunId, kTargetEventId, kTargetTrackId, &evP, &trajP, &fP)) {
        std::cout << "ERROR: track not found in patched file.\n";
        return 1;
    }
    double patchedLen;
    std::vector<double> patchedMomenta;
    printTrajectoryTable(*trajP, patchedLen, patchedMomenta);

    std::cout << "\nSecondaries (ParentId == " << kTargetTrackId << "):\n";
    G4ReweightTraj rwTraj(kTargetTrackId, 2212, -1, kTargetRunId, {0,0});
    for (const auto& sec : evP->Trajectories) {
        if (sec.GetParentId() != kTargetTrackId) continue;
        int pdg = sec.GetPDGCode();
        double p0 = sec.Points.empty() ? -1 : sec.Points.front().GetMomentum().Mag();
        std::cout << "  TrackId=" << sec.GetTrackId() << "  PDG=" << pdg
                   << "  initial |p|=" << p0 << " MeV\n";
        rwTraj.AddChild(new G4ReweightTraj(sec.GetTrackId(), pdg, kTargetTrackId, kTargetRunId, {0,0}));
    }

    // UNCONFIRMED: assumes GetInteractionSubtype() exists on G4ProtonReweighter
    // with the same signature/behavior as G4PiPlusReweighter's.
    std::string channel = reweighter->GetInteractionSubtype(rwTraj);
    std::cout << "\nClassification (from reweighter's own GetInteractionSubtype()): \"" << channel << "\"\n";

    // ---- Sigma_I(p) at every PATCHED step momentum ----
    std::cout << "\n================================================================\n";
    std::cout << " Sigma_I_nominal(p) -- read directly via reweighter->GetInelasticXSec()\n";
    std::cout << "================================================================\n";

    std::vector<double> patchedStepLengths;
    std::vector<double> patchedPreStepMomenta;
    for (size_t j = 1; j < trajP->Points.size(); ++j) {
        const auto& prev = trajP->Points[j-1];
        const auto& pt   = trajP->Points[j];
        auto d = pt.GetPosition() - prev.GetPosition();
        patchedStepLengths.push_back(d.Vect().Mag() / 10.0);
        patchedPreStepMomenta.push_back(prev.GetMomentum().Mag());
    }

    std::map<double, double> xsecCache;
    for (double p : patchedPreStepMomenta) xsecCache[p] = reweighter->GetInelasticXSec(p, false);
    for (auto& kv : xsecCache) {
        std::cout << "  p=" << std::setprecision(6) << kv.first
                   << " MeV  ->  Sigma_I_nom=" << kv.second << " 1/cm\n";
    }

    // ---- Predicted weight (reduced formula) ----
    std::vector<double> scaleFactors = {0.34, 0.67, 1.33, 1.66, 2.00};

    auto computeWeight = [&](const std::vector<double>& momenta,
                              const std::vector<double>& lengths, double f) {
        double S = 0.;
        for (size_t i = 0; i < momenta.size(); ++i) S += lengths[i] * xsecCache[momenta[i]];
        return std::exp((1.0 - f) * S) * f;
    };

    std::cout << "\n================================================================\n";
    std::cout << " PREDICTED WEIGHTS (reduced 'reac' formula, using real Sigma_I)\n";
    std::cout << "================================================================\n";
    std::cout << std::setw(10) << "f" << std::setw(18) << "patched w\n";
    for (double f : scaleFactors) {
        double wP = computeWeight(patchedPreStepMomenta, patchedStepLengths, f);
        std::cout << std::setw(10) << f << std::setw(18) << wP << "\n";
    }

    // ---- Per-step and running weight vs. cumulative length ----
    std::cout << "\n================================================================\n";
    std::cout << " STEP + RUNNING WEIGHT vs CUMULATIVE LENGTH\n";
    std::cout << "================================================================\n";

    TFile* fOutGraphs = new TFile("validate_proton_example_5.root", "RECREATE");

    for (double f : scaleFactors) {
        std::vector<double> xVals, yVals, stepWVals;
        double S = 0., cumLen = 0.;
        std::cout << "\n-- f = " << f << " --\n";
        std::cout << std::setw(14) << "cum.len[cm]" << std::setw(14) << "step w" << std::setw(14) << "running w\n";
        for (size_t i = 0; i < patchedPreStepMomenta.size(); ++i) {
            bool isLastStep = (i == patchedPreStepMomenta.size() - 1);
            double stepW = std::exp((1.0 - f) * patchedStepLengths[i] * xsecCache[patchedPreStepMomenta[i]]);
            if (isLastStep) stepW *= f;

            S += patchedStepLengths[i] * xsecCache[patchedPreStepMomenta[i]];
            cumLen += patchedStepLengths[i];
            double w = isLastStep ? std::exp((1.0 - f) * S) * f : std::exp((1.0 - f) * S);
            xVals.push_back(cumLen);
            yVals.push_back(w);
            stepWVals.push_back(stepW);
            std::cout << std::setw(14) << cumLen << std::setw(14) << stepW << std::setw(14) << w << "\n";
        }
        std::string label = "f" + std::to_string(f).substr(0, std::to_string(f).find('.') + 3);
        TGraph* g = new TGraph((int)xVals.size(), xVals.data(), yVals.data());
        g->SetName(("running_weight_patched_" + label).c_str());
        g->GetXaxis()->SetTitle("Cumulative Track Length [cm]");
        g->GetYaxis()->SetTitle("G4RW Weight");
        g->Write();

        TGraph* gStep = new TGraph((int)xVals.size(), xVals.data(), stepWVals.data());
        gStep->SetName(("step_weight_patched_" + label).c_str());
        gStep->GetXaxis()->SetTitle("Cumulative Track Length [cm]");
        gStep->GetYaxis()->SetTitle("G4RW Step Weight");
        gStep->Write();
    }
    fOutGraphs->Close();
    std::cout << "\nWrote validate_proton_example_5.root (running_weight_*, step_weight_* TGraphs, patched)\n";

    // ---- Compare to saved production weight ----
    std::cout << "\n================================================================\n";
    std::cout << " SAVED WEIGHTS IN PRODUCTION OUTPUT FILE\n";
    std::cout << "================================================================\n";
    TFile* fW = TFile::Open(kWeightFile.c_str());
    TTree* wtree = (TTree*)fW->Get("g4rwEvt");
    int file_idx, run_id, event_id, track_id, pdg;
    double weight_1p33, weight_0p67, weight_1p66, weight_0p34, weight_2p00, len;
    wtree->SetBranchAddress("file_idx", &file_idx);
    wtree->SetBranchAddress("run_id",   &run_id);
    wtree->SetBranchAddress("event_id", &event_id);
    wtree->SetBranchAddress("track_id", &track_id);
    wtree->SetBranchAddress("pdg",      &pdg);
    wtree->SetBranchAddress("weight_1p33", &weight_1p33);
    wtree->SetBranchAddress("weight_0p67", &weight_0p67);
    wtree->SetBranchAddress("weight_1p66", &weight_1p66);
    wtree->SetBranchAddress("weight_0p34", &weight_0p34);
    wtree->SetBranchAddress("weight_2p00", &weight_2p00);
    wtree->SetBranchAddress("len",      &len);

    // event_id is written as a 32-bit int = run_id*1e6 + EventId in
    // do_weighting_v2_protons.c__ (same convention/overflow risk as pion
    // -- decode using run_id, same as caf_edepsim_weights_matcher.C).
    int targetEventIdCombined = kTargetRunId * 1000000 + kTargetEventId;
    bool foundWeight = false;
    for (Long64_t i = 0; i < wtree->GetEntries(); ++i) {
        wtree->GetEntry(i);
        if (event_id == targetEventIdCombined && track_id == kTargetTrackId) {
            std::cout << "  weight_0p34 = " << weight_0p34 << "\n";
            std::cout << "  weight_0p67 = " << weight_0p67 << "\n";
            std::cout << "  weight_1p33 = " << weight_1p33 << "\n";
            std::cout << "  weight_1p66 = " << weight_1p66 << "\n";
            std::cout << "  weight_2p00 = " << weight_2p00 << "\n";
            std::cout << "  len (saved) = " << len << " cm\n";
            foundWeight = true;
            break;
        }
    }
    if (!foundWeight) std::cout << "WARNING: no matching saved weight found.\n";

    std::cout << "\nDone.\n";
    return 0;
}
