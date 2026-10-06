// validate_pion_example_6.c++
//
// Template incorporating all edits agreed in the g4rw_formula_trace
// validation thread:
//   - per-step (non-cumulative) weight printed alongside running weight
//   - step_weight_* TGraphs saved (in addition to running_weight_*)
//   - Eq. 3.10 exact interaction weight vs. code's W_I=f approximation,
//     both at the last step and propagated to the full-track weight
//
// Target track and file identifiers are set by the block below -- change
// ONLY that block to adapt this template to a different track.
//
// Compile:
//   g++ -o validate_pion_example_6 validate_pion_example_6.c++ \
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
//   ./validate_pion_example_6 reweight_PiPlus.fcl > validate_pion_example_6_output.txt 2>&1

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
#include "geant4reweight/src/ReweightBase/G4PiPlusReweighter.hh"
#include "geant4reweight/src/ReweightBase/G4ReweightTraj.hh"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <map>

// ---------------------------------------------------------------------
// TARGET TRACK -- the only block that changes between example scripts.
// ---------------------------------------------------------------------
const int    kTargetRunId   = 9256;
const int    kTargetEventId = 2;
const int    kTargetTrackId = 3;
const std::string kPatchedFile =
    "/global/cfs/cdirs/dune/users/mazam/EdepSim_G4RW_nu/patched_withRock_genieNu/MiniRun5_1E19_RHC.nu.0000925.G4RW.EDEPSIM.root";
const std::string kUnpatchedFile =
    "/global/cfs/cdirs/dune/users/mazam/EdepSim_G4RW_nu/unpatched_withRock_genieNu/MiniRun5_1E19_RHC.nu.0000925.G4RW.EDEPSIM.root";
const std::string kWeightFile =
    "/global/cfs/cdirs/dune/users/mazam/G4RW/output_weights/patched_withRock_genieNu/weights_1000_all_LAr_PiPlus_cascade_0_9000_5weights.root";

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
            throw cet::exception("validate_pion_example") << "Could not find " << filename;
        theFile = new TFile(found_filename.c_str());
    }
    if (!theFile || theFile->IsZombie() || !theFile->IsOpen())
        throw cet::exception("validate_pion_example") << "Could not open " << filename;
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
            if (traj.GetPDGCode() != 211) continue;
            *outEvent = eventPtr;
            *outTraj = &traj;
            *outFile = f;
            return true;
        }
    }
    return false;
}

int main(int argc, char** argv) {

    std::string fclFile = (argc > 1) ? argv[1] : "reweight_PiPlus.fcl";

    std::cout << "================================================================\n";
    std::cout << " VALIDATE PION EXAMPLE\n";
    std::cout << " Target: File 0000925, RunId=" << kTargetRunId
               << ", EventId=" << kTargetEventId
               << ", TrackId=" << kTargetTrackId << ", PDG=211\n";
    std::cout << "================================================================\n";

    auto fclPars   = MakeFCLPars(fclFile);
    auto material  = fclPars.get<fhicl::ParameterSet>("Material");
    auto fracsName = fclPars.get<std::string>("FracsFile");
    TFile* fracsFile = OpenFile(fracsName);

    G4ReweightManager* rwManager = new G4ReweightManager({material});

    std::map<std::string, TH1D*> dummyScales;
    for (auto name : {"inel", "cex", "abs", "dcex", "prod"}) {
        TH1D* h = new TH1D(name, "", 1, 0., 5000.);
        h->SetBinContent(1, 1.0);
        dummyScales[name] = h;
    }
    G4PiPlusReweighter* reweighter = new G4PiPlusReweighter(
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
    G4ReweightTraj rwTraj(kTargetTrackId, 211, -1, kTargetRunId, {0,0});
    for (const auto& sec : evP->Trajectories) {
        if (sec.GetParentId() != kTargetTrackId) continue;
        int pdg = sec.GetPDGCode();
        double p0 = sec.Points.empty() ? -1 : sec.Points.front().GetMomentum().Mag();
        std::cout << "  TrackId=" << sec.GetTrackId() << "  PDG=" << pdg
                   << "  initial |p|=" << p0 << " MeV\n";
        rwTraj.AddChild(new G4ReweightTraj(sec.GetTrackId(), pdg, kTargetTrackId, kTargetRunId, {0,0}));
    }

    std::string channel = reweighter->GetInteractionSubtype(rwTraj);
    std::cout << "\nClassification (from G4RW's own GetInteractionSubtype()): \"" << channel << "\"\n";

    // ---- UNPATCHED trajectory ----
    std::cout << "\n================================================================\n";
    std::cout << " UNPATCHED trajectory (same RunId/EventId/TrackId)\n";
    std::cout << "================================================================\n";

    TG4Event* evU = nullptr; const TG4Trajectory* trajU = nullptr; TFile* fU = nullptr;
    std::vector<double> unpatchedMomenta;
    std::vector<double> unpatchedLengths;
    if (!findTrack(kUnpatchedFile, kTargetRunId, kTargetEventId, kTargetTrackId, &evU, &trajU, &fU)) {
        std::cout << "NOTE: matching track not found in unpatched sample.\n";
    } else {
        double dummyLen; std::vector<double> dummyMom;
        printTrajectoryTable(*trajU, dummyLen, dummyMom);
        for (size_t j = 1; j < trajU->Points.size(); ++j) {
            const auto& prev = trajU->Points[j-1];
            const auto& pt   = trajU->Points[j];
            auto d = pt.GetPosition() - prev.GetPosition();
            unpatchedLengths.push_back(d.Vect().Mag() / 10.0);
            unpatchedMomenta.push_back(prev.GetMomentum().Mag());
        }
    }

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
    for (double p : unpatchedMomenta)      xsecCache[p] = reweighter->GetInelasticXSec(p, false);
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
    std::cout << std::setw(10) << "f" << std::setw(18) << "patched w" << std::setw(18) << "unpatched w\n";
    for (double f : scaleFactors) {
        double wP = computeWeight(patchedPreStepMomenta, patchedStepLengths, f);
        double wU = unpatchedMomenta.empty() ? -1 : computeWeight(unpatchedMomenta, unpatchedLengths, f);
        std::cout << std::setw(10) << f << std::setw(18) << wP << std::setw(18) << wU << "\n";
    }

    // ---- Running + per-step weight vs. cumulative length, saved as TGraphs ----
    std::cout << "\n================================================================\n";
    std::cout << " RUNNING + STEP WEIGHT vs CUMULATIVE LENGTH (patched and unpatched)\n";
    std::cout << "================================================================\n";

    TFile* fOutGraphs = new TFile("validate_pion_example_6.root", "RECREATE");

    auto saveRunningWeightGraphs = [&](const std::vector<double>& momenta,
                                        const std::vector<double>& lengths,
                                        const std::string& sampleTag) {
        std::cout << "\n--- " << sampleTag << " (" << momenta.size() << " step(s)) ---\n";
        for (double f : scaleFactors) {
            std::vector<double> xVals, yVals, stepWVals;
            double S = 0., cumLen = 0.;
            std::cout << "\n-- f = " << f << " --\n";
            std::cout << std::setw(14) << "cum.len[cm]" << std::setw(14) << "step w" << std::setw(14) << "running w\n";
            for (size_t i = 0; i < momenta.size(); ++i) {
                bool isLastStep = (i == momenta.size() - 1);
                double stepW = std::exp((1.0 - f) * lengths[i] * xsecCache[momenta[i]]);
                if (isLastStep) stepW *= f;   // W_I = f, the code's approximation (see Eq. 3.10 block below)

                S += lengths[i] * xsecCache[momenta[i]];
                cumLen += lengths[i];
                double w = isLastStep ? std::exp((1.0 - f) * S) * f : std::exp((1.0 - f) * S);
                xVals.push_back(cumLen);
                yVals.push_back(w);
                stepWVals.push_back(stepW);
                std::cout << std::setw(14) << cumLen << std::setw(14) << stepW << std::setw(14) << w << "\n";
            }
            TGraph* g = new TGraph((int)xVals.size(), xVals.data(), yVals.data());
            std::string label = "f" + std::to_string(f).substr(0, std::to_string(f).find('.') + 3);
            g->SetName(("running_weight_" + sampleTag + "_" + label).c_str());
            g->SetTitle((sampleTag + " running weight vs cumulative length, f=" + std::to_string(f)).c_str());
            g->GetXaxis()->SetTitle("Cumulative Track Length [cm]");
            g->GetYaxis()->SetTitle("G4RW Weight");
            g->Write();

            TGraph* gStep = new TGraph((int)xVals.size(), xVals.data(), stepWVals.data());
            gStep->SetName(("step_weight_" + sampleTag + "_" + label).c_str());
            gStep->SetTitle((sampleTag + " individual step weight vs cumulative length, f=" + std::to_string(f)).c_str());
            gStep->GetXaxis()->SetTitle("Cumulative Track Length [cm]");
            gStep->GetYaxis()->SetTitle("G4RW Step Weight");
            gStep->Write();
        }
    };

    saveRunningWeightGraphs(patchedPreStepMomenta, patchedStepLengths, "patched");
    if (!unpatchedMomenta.empty()) {
        saveRunningWeightGraphs(unpatchedMomenta, unpatchedLengths, "unpatched");
    } else {
        std::cout << "\nNOTE: no unpatched momenta available -- skipping unpatched running-weight graphs.\n";
    }

    fOutGraphs->Close();
    std::cout << "\nWrote validate_pion_example_6.root (running_weight_*, step_weight_* TGraphs, patched + unpatched)\n";

    // ================================================================
    // EQ. 3.10 EXACT INTERACTION WEIGHT vs CODE'S W_I = f APPROXIMATION
    // (patched track, interacting = last step)
    // ================================================================
    std::cout << "\n================================================================\n";
    std::cout << " EQ. 3.10 EXACT INTERACTION WEIGHT vs CODE'S W_I = f APPROXIMATION\n";
    std::cout << " (patched track, interacting = last step)\n";
    std::cout << "================================================================\n";

    size_t lastIdx = patchedPreStepMomenta.size() - 1;
    double deltaL_last = patchedStepLengths[lastIdx];       // <--- Delta L_i for Eq. 3.10
    double p_last      = patchedPreStepMomenta[lastIdx];    // <--- momentum at which sigma_I,i, sigma_E,i are evaluated

    double sigma_I_nom = xsecCache[p_last];                          // <--- sigma_I,i for Eq. 3.10
    double sigma_E_nom = reweighter->GetElasticXSec(p_last, false);  // <--- sigma_E,i for Eq. 3.10
    double sigma_X_nom = 0.0;                                        // <--- sigma_X,i for Eq. 3.10 -- NOT AVAILABLE in G4RW

    double sigma_nom = sigma_I_nom + sigma_E_nom + sigma_X_nom;      // <--- sigma_i (total, denominator of Eq. 3.10)

    std::cout << "  Delta L_i (last step length)       = " << deltaL_last << " cm   (<--- Delta L_i for Eq. 3.10)\n";
    std::cout << "  p_i (pre-step momentum, last step) = " << p_last << " MeV  (<--- momentum at which sigma_I,i, sigma_E,i are evaluated)\n";
    std::cout << "  sigma_I,i (nominal inelastic)       = " << sigma_I_nom << " /cm  (<--- sigma_I,i for Eq. 3.10)\n";
    std::cout << "  sigma_E,i (nominal elastic)         = " << sigma_E_nom << " /cm  (<--- sigma_E,i for Eq. 3.10)\n";
    std::cout << "  sigma_X,i (Coulomb+decay)           = " << sigma_X_nom << " /cm  (<--- sigma_X,i for Eq. 3.10 -- NOT MODELED IN G4RW, set to 0)\n";
    std::cout << "  sigma_i (= sigma_I+sigma_E+sigma_X)  = " << sigma_nom << " /cm  (<--- sigma_i, denominator of Eq. 3.10)\n";

    std::cout << "\n" << std::setw(6) << "f" << std::setw(16) << "sigma'_I,i" << std::setw(16) << "sigma'_i"
               << std::setw(18) << "W_I approx (=f)" << std::setw(18) << "W_I exact (3.10)"
               << std::setw(14) << "% diff\n";

    for (double f : scaleFactors) {
        double sigma_I_prime = f * sigma_I_nom;
        double sigma_E_prime = sigma_E_nom;
        double sigma_prime   = sigma_I_prime + sigma_E_prime + sigma_X_nom;

        double W_I_approx = f;
        double W_I_exact  = f * (sigma_nom / sigma_prime)
                           * (1.0 - std::exp(-deltaL_last * sigma_prime))
                           / (1.0 - std::exp(-deltaL_last * sigma_nom));

        double pctDiff = 100.0 * (W_I_exact - W_I_approx) / W_I_approx;

        std::cout << std::setw(6) << f << std::setw(16) << sigma_I_prime << std::setw(16) << sigma_prime
                   << std::setw(18) << W_I_approx << std::setw(18) << W_I_exact
                   << std::setw(13) << pctDiff << "%\n";
    }

    std::cout << "\n-- Full-track weight: code's w vs. code's w with Eq. 3.10 swapped in --\n";
    std::cout << std::setw(6) << "f" << std::setw(18) << "w_code (=f approx)" << std::setw(20) << "w_with_eq310"
               << std::setw(14) << "% diff\n";
    for (double f : scaleFactors) {
        double S_total = 0.;
        for (size_t i = 0; i < patchedPreStepMomenta.size(); ++i)
            S_total += patchedStepLengths[i] * xsecCache[patchedPreStepMomenta[i]];

        double sigma_I_prime = f * sigma_I_nom;
        double sigma_prime   = sigma_I_prime + sigma_E_nom + sigma_X_nom;
        double W_I_exact = f * (sigma_nom / sigma_prime)
                          * (1.0 - std::exp(-deltaL_last * sigma_prime))
                          / (1.0 - std::exp(-deltaL_last * sigma_nom));

        double w_code       = std::exp((1.0 - f) * S_total) * f;
        double w_with_eq310 = std::exp((1.0 - f) * S_total) * W_I_exact;
        double pctDiff = 100.0 * (w_with_eq310 - w_code) / w_code;

        std::cout << std::setw(6) << f << std::setw(18) << w_code << std::setw(20) << w_with_eq310
                   << std::setw(13) << pctDiff << "%\n";
    }

    // ---- Compare to saved production weight ----
    std::cout << "\n================================================================\n";
    std::cout << " SAVED WEIGHTS IN PRODUCTION OUTPUT FILE\n";
    std::cout << "================================================================\n";
    TFile* fW = TFile::Open(kWeightFile.c_str());
    TTree* wtree = (TTree*)fW->Get("g4rwEvt");
    int file_idx, event_id, track_id, pdg;
    double weight_1p33, weight_0p67, weight_1p66, weight_0p34, weight_2p00, len;
    wtree->SetBranchAddress("file_idx", &file_idx);
    wtree->SetBranchAddress("event_id", &event_id);
    wtree->SetBranchAddress("track_id", &track_id);
    wtree->SetBranchAddress("pdg", &pdg);
    wtree->SetBranchAddress("weight_1p33", &weight_1p33);
    wtree->SetBranchAddress("weight_0p67", &weight_0p67);
    wtree->SetBranchAddress("weight_1p66", &weight_1p66);
    wtree->SetBranchAddress("weight_0p34", &weight_0p34);
    wtree->SetBranchAddress("weight_2p00", &weight_2p00);
    wtree->SetBranchAddress("len", &len);

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
