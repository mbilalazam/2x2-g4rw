///////////////////////////////////////////////////////////////////////////////
// plot_caf_weight_vs_length_2d.C
//
// 2D histogram: G4RW weight vs. pion track length, one plot per scale
// factor, for BOTH truth length and reco length, using caf_pion_list.root
// (produced by dump_caf_pion_list.C, joined against syst_g4rw_pion_reco).
// Every point plotted is a CAF-matched pi+ with a guaranteed-traceable
// weight.
//
// Motivation: recent findings show a subset of tracks (long length, pion
// escapes the detector before any point is tagged proc=4/subproc=121)
// get NO interaction-kick factor applied at all -- their weight is pure
// e^{(1-f)S} survival. The truth-length plot shows that population as a
// visually distinct low-weight band at large length. The reco-length
// plot is a population-level test of WHY: if that band corresponds to
// pions leaving the (re)constructible volume, reco length should be
// truncated relative to truth length for those same tracks, and the
// band should compress or thin out here relative to the truth plot.
// Reco-length points require hasReco==1 (dump_caf_pion_list.C sets this
// via the nu_id/cafFileIdx/G4ID join against syst_g4rw_pion_reco).
//
// Output (10 plots total, 5 scale factors x {truth, reco}):
//   plots/caf_weight_vs_length_2d_<f>.png       (truth length)
//   plots/caf_weight_vs_recolength_2d_<f>.png   (reco length)
//
// Compile and run:
//   g++ -std=c++17 -O2 -g plot_caf_weight_vs_length_2d.C -o plot_caf_weight_vs_length_2d `root-config --cflags --libs`
//   ./plot_caf_weight_vs_length_2d
///////////////////////////////////////////////////////////////////////////////

#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"

#include <array>
#include <string>
#include <iostream>
#include <cmath>

const std::string kInFile  = "caf_pion_list_eq3p10.root";
const std::string kOutDir  = "2d_plots_eq3p10";

const std::array<std::string,5> kBranchNames = {"w_0p34", "w_0p67", "w_1p33", "w_1p66", "w_2p00"};
const std::array<std::string,5> kLabels      = {"f = 0.34", "f = 0.67", "f = 1.33", "f = 1.66", "f = 2.00"};
const std::array<std::string,5> kFileTags    = {"0p34", "0p67", "1p33", "1p66", "2p00"};

const int nBinsX = 80, nBinsY = 100;
const double xMin = 0.0, xMax = 400.0;
const double yMin = 0.0, yMax = 2.0;

// Same +-55 cm fiducial cube used in purity_dlp_multiplicity.cc's
// contained_geom cut and in plot_caf_weight_vs_endpos.C.
const double kFidBound = 55.0;

void DrawDuneLabel() {
    TLatex tL; tL.SetNDC();
    tL.DrawLatex(0.13, 0.94, "#bf{DUNE:ND-LAr 2x2}");
}

// ---------------------------------------------------------------------
// Fills and draws all 5 scale-factor histograms for one length variable.
// lengthKind is used only for titles/filenames -- "Truth" or "Reco".
// containSel: -1 = no containment cut, 0 = non-contained only,
//             1 = contained only (contained = |endX|,|endY|,|endZ| < kFidBound,
//             using TRUTH end position regardless of useReco -- containment
//             is a property of the pion's actual trajectory, not of what
//             reco happened to measure).
// ---------------------------------------------------------------------
void MakeWeightVsLengthPlots(TTree* t, bool useReco, const std::string& lengthKind,
                              const std::string& xAxisTitle, const std::string& fileTag,
                              int containSel = -1) {

    int pdg, hasReco = 1;   // default 1 so the (!useReco) branch never gates anything
    double length = -999, recoLength = -999;
    double endX = -9999, endY = -9999, endZ = -9999;
    double w[5];
    t->SetBranchAddress("pdg",    &pdg);
    t->SetBranchAddress("length", &length);
    for (int k = 0; k < 5; ++k) t->SetBranchAddress(kBranchNames[k].c_str(), &w[k]);
    if (useReco) {
        t->SetBranchAddress("recoLength", &recoLength);
        t->SetBranchAddress("hasReco",    &hasReco);
    }
    bool hasEndXYZ = (containSel != -1) &&
                     t->GetBranch("endX") && t->GetBranch("endY") && t->GetBranch("endZ");
    if (containSel != -1) {
        if (hasEndXYZ) {
            t->SetBranchAddress("endX", &endX);
            t->SetBranchAddress("endY", &endY);
            t->SetBranchAddress("endZ", &endZ);
        } else {
            std::cerr << "WARNING: containment cut requested but endX/endY/endZ branches not found --\n"
                       << "skipping " << lengthKind << " (containSel=" << containSel << ") plots.\n";
            return;
        }
    }

    std::array<TH2D*,5> hists;
    for (int k = 0; k < 5; ++k) {
        std::string name = "h_weight_vs_" + fileTag + "_" + kFileTags[k];
        hists[k] = new TH2D(name.c_str(),
            // (kLabels[k] + ";" + xAxisTitle + ";G4RW Weight").c_str(),
            (";" + xAxisTitle + ";G4RW Weight").c_str(),
            nBinsX, xMin, xMax, nBinsY, yMin, yMax);
    }

    Long64_t N = t->GetEntries();
    long nFilled = 0, nSkippedNoReco = 0, nSkippedPdg = 0, nSkippedContain = 0;
    for (Long64_t i = 0; i < N; ++i) {
        t->GetEntry(i);
        if (pdg != 211) { nSkippedPdg++; continue; }
        if (useReco && hasReco != 1) { nSkippedNoReco++; continue; }

        if (containSel != -1) {
            bool contained = (std::fabs(endX) < kFidBound) && (std::fabs(endY) < kFidBound) && (std::fabs(endZ) < kFidBound);
            if ((containSel == 1 && !contained) || (containSel == 0 && contained)) { nSkippedContain++; continue; }
        }

        double xVal = useReco ? recoLength : length;
        for (int k = 0; k < 5; ++k) hists[k]->Fill(xVal, w[k]);
        nFilled++;
    }
    std::cout << "[" << lengthKind << ", containSel=" << containSel << "] Filled from " << nFilled
               << " pi+ entries (of " << N << " total rows; skipped " << nSkippedPdg << " non-pi+";
    if (useReco) std::cout << ", " << nSkippedNoReco << " with no reco match";
    if (containSel != -1) std::cout << ", " << nSkippedContain << " failing containment selection";
    std::cout << ").\n";

    for (int k = 0; k < 5; ++k) {
        TCanvas* c = new TCanvas(("c_" + fileTag + "_" + kFileTags[k]).c_str(), "", 800, 650);
        c->SetLeftMargin(0.14);
        c->SetRightMargin(0.14);
        // c->SetLogz();  // linear z-scale (color = raw bin counts, not log10)

        hists[k]->GetYaxis()->SetTitleOffset(1.3);
        hists[k]->Draw("COLZ");

        DrawDuneLabel();
        TLatex tR; tR.SetNDC(); tR.SetTextSize(0.035);
        std::string niceLabel = lengthKind;
        for (auto& ch : niceLabel) if (ch == '_') ch = ' ';
        tR.DrawLatex(0.55, 0.94, (niceLabel + " #pi^{+} (CAF), " + kLabels[k]).c_str());

        std::string outPng = kOutDir + "/caf_weight_vs_" + fileTag + "_2d_" + kFileTags[k] + ".png";
        c->SaveAs(outPng.c_str());
        delete c;
    }
}

void plot_caf_weight_vs_length_2d() {

    gROOT->LoadMacro("protoDUNEStyle.C");
    gROOT->SetStyle("protoDUNEStyle");
    gROOT->ForceStyle();
    gSystem->mkdir(kOutDir.c_str(), true);

    TFile* f = TFile::Open(kInFile.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: cannot open " << kInFile << " -- run dump_caf_pion_list.C first.\n";
        return;
    }
    TTree* t = (TTree*)f->Get("caf_pion_list");
    if (!t) {
        std::cerr << "ERROR: no caf_pion_list tree in " << kInFile << "\n";
        return;
    }

    bool hasRecoBranches = (t->GetBranch("recoLength") != nullptr) && (t->GetBranch("hasReco") != nullptr);
    bool hasEndBranches  = (t->GetBranch("endX") != nullptr) && (t->GetBranch("endY") != nullptr) && (t->GetBranch("endZ") != nullptr);

    MakeWeightVsLengthPlots(t, /*useReco=*/false, "Truth", "Track Length [cm]", "length");

    if (hasEndBranches) {
        MakeWeightVsLengthPlots(t, /*useReco=*/false, "Truth_Contained",    "Track Length [cm]", "length_contained",     /*containSel=*/1);
        MakeWeightVsLengthPlots(t, /*useReco=*/false, "Truth_NonContained", "Track Length [cm]", "length_noncontained",  /*containSel=*/0);
    } else {
        std::cerr << "WARNING: 'endX'/'endY'/'endZ' branches not found in caf_pion_list.root --\n"
                   << "skipping containment-split plots. Re-run dump_caf_pion_list.C with the\n"
                   << "endpos extension to produce them.\n";
    }

    if (hasRecoBranches) {
        MakeWeightVsLengthPlots(t, /*useReco=*/true, "Reco", "Reco Track Length [cm]", "recolength");
    } else {
        std::cerr << "WARNING: 'recoLength'/'hasReco' branches not found in caf_pion_list.root --\n"
                   << "skipping reco-length plots. Re-run the updated dump_caf_pion_list.C\n"
                   << "(with the syst_g4rw_pion_reco join) to produce them.\n";
    }

    std::cout << "\nDone. Wrote plots to " << kOutDir << "/\n";
    f->Close();
}

int main() {
    plot_caf_weight_vs_length_2d();
    return 0;
}