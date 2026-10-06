///////////////////////////////////////////////////////////////////////////////
// plot_truth_pion_weights_patched.C
//
// Truth pion (matched interactions) G4RW weight distribution, patched sample
// only, all five cross-section scale factors overlaid on one plot.
//
// Unlike g4rw_weight_compare_patched_vs_unpatched_pions.C, no join against
// the raw g4rwEvt weight files is needed here: the weights are already
// stored directly in syst_g4rw_pion_truth (branch w_patched[5]), so this
// macro just reads the tree and fills histograms.
//
// Compile:
// g++ -std=c++17 -O2 -g plot_truth_pion_weights_patched.C -o plot_truth_pion_weights_patched `root-config --cflags --libs`
// Run:
//   ./plot_truth_pion_weights_patched
///////////////////////////////////////////////////////////////////////////////

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStyle.h"

#include <array>
#include <string>
#include <iostream>
#include <cmath>

// =============================================================
// Paths
// =============================================================
const std::string CAF_FILE =
    "/global/u1/m/mazam/Geant4ReWeighting/2x2_G4RW/test/on_CAFS/31Aug2026/output.root";

const std::string OUT_DIR =
    "/global/u1/m/mazam/Geant4ReWeighting/2x2_G4RW/test/on_CAFS/31Aug2026/truth_pion_weights";

// Fixed array convention used when w_patched[5] was filled:
// [0]=0.34, [1]=0.67, [2]=1.33, [3]=1.66, [4]=2.00
const std::array<double,5>      kScaleFactors = {0.34, 0.67, 1.33, 1.66, 2.00};
const std::array<std::string,5> kLabels       = {"f = 0.34", "f = 0.67", "f = 1.33", "f = 1.66", "f = 2.00"};
const std::array<int,5>         kColors       = {kCyan+2, kMagenta+1, kOrange+1, kRed, kGreen+2};

void DrawDuneLabel() {
    TLatex tL; tL.SetNDC();
    tL.DrawLatex(0.20, 0.94, "#bf{DUNE:ND-LAr 2x2}");
}

// =============================================================
// Main
// =============================================================
void plot_truth_pion_weights_patched()
{
    gROOT->LoadMacro("protoDUNEStyle.C");
    gROOT->SetStyle("protoDUNEStyle");
    gROOT->ForceStyle();
    gStyle->SetTitleX(0.35);
    gSystem->mkdir(OUT_DIR.c_str(), true);

    TFile* fCAF = TFile::Open(CAF_FILE.c_str());
    if (!fCAF || fCAF->IsZombie()) {
        std::cerr << "ERROR: Cannot open CAF file: " << CAF_FILE << "\n";
        return;
    }
    TTree* tTruth = (TTree*)fCAF->Get("syst_g4rw_pion_truth");
    if (!tTruth) {
        std::cerr << "ERROR: Cannot find syst_g4rw_pion_truth\n";
        return;
    }

    double w_patched[5];
    tTruth->SetBranchAddress("w_patched", w_patched);
    int pdg;
    tTruth->SetBranchAddress("pdg", &pdg);

    // ── Histograms, one per scale factor ────────────────────────────────────
    std::array<TH1D*,5> hists;
    for (int k = 0; k < 5; k++) {
        std::string name = "h_wpatched_" + std::to_string(k);
        hists[k] = new TH1D(name.c_str(),
            "Truth Pions (Matched Interactions);G4RW Weight (Patched EdepSim);Number of Pions",
            100, 0, 2.0);
        hists[k]->SetLineColor(kColors[k]);
        hists[k]->SetLineWidth(3);
    }

    // ── Loop and fill, plus a simple w==1 (no reweighting) sanity count ─────
    long N = tTruth->GetEntries();
    std::array<long,5> n_eq1 = {0,0,0,0,0};

    for (long i = 0; i < N; i++) {
        tTruth->GetEntry(i);
        if (pdg != 211) continue;   // exclude pi- (contaminated, no real weights yet)
        for (int k = 0; k < 5; k++) {
            hists[k]->Fill(w_patched[k]);
            if (std::abs(w_patched[k] - 1.0) < 1e-6) n_eq1[k]++;
        }
    }

    std::cout << "\nTotal truth pions (matched interactions): " << N << "\n";
    std::cout << "Fraction with w_patched == 1 (no reweighting), per scale factor:\n";
    for (int k = 0; k < 5; k++) {
        std::cout << "  " << kLabels[k] << " : "
                  << n_eq1[k] << " / " << N
                  << " (" << 100.0 * n_eq1[k] / N << "%)\n";
    }

    // ── Plot: all five overlaid on one canvas ───────────────────────────────
    TCanvas* c = new TCanvas("c_wpatched", "", 900, 650);
    c->SetLeftMargin(0.14);

    double ymax = 0;
    for (int k = 0; k < 5; k++) ymax = std::max(ymax, hists[k]->GetMaximum());

    hists[0]->GetYaxis()->SetRangeUser(0, ymax * 1.4);
    hists[0]->GetYaxis()->SetTitleOffset(1.3);
    hists[0]->Draw("HIST");
    for (int k = 1; k < 5; k++) hists[k]->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.62, 0.60, 0.88, 0.88);
    leg->SetTextFont(43); leg->SetTextSize(18);
    for (int k = 0; k < 5; k++) leg->AddEntry(hists[k], kLabels[k].c_str(), "l");
    leg->Draw();

    DrawDuneLabel();

    std::string out = OUT_DIR + "/truth_pion_weight_patched_all_scales.png";
    c->SaveAs(out.c_str());

    std::cout << "\nDone. Plot saved to:\n" << out << "\n";
    fCAF->Close();
}

int main()
{
    plot_truth_pion_weights_patched();
    return 0;
}