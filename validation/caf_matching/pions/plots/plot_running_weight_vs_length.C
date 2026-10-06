///////////////////////////////////////////////////////////////////////////////
// plot_running_weight_vs_length.C
//
// Answers "what weight is assigned at each step, across the track length"
// for a single track -- plots the running (cumulative) G4RW weight as a
// function of cumulative path length, for Example 1 and Example 2, using
// the TGraphs saved by validate_pion_example_1.c++ / _2.c++.
//
// One plot per example, one curve per scale factor (all 5 overlaid).
//
// One plot per (example, sample) combination, one curve per scale factor:
//   running_weight_vs_length_example_1_patched.png
//   running_weight_vs_length_example_1_unpatched.png
//   running_weight_vs_length_example_2_patched.png
//   running_weight_vs_length_example_2_unpatched.png
//
// Note: unpatched graphs have far fewer points (Example 1: 1 step;
// Example 2: 4 steps) than patched (45/6 steps) -- this is expected and
// is itself part of the visual story: the unpatched trajectory has much
// coarser granularity, which is exactly what the patch fixes.
//
// Compile and run:
//   g++ -std=c++17 -O2 -g plot_running_weight_vs_length.C -o plot_running_weight_vs_length `root-config --cflags --libs`
//   ./plot_running_weight_vs_length
///////////////////////////////////////////////////////////////////////////////

#include "TFile.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TMultiGraph.h"

#include <string>
#include <vector>
#include <iostream>

const std::string IN_DIR  = "/global/cfs/cdirs/dune/users/mazam/softwares/patches/2x2_g4rw/validation";
const std::string OUT_DIR = IN_DIR + "/plots_eq3p10";

// scale-factor labels as saved by validate_pion_example_*.c++, and their
// display colors/values.
struct ScaleInfo { std::string graphSuffix; std::string display; int color; int style; };
const std::vector<ScaleInfo> kScales = {
    {"f0.34", "0.34", kBlue+2,   1},
    {"f0.67", "0.67", kAzure+1,  1},
    {"f1.33", "1.33", kOrange+1, 1},
    {"f1.66", "1.66", kRed,      1},
    {"f2.00", "2.00", kRed+3,    1}
};

void DrawRunningWeight(const std::string& inFile, const std::string& exampleLabel,
                        const std::string& sampleTag, const std::string& sampleDisplay,
                        const std::string& outPng,
                        const std::string& graphPrefix = "running_weight_",
                        const std::string& yAxisTitle  = "G4RW Weight")
{
    TFile* f = TFile::Open(inFile.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: cannot open " << inFile << "\n";
        return;
    }

    TCanvas* c = new TCanvas("c_running", "c_running", 800, 600);
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.08);

    TMultiGraph* mg = new TMultiGraph();
    TLegend* leg = new TLegend(0.55, 0.15, 0.90, 0.45);
    leg->SetTextFont(133);
    leg->SetTextSize(16);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);

    bool anyFound = false;
    for (const auto& s : kScales) {
        std::string gname = graphPrefix + sampleTag + "_" + s.graphSuffix;
        TGraph* g = (TGraph*)f->Get(gname.c_str());
        if (!g) {
            std::cerr << "WARNING: graph not found: " << gname << " in " << inFile << "\n";
            continue;
        }
        anyFound = true;
        g->SetLineColor(s.color);
        g->SetLineWidth(3);
        g->SetLineStyle(s.style);
        g->SetMarkerColor(s.color);
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.1); // slightly larger, since unpatched graphs can have very few points
        mg->Add(g, "LP");
        leg->AddEntry(g, Form("f = %s", s.display.c_str()), "l");
    }
    if (!anyFound) {
        std::cerr << "ERROR: no running-weight graphs found for sample '" << sampleTag
                   << "' in " << inFile << " -- skipping this plot.\n";
        f->Close();
        return;
    }

    mg->SetTitle("");
    mg->Draw("A");
    mg->GetXaxis()->SetTitle("Cumulative Track Length [cm]");
    mg->GetYaxis()->SetTitle(yAxisTitle.c_str());
    mg->GetYaxis()->SetTitleOffset(1.1);

    // horizontal reference line at w=1 (the f=1.00 sanity-check level)
    double xmin = mg->GetXaxis()->GetXmin();
    double xmax = mg->GetXaxis()->GetXmax();
    TLine* refLine = new TLine(xmin, 1.0, xmax, 1.0);
    refLine->SetLineColor(kGray+1);
    refLine->SetLineStyle(2);
    refLine->Draw("SAME");

    leg->Draw();

    TLatex tL; tL.SetNDC();
    tL.DrawLatex(0.15, 0.94, "#bf{DUNE:ND-LAr 2x2}");
    TLatex tR; tR.SetNDC(); tR.SetTextSize(0.035);
    tR.DrawLatex(0.55, 0.94, Form("%s (%s)", exampleLabel.c_str(), sampleDisplay.c_str()));

    c->SaveAs(outPng.c_str());
    delete c;
    f->Close();
}

void plot_running_weight_vs_length()
{
    gROOT->LoadMacro("protoDUNEStyle.C");
    gROOT->SetStyle("protoDUNEStyle");
    gROOT->ForceStyle();
    gStyle->SetTitleX(0.25);
    gSystem->mkdir(OUT_DIR.c_str(), true);

    DrawRunningWeight(IN_DIR + "/validate_pion_example_1.root",
                       "Example 1: #pi^{+}, 6.83 cm (quasielastic)",
                       "patched", "Patched",
                       OUT_DIR + "/running_weight_vs_length_example_1_patched.png");
    DrawRunningWeight(IN_DIR + "/validate_pion_example_1.root",
                       "Example 1: #pi^{+}, 6.83 cm (quasielastic)",
                       "unpatched", "Unpatched",
                       OUT_DIR + "/running_weight_vs_length_example_1_unpatched.png");

    DrawRunningWeight(IN_DIR + "/validate_pion_example_2_eq3p10.root",
                       "Example 2: #pi^{+}, 14.4 cm (absorption)",
                       "patched", "Patched",
                       OUT_DIR + "/running_weight_vs_length_example_2_patched.png");
    DrawRunningWeight(IN_DIR + "/validate_pion_example_2_eq3p10.root",
                       "Example 2: #pi^{+}, 14.4 cm (absorption)",
                       "unpatched", "Unpatched",
                       OUT_DIR + "/running_weight_vs_length_example_2_unpatched.png");

    // DrawRunningWeight(IN_DIR + "/validate_pion_example_1.root",
    //                    "Example 1: #pi^{+}, 6.83 cm (quasielastic)",
    //                    "patched", "Patched",
    //                    OUT_DIR + "/step_weight_vs_length_example_1_patched.png",
    //                    "step_weight_", "G4RW Step Weight");
    // DrawRunningWeight(IN_DIR + "/validate_pion_example_1.root",
    //                    "Example 1: #pi^{+}, 6.83 cm (quasielastic)",
    //                    "unpatched", "Unpatched",
    //                    OUT_DIR + "/step_weight_vs_length_example_1_unpatched.png",
    //                    "step_weight_", "G4RW Step Weight");

    DrawRunningWeight(IN_DIR + "/validate_pion_example_2_eq3p10.root",
                       "Example 2: #pi^{+}, 14.4 cm (absorption)",
                       "patched", "Patched",
                       OUT_DIR + "/step_weight_vs_length_example_2_patched.png",
                       "step_weight_", "G4RW Step Weight");
    DrawRunningWeight(IN_DIR + "/validate_pion_example_2_eq3p10.root",
                       "Example 2: #pi^{+}, 14.4 cm (absorption)",
                       "unpatched", "Unpatched",
                       OUT_DIR + "/step_weight_vs_length_example_2_unpatched.png",
                       "step_weight_", "G4RW Step Weight");


    // Example 3 (Band A) and Example 5 (Band B anomaly): the pion ranges
    // out to p=0 via ionization BEFORE any recorded point carries the
    // proc=4/subproc=121 tag do_weighting_v2_PiPlus.c__ checks for --
    // confirmed by inspecting each track's own trajectory printout.
    // last_is_inel is therefore always false for these two, so
    // GetExclusiveFactor(f) is never applied: the saved weight is pure
    // e^{(1-f)S} survival, with NO interaction-kick term. Labeled
    // explicitly below rather than presented as an ordinary example.
    DrawRunningWeight(IN_DIR + "/validate_pion_example_3.root",
                       "Example 3: #pi^{+}, 280.2 cm (NO inelastic tag found -- pure survival, no kick)",
                       "patched", "Patched",
                       OUT_DIR + "/running_weight_vs_length_example_3_patched.png");
    DrawRunningWeight(IN_DIR + "/validate_pion_example_3.root",
                       "Example 3: #pi^{+}, 280.2 cm (NO inelastic tag found -- pure survival, no kick)",
                       "patched", "Patched",
                       OUT_DIR + "/step_weight_vs_length_example_3_patched.png",
                       "step_weight_", "G4RW Step Weight");

    DrawRunningWeight(IN_DIR + "/validate_pion_example_4.root",
                       "Example 4: #pi^{+}, 234.2 cm (clean -- inelastic tag correctly found)",
                       "patched", "Patched",
                       OUT_DIR + "/running_weight_vs_length_example_4_patched.png");
    DrawRunningWeight(IN_DIR + "/validate_pion_example_4.root",
                       "Example 4: #pi^{+}, 234.2 cm (clean -- inelastic tag correctly found)",
                       "patched", "Patched",
                       OUT_DIR + "/step_weight_vs_length_example_4_patched.png",
                       "step_weight_", "G4RW Step Weight");

    DrawRunningWeight(IN_DIR + "/validate_pion_example_5.root",
                       "Example 5: #pi^{+}, 145.1 cm (NO inelastic tag found -- pure survival, no kick)",
                       "patched", "Patched",
                       OUT_DIR + "/running_weight_vs_length_example_5_patched.png");
    DrawRunningWeight(IN_DIR + "/validate_pion_example_5.root",
                       "Example 5: #pi^{+}, 145.1 cm (NO inelastic tag found -- pure survival, no kick)",
                       "patched", "Patched",
                       OUT_DIR + "/step_weight_vs_length_example_5_patched.png",
                       "step_weight_", "G4RW Step Weight");

    DrawRunningWeight(IN_DIR + "/validate_pion_example_6.root",
                       "Example 6: #pi^{+}, 70.7 cm (clean -- inelastic tag correctly found)",
                       "patched", "Patched",
                       OUT_DIR + "/running_weight_vs_length_example_6_patched.png");
    DrawRunningWeight(IN_DIR + "/validate_pion_example_6.root",
                       "Example 6: #pi^{+}, 70.7 cm (clean -- inelastic tag correctly found)",
                       "patched", "Patched",
                       OUT_DIR + "/step_weight_vs_length_example_6_patched.png",
                       "step_weight_", "G4RW Step Weight");

    DrawRunningWeight(IN_DIR + "/validate_pion_example_3.root",
                       "Example 3: #pi^{+}, 280.2 cm (NO inelastic tag found -- pure survival, no kick)",
                       "unpatched", "Unpatched",
                       OUT_DIR + "/running_weight_vs_length_example_3_unpatched.png");
    DrawRunningWeight(IN_DIR + "/validate_pion_example_3.root",
                       "Example 3: #pi^{+}, 280.2 cm (NO inelastic tag found -- pure survival, no kick)",
                       "unpatched", "Unpatched",
                       OUT_DIR + "/step_weight_vs_length_example_3_unpatched.png",
                       "step_weight_", "G4RW Step Weight");

    DrawRunningWeight(IN_DIR + "/validate_pion_example_4.root",
                       "Example 4: #pi^{+}, 234.2 cm (clean -- inelastic tag correctly found)",
                       "unpatched", "Unpatched",
                       OUT_DIR + "/running_weight_vs_length_example_4_unpatched.png");
    DrawRunningWeight(IN_DIR + "/validate_pion_example_4.root",
                       "Example 4: #pi^{+}, 234.2 cm (clean -- inelastic tag correctly found)",
                       "unpatched", "Unpatched",
                       OUT_DIR + "/step_weight_vs_length_example_4_unpatched.png",
                       "step_weight_", "G4RW Step Weight");

    DrawRunningWeight(IN_DIR + "/validate_pion_example_5.root",
                       "Example 5: #pi^{+}, 145.1 cm (NO inelastic tag found -- pure survival, no kick)",
                       "unpatched", "Unpatched",
                       OUT_DIR + "/running_weight_vs_length_example_5_unpatched.png");
    DrawRunningWeight(IN_DIR + "/validate_pion_example_5.root",
                       "Example 5: #pi^{+}, 145.1 cm (NO inelastic tag found -- pure survival, no kick)",
                       "unpatched", "Unpatched",
                       OUT_DIR + "/step_weight_vs_length_example_5_unpatched.png",
                       "step_weight_", "G4RW Step Weight");

    DrawRunningWeight(IN_DIR + "/validate_pion_example_6.root",
                       "Example 6: #pi^{+}, 70.7 cm (clean -- inelastic tag correctly found)",
                       "unpatched", "Unpatched",
                       OUT_DIR + "/running_weight_vs_length_example_6_unpatched.png");
    DrawRunningWeight(IN_DIR + "/validate_pion_example_6.root",
                       "Example 6: #pi^{+}, 70.7 cm (clean -- inelastic tag correctly found)",
                       "unpatched", "Unpatched",
                       OUT_DIR + "/step_weight_vs_length_example_6_unpatched.png",
                       "step_weight_", "G4RW Step Weight");
    
    std::cout << "\nDone.\n";
}

int main()
{
    plot_running_weight_vs_length();
    return 0;
}
