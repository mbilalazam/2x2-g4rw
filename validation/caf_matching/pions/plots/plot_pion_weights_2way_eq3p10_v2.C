///////////////////////////////////////////////////////////////////////////////
// plot_pion_weights_2way_eq3p10_v2.C
//
// Truth pion (CAF-matched, from caf_pion_list_v2.root) G4RW weight_1p33
// (f=1.33) distribution, comparing:
//
//   1) Patched edepsim, old formula (W_I = f)     -- WFILE_OLD, looked up
//   2) Patched edepsim, Eq. 3.10 formula          -- WFILE_EQ310, looked up
//
// FIX vs the v1 macro: caf_pion_list_v2.root's embedded w_1p33 comes from
// a DIFFERENT, older weight production (old_files/..._addChild_range_...)
// than either WFILE_OLD or WFILE_EQ310 -- not a fair comparison. This
// version uses caf_pion_list_v2 purely as an IDENTITY list (file_idx,
// event_id, track_id, pdg for genuinely CAF-matched pi+ tracks) and looks
// up fresh weights from BOTH current weight files by that same key, so
// the two histograms are apples-to-apples.
//
// Matching key: (file_idx, event_id, track_id) -- this is the REAL
// convention purity_dlp_multiplicity.cc itself uses (WKey5), confirmed by
// inspecting the actual weight-lookup call site. NOT the invented
// RunId/EventId decode from the original dump_caf_pion_list.C.
//
// Compile:
//   g++ -std=c++17 -O2 -g plot_pion_weights_2way_eq3p10_v2.C -o plot_pion_weights_2way_eq3p10_v2 `root-config --cflags --libs`
// Run:
//   ./plot_pion_weights_2way_eq3p10_v2
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

#include <map>
#include <tuple>
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>

// =============================================================
// Paths
// =============================================================
const std::string CAF_PION_LIST_V2 =
    "/global/cfs/cdirs/dune/users/mazam/softwares/patches/2x2_g4rw/validation/caf_pion_list_v2.root";

const std::string WFILE_OLD =
    "/global/cfs/cdirs/dune/users/mazam/G4RW/output_weights/patched_withRock_genieNu/"
    "weights_1000_all_LAr_PiPlus_cascade_0_9000_5weights.root";

const std::string WFILE_EQ310 =
    "/global/cfs/cdirs/dune/users/mazam/G4RW/output_weights/eq310_patched_withRock_genieNu/"
    "weights_1000_all_LAr_PiPlus_cascade_0_9000_5weights.root";

const std::string OUT_DIR =
    "/global/cfs/cdirs/dune/users/mazam/softwares/patches/2x2_g4rw/validation_pion_Eq3p10";

// =============================================================
// Matching key: (file_idx, event_id, track_id) -- the REAL convention.
// =============================================================
using Key = std::tuple<int,int,int>;

std::map<Key,double> LoadWeightMap(const std::string& file,
                                    const std::string& branchName = "weight_1p33")
{
    std::map<Key,double> wmap;

    TFile* f = TFile::Open(file.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: Cannot open weight file: " << file << "\n";
        return wmap;
    }
    TTree* t = (TTree*)f->Get("g4rwEvt");
    if (!t) {
        std::cerr << "ERROR: Cannot find g4rwEvt in " << file << "\n";
        f->Close();
        return wmap;
    }

    int file_idx, event_id, track_id, pdg;
    double weight;

    t->SetBranchAddress("file_idx",  &file_idx);
    t->SetBranchAddress("event_id",  &event_id);
    t->SetBranchAddress("track_id",  &track_id);
    t->SetBranchAddress("pdg",       &pdg);
    t->SetBranchAddress(branchName.c_str(), &weight);

    long N = t->GetEntries();
    for (long i = 0; i < N; i++) {
        t->GetEntry(i);
        if (pdg != 211) continue;   // pi+ only
        Key key = std::make_tuple(file_idx, event_id, track_id);
        wmap[key] = weight;
    }

    std::cout << "Loaded " << wmap.size() << " pi+ weight entries from " << file << "\n";

    f->Close();
    return wmap;
}

void DrawDuneLabel() {
    TLatex tL; tL.SetNDC();
    tL.DrawLatex(0.20, 0.94, "#bf{DUNE:ND-LAr 2x2}");
}

// =============================================================
// Main
// =============================================================
void plot_pion_weights_2way_eq3p10_v2()
{
    gROOT->LoadMacro("protoDUNEStyle.C");
    gROOT->SetStyle("protoDUNEStyle");
    gROOT->ForceStyle();
    gStyle->SetTitleX(0.35);
    gSystem->mkdir(OUT_DIR.c_str(), true);

    // ── Load both current-weight-file lookup maps ───────────────────────────
    std::map<Key,double> wmap_old   = LoadWeightMap(WFILE_OLD);
    std::map<Key,double> wmap_eq310 = LoadWeightMap(WFILE_EQ310);

    // ── Open the CAF-matched pion identity list ─────────────────────────────
    TFile* fList = TFile::Open(CAF_PION_LIST_V2.c_str());
    if (!fList || fList->IsZombie()) {
        std::cerr << "ERROR: Cannot open caf_pion_list_v2: " << CAF_PION_LIST_V2 << "\n";
        return;
    }
    TTree* tList = (TTree*)fList->Get("caf_pion_list");
    if (!tList) {
        std::cerr << "ERROR: Cannot find caf_pion_list tree\n";
        return;
    }

    int file_idx, event_id, track_id, pdg;
    tList->SetBranchAddress("file_idx", &file_idx);
    tList->SetBranchAddress("event_id", &event_id);
    tList->SetBranchAddress("track_id", &track_id);
    tList->SetBranchAddress("pdg",      &pdg);

    // ── Histograms ───────────────────────────────────────────────────────────
    TH1D* h_old = new TH1D("h_old",
        "Truth #pi^{+} (CAF-matched) f=1.33;G4RW Weight;Number of #pi^{+}", 100, 0, 2.0);
    TH1D* h_eq310 = new TH1D("h_eq310",
        "Truth #pi^{+} (CAF-matched) f=1.33;G4RW Weight;Number of #pi^{+}", 100, 0, 2.0);

    h_old->SetLineColor(kRed+1);   h_old->SetLineWidth(3);
    h_eq310->SetLineColor(kGreen+2);    h_eq310->SetLineWidth(3);
    h_eq310->SetLineStyle(2);   // dashed -- the two curves sit very close together

    long N = tList->GetEntries();
    long n_matched_old = 0, n_matched_eq310 = 0, n_pi_plus = 0;

    for (long i = 0; i < N; i++) {
        tList->GetEntry(i);
        if (pdg != 211) continue;   // pi+ only
        n_pi_plus++;

        Key key = std::make_tuple(file_idx, event_id, track_id);

        auto itOld = wmap_old.find(key);
        if (itOld != wmap_old.end()) {
            h_old->Fill(itOld->second);
            n_matched_old++;
        }

        auto itEq = wmap_eq310.find(key);
        if (itEq != wmap_eq310.end()) {
            h_eq310->Fill(itEq->second);
            n_matched_eq310++;
        }
    }

    std::cout << "\nTotal pi+ in caf_pion_list_v2: " << n_pi_plus << "\n";
    std::cout << "Matched into old (production) sample: " << n_matched_old
              << " / " << n_pi_plus << " (" << 100.0*n_matched_old/n_pi_plus << "%)\n";
    std::cout << "Matched into eq3.10 sample:           " << n_matched_eq310
              << " / " << n_pi_plus << " (" << 100.0*n_matched_eq310/n_pi_plus << "%)\n";

    // ── Plot ──────────────────────────────────────────────────────────────────
    TCanvas* c = new TCanvas("c_pion2way", "", 900, 650);
    c->SetLeftMargin(0.14);

    double ymax = std::max(h_old->GetMaximum(), h_eq310->GetMaximum());

    h_old->GetYaxis()->SetRangeUser(0, ymax * 1.4);
    h_old->GetYaxis()->SetTitleOffset(1.3);
    h_old->Draw("HIST");
    h_eq310->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.55, 0.72, 0.88, 0.88);
    leg->SetTextFont(43); leg->SetTextSize(18);
    leg->AddEntry(h_old,   "W_{I} = f_{I}", "l");
    leg->AddEntry(h_eq310, "W_{I} = f_{I} #frac{#sigma_{i}}{#sigma_{i}^{'}} #frac{1-e^{-#Delta L_{i}#sigma_{i}^{'}}}{1-e^{-#Delta L_{i}#sigma_{i}}}", "l");
    leg->Draw();

    DrawDuneLabel();

    std::string out = OUT_DIR + "/pion_weight_1p33_2way_comparison_v2.png";
    c->SaveAs(out.c_str());

    std::cout << "\nMeans -- old: " << h_old->GetMean()
              << "  eq3.10: " << h_eq310->GetMean() << "\n";

    std::cout << "\nDone. Plot saved to:\n" << out << "\n";

    fList->Close();
}

int main()
{
    plot_pion_weights_2way_eq3p10_v2();
    return 0;
}