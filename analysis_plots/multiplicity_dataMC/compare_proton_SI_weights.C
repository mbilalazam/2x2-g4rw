///////////////////////////////////////////////////////////////////////////////
// compare_proton_SI_weights.C
//
// Sanity check: plot the distribution of matched G4RW/SI weights
// (weight_1p33, weight_0p67) for the species going into the multiplicity
// distribution.
//
// Two comparisons are made:
//
//   A) Protons, TWO different reco trees against the SAME proton weight
//      file:
//        1) MC_FILE / syst_g4rw_mult_reco  -- the multiplicity study's
//           tree (nominal_mc.cc, L>2cm, no containment cut, all species
//           mixed in one tree, filtered here to truePDG==2212)
//        2) PROTON_CAF_FILE / syst_g4rw_proton_reco -- the proton study's
//           tree (proton_technote.cc-derived code, L>3cm, containment ON,
//           protons only by construction)
//      Both trees are keyed the same way -- (cafFileIdx, nu_id, G4ID) --
//      so any large difference in the resulting weight distributions
//      would point to a real inconsistency between the two
//      selections/samples, not a lookup bug.
//
//   B) pi+, pi-, and (for reference) proton, ALL from MC_FILE /
//      syst_g4rw_mult_reco only -- i.e. the actual species mix that goes
//      into the multiplicity distribution, each looked up against its own
//      weight file (WFILE_PIPLUS/WFILE_PIMINUS/WFILE_PROTON).
//
// All weights are plotted RAW (no San() [0.01,100] clamp), linear x-axis
// 0 to 2, raw counts (no normalization) -- matching the style of the
// existing "Truth Pions (CAF-matched, f=1.33)" reference plot.
//
// Compile:
//   g++ -std=c++17 -O2 -g compare_proton_SI_weights.C \
//       -o compare_proton_SI_weights `root-config --cflags --libs`
//
// Run:
//   ./compare_proton_SI_weights
///////////////////////////////////////////////////////////////////////////////

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"

#include <map>
#include <tuple>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <algorithm>

// =============================================================================
// Paths
// =============================================================================
const std::string MC_FILE =
    "/global/cfs/cdirs/dune/users/mazam/talks/workshop_Sep2026/nominal_mc.root";

const std::string PROTON_CAF_FILE =
    "/global/cfs/cdirs/dune/users/mazam/softwares/patches/2x2_g4rw/validation_proton/codes_for_slides/proton_plots/output.root";

const std::string WFILE_PROTON =
    "/global/cfs/cdirs/dune/users/mazam/G4RW/output_weights/patched_withRock_genieNu_eq3p10/weights_1000_all_LAr_P_cascade_0_9000_5weights.root";

const std::string WDIR_PAT =
    "/global/cfs/cdirs/dune/users/mazam/G4RW/output_weights/patched_withRock_genieNu_eq3p10/";
const std::string WFILE_PIPLUS  = WDIR_PAT + "weights_1000_all_LAr_PiPlus_cascade_0_9000_5weights.root";
const std::string WFILE_PIMINUS = WDIR_PAT + "weights_1000_all_LAr_PiMinus_cascade_0_9000_5weights.root";

const std::string OUT_DIR =
    "/global/cfs/cdirs/dune/users/mazam/talks/workshop_Sep2026/multiplicity_systematics";

const char* LAB_DUNE = "#bf{DUNE:ND-LAr 2x2}";

// =============================================================================
// Weight map -- RAW weight_1p33 and weight_0p67, no San() clamp
// =============================================================================
using WKey = std::tuple<int,int,int>;  // (cafFileIdx, nu_id, G4ID)
struct W2 { double w1p33, w0p67; };

void LoadWeightsRaw(const std::string& file, std::map<WKey,W2>& wmap)
{
    TFile* f = TFile::Open(file.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: Cannot open: " << file << "\n";
        return;
    }
    TTree* t = (TTree*)f->Get("g4rwEvt");
    if (!t) {
        std::cerr << "ERROR: No g4rwEvt tree in: " << file << "\n";
        f->Close();
        return;
    }

    int file_idx, event_id, track_id;
    double w_1p33, w_0p67;
    t->SetBranchAddress("file_idx",    &file_idx);
    t->SetBranchAddress("event_id",    &event_id);
    t->SetBranchAddress("track_id",    &track_id);
    t->SetBranchAddress("weight_1p33", &w_1p33);
    t->SetBranchAddress("weight_0p67", &w_0p67);

    long N = t->GetEntries();
    for (long i = 0; i < N; ++i) {
        t->GetEntry(i);
        wmap[{file_idx, event_id, track_id}] = {w_1p33, w_0p67};
    }
    f->Close();
    std::cout << "Loaded " << N << " entries from: " << file << "\n";
}

// =============================================================================
// Fill up/down histograms with RAW matched proton weight values from a
// reco tree. Both trees (syst_g4rw_mult_reco, syst_g4rw_proton_reco) share
// the same branch names: nu_id, cafFileIdx, G4ID, truePDG.
// =============================================================================
struct FillStats {
    long nTotal = 0, nProtonTrue = 0, nMatched = 0;
    double minUp = 1e300, maxUp = -1e300, minDn = 1e300, maxDn = -1e300;
};

FillStats FillProtonWeightHists(TTree* tree,
                                 const std::map<WKey,W2>& wmap,
                                 TH1D* hUp, TH1D* hDn)
{
    FillStats st;

    long long nu_id;
    int cafFileIdx, G4ID, truePDG;
    tree->SetBranchAddress("nu_id",      &nu_id);
    tree->SetBranchAddress("cafFileIdx", &cafFileIdx);
    tree->SetBranchAddress("G4ID",       &G4ID);
    tree->SetBranchAddress("truePDG",    &truePDG);

    long N = tree->GetEntries();
    for (long i = 0; i < N; ++i) {
        tree->GetEntry(i);
        st.nTotal++;
        if (abs(truePDG) != 2212) continue;
        st.nProtonTrue++;

        WKey key = {cafFileIdx, (int)nu_id, G4ID};
        auto it = wmap.find(key);
        if (it == wmap.end()) continue;
        st.nMatched++;

        double wu = it->second.w1p33;
        double wd = it->second.w0p67;
        hUp->Fill(wu);
        hDn->Fill(wd);

        st.minUp = std::min(st.minUp, wu); st.maxUp = std::max(st.maxUp, wu);
        st.minDn = std::min(st.minDn, wd); st.maxDn = std::max(st.maxDn, wd);
    }

    tree->ResetBranchAddresses();
    return st;
}

void PrintStats(const char* label, const FillStats& st)
{
    std::cout << "\n--- " << label << " ---\n";
    std::cout << "  Total tree entries:        " << st.nTotal << "\n";
    std::cout << "  truePDG == proton:         " << st.nProtonTrue << "\n";
    std::cout << "  Matched in weight file:    " << st.nMatched << "\n";
    if (st.nMatched > 0) {
        std::cout << "  w_1p33 range: [" << st.minUp << ", " << st.maxUp << "]\n";
        std::cout << "  w_0p67 range: [" << st.minDn << ", " << st.maxDn << "]\n";
    }
}

// =============================================================================
// Generic species fill -- exact truePDG match (not abs()), so pi+ (211) and
// pi- (-211) are kept distinct. Same tree schema (nu_id, cafFileIdx, G4ID,
// truePDG) as FillProtonWeightHists above.
// =============================================================================
FillStats FillSpeciesWeightHists(TTree* tree,
                                  const std::map<WKey,W2>& wmap,
                                  int targetTruePDG,
                                  TH1D* hUp, TH1D* hDn)
{
    FillStats st;

    long long nu_id;
    int cafFileIdx, G4ID, truePDG;
    tree->SetBranchAddress("nu_id",      &nu_id);
    tree->SetBranchAddress("cafFileIdx", &cafFileIdx);
    tree->SetBranchAddress("G4ID",       &G4ID);
    tree->SetBranchAddress("truePDG",    &truePDG);

    long N = tree->GetEntries();
    for (long i = 0; i < N; ++i) {
        tree->GetEntry(i);
        st.nTotal++;
        if (truePDG != targetTruePDG) continue;
        st.nProtonTrue++;   // reused field: count of entries matching this species

        WKey key = {cafFileIdx, (int)nu_id, G4ID};
        auto it = wmap.find(key);
        if (it == wmap.end()) continue;
        st.nMatched++;

        double wu = it->second.w1p33;
        double wd = it->second.w0p67;
        hUp->Fill(wu);
        hDn->Fill(wd);

        st.minUp = std::min(st.minUp, wu); st.maxUp = std::max(st.maxUp, wu);
        st.minDn = std::min(st.minDn, wd); st.maxDn = std::max(st.maxDn, wd);
    }

    tree->ResetBranchAddresses();
    return st;
}

void PrintSpeciesStats(const char* label, const char* speciesName, const FillStats& st)
{
    std::cout << "\n--- " << label << " ---\n";
    std::cout << "  Total tree entries:        " << st.nTotal << "\n";
    std::cout << "  truePDG == " << speciesName << ":            " << st.nProtonTrue << "\n";
    std::cout << "  Matched in weight file:    " << st.nMatched << "\n";
    if (st.nMatched > 0) {
        std::cout << "  w_1p33 range: [" << st.minUp << ", " << st.maxUp << "]\n";
        std::cout << "  w_0p67 range: [" << st.minDn << ", " << st.maxDn << "]\n";
    }
}

void DrawComparison(TH1D* hA, TH1D* hB,
                     const char* labelA, const char* labelB,
                     const char* title, const std::string& outPng)
{
    TCanvas* c = new TCanvas("c_cmp", "", 800, 600);
    c->SetLeftMargin(0.14);

    hA->SetLineColor(kGreen+2);
    hA->SetLineWidth(2);
    hB->SetLineColor(kRed+1);
    hB->SetLineWidth(2);
    hB->SetLineStyle(2);

    hA->SetTitle(title);
    hA->GetXaxis()->SetTitle("G4RW Weight");
    hA->GetXaxis()->SetRangeUser(0, 2);
    hA->GetYaxis()->SetTitle("Number of Protons");

    double ymax = 1.3 * std::max(hA->GetMaximum(), hB->GetMaximum());
    hA->GetYaxis()->SetRangeUser(0, ymax);

    hA->Draw("HIST");
    hB->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.45, 0.75, 0.90, 0.88);
    leg->SetTextFont(43); leg->SetTextSize(16);
    leg->SetFillStyle(0); leg->SetBorderSize(0);
    leg->AddEntry(hA, labelA, "l");
    leg->AddEntry(hB, labelB, "l");
    leg->Draw();

    TLatex tL; tL.SetNDC();
    tL.DrawLatex(0.15, 0.94, LAB_DUNE);

    c->SaveAs(outPng.c_str());
    delete c;
}

// =============================================================================
// Generalized overlay of N histograms (used for the pi+/pi- comparison and
// the proton+pi+/pi- combined comparison) -- same style (linear 0-2, raw
// counts) as DrawComparison above.
// =============================================================================
struct HistSpec { TH1D* h; std::string label; int color; int style; };

void DrawMultiSpecies(const std::vector<HistSpec>& specs,
                      const char* title, const std::string& outPng)
{
    TCanvas* c = new TCanvas("c_multi", "", 800, 600);
    c->SetLeftMargin(0.14);

    double ymax = 0.0;
    for (const auto& s : specs) {
        s.h->SetLineColor(s.color);
        s.h->SetLineWidth(2);
        s.h->SetLineStyle(s.style);
        ymax = std::max(ymax, s.h->GetMaximum());
    }
    ymax *= 1.3;

    specs[0].h->SetTitle(title);
    specs[0].h->GetXaxis()->SetTitle("G4RW Weight");
    specs[0].h->GetXaxis()->SetRangeUser(0, 2);
    specs[0].h->GetYaxis()->SetTitle("Number of Tracks");
    specs[0].h->GetYaxis()->SetRangeUser(0, ymax);
    specs[0].h->Draw("HIST");
    for (size_t i = 1; i < specs.size(); ++i)
        specs[i].h->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.45, 0.70, 0.90, 0.88);
    leg->SetTextFont(43); leg->SetTextSize(16);
    leg->SetFillStyle(0); leg->SetBorderSize(0);
    for (const auto& s : specs) leg->AddEntry(s.h, s.label.c_str(), "l");
    leg->Draw();

    TLatex tL; tL.SetNDC();
    tL.DrawLatex(0.15, 0.94, LAB_DUNE);

    c->SaveAs(outPng.c_str());
    delete c;
}

// =============================================================================
// Main
// =============================================================================
void compare_proton_SI_weights()
{
    gROOT->LoadMacro("protoDUNEStyle.C");
    gROOT->SetStyle("protoDUNEStyle");
    gROOT->ForceStyle();
    TH1::AddDirectory(false);
    gStyle->SetTitleX(0.35);
    gSystem->mkdir(OUT_DIR.c_str(), true);

    // ── Load weight files once (raw, no clamp) ───────────────────────────────
    std::map<WKey,W2> wmap_proton, wmap_piplus, wmap_piminus;
    LoadWeightsRaw(WFILE_PROTON,  wmap_proton);
    LoadWeightsRaw(WFILE_PIPLUS,  wmap_piplus);
    LoadWeightsRaw(WFILE_PIMINUS, wmap_piminus);

    // ── Multiplicity study tree ──────────────────────────────────────────────
    TFile* fMC = TFile::Open(MC_FILE.c_str());
    if (!fMC || fMC->IsZombie()) { std::cerr << "ERROR: Cannot open " << MC_FILE << "\n"; return; }
    TTree* tMult = (TTree*)fMC->Get("syst_g4rw_mult_reco");
    if (!tMult) { std::cerr << "ERROR: No syst_g4rw_mult_reco in " << MC_FILE << "\n"; return; }

    // ── Proton study tree ─────────────────────────────────────────────────
    TFile* fProton = TFile::Open(PROTON_CAF_FILE.c_str());
    if (!fProton || fProton->IsZombie()) { std::cerr << "ERROR: Cannot open " << PROTON_CAF_FILE << "\n"; return; }
    TTree* tProton = (TTree*)fProton->Get("syst_g4rw_proton_reco");
    if (!tProton) { std::cerr << "ERROR: No syst_g4rw_proton_reco in " << PROTON_CAF_FILE << "\n"; return; }

    // ── Build histograms: linear, 0 to 2, 100 bins ──────────────────────────
    TH1D* hUp_mult    = new TH1D("hUp_mult",    "", 100, 0, 2);
    TH1D* hDn_mult    = new TH1D("hDn_mult",    "", 100, 0, 2);
    TH1D* hUp_proton  = new TH1D("hUp_proton",  "", 100, 0, 2);
    TH1D* hDn_proton  = new TH1D("hDn_proton",  "", 100, 0, 2);
    TH1D* hUp_piplus  = new TH1D("hUp_piplus",  "", 100, 0, 2);
    TH1D* hDn_piplus  = new TH1D("hDn_piplus",  "", 100, 0, 2);
    TH1D* hUp_piminus = new TH1D("hUp_piminus", "", 100, 0, 2);
    TH1D* hDn_piminus = new TH1D("hDn_piminus", "", 100, 0, 2);

    // All reads from tMult (proton, pi+, pi-) happen BEFORE fMC->Close() --
    // otherwise the tree would be gone for the later fills.
    FillStats stMult    = FillProtonWeightHists(tMult, wmap_proton, hUp_mult, hDn_mult);
    FillStats stPiPlus  = FillSpeciesWeightHists(tMult, wmap_piplus,   211, hUp_piplus,  hDn_piplus);
    FillStats stPiMinus = FillSpeciesWeightHists(tMult, wmap_piminus, -211, hUp_piminus, hDn_piminus);
    FillStats stProton  = FillProtonWeightHists(tProton, wmap_proton, hUp_proton, hDn_proton);

    fMC->Close();
    fProton->Close();

    PrintStats("Multiplicity study (syst_g4rw_mult_reco)", stMult);
    PrintStats("Proton study (syst_g4rw_proton_reco)",     stProton);
    PrintSpeciesStats("Multiplicity study (syst_g4rw_mult_reco)", "pi+", stPiPlus);
    PrintSpeciesStats("Multiplicity study (syst_g4rw_mult_reco)", "pi-", stPiMinus);

    // ── Draw: proton, mult study vs proton study ─────────────────────────
    DrawComparison(hUp_mult, hUp_proton,
                   "Nominal Multiplicity Protons", "Proton Study Protons",
                   "Protons from CAFs (f=1.33)",
                   OUT_DIR + "/compare_protonSIweights_w1p33.png");

    DrawComparison(hDn_mult, hDn_proton,
                   "Nominal Multiplicity Protons", "Proton Study Protons",
                   "Protons from CAFs (f=0.67)",
                   OUT_DIR + "/compare_protonSIweights_w0p67.png");

    // ── Draw: pi+ vs pi-, from nominal_mc.root only ──────────────────────
    DrawMultiSpecies({
            {hUp_piplus,  "pi+", kBlue+1,    1},
            {hUp_piminus, "pi-", kMagenta+2, 2},
        },
        "Pions from CAFs (f=1.33)",
        OUT_DIR + "/compare_pionSIweights_w1p33.png");

    DrawMultiSpecies({
            {hDn_piplus,  "pi+", kBlue+1,    1},
            {hDn_piminus, "pi-", kMagenta+2, 2},
        },
        "Pions from CAFs (f=0.67)",
        OUT_DIR + "/compare_pionSIweights_w0p67.png");

    // ── Draw: proton + pi+ + pi-, all from nominal_mc.root ───────────────
    DrawMultiSpecies({
            {hUp_mult,    "proton", kGreen+2,   1},
            {hUp_piplus,  "pi+",    kBlue+1,    1},
            {hUp_piminus, "pi-",    kMagenta+2, 2},
        },
        "Protons + Pions from CAFs (f=1.33)",
        OUT_DIR + "/compare_allSpeciesSIweights_w1p33.png");

    DrawMultiSpecies({
            {hDn_mult,    "proton", kGreen+2,   1},
            {hDn_piplus,  "pi+",    kBlue+1,    1},
            {hDn_piminus, "pi-",    kMagenta+2, 2},
        },
        "Protons + Pions from CAFs (f=0.67)",
        OUT_DIR + "/compare_allSpeciesSIweights_w0p67.png");

    std::cout << "\nDone. Plots saved in:\n" << OUT_DIR << "\n";
}

int main()
{
    compare_proton_SI_weights();
    return 0;
}