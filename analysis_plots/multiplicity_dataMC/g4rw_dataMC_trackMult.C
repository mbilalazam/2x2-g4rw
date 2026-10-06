///////////////////////////////////////////////////////////////////////////////
// g4rw_dataMC_trackMult.C
//
// Track multiplicity SI band, built the SAME WAY as the proton-only
// length/cosTheta plots (plot_dataMC_recoProton_truthSplit_all.C): each
// hadron track contributes its OWN weight independently, with NO
// event-level combination -- no multiplying several tracks' weights
// together, no quadrature combination either.
//
//   - h_nom_events: ONE entry per interaction (unweighted), at that
//     interaction's trackMult -- this is the displayed "Number of
//     Interactions" curve.
//   - h_nom_tracks / h_up_tracks / h_dn_tracks: ONE entry PER HADRON
//     TRACK (proton, pi+, pi- -- muon/kaon/other default to w=1, since
//     they have no G4RW/SI weight file), filled at the MULTIPLICITY OF
//     THE EVENT that track belongs to, weighted by THAT TRACK'S OWN
//     weight. No combination across tracks in the same event at all --
//     exactly mirroring how the proton script fills
//     h_up_len->Fill(length, w) once per proton, independently,
//     regardless of what else is in that event.
//   - Fractions: symmetric envelope, sigma_i = 0.5*|sum_up_i -
//     sum_dn_i| / nTracks_i, computed from these TRACK-level sums,
//     pooling all hadron species together. (The proton script's own
//     asymmetric max(0, up-nom)/max(0, nom-down) formula was tried
//     first and collapsed to zero in every bin here: it assumes
//     meanW(1.33) > 1 and meanW(0.67) < 1, but for this pooled hadron
//     population BOTH means point the other way in every bin, so
//     max(0, ...) clamped both directions to zero. The symmetric
//     envelope doesn't require knowing which universe is numerically
//     larger, so it doesn't have that failure mode.)
//   - Band applied to h_nom_events (the actual interaction count).
//
// mult=1 note (mostly lone muons -- see the "first bin error=0"
// question): muons have no G4RW/SI weight file, so a bin dominated by
// single-muon interactions naturally gets essentially zero systematic
// here, since there's nothing being reweighted in it.
//
// Produces THREE plots:
//   1. g4rw_MConly_trackMult.png       (MC only)
//   2. g4rw_dataMC_trackMult_band.png  (data+MC+band, no frac panel)
//   3. g4rw_dataMC_trackMult_splitPad.png (data+MC+band+frac panel)
//
// Compile:
//   g++ -std=c++17 -O2 -g g4rw_dataMC_trackMult.C -o g4rw_dataMC_trackMult `root-config --cflags --libs`
// Run:
//   ./g4rw_dataMC_trackMult
///////////////////////////////////////////////////////////////////////////////

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"

#include <map>
#include <tuple>
#include <array>
#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>

// =============================================================================
// Paths
// =============================================================================
const std::string DATA_FILE =
    "/global/cfs/cdirs/dune/users/mazam/talks/workshop_Sep2026/nominal_data.root";
const std::string MC_FILE =
    "/global/cfs/cdirs/dune/users/mazam/talks/workshop_Sep2026/nominal_mc.root";

const std::string WDIR_PAT =
    "/global/cfs/cdirs/dune/users/mazam/G4RW/output_weights/patched_withRock_genieNu_eq3p10/";
const std::string WFILE_PIPLUS  = WDIR_PAT + "weights_1000_all_LAr_PiPlus_cascade_0_9000_5weights.root";
const std::string WFILE_PIMINUS = WDIR_PAT + "weights_1000_all_LAr_PiMinus_cascade_0_9000_5weights.root";
const std::string WFILE_PROTON  = WDIR_PAT + "weights_1000_all_LAr_P_cascade_0_9000_5weights.root";

const std::string OUT_DIR =
    "/global/cfs/cdirs/dune/users/mazam/talks/workshop_Sep2026/multiplicity_systematics";

const char* LAB_DUNE = "#bf{DUNE:ND-LAr 2x2}";

// Locked cosmetics -- same as the rest of the trackMult family
// (detector_dataMC_trackMult.C, nusyst_dataMC_trackMult.C).
const double Y_MAX = 80.0;
const double LEG_X1 = 0.55, LEG_Y1 = 0.60, LEG_X2 = 0.88, LEG_Y2 = 0.88;

const int    nBinsM = 15;
const double xMinM  = 0.0, xMaxM = 15.0;

// =============================================================================
// Weight map -- San()-clamped, same clamp used throughout this workstream.
// =============================================================================
using WKey = std::tuple<int,int,int>;
using W2   = std::array<double,2>;  // [0]=w_1p33, [1]=w_0p67

double San(double w) { return (w < 0.01 || w > 100.0) ? 1.0 : w; }

void LoadWeights(const std::string& file, std::map<WKey,W2>& wmap)
{
    TFile* f = TFile::Open(file.c_str());
    if (!f || f->IsZombie()) { std::cerr << "ERROR: Cannot open: " << file << "\n"; return; }
    TTree* t = (TTree*)f->Get("g4rwEvt");
    if (!t) { std::cerr << "ERROR: No g4rwEvt in: " << file << "\n"; f->Close(); return; }

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
        wmap[{file_idx, event_id, track_id}] = {San(w_1p33), San(w_0p67)};
    }
    f->Close();
    std::cout << "Loaded " << N << " entries from: " << file << "\n";
}

// =============================================================================
// Fractions -- symmetric envelope, from TRACK-level sums pooling all
// hadron species. See header comment for why this replaced the proton
// script's original asymmetric max(0, up-nom)/max(0, nom-down) formula.
// =============================================================================
struct Fractions { std::vector<double> up, dn; };

Fractions ComputeFractions(TH1D* h_nom, TH1D* h_up, TH1D* h_dn)
{
    Fractions frac;
    int nB = h_nom->GetNbinsX();
    for (int i = 1; i <= nB; i++) {
        double nTr = h_nom->GetBinContent(i);
        double sumUp = h_up->GetBinContent(i);
        double sumDn = h_dn->GetBinContent(i);
        double sigma = (nTr > 0.0) ? 0.5 * std::abs(sumUp - sumDn) / nTr : 0.0;
        // Symmetric by construction -- up and down get the same value,
        // so the band and the bottom panel are both symmetric too.
        frac.up.push_back(sigma);
        frac.dn.push_back(sigma);
    }
    return frac;
}

TGraphAsymmErrors* BuildBand(TH1D* h_mc_display, const Fractions& frac, const char* name)
{
    TGraphAsymmErrors* gr = new TGraphAsymmErrors();
    gr->SetName(name);
    int nB = h_mc_display->GetNbinsX();
    for (int i = 1; i <= nB; i++) {
        double x   = h_mc_display->GetBinCenter(i);
        double yMC = h_mc_display->GetBinContent(i);
        double ex  = 0.5 * h_mc_display->GetBinWidth(i);
        gr->SetPoint(i-1, x, yMC);
        gr->SetPointError(i-1, ex, ex, frac.dn[i-1]*yMC, frac.up[i-1]*yMC);
    }
    gr->SetFillColorAlpha(kOrange+1, 0.30);
    gr->SetLineColor(kOrange+1);
    return gr;
}

TH1D* MakeFracHist(TH1D* h_ref, const std::vector<double>& frac, const char* name)
{
    TH1D* h = (TH1D*)h_ref->Clone(name);
    h->Reset();
    for (int i = 1; i <= h->GetNbinsX(); i++) h->SetBinContent(i, frac[i-1]);
    return h;
}

void DrawDuneLabel(double x) { TLatex tL; tL.SetNDC(); tL.DrawLatex(x, 0.94, LAB_DUNE); }

// =============================================================================
// Plot 1: MC-only (unscaled)
// =============================================================================
void PlotMCOnly(TH1D* hMC_unscaled, TGraphAsymmErrors* grBand, const std::string& outPng)
{
    TCanvas* c = new TCanvas("c_mconly", "", 800, 600);
    c->SetLeftMargin(0.14);

    hMC_unscaled->SetLineColor(kRed); hMC_unscaled->SetLineWidth(3);
    hMC_unscaled->SetTitle("Track Multiplicity");
    hMC_unscaled->GetXaxis()->SetTitle("Number of Reconstructed Tracks");
    hMC_unscaled->GetXaxis()->SetRangeUser(0, 10);
    hMC_unscaled->GetYaxis()->SetTitle("Number of Interactions");
    hMC_unscaled->GetYaxis()->SetTitleOffset(1.1);
    double ymax = 1.35 * hMC_unscaled->GetMaximum();
    hMC_unscaled->GetYaxis()->SetRangeUser(0, ymax);

    hMC_unscaled->Draw("HIST");
    grBand->Draw("2 SAME");
    hMC_unscaled->Draw("HIST SAME");

    TLegend* leg = new TLegend(LEG_X1, LEG_Y1, LEG_X2, LEG_Y2);
    leg->SetTextFont(43); leg->SetTextSize(18);
    leg->SetFillStyle(0); leg->SetBorderSize(0);
    leg->AddEntry(hMC_unscaled, "MC (unscaled)", "l");
    leg->AddEntry(grBand, "SI Systematics", "f");
    leg->Draw();
    DrawDuneLabel(0.15);

    c->SaveAs(outPng.c_str());
    delete c;
}

// =============================================================================
// Plot 2: Data + MC + band, no fractional panel
// =============================================================================
void PlotBandOnly(TH1D* hData, TH1D* hMC, TGraphAsymmErrors* grBand, const std::string& outPng)
{
    TCanvas* c = new TCanvas("c_band", "", 800, 600);
    c->SetLeftMargin(0.14);

    hData->SetMarkerStyle(20); hData->SetMarkerSize(1.2);
    hData->SetLineColor(kBlack); hData->SetMarkerColor(kBlack);
    hData->SetTitle("Track Multiplicity");
    hData->GetXaxis()->SetTitle("Number of Reconstructed Tracks");
    hData->GetXaxis()->SetRangeUser(0, 10);
    hData->GetYaxis()->SetTitle("Number of Interactions");
    hData->GetYaxis()->SetTitleOffset(1.1);
    hData->GetYaxis()->SetRangeUser(0, Y_MAX);

    hMC->SetLineColor(kRed); hMC->SetLineWidth(3);

    hData->Draw("E1");
    grBand->Draw("2 SAME");
    hMC->Draw("HIST SAME");
    hData->Draw("E1 SAME");

    TLegend* leg = new TLegend(LEG_X1, LEG_Y1, LEG_X2, LEG_Y2);
    leg->SetTextFont(43); leg->SetTextSize(18);
    leg->SetFillStyle(0); leg->SetBorderSize(0);
    leg->AddEntry(hData, "Data", "lep");
    leg->AddEntry(hMC,   "MC (norm. to data)", "l");
    leg->AddEntry(grBand,"SI Systematics", "f");
    leg->Draw();
    DrawDuneLabel(0.18);

    c->SaveAs(outPng.c_str());
    delete c;
}

// =============================================================================
// Plot 3: Data + MC + band + fractional panel -- single filled orange
// panel (frac.up == frac.dn, symmetric by construction), matching the
// nusyst_dataMC_trackMult.C / GENIE-systematics splitpad style.
// =============================================================================
void PlotBandWithFrac(TH1D* hData, TH1D* hMC, TGraphAsymmErrors* grBand,
                       TH1D* hFrac, const std::string& outPng)
{
    TCanvas* c = new TCanvas("c_bandfrac", "", 800, 750);

    TPad* p1 = new TPad("p1","",0,0.3,1,1.0);
    p1->SetLeftMargin(0.14); p1->SetRightMargin(0.05);
    p1->SetBottomMargin(0.02); p1->SetTopMargin(0.08);
    p1->Draw(); p1->cd();

    hData->SetMarkerStyle(20); hData->SetMarkerSize(1.2);
    hData->SetLineColor(kBlack); hData->SetMarkerColor(kBlack);
    hData->SetTitle("Track Multiplicity");
    hData->GetXaxis()->SetLabelSize(0); hData->GetXaxis()->SetTitleSize(0);
    hData->GetYaxis()->SetTitle("Number of Interactions");
    hData->GetYaxis()->SetTitleOffset(1.1);
    hData->GetYaxis()->SetRangeUser(0, Y_MAX);

    hMC->SetLineColor(kRed); hMC->SetLineWidth(3);

    hData->Draw("E1");
    grBand->Draw("2 SAME");
    hMC->Draw("HIST SAME");
    hData->Draw("E1 SAME");

    TLegend* leg = new TLegend(LEG_X1, LEG_Y1, LEG_X2, LEG_Y2);
    leg->SetTextFont(43); leg->SetTextSize(18);
    leg->SetFillStyle(0); leg->SetBorderSize(0);
    leg->AddEntry(hData, "Data", "lep");
    leg->AddEntry(hMC,   "MC (norm. to data)", "l");
    leg->AddEntry(grBand,"SI Systematics", "f");
    leg->Draw();
    DrawDuneLabel(0.15);

    c->cd();
    TPad* p2 = new TPad("p2","",0,0.0,1,0.3);
    p2->SetLeftMargin(0.14); p2->SetRightMargin(0.05);
    p2->SetTopMargin(0.02); p2->SetBottomMargin(0.35);
    p2->Draw(); p2->cd();

    hFrac->SetLineColor(kOrange+1);
    hFrac->SetFillColorAlpha(kOrange+1, 0.30);
    hFrac->SetLineWidth(2);
    hFrac->SetTitle("");
    hFrac->GetXaxis()->SetTitle("Number of Reconstructed Tracks");
    hFrac->GetYaxis()->SetTitle("#sigma/N");
    hFrac->GetXaxis()->SetRangeUser(0, 10);
    hFrac->GetXaxis()->SetLabelSize(0.12);
    hFrac->GetXaxis()->SetTitleSize(0.14);
    hFrac->GetXaxis()->SetTitleOffset(1.0);
    hFrac->GetYaxis()->SetLabelSize(0.10);
    hFrac->GetYaxis()->SetTitleSize(0.12);
    hFrac->GetYaxis()->SetTitleOffset(0.55);
    hFrac->GetYaxis()->SetNdivisions(505);
    double fmax = hFrac->GetMaximum();
    hFrac->GetYaxis()->SetRangeUser(0, 1.3*(fmax>0?fmax:0.3));
    hFrac->Draw("HIST");

    c->SaveAs(outPng.c_str());
    delete c;
}

// =============================================================================
// Main
// =============================================================================
void g4rw_dataMC_trackMult()
{
    gROOT->LoadMacro("protoDUNEStyle.C");
    gROOT->SetStyle("protoDUNEStyle");
    gROOT->ForceStyle();
    gStyle->SetTitleX(0.35);
    TH1::AddDirectory(false);
    gSystem->mkdir(OUT_DIR.c_str(), true);

    // ── Load weight maps ──────────────────────────────────────────────────
    std::map<WKey,W2> wmap_piplus, wmap_piminus, wmap_proton;
    LoadWeights(WFILE_PIPLUS,  wmap_piplus);
    LoadWeights(WFILE_PIMINUS, wmap_piminus);
    LoadWeights(WFILE_PROTON,  wmap_proton);

    // ── Data ──────────────────────────────────────────────────────────────
    TFile* fData = TFile::Open(DATA_FILE.c_str());
    if (!fData || fData->IsZombie()) { std::cerr << "ERROR: Cannot open data file\n"; return; }
    TH1D* hData = (TH1D*)fData->Get("track_mult");
    if (!hData) { std::cerr << "ERROR: track_mult not found\n"; return; }
    hData->SetDirectory(0);
    fData->Close();

    // ── MC: syst_g4rw_mult_reco ──────────────────────────────────────────
    TFile* fMC = TFile::Open(MC_FILE.c_str());
    if (!fMC || fMC->IsZombie()) { std::cerr << "ERROR: Cannot open MC file\n"; return; }
    TTree* tReco = (TTree*)fMC->Get("syst_g4rw_mult_reco");
    if (!tReco) { std::cerr << "ERROR: syst_g4rw_mult_reco not found\n"; return; }

    long long nu_id; int cafFileIdx, G4ID, truePDG;
    tReco->SetBranchAddress("nu_id",      &nu_id);
    tReco->SetBranchAddress("cafFileIdx", &cafFileIdx);
    tReco->SetBranchAddress("G4ID",       &G4ID);
    tReco->SetBranchAddress("truePDG",    &truePDG);

    // First pass: count tracks per interaction -> that count IS the
    // interaction's trackMult (one row per qualifying track, same
    // selection as the trackMult counting loop in nominal_mc.cc).
    using EventKey = std::pair<int,long long>;
    std::map<EventKey,int> eventMult;
    {
        long N = tReco->GetEntries();
        for (long i = 0; i < N; ++i) {
            tReco->GetEntry(i);
            eventMult[{cafFileIdx, nu_id}]++;
        }
    }
    std::cout << "Number of interactions: " << eventMult.size() << "\n";

    // ── Histograms ────────────────────────────────────────────────────────
    TH1D* h_nom_events = new TH1D("h_nom_events", "", nBinsM, xMinM, xMaxM); // one entry per interaction
    TH1D* h_nom_tracks = new TH1D("h_nom_tracks", "", nBinsM, xMinM, xMaxM); // one entry per hadron track
    TH1D* h_up_tracks  = new TH1D("h_up_tracks",  "", nBinsM, xMinM, xMaxM);
    TH1D* h_dn_tracks  = new TH1D("h_dn_tracks",  "", nBinsM, xMinM, xMaxM);

    for (const auto& kv : eventMult) h_nom_events->Fill(kv.second);

    // Second pass: one entry per TRACK, at the mult of the event it
    // belongs to, weighted by THAT TRACK'S OWN weight -- no combination
    // across tracks in the same event, mirroring the proton script's
    // h_up_len->Fill(length, w) treatment exactly.
    long N = tReco->GetEntries();
    long nTotalTracks = 0, nMatched = 0;
    for (long i = 0; i < N; ++i) {
        tReco->GetEntry(i);
        int mult = eventMult[{cafFileIdx, nu_id}];

        WKey key = {cafFileIdx, (int)nu_id, G4ID};
        double w_up = 1.0, w_dn = 1.0;
        nTotalTracks++;

        if (truePDG == 211) {
            auto it = wmap_piplus.find(key);
            if (it != wmap_piplus.end()) { w_up = it->second[0]; w_dn = it->second[1]; nMatched++; }
        } else if (truePDG == -211) {
            auto it = wmap_piminus.find(key);
            if (it != wmap_piminus.end()) { w_up = it->second[0]; w_dn = it->second[1]; nMatched++; }
        } else if (abs(truePDG) == 2212) {
            auto it = wmap_proton.find(key);
            if (it != wmap_proton.end()) { w_up = it->second[0]; w_dn = it->second[1]; nMatched++; }
        }
        // muon/kaon/other -> w stays 1.0, no weight file, but still
        // contributes a (mult, w=1) entry, same as everywhere else in
        // this workstream.

        h_nom_tracks->Fill(mult);
        h_up_tracks->Fill(mult, w_up);
        h_dn_tracks->Fill(mult, w_dn);
    }
    fMC->Close();

    std::cout << "Total hadron tracks: " << nTotalTracks
               << "  (matched to a weight file: " << nMatched << ")\n";

    // ── Fractions, from TRACK-level sums pooling all hadron species. ─────
    Fractions frac = ComputeFractions(h_nom_tracks, h_up_tracks, h_dn_tracks);

    std::cout << "\n=== Per-bin: independent-track mean weight vs. band fraction (symmetric envelope) ===\n";
    for (int i = 1; i <= nBinsM; i++) {
        double nTr = h_nom_tracks->GetBinContent(i);
        if (nTr <= 0) continue;
        double meanUp = h_up_tracks->GetBinContent(i) / nTr;
        double meanDn = h_dn_tracks->GetBinContent(i) / nTr;
        std::cout << "  mult=" << (i-1) << " : nTracks=" << nTr
                   << "  meanW(1.33)=" << meanUp << "  meanW(0.67)=" << meanDn
                   << "  frac=" << frac.up[i-1] << " (" << (100.0*frac.up[i-1]) << "%)\n";
    }

    // ── Scale to data ─────────────────────────────────────────────────────
    double scale = (h_nom_events->Integral() > 0.0) ? hData->Integral() / h_nom_events->Integral() : 1.0;
    TH1D* hMC_scaled = (TH1D*)h_nom_events->Clone("hMC_scaled");
    hMC_scaled->Scale(scale);
    std::cout << "\nNormalization scale factor: " << scale << "\n";

    // ── Band + fraction histograms ───────────────────────────────────────
    auto* grBand_scaled   = BuildBand(hMC_scaled,   frac, "gr_band_scaled");
    auto* grBand_unscaled = BuildBand(h_nom_events,  frac, "gr_band_unscaled");
    TH1D* hFrac = MakeFracHist(h_nom_events, frac.up, "hFrac");

    // ── Three plots ───────────────────────────────────────────────────────
    PlotMCOnly(h_nom_events, grBand_unscaled,
               OUT_DIR + "/g4rw_MConly_trackMult.png");

    PlotBandOnly(hData, hMC_scaled, grBand_scaled,
                 OUT_DIR + "/g4rw_dataMC_trackMult_band.png");

    PlotBandWithFrac(hData, hMC_scaled, grBand_scaled, hFrac,
                      OUT_DIR + "/g4rw_dataMC_trackMult_splitPad.png");

    std::cout << "\n========================================================\n";
    std::cout << "g4rw_dataMC_trackMult complete.\n";
    std::cout << "  Plot: " << OUT_DIR << "/g4rw_MConly_trackMult.png\n";
    std::cout << "  Plot: " << OUT_DIR << "/g4rw_dataMC_trackMult_band.png\n";
    std::cout << "  Plot: " << OUT_DIR << "/g4rw_dataMC_trackMult_splitPad.png\n";
    std::cout << "========================================================\n";
}

int main()
{
    g4rw_dataMC_trackMult();
    return 0;
}