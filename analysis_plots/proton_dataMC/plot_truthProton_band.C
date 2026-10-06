///////////////////////////////////////////////////////////////////////////////
// plot_truthProton_band.C
//
// Truth-level counterpart to plot_dataMC_recoProton_truthSplit_all.C.
// Reads syst_g4rw_proton (truth-matched backtracked protons: nu_id,
// cafFileIdx, G4ID, cosTheta, length -- already proton-only, no truePDG
// branch needed) instead of syst_g4rw_proton_reco.
//
// No data overlay (truth has no data equivalent) and no data-normalization
// scale factor -- hMC_cos/hMC_len are the RAW truth-level nominal counts.
// Otherwise identical band-building logic (ComputeFractions / BuildG4RWBand)
// to the reco version, using the SAME weight map (weights are intrinsically
// truth-level, keyed by cafFileIdx/nu_id/G4ID).
//
// NEW: PrintFractionTable prints the per-bin frac.up/frac.dn numbers for
// both cosTheta and length, plus a weight-lookup hit/miss counter, so the
// band's actual size can be checked numerically instead of judged by eye
// off the plot (important for cosTheta, where the band is visually small
// in most bins).
//
// NEW: also produces a RECO-level counterpart plot, using the reco cosTheta/
// length of the SAME accepted, correctly-backtracked proton candidates
// (syst_g4rw_proton_reco, filtered to truePDG==2212 -- the 2257-entry
// population, same as the truth tree). Raw MC counts, no normalization to
// data (matching the truth plot's convention, and making the two directly
// comparable: same underlying tracks, same SI weights, just reco vs. true
// kinematics on the x-axis).
//
// NEW: NON-PROTON reco proton candidates (matched true pi-, pi+, other
// species), track length only, raw counts, same band construction as the
// proton plots above:
//   reco  : syst_g4rw_proton_reco, truePDG != 2212 && truePDG != 0
//   truth : syst_g4rw_nonproton (same candidates, true track length)
// Weights are by TRUE species (pi+ -> PiPlus file, pi- -> PiMinus file;
// other matched species keep w = 1, no weight file). Candidates with
// truePDG == 0 (no truth match) have no truth counterpart, so they are
// excluded from both non-proton plots; both plots show the same tracks.
//
// Compile and run:
//   g++ -std=c++17 -O2 -g plot_truthProton_band.C -o plot_truthProton_band `root-config --cflags --libs`
//   ./plot_truthProton_band
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
#include "TGraphAsymmErrors.h"

#include <map>
#include <tuple>
#include <array>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <string>

// =============================================================
// Paths -- same CAF file and weight file as the reco version
// =============================================================
const std::string CAF_FILE =
    "/global/cfs/cdirs/dune/users/mazam/softwares/patches/2x2_g4rw/validation_proton/codes_for_slides/proton_plots/output.root";

const std::string WFILE_PROTON =
    "/global/cfs/cdirs/dune/users/mazam/G4RW/output_weights/patched_withRock_genieNu_eq3p10/weights_1000_all_LAr_P_cascade_0_9000_5weights.root";

// NEW: pion weight files, used only by the non-proton plots at the end
const std::string WFILE_PIPLUS =
    "/global/cfs/cdirs/dune/users/mazam/G4RW/output_weights/patched_withRock_genieNu_eq3p10/weights_1000_all_LAr_PiPlus_cascade_0_9000_5weights.root";

const std::string WFILE_PIMINUS =
    "/global/cfs/cdirs/dune/users/mazam/G4RW/output_weights/patched_withRock_genieNu_eq3p10/weights_1000_all_LAr_PiMinus_cascade_0_9000_5weights.root";

const std::string OUT_DIR =
    "/global/cfs/cdirs/dune/users/mazam/softwares/patches/2x2_g4rw/validation_proton/codes_for_slides/total_systematics";

// =============================================================
// Binning -- same as reco version, for direct visual comparability
// =============================================================
const int    kNBinsC = 15;
const double kXMinC  = -1.0;
const double kXMaxC  =  1.0;

const int    kNBinsL = 15;
const double kXMinL  =  3.0;
const double kXMaxL  = 83.0;

using WKey = std::tuple<int,int,int>;
using W2   = std::array<double,2>;

void LoadWeights(const std::string& file, std::map<WKey,W2>& wmap)
{
    TFile* f = TFile::Open(file.c_str());
    if (!f || f->IsZombie()) { std::cerr << "ERROR: Cannot open: " << file << "\n"; return; }
    TTree* t = (TTree*)f->Get("g4rwEvt");
    if (!t) { std::cerr << "ERROR: No g4rwEvt tree in: " << file << "\n"; f->Close(); return; }
    int file_idx, event_id, track_id;
    double w_1p33, w_0p67;
    t->SetBranchAddress("file_idx",    &file_idx);
    t->SetBranchAddress("event_id",    &event_id);
    t->SetBranchAddress("track_id",    &track_id);
    t->SetBranchAddress("weight_1p33", &w_1p33);
    t->SetBranchAddress("weight_0p67", &w_0p67);
    long N = t->GetEntries();
    for (long i = 0; i < N; i++) {
        t->GetEntry(i);
        wmap[{file_idx, event_id, track_id}] = {w_1p33, w_0p67};
    }
    f->Close();
    std::cout << "Loaded " << N << " entries from: " << file << "\n";
}

// NEW: weight by TRUE species. pi+ -> PiPlus file, pi- -> PiMinus file,
// any other matched species keeps w = 1 (no weight file).
W2 SpeciesWeight(int pdg, const WKey& key,
                 const std::map<WKey,W2>& mPip, const std::map<WKey,W2>& mPim,
                 long& nPiFound, long& nPiMissing)
{
    W2 w = {1.0, 1.0};
    const std::map<WKey,W2>* m = nullptr;
    if (pdg == 211) m = &mPip; else if (pdg == -211) m = &mPim;
    if (m) {
        auto it = m->find(key);
        if (it != m->end()) { w = it->second; nPiFound++; } else nPiMissing++;
    }
    return w;
}

struct Fractions { std::vector<double> up, dn; };

Fractions ComputeFractions(TH1D* h_nom, TH1D* h_up, TH1D* h_dn)
{
    Fractions frac;
    int nB = h_nom->GetNbinsX();
    for (int i = 1; i <= nB; i++) {
        double yNom = h_nom->GetBinContent(i);
        double eup = std::max(0.0, h_up->GetBinContent(i) - yNom);
        double edn = std::max(0.0, yNom - h_dn->GetBinContent(i));
        frac.up.push_back((yNom > 0.0) ? eup / yNom : 0.0);
        frac.dn.push_back((yNom > 0.0) ? edn / yNom : 0.0);
    }
    return frac;
}

// NEW: prints the per-bin nominal count, summed weights, and resulting
// fractional up/down deviation -- lets you check the actual numbers behind
// the band instead of judging its width by eye off the plot.
void PrintFractionTable(const char* label, TH1D* h_nom, TH1D* h_up, TH1D* h_dn,
                         const Fractions& frac)
{
    std::cout << "\n=== Truth-level per-bin fractions: " << label << " ===\n";
    std::cout << std::setw(5)  << "Bin"
               << std::setw(12) << "N(nom)"
               << std::setw(12) << "sum(w_up)"
               << std::setw(12) << "sum(w_dn)"
               << std::setw(12) << "frac.up"
               << std::setw(12) << "frac.dn" << "\n";
    int nB = h_nom->GetNbinsX();
    for (int i = 1; i <= nB; i++) {
        std::cout << std::setw(5) << i
                   << std::setw(12) << h_nom->GetBinContent(i)
                   << std::setw(12) << h_up->GetBinContent(i)
                   << std::setw(12) << h_dn->GetBinContent(i)
                   << std::setw(12) << frac.up[i-1]
                   << std::setw(12) << frac.dn[i-1] << "\n";
    }
}

TGraphAsymmErrors* BuildG4RWBand(TH1D* h_nom, const Fractions& frac)
{
    TGraphAsymmErrors* gr = new TGraphAsymmErrors();
    int nB = h_nom->GetNbinsX();
    for (int i = 1; i <= nB; i++) {
        double x   = h_nom->GetBinCenter(i);
        double y   = h_nom->GetBinContent(i);
        double ex  = 0.5 * h_nom->GetBinWidth(i);
        gr->SetPoint(i-1, x, y);
        gr->SetPointError(i-1, ex, ex, frac.dn[i-1] * y, frac.up[i-1] * y);
    }
    gr->SetFillColorAlpha(kOrange+1, 0.50);
    gr->SetLineColor(kOrange+1);
    return gr;
}

// legendLabel lets the truth and reco calls use different nominal-line
// wording ("MC Truth (backtracked)" vs "MC Reco") while sharing the same
// drawing code; defaults to the truth wording so existing calls need no change.
void DrawTruthBand(TH1D* hNom, TGraphAsymmErrors* grBand,
                    const char* xTitle, const char* title,
                    const std::string& outPng, double ymax,
                    const char* legendLabel = "MC Truth (backtracked)")
{
    TCanvas* c = new TCanvas("c_truthband", "c_truthband", 800, 600);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.15);

    hNom->SetLineColor(kRed); hNom->SetLineWidth(3);
    hNom->SetTitle(title);
    hNom->GetXaxis()->SetTitle(xTitle);
    hNom->GetYaxis()->SetTitle("Counts");
    hNom->GetYaxis()->SetRangeUser(0, ymax);

    hNom->Draw("HIST");
    grBand->Draw("2 SAME");
    hNom->Draw("HIST SAME");   // redraw line on top of the band fill

    TLegend* leg = new TLegend(0.50, 0.75, 0.80, 0.93);
    leg->SetTextFont(133);
    leg->SetTextSize(18);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hNom,   legendLabel, "l");
    leg->AddEntry(grBand, "SI Systematics (#pm1#sigma)", "f");
    leg->Draw();

    TLatex tL; tL.SetNDC();
    tL.DrawLatex(0.20, 0.94, "#bf{DUNE:ND-LAr 2x2}");

    c->SaveAs(outPng.c_str());
    delete c;
}

void plot_truthProton_band()
{
    gROOT->LoadMacro("protoDUNEStyle.C");
    gROOT->SetStyle("protoDUNEStyle");
    gROOT->ForceStyle();
    gStyle->SetTitleX(0.25);
    gSystem->mkdir(OUT_DIR.c_str(), true);

    std::map<WKey,W2> wmap_proton;
    LoadWeights(WFILE_PROTON, wmap_proton);

    std::map<WKey,W2> wmap_piplus, wmap_piminus;   // NEW
    LoadWeights(WFILE_PIPLUS,  wmap_piplus);
    LoadWeights(WFILE_PIMINUS, wmap_piminus);

    TFile* fCAF = TFile::Open(CAF_FILE.c_str());
    if (!fCAF || fCAF->IsZombie()) { std::cerr << "ERROR: Cannot open CAF file\n"; return; }

    // =============================================================
    // TRUTH-LEVEL (unchanged)
    // =============================================================
    TTree* tTruth = (TTree*)fCAF->Get("syst_g4rw_proton");
    if (!tTruth) { std::cerr << "ERROR: Cannot find syst_g4rw_proton tree\n"; return; }

    TH1D* h_nom_cos = new TH1D("h_nom_truth_cos", "", kNBinsC, kXMinC, kXMaxC);
    TH1D* h_up_cos  = new TH1D("h_up_truth_cos",  "", kNBinsC, kXMinC, kXMaxC);
    TH1D* h_dn_cos  = new TH1D("h_dn_truth_cos",  "", kNBinsC, kXMinC, kXMaxC);

    TH1D* h_nom_len = new TH1D("h_nom_truth_len", "", kNBinsL, kXMinL, kXMaxL);
    TH1D* h_up_len  = new TH1D("h_up_truth_len",  "", kNBinsL, kXMinL, kXMaxL);
    TH1D* h_dn_len  = new TH1D("h_dn_truth_len",  "", kNBinsL, kXMinL, kXMaxL);

    int cafFileIdx, G4ID;
    long long nu_id;
    double cosT, length;

    tTruth->SetBranchAddress("cafFileIdx", &cafFileIdx);
    tTruth->SetBranchAddress("nu_id",      &nu_id);
    tTruth->SetBranchAddress("G4ID",       &G4ID);
    tTruth->SetBranchAddress("cosTheta",   &cosT);
    tTruth->SetBranchAddress("length",     &length);

    long N = tTruth->GetEntries();
    std::cout << "Truth-matched proton entries: " << N << "\n";

    long nFound = 0, nMissing = 0;   // NEW: track weight-lookup hit rate

    for (long i = 0; i < N; i++) {
        tTruth->GetEntry(i);
        WKey key = {cafFileIdx, (int)nu_id, G4ID};

        W2 w = {1.0, 1.0};
        auto it = wmap_proton.find(key);
        if (it != wmap_proton.end()) { w = it->second; nFound++; }   // NEW
        else { nMissing++; }                                         // NEW

        h_nom_cos->Fill(cosT);
        h_up_cos->Fill(cosT, w[0]);  h_dn_cos->Fill(cosT, w[1]);

        h_nom_len->Fill(length);
        h_up_len->Fill(length, w[0]); h_dn_len->Fill(length, w[1]);
    }

    std::cout << "Weight lookup: " << nFound << " found, " << nMissing   // NEW
              << " missing (of " << N << ")\n";                          // NEW

    Fractions frac_cos = ComputeFractions(h_nom_cos, h_up_cos, h_dn_cos);
    Fractions frac_len = ComputeFractions(h_nom_len, h_up_len, h_dn_len);

    PrintFractionTable("truth cosTheta", h_nom_cos, h_up_cos, h_dn_cos, frac_cos);   // NEW
    PrintFractionTable("truth length",   h_nom_len, h_up_len, h_dn_len, frac_len);   // NEW

    auto* gr_cos = BuildG4RWBand(h_nom_cos, frac_cos);
    auto* gr_len = BuildG4RWBand(h_nom_len, frac_len);

    double ymaxCos = 1.3 * h_nom_cos->GetMaximum();
    double ymaxLen = 1.3 * h_nom_len->GetMaximum();

    DrawTruthBand(h_nom_cos, gr_cos, "cos(#theta)", "Truth Protons",
                  OUT_DIR + "/band_truthProtonCosTheta.png", ymaxCos);
    DrawTruthBand(h_nom_len, gr_len, "Track Length [cm]", "Truth Protons",
                  OUT_DIR + "/band_truthProtonLength.png", ymaxLen);

    // =============================================================
    // NEW: RECO-LEVEL, raw counts, no data normalization -- SAME accepted,
    // correctly-backtracked candidates as the truth plot above (truePDG==2212
    // filter gives the identical 2257-entry population), just plotted against
    // their reco cosTheta/length instead of true cosTheta/length. Directly
    // comparable to the truth plot: same tracks, same SI weights.
    // =============================================================
    TTree* tRecoAll = (TTree*)fCAF->Get("syst_g4rw_proton_reco");
    if (!tRecoAll) { std::cerr << "ERROR: Cannot find syst_g4rw_proton_reco tree\n"; return; }

    TH1D* h_nom_reco_cos = new TH1D("h_nom_reco_cos", "", kNBinsC, kXMinC, kXMaxC);
    TH1D* h_up_reco_cos  = new TH1D("h_up_reco_cos",  "", kNBinsC, kXMinC, kXMaxC);
    TH1D* h_dn_reco_cos  = new TH1D("h_dn_reco_cos",  "", kNBinsC, kXMinC, kXMaxC);

    TH1D* h_nom_reco_len = new TH1D("h_nom_reco_len", "", kNBinsL, kXMinL, kXMaxL);
    TH1D* h_up_reco_len  = new TH1D("h_up_reco_len",  "", kNBinsL, kXMinL, kXMaxL);
    TH1D* h_dn_reco_len  = new TH1D("h_dn_reco_len",  "", kNBinsL, kXMinL, kXMaxL);

    int recoCafFileIdx, recoG4ID, truePDG;
    long long recoNuId;
    double recoCosT, recoLength;

    tRecoAll->SetBranchAddress("cafFileIdx", &recoCafFileIdx);
    tRecoAll->SetBranchAddress("nu_id",      &recoNuId);
    tRecoAll->SetBranchAddress("G4ID",       &recoG4ID);
    tRecoAll->SetBranchAddress("truePDG",    &truePDG);
    tRecoAll->SetBranchAddress("cosTheta",   &recoCosT);
    tRecoAll->SetBranchAddress("length",     &recoLength);

    long nRecoTotal = tRecoAll->GetEntries();
    long nRecoProtonMatched = 0;
    long nRecoFound = 0, nRecoMissing = 0;

    for (long i = 0; i < nRecoTotal; i++) {
        tRecoAll->GetEntry(i);
        if (truePDG != 2212) continue;   // same population as the truth tree
        nRecoProtonMatched++;

        WKey key = {recoCafFileIdx, (int)recoNuId, recoG4ID};

        W2 w = {1.0, 1.0};
        auto it = wmap_proton.find(key);
        if (it != wmap_proton.end()) { w = it->second; nRecoFound++; }
        else { nRecoMissing++; }

        h_nom_reco_cos->Fill(recoCosT);
        h_up_reco_cos->Fill(recoCosT, w[0]);  h_dn_reco_cos->Fill(recoCosT, w[1]);

        h_nom_reco_len->Fill(recoLength);
        h_up_reco_len->Fill(recoLength, w[0]); h_dn_reco_len->Fill(recoLength, w[1]);
    }

    std::cout << "\nReco proton candidates (truePDG==2212): " << nRecoProtonMatched
              << " / " << nRecoTotal << " total accepted candidates\n";
    std::cout << "Weight lookup (reco): " << nRecoFound << " found, " << nRecoMissing
              << " missing (of " << nRecoProtonMatched << ")\n";

    Fractions frac_reco_cos = ComputeFractions(h_nom_reco_cos, h_up_reco_cos, h_dn_reco_cos);
    Fractions frac_reco_len = ComputeFractions(h_nom_reco_len, h_up_reco_len, h_dn_reco_len);

    PrintFractionTable("reco cosTheta", h_nom_reco_cos, h_up_reco_cos, h_dn_reco_cos, frac_reco_cos);
    PrintFractionTable("reco length",   h_nom_reco_len, h_up_reco_len, h_dn_reco_len, frac_reco_len);

    auto* gr_reco_cos = BuildG4RWBand(h_nom_reco_cos, frac_reco_cos);
    auto* gr_reco_len = BuildG4RWBand(h_nom_reco_len, frac_reco_len);

    double ymaxRecoCos = 1.3 * h_nom_reco_cos->GetMaximum();
    double ymaxRecoLen = 1.3 * h_nom_reco_len->GetMaximum();

    DrawTruthBand(h_nom_reco_cos, gr_reco_cos, "cos(#theta)", "Reco Protons",
                  OUT_DIR + "/band_recoProtonCosTheta_raw.png", ymaxRecoCos, "MC Reco");
    DrawTruthBand(h_nom_reco_len, gr_reco_len, "Track Length [cm]", "Reco Protons",
                  OUT_DIR + "/band_recoProtonLength_raw.png", ymaxRecoLen, "MC Reco");

    // =============================================================
    // NEW: NON-PROTON reco proton candidates (matched true pi-, pi+, other
    // species), track length only, raw counts, same band construction as the
    // proton plots above.
    //   reco  : syst_g4rw_proton_reco, truePDG != 2212 && truePDG != 0
    //   truth : syst_g4rw_nonproton (same candidates, true track length)
    // Candidates with truePDG == 0 (no truth match) have no truth counterpart,
    // so they are excluded from both plots and both show the same tracks.
    // =============================================================
    TTree* tNP = (TTree*)fCAF->Get("syst_g4rw_nonproton");
    if (!tNP) { std::cerr << "ERROR: Cannot find syst_g4rw_nonproton tree\n"; return; }

    TH1D* h_nom_np_reco_len  = new TH1D("h_nom_np_reco_len",  "", kNBinsL, kXMinL, kXMaxL);
    TH1D* h_up_np_reco_len   = new TH1D("h_up_np_reco_len",   "", kNBinsL, kXMinL, kXMaxL);
    TH1D* h_dn_np_reco_len   = new TH1D("h_dn_np_reco_len",   "", kNBinsL, kXMinL, kXMaxL);
    TH1D* h_nom_np_truth_len = new TH1D("h_nom_np_truth_len", "", kNBinsL, kXMinL, kXMaxL);
    TH1D* h_up_np_truth_len  = new TH1D("h_up_np_truth_len",  "", kNBinsL, kXMinL, kXMaxL);
    TH1D* h_dn_np_truth_len  = new TH1D("h_dn_np_truth_len",  "", kNBinsL, kXMinL, kXMaxL);

    // ---- reco pass (tRecoAll branches are still bound from the proton section) ----
    long npProton = 0, npUnmatched = 0, npMatched = 0, npPim = 0, npPip = 0;
    long npRecoPiFound = 0, npRecoPiMissing = 0;
    for (long i = 0; i < nRecoTotal; i++) {
        tRecoAll->GetEntry(i);
        if (truePDG == 2212) { npProton++;    continue; }
        if (truePDG == 0) npUnmatched++;   // no truth match: kept in the reco plot, w = 1
        else              npMatched++;
        npMatched++;
        if (truePDG == -211) npPim++; else if (truePDG == 211) npPip++;

        WKey key = {recoCafFileIdx, (int)recoNuId, recoG4ID};
        W2 w = SpeciesWeight(truePDG, key, wmap_piplus, wmap_piminus, npRecoPiFound, npRecoPiMissing);

        h_nom_np_reco_len->Fill(recoLength);
        h_up_np_reco_len->Fill(recoLength, w[0]);
        h_dn_np_reco_len->Fill(recoLength, w[1]);
    }

    // ---- truth pass ----
    int npCafFileIdx, npG4ID, npTruePDG;
    long long npNuId;
    double npCosT, npLength;
    tNP->SetBranchAddress("cafFileIdx", &npCafFileIdx);
    tNP->SetBranchAddress("nu_id",      &npNuId);
    tNP->SetBranchAddress("G4ID",       &npG4ID);
    tNP->SetBranchAddress("truePDG",    &npTruePDG);
    tNP->SetBranchAddress("cosTheta",   &npCosT);
    tNP->SetBranchAddress("length",     &npLength);

    long npTruth = tNP->GetEntries();
    long npTruthPiFound = 0, npTruthPiMissing = 0;
    for (long i = 0; i < npTruth; i++) {
        tNP->GetEntry(i);
        WKey key = {npCafFileIdx, (int)npNuId, npG4ID};
        W2 w = SpeciesWeight(npTruePDG, key, wmap_piplus, wmap_piminus, npTruthPiFound, npTruthPiMissing);

        h_nom_np_truth_len->Fill(npLength);
        h_up_np_truth_len->Fill(npLength, w[0]);
        h_dn_np_truth_len->Fill(npLength, w[1]);
    }

    // ---- checks ----
    long npOther = npMatched - npPim - npPip;
    std::cout << "\n=== Composition of the " << nRecoTotal << " reco proton candidates ===\n"
              << "  proton         : " << npProton    << "\n"
              << "  pi-            : " << npPim       << "\n"
              << "  pi+            : " << npPip       << "\n"
              << "  other species  : " << npOther     << "  (w = 1, no weight file)\n"
              << "  no truth match : " << npUnmatched << "  (w = 1, excluded from the non-proton plots)\n";
    std::cout << "Matched non-proton candidates: reco tree " << npMatched << ", truth tree " << npTruth
              << (npMatched == npTruth ? "  (consistent)" : "  ** MISMATCH **") << "\n";
    std::cout << "Pion weight lookup (reco):  " << npRecoPiFound  << " found, " << npRecoPiMissing  << " missing\n";
    std::cout << "Pion weight lookup (truth): " << npTruthPiFound << " found, " << npTruthPiMissing << " missing\n";
    std::cout << "Outside [" << kXMinL << ", " << kXMaxL << "] cm (under+overflow): reco "
              << h_nom_np_reco_len->GetBinContent(0)  + h_nom_np_reco_len->GetBinContent(kNBinsL + 1)
              << ", truth "
              << h_nom_np_truth_len->GetBinContent(0) + h_nom_np_truth_len->GetBinContent(kNBinsL + 1) << "\n";

    // ---- band, table, plots (same functions as the proton plots) ----
    Fractions frac_np_reco_len  = ComputeFractions(h_nom_np_reco_len,  h_up_np_reco_len,  h_dn_np_reco_len);
    Fractions frac_np_truth_len = ComputeFractions(h_nom_np_truth_len, h_up_np_truth_len, h_dn_np_truth_len);

    PrintFractionTable("non-proton reco length",  h_nom_np_reco_len,  h_up_np_reco_len,  h_dn_np_reco_len,  frac_np_reco_len);
    PrintFractionTable("non-proton truth length", h_nom_np_truth_len, h_up_np_truth_len, h_dn_np_truth_len, frac_np_truth_len);

    auto* gr_np_reco_len  = BuildG4RWBand(h_nom_np_reco_len,  frac_np_reco_len);
    auto* gr_np_truth_len = BuildG4RWBand(h_nom_np_truth_len, frac_np_truth_len);

    double ymaxNPReco  = 1.3 * h_nom_np_reco_len->GetMaximum();
    double ymaxNPTruth = 1.3 * h_nom_np_truth_len->GetMaximum();

    DrawTruthBand(h_nom_np_reco_len, gr_np_reco_len, "Track Length [cm]", "Reco Non-Protons",
                  OUT_DIR + "/band_recoNonProtonLength_raw.png", ymaxNPReco, "MC Reco");
    DrawTruthBand(h_nom_np_truth_len, gr_np_truth_len, "Track Length [cm]", "Truth Non-Protons",
                  OUT_DIR + "/band_truthNonProtonLength.png", ymaxNPTruth);

    std::cout << "\nDone. Truth- and reco-level band plots (protons and non-protons) saved in:\n" << OUT_DIR << "\n";
    fCAF->Close();
}

int main()
{
    plot_truthProton_band();
    return 0;
}