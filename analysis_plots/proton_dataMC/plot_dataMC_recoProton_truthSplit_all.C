///////////////////////////////////////////////////////////////////////////////
// plot_dataMC_recoProton_truthSplit_all.C
//
// USES THE ORIGINAL BINNING APPROACH: nominal MC (hMC_cos/hMC_len, used
// as both the plotted red line AND the scale target for the band/lines)
// comes from CAF's own pre-binned (300-bin) reco_matched_proton_cosL/
// _length histograms, rebinned down to 15 bins -- same method as the
// original proton-study scripts. Kept this way deliberately for direct
// comparability with existing plots.
//
// h_nom_cos/h_nom_len are STILL filled locally (needed as the reference
// point for ComputeFractions -- the per-bin fraction is computed relative
// to the unweighted LOCAL count, same as always), but they are NOT used
// to build hMC_cos/hMC_len in this version.
//
// Band and lines are built from the SAME shared Fractions object either
// way, so lines still coincide with band edges by construction --
// that part of the earlier fix is independent of which histogram is
// used as the scale target, and is kept.
//
// Produces plots per variable (cosTheta, length):
//   1. Lines only:  Data, MC nominal, f=1.33 line (green+2), f=0.67 line (blue+1)
//   2. Band only:   Data, MC nominal, +-1sigma SI Systematics band
//   3. Band+lines:  Data, MC nominal, band, AND the f=1.33/f=0.67 lines together
//   4. NEW: Mean G4RW weight per bin (meanUp/meanDn vs. bin, with a
//      reference line at w=1) -- see PrintMeanWeightTable/DrawMeanWeightVsBin
//      below. Added to answer directly: "is the systematic just the mean
//      of the G4RW weight?" -- yes, PER BIN (not globally): frac.up ==
//      max(0, meanWeight_up - 1) and frac.dn == max(0, 1 - meanWeight_dn),
//      by construction from ComputeFractions. The table prints this
//      exact relationship bin-by-bin with a pass/fail check; the plot
//      shows meanUp/meanDn as curves so it's visible directly against
//      the band's shape (wider where local mean weight moves further
//      from 1, narrower where it sits close to 1).
//
// Compile and run:
//   g++ -std=c++17 -O2 -g plot_dataMC_recoProton_truthSplit_all.C -o plot_dataMC_recoProton_truthSplit_all `root-config --cflags --libs`
//   ./plot_dataMC_recoProton_truthSplit_all
///////////////////////////////////////////////////////////////////////////////

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TGraphAsymmErrors.h"

#include <map>
#include <tuple>
#include <array>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>

// =============================================================
// Paths
// =============================================================
const std::string CAF_FILE =
    "/global/cfs/cdirs/dune/users/mazam/softwares/patches/2x2_g4rw/validation_proton/codes_for_slides/proton_plots/output.root";

const std::string DATA_FILE =
    "/global/u1/m/mazam/sandbox/v11/sandbox_v11_minervaOff_geomContainRecoOnly_Lgt3cm_300cm_endZneg0to5cmVeto_cosZgtNeg0p9.root";

const std::string WFILE_PROTON =
    "/global/cfs/cdirs/dune/users/mazam/G4RW/output_weights/patched_withRock_genieNu_eq3p10/weights_1000_all_LAr_P_cascade_0_9000_5weights.root";

const std::string WFILE_PIPLUS =
    "/global/cfs/cdirs/dune/users/mazam/G4RW/output_weights/patched_withRock_genieNu_eq3p10/weights_1000_all_LAr_PiPlus_cascade_0_9000_5weights.root";

const std::string WFILE_PIMINUS =
    "/global/cfs/cdirs/dune/users/mazam/G4RW/output_weights/patched_withRock_genieNu_eq3p10/weights_1000_all_LAr_PiMinus_cascade_0_9000_5weights.root";

const std::string OUT_DIR =
    "/global/cfs/cdirs/dune/users/mazam/softwares/patches/2x2_g4rw/validation_proton/codes_for_slides/total_systematics";

// =============================================================
// Binning
// =============================================================
const int    kNBinsC = 15;
const double kXMinC  = -1.0;
const double kXMaxC  =  1.0;

const int    kNBinsL = 15;
const double kXMinL  =  3.0;
const double kXMaxL  = 83.0;

// =============================================================
// Weight map -- [0]=w_1p33, [1]=w_0p67
// =============================================================
using WKey = std::tuple<int,int,int>;
using W2   = std::array<double,2>;

void LoadWeights(const std::string& file, std::map<WKey,W2>& wmap)
{
    TFile* f = TFile::Open(file.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: Cannot open: " << file << "\n"; return;
    }
    TTree* t = (TTree*)f->Get("g4rwEvt");
    if (!t) {
        std::cerr << "ERROR: No g4rwEvt tree in: " << file << "\n";
        f->Close(); return;
    }
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

// =============================================================
// Rebin preserving errors
// =============================================================
TH1D* RebinPreserveErrors(TH1D* h, int nbins, double xmin, double xmax,
                           const char* name)
{
    TH1D* hR = new TH1D(name, "", nbins, xmin, xmax);
    hR->Sumw2();
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        double x = h->GetBinCenter(i);
        if (x < xmin || x >= xmax) continue;
        int b = hR->FindBin(x);
        hR->SetBinContent(b, hR->GetBinContent(b) + h->GetBinContent(i));
        hR->SetBinError(b, std::sqrt(hR->GetBinError(b)*hR->GetBinError(b)
                                   + h->GetBinError(i)*h->GetBinError(i)));
    }
    return hR;
}

// =============================================================
// Compute per-bin fractional up/down deviation from h_nom/h_up/h_dn
// (all UNSCALED, all from the SAME local fill loop). Single source of
// truth for BOTH the band and the lines below.
// =============================================================
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

// =============================================================
// Per-bin local mean weight, and its relationship to the band's
// per-bin fraction (Fractions.up/dn) -- printed as a table, and also
// plotted as its own curve, to directly answer:
//   "is the systematic just the mean of the G4RW weight?"
// Answer: YES, per bin (not globally) -- frac.up = max(0, meanUp-1),
// frac.dn = max(0, 1-meanDn), by construction from ComputeFractions
// above (sum_of_weights / count = mean weight of that bin).
// =============================================================
void PrintMeanWeightTable(const char* label, TH1D* h_nom, TH1D* h_up, TH1D* h_dn,
                           const Fractions& frac)
{
    std::cout << "\n=== Per-bin mean G4RW weight vs. band fraction: " << label << " ===\n";
    std::cout << std::setw(5)  << "Bin"
               << std::setw(15) << "N(unweighted)"
               << std::setw(15) << "MeanW(1.33)"
               << std::setw(15) << "MeanW(0.67)"
               << std::setw(12) << "frac.up"
               << std::setw(12) << "frac.dn"
               << std::setw(10) << "Match?" << "\n";

    int nB = h_nom->GetNbinsX();
    for (int i = 1; i <= nB; i++) {
        double N      = h_nom->GetBinContent(i);
        double meanUp = (N > 0.0) ? h_up->GetBinContent(i) / N : 0.0;
        double meanDn = (N > 0.0) ? h_dn->GetBinContent(i) / N : 0.0;

        // What frac.up/frac.dn SHOULD be, derived independently from the
        // mean weight -- compared against the actual band fraction below.
        double expectUp = std::max(0.0, meanUp - 1.0);
        double expectDn = std::max(0.0, 1.0 - meanDn);
        bool matchUp = std::abs(expectUp - frac.up[i-1]) < 1e-9;
        bool matchDn = std::abs(expectDn - frac.dn[i-1]) < 1e-9;

        std::cout << std::setw(5)  << i
                   << std::setw(15) << N
                   << std::setw(15) << meanUp
                   << std::setw(15) << meanDn
                   << std::setw(12) << frac.up[i-1]
                   << std::setw(12) << frac.dn[i-1]
                   << std::setw(10) << ((matchUp && matchDn) ? "yes" : "NO") << "\n";
    }
}

// Plot the local mean weight (meanUp = h_up/h_nom, meanDn = h_dn/h_nom)
// as a function of bin, with a horizontal reference line at w=1 (no
// reweighting). Directly comparable, bin-for-bin, to the SI band's width
// in the band_*_truthSplit.png plots -- wider band <-> mean further
// from 1, by the exact relationship printed in PrintMeanWeightTable.
void DrawMeanWeightVsBin(TH1D* h_nom, TH1D* h_up, TH1D* h_dn,
                          const char* xTitle, const char* title,
                          const std::string& outPng)
{
    int nB = h_nom->GetNbinsX();
    TH1D* hMeanUp = (TH1D*)h_nom->Clone("hMeanUp"); hMeanUp->Reset();
    TH1D* hMeanDn = (TH1D*)h_nom->Clone("hMeanDn"); hMeanDn->Reset();

    for (int i = 1; i <= nB; i++) {
        double N = h_nom->GetBinContent(i);
        hMeanUp->SetBinContent(i, (N > 0.0) ? h_up->GetBinContent(i) / N : 0.0);
        hMeanDn->SetBinContent(i, (N > 0.0) ? h_dn->GetBinContent(i) / N : 0.0);
    }

    TCanvas* c = new TCanvas("c_meanw", "c_meanw", 800, 600);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.15);

    hMeanUp->SetLineColor(kGreen+2); hMeanUp->SetLineWidth(3);
    hMeanDn->SetLineColor(kBlue+1);  hMeanDn->SetLineWidth(3);

    hMeanUp->SetTitle(title);
    hMeanUp->GetXaxis()->SetTitle(xTitle);
    hMeanUp->GetYaxis()->SetTitle("Local mean G4RW weight");

    double hiUp = hMeanUp->GetMaximum(), hiDn = hMeanDn->GetMaximum();
    double loUp = hMeanUp->GetMinimum(0.0), loDn = hMeanDn->GetMinimum(0.0);
    double ymax = 1.3 * std::max(hiUp, hiDn);
    double ymin = 0.7 * std::min(loUp > 0 ? loUp : 1.0, loDn > 0 ? loDn : 1.0);
    hMeanUp->GetYaxis()->SetRangeUser(std::min(0.0, ymin), ymax);

    hMeanUp->Draw("HIST");
    hMeanDn->Draw("HIST SAME");

    TLine* lOne = new TLine(h_nom->GetXaxis()->GetXmin(), 1.0,
                             h_nom->GetXaxis()->GetXmax(), 1.0);
    lOne->SetLineColor(kBlack); lOne->SetLineStyle(3); lOne->SetLineWidth(2);
    lOne->Draw("SAME");

    TLegend* leg = new TLegend(0.48, 0.68, 0.80, 0.93);
    leg->SetTextFont(133); leg->SetTextSize(17);
    leg->SetFillStyle(0); leg->SetBorderSize(0);
    leg->AddEntry(hMeanUp, "Mean weight, f=1.33", "l");
    leg->AddEntry(hMeanDn, "Mean weight, f=0.67", "l");
    leg->AddEntry(lOne,    "No reweighting (w=1)", "l");
    leg->Draw();

    TLatex tL; tL.SetNDC();
    tL.DrawLatex(0.20, 0.94, "#bf{DUNE:ND-LAr 2x2}");

    c->SaveAs(outPng.c_str());
    delete c;
}

// =============================================================
// Band, built from hMC_scaled and the shared fractions.
// =============================================================
TGraphAsymmErrors* BuildG4RWBand(TH1D* h_mc_scaled, const Fractions& frac)
{
    TGraphAsymmErrors* gr = new TGraphAsymmErrors();
    int nB = h_mc_scaled->GetNbinsX();
    for (int i = 1; i <= nB; i++) {
        double x   = h_mc_scaled->GetBinCenter(i);
        double yMC = h_mc_scaled->GetBinContent(i);
        double ex  = 0.5 * h_mc_scaled->GetBinWidth(i);
        gr->SetPoint(i-1, x, yMC);
        gr->SetPointError(i-1, ex, ex, frac.dn[i-1] * yMC, frac.up[i-1] * yMC);
    }
    gr->SetFillColorAlpha(kOrange+1, 0.50);
    gr->SetLineColor(kOrange+1);
    return gr;
}

// =============================================================
// Lines, built from hMC_scaled and the SAME shared fractions -- these
// are the band's edges, expressed as histograms. Guaranteed to coincide
// with the band by construction.
// =============================================================
TH1D* BuildLineFromFraction(TH1D* h_mc_scaled, const std::vector<double>& frac,
                             bool isUpper, const char* name)
{
    TH1D* h = (TH1D*)h_mc_scaled->Clone(name);
    int nB = h->GetNbinsX();
    for (int i = 1; i <= nB; i++) {
        double yMC = h_mc_scaled->GetBinContent(i);
        double y = isUpper ? yMC * (1.0 + frac[i-1]) : yMC * (1.0 - frac[i-1]);
        h->SetBinContent(i, y);
    }
    return h;
}

// =============================================================
// Plot 1: Lines only -- Data, MC nominal, f=1.33 (green+2), f=0.67 (blue+1)
// =============================================================
void DrawLinesOnly(TH1D* hData, TH1D* hMC, TH1D* hUp, TH1D* hDn,
                    const char* xTitle, const char* title,
                    const std::string& outPng, double ymax)
{
    TCanvas* c = new TCanvas("c_lines", "c_lines", 800, 600);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.15);

    hData->SetLineColor(kBlack);      hData->SetLineWidth(3);
    hMC->SetLineColor(kRed);          hMC->SetLineWidth(3);
    hUp->SetLineColor(kGreen+2);      hUp->SetLineWidth(3); hUp->SetLineStyle(2);
    hDn->SetLineColor(kBlue+1);       hDn->SetLineWidth(3); hDn->SetLineStyle(2);

    hData->SetTitle(title);
    hData->GetXaxis()->SetTitle(xTitle);
    hData->GetYaxis()->SetTitle("Counts");
    hData->GetYaxis()->SetRangeUser(0, ymax);

    hData->Draw("HIST E");
    hMC->Draw("HIST SAME");
    hUp->Draw("HIST SAME");
    hDn->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.48, 0.68, 0.80, 0.93);
    leg->SetTextFont(133);
    leg->SetTextSize(17);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hData, "Data Reco",                 "l");
    leg->AddEntry(hMC,   "MC Reco (norm. to Data)",   "l");
    leg->AddEntry(hUp,   "G4RW, f=1.33",               "l");
    leg->AddEntry(hDn,   "G4RW, f=0.67",               "l");
    leg->Draw();

    TLatex tL; tL.SetNDC();
    tL.DrawLatex(0.20, 0.94, "#bf{DUNE:ND-LAr 2x2}");

    c->SaveAs(outPng.c_str());
    delete c;
}

// =============================================================
// Plot 2: Band only -- Data, MC nominal, SI Systematics band
// =============================================================
void DrawBandOnly(TH1D* hData, TH1D* hMC, TGraphAsymmErrors* grBand,
                   const char* xTitle, const char* title,
                   const std::string& outPng, double ymax)
{
    TCanvas* c = new TCanvas("c_band", "c_band", 800, 600);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.15);

    hData->SetLineColor(kBlack); hData->SetLineWidth(3);
    hMC->SetLineColor(kRed);     hMC->SetLineWidth(3);

    hData->SetTitle(title);
    hData->GetXaxis()->SetTitle(xTitle);
    hData->GetYaxis()->SetTitle("Counts");
    hData->GetYaxis()->SetRangeUser(0, ymax);

    hData->Draw("HIST E");
    grBand->Draw("2 SAME");
    hMC->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.50, 0.75, 0.80, 0.93);
    leg->SetTextFont(133);
    leg->SetTextSize(18);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hData,  "Data Reco",                        "l");
    leg->AddEntry(hMC,    "MC Reco (norm. to Data)",          "l");
    leg->AddEntry(grBand, "SI Systematics (#pm1#sigma)",       "f");
    leg->Draw();

    TLatex tL; tL.SetNDC();
    tL.DrawLatex(0.20, 0.94, "#bf{DUNE:ND-LAr 2x2}");

    c->SaveAs(outPng.c_str());
    delete c;
}

// =============================================================
// Plot 3: Band + lines together
// =============================================================
void DrawBandPlusLines(TH1D* hData, TH1D* hMC, TGraphAsymmErrors* grBand,
                        TH1D* hUp, TH1D* hDn,
                        const char* xTitle, const char* title,
                        const std::string& outPng, double ymax)
{
    TCanvas* c = new TCanvas("c_bandlines", "c_bandlines", 800, 600);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.15);

    hData->SetLineColor(kBlack); hData->SetLineWidth(3);
    hMC->SetLineColor(kRed);     hMC->SetLineWidth(3);
    hUp->SetLineColor(kGreen+2); hUp->SetLineWidth(3); hUp->SetLineStyle(2);
    hDn->SetLineColor(kBlue+1);  hDn->SetLineWidth(3); hDn->SetLineStyle(2);

    hData->SetTitle(title);
    hData->GetXaxis()->SetTitle(xTitle);
    hData->GetYaxis()->SetTitle("Counts");
    hData->GetYaxis()->SetRangeUser(0, ymax);

    hData->Draw("HIST E");
    grBand->Draw("2 SAME");
    hMC->Draw("HIST SAME");
    hUp->Draw("HIST SAME");
    hDn->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.46, 0.60, 0.80, 0.93);
    leg->SetTextFont(133);
    leg->SetTextSize(16);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hData,  "Data Reco",                        "l");
    leg->AddEntry(hMC,    "MC Reco (norm. to Data)",          "l");
    leg->AddEntry(grBand, "SI Systematics (#pm1#sigma)",       "f");
    leg->AddEntry(hUp,    "G4RW, f=1.33",                      "l");
    leg->AddEntry(hDn,    "G4RW, f=0.67",                      "l");
    leg->Draw();

    TLatex tL; tL.SetNDC();
    tL.DrawLatex(0.20, 0.94, "#bf{DUNE:ND-LAr 2x2}");

    c->SaveAs(outPng.c_str());
    delete c;
}

// =============================================================
// Main
// =============================================================
void plot_dataMC_recoProton_truthSplit_all()
{
    gROOT->LoadMacro("protoDUNEStyle.C");
    gROOT->SetStyle("protoDUNEStyle");
    gROOT->ForceStyle();
    gStyle->SetTitleX(0.25);
    gSystem->mkdir(OUT_DIR.c_str(), true);

    // ── Load the three weight maps ───────────────────────────────────────
    std::map<WKey,W2> wmap_proton, wmap_piplus, wmap_piminus;
    LoadWeights(WFILE_PROTON,  wmap_proton);
    LoadWeights(WFILE_PIPLUS,  wmap_piplus);
    LoadWeights(WFILE_PIMINUS, wmap_piminus);

    // ── Open CAF and data files ──────────────────────────────────────────
    TFile* fCAF = TFile::Open(CAF_FILE.c_str());
    if (!fCAF || fCAF->IsZombie()) { std::cerr << "ERROR: Cannot open CAF file\n"; return; }
    TTree* tReco = (TTree*)fCAF->Get("syst_g4rw_proton_reco");
    if (!tReco) { std::cerr << "ERROR: Cannot find syst_g4rw_proton_reco tree\n"; return; }

    TFile* fData = TFile::Open(DATA_FILE.c_str());
    if (!fData || fData->IsZombie()) { std::cerr << "ERROR: Cannot open data file\n"; return; }
    TH1D* hDataCos_raw = (TH1D*)fData->Get("recoProtonCosLAll");
    TH1D* hDataLen_raw = (TH1D*)fData->Get("recoProtonLengthAll");
    if (!hDataCos_raw || !hDataLen_raw) { std::cerr << "ERROR: Missing data histograms\n"; return; }
    TH1D* hDataCos = RebinPreserveErrors(hDataCos_raw, kNBinsC, kXMinC, kXMaxC, "hDataCos");
    TH1D* hDataLen = RebinPreserveErrors(hDataLen_raw, kNBinsL, kXMinL, kXMaxL, "hDataLen");

    // ── ORIGINAL BINNING APPROACH: CAF's own pre-binned histograms,
    //    rebinned down -- this is what hMC_cos/hMC_len will be built
    //    from below, NOT h_nom. ──────────────────────────────────────────
    TH1D* hMC_cos_raw = RebinPreserveErrors((TH1D*)fCAF->Get("reco_matched_proton_cosL"), kNBinsC, kXMinC, kXMaxC, "hMC_cos_raw");
    TH1D* hMC_len_raw = RebinPreserveErrors((TH1D*)fCAF->Get("reco_matched_proton_length"), kNBinsL, kXMinL, kXMaxL, "hMC_len_raw");

    // ── h_nom / h_up / h_dn (unscaled), ALL filled locally in the SAME loop
    //    -- h_nom is used ONLY as the reference point for ComputeFractions
    //    below (the fraction is still computed relative to the local,
    //    unweighted count), NOT as the source for hMC_cos/hMC_len. ──────
    TH1D* h_nom_cos = new TH1D("h_nom_cos", "", kNBinsC, kXMinC, kXMaxC);
    TH1D* h_up_cos  = new TH1D("h_up_cos",  "", kNBinsC, kXMinC, kXMaxC);
    TH1D* h_dn_cos  = new TH1D("h_dn_cos",  "", kNBinsC, kXMinC, kXMaxC);

    TH1D* h_nom_len = new TH1D("h_nom_len", "", kNBinsL, kXMinL, kXMaxL);
    TH1D* h_up_len  = new TH1D("h_up_len",  "", kNBinsL, kXMinL, kXMaxL);
    TH1D* h_dn_len  = new TH1D("h_dn_len",  "", kNBinsL, kXMinL, kXMaxL);

    int cafFileIdx, G4ID, truePDG;
    long long nu_id;
    double cosT, length;

    tReco->SetBranchAddress("cafFileIdx", &cafFileIdx);
    tReco->SetBranchAddress("nu_id",      &nu_id);
    tReco->SetBranchAddress("G4ID",       &G4ID);
    tReco->SetBranchAddress("truePDG",    &truePDG);
    tReco->SetBranchAddress("cosTheta",   &cosT);
    tReco->SetBranchAddress("length",     &length);

    long N = tReco->GetEntries();
    long nProton = 0, nPiPlus = 0, nPiMinus = 0, nOther = 0;
    std::cout << "Reco entries: " << N << "\n";

    // ── Capture every entry landing in LENGTH BIN 1, for a worked-example
    //    printout of the band/line math below. Bin edges computed from
    //    kXMinL/kXMaxL/kNBinsL, same convention TH1D itself uses. ────────
    double binWidthL = (kXMaxL - kXMinL) / kNBinsL;
    double bin1LoL = kXMinL, bin1HiL = kXMinL + binWidthL;
    struct Bin1Entry { int truePDG; double length, w_up, w_dn; std::string species; };
    std::vector<Bin1Entry> bin1Entries;

    for (long i = 0; i < N; i++) {
        tReco->GetEntry(i);
        WKey key = {cafFileIdx, (int)nu_id, G4ID};

        W2 w = {1.0, 1.0};
        std::string species = "other/muon";
        if (truePDG == 2212) {
            if (wmap_proton.count(key)) w = wmap_proton[key];
            nProton++; species = "proton";
        } else if (truePDG == 211) {
            if (wmap_piplus.count(key)) w = wmap_piplus[key];
            nPiPlus++; species = "pi+";
        } else if (truePDG == -211) {
            if (wmap_piminus.count(key)) w = wmap_piminus[key];
            nPiMinus++; species = "pi-";
        } else {
            nOther++;
        }

        h_nom_cos->Fill(cosT);
        h_up_cos->Fill(cosT, w[0]);  h_dn_cos->Fill(cosT, w[1]);

        h_nom_len->Fill(length);
        h_up_len->Fill(length, w[0]); h_dn_len->Fill(length, w[1]);

        if (length >= bin1LoL && length < bin1HiL) {
            bin1Entries.push_back({truePDG, length, w[0], w[1], species});
        }
    }

    std::cout << "\n--- Truth-species breakdown ---\n";
    std::cout << "Proton (2212): " << nProton  << " / " << N << "\n";
    std::cout << "pi+     (211): " << nPiPlus  << " / " << N << "\n";
    std::cout << "pi-    (-211): " << nPiMinus << " / " << N << "\n";
    std::cout << "Muon/Other (w=1, no G4RW weight applied): " << nOther << " / " << N << "\n";

    // ── Worked example: length bin 1 -- print only the first few entries,
    //    but sum over ALL of them so the downstream math stays correct. ──
    std::cout << "\n=== Worked example: LENGTH BIN 1 (" << bin1LoL << " - " << bin1HiL << " cm) ===\n";
    std::cout << bin1Entries.size() << " entries total in this bin (showing first 5):\n";
    std::cout << std::setw(8) << "truePDG" << std::setw(10) << "species" << std::setw(10) << "length"
               << std::setw(12) << "w_1.33" << std::setw(12) << "w_0.67" << "\n";
    double sumUp = 0, sumDn = 0;
    int nPrinted = 0;
    for (auto& e : bin1Entries) {
        if (nPrinted < 5) {
            std::cout << std::setw(8) << e.truePDG << std::setw(10) << e.species << std::setw(10) << e.length
                       << std::setw(12) << e.w_up << std::setw(12) << e.w_dn << "\n";
            nPrinted++;
        }
        sumUp += e.w_up; sumDn += e.w_dn;
    }
    if ((int)bin1Entries.size() > 5) std::cout << "  ... (" << (bin1Entries.size() - 5) << " more not shown, but included in the sums below)\n";
    std::cout << "\nSum of w_1.33 over this bin (= h_up_len bin 1) = " << sumUp << "\n";
    std::cout << "Sum of w_0.67 over this bin (= h_dn_len bin 1) = " << sumDn << "\n";
    std::cout << "Count in this bin, unweighted (= h_nom_len bin 1) = " << bin1Entries.size() << "\n";

    // ── Nominal MC used for the plotted red line and as the scale target.
    //    ORIGINAL BINNING APPROACH: uses CAF's own pre-binned histogram
    //    (hMC_cos_raw/hMC_len_raw), NOT h_nom. ──────────────────────────
    double scaleCos = (hMC_cos_raw->Integral("width") > 0.0)
                     ? hDataCos->Integral("width") / hMC_cos_raw->Integral("width")
                     : 1.0;
    TH1D* hMC_cos = (TH1D*)hMC_cos_raw->Clone("hMC_cos");
    hMC_cos->Scale(scaleCos);

    double scaleLen = (hMC_len_raw->Integral("width") > 0.0)
                     ? hDataLen->Integral("width") / hMC_len_raw->Integral("width")
                     : 1.0;
    TH1D* hMC_len = (TH1D*)hMC_len_raw->Clone("hMC_len");
    hMC_len->Scale(scaleLen);

    // ── Verification print: h_nom (local, unweighted) vs hMC_len_raw
    //    (CAF's own histogram) -- these are NOT expected to match exactly
    //    in this version (that's the known, accepted binning-method
    //    difference), shown here just so the two are visible side by side. ─
    std::cout << "\n--- h_nom (local, unweighted) vs hMC_len_raw (CAF's own histogram), per bin ---\n";
    std::cout << "cosTheta bins:\n";
    for (int i = 1; i <= kNBinsC; i++)
        std::cout << "  bin " << i << ": h_nom=" << h_nom_cos->GetBinContent(i)
                   << "  CAF hist=" << hMC_cos_raw->GetBinContent(i) << "\n";
    std::cout << "length bins:\n";
    for (int i = 1; i <= kNBinsL; i++)
        std::cout << "  bin " << i << ": h_nom=" << h_nom_len->GetBinContent(i)
                   << "  CAF hist=" << hMC_len_raw->GetBinContent(i) << "\n";

    std::cout << "\nscaleCos = " << scaleCos << ", scaleLen = " << scaleLen << "\n";
    std::cout << "hDataCos total (width-integral): " << hDataCos->Integral("width") << "\n";
    std::cout << "hMC_cos total (scaled, width-integral): " << hMC_cos->Integral("width") << "\n";
    std::cout << "hDataLen total (width-integral): " << hDataLen->Integral("width") << "\n";
    std::cout << "hMC_len total (scaled, width-integral): " << hMC_len->Integral("width") << "\n";

    // ── Fractions computed ONCE, shared by band and lines -- STILL uses
    //    h_nom as the fraction's reference point (unchanged). ─────────────
    Fractions frac_cos = ComputeFractions(h_nom_cos, h_up_cos, h_dn_cos);
    Fractions frac_len = ComputeFractions(h_nom_len, h_up_len, h_dn_len);

    // ── NEW: per-bin mean weight vs. band fraction, printed table + plot.
    //    Answers directly: is the systematic just the mean G4RW weight?
    //    Yes, PER BIN -- see function comments above. ─────────────────────
    PrintMeanWeightTable("cosTheta", h_nom_cos, h_up_cos, h_dn_cos, frac_cos);
    PrintMeanWeightTable("length",   h_nom_len, h_up_len, h_dn_len, frac_len);

    DrawMeanWeightVsBin(h_nom_cos, h_up_cos, h_dn_cos,
                        "cos(#theta)", "Protons: Mean G4RW Weight per Bin",
                        OUT_DIR + "/meanWeight_cosTheta_truthSplit.png");
    DrawMeanWeightVsBin(h_nom_len, h_up_len, h_dn_len,
                        "Track Length [cm]", "Protons: Mean G4RW Weight per Bin",
                        OUT_DIR + "/meanWeight_length_truthSplit.png");

    // ── Finish the worked example: show exactly how bin 1's sums above
    //    turn into the plotted band edges / lines. NOTE: hMC_len bin 1
    //    now comes from hMC_len_raw (CAF's histogram) x scaleLen, NOT
    //    from h_nom_len x scaleLen -- updated below to reflect that. ─────
    std::cout << "\n--- Continuing worked example: bin 1 math ---\n";
    std::cout << "frac_up = (sum_w1.33 - h_nom) / h_nom = ("
               << h_up_len->GetBinContent(1) << " - " << h_nom_len->GetBinContent(1) << ") / "
               << h_nom_len->GetBinContent(1) << " = " << frac_len.up[0] << "\n";
    std::cout << "frac_dn = (h_nom - sum_w0.67) / h_nom = ("
               << h_nom_len->GetBinContent(1) << " - " << h_dn_len->GetBinContent(1) << ") / "
               << h_nom_len->GetBinContent(1) << " = " << frac_len.dn[0] << "\n";
    std::cout << "scaleLen (data/MC overall normalization) = " << scaleLen << "\n";
    std::cout << "hMC_len bin 1 (= hMC_len_raw[CAF hist] bin 1 x scaleLen) = "
               << hMC_len_raw->GetBinContent(1) << " x " << scaleLen << " = " << hMC_len->GetBinContent(1) << "\n";
    std::cout << "\nBand/line upper edge = hMC_len x (1+frac_up) = "
               << hMC_len->GetBinContent(1) << " x (1+" << frac_len.up[0] << ") = "
               << hMC_len->GetBinContent(1) * (1.0 + frac_len.up[0]) << "  <- this is the green f=1.33 line's height\n";
    std::cout << "Band/line lower edge = hMC_len x (1-frac_dn) = "
               << hMC_len->GetBinContent(1) << " x (1-" << frac_len.dn[0] << ") = "
               << hMC_len->GetBinContent(1) * (1.0 - frac_len.dn[0]) << "  <- this is the blue f=0.67 line's height\n";

    auto* gr_cos = BuildG4RWBand(hMC_cos, frac_cos);
    auto* gr_len = BuildG4RWBand(hMC_len, frac_len);

    TH1D* h_up_cos_line = BuildLineFromFraction(hMC_cos, frac_cos.up, true,  "h_up_cos_line");
    TH1D* h_dn_cos_line = BuildLineFromFraction(hMC_cos, frac_cos.dn, false, "h_dn_cos_line");
    TH1D* h_up_len_line = BuildLineFromFraction(hMC_len, frac_len.up, true,  "h_up_len_line");
    TH1D* h_dn_len_line = BuildLineFromFraction(hMC_len, frac_len.dn, false, "h_dn_len_line");

    // Separate clones for the lines-only plot vs. the band+lines plot,
    // since ROOT histograms carry draw styling as object state.
    TH1D* h_up_cos_l1 = (TH1D*)h_up_cos_line->Clone("h_up_cos_l1");
    TH1D* h_dn_cos_l1 = (TH1D*)h_dn_cos_line->Clone("h_dn_cos_l1");
    TH1D* h_up_len_l1 = (TH1D*)h_up_len_line->Clone("h_up_len_l1");
    TH1D* h_dn_len_l1 = (TH1D*)h_dn_len_line->Clone("h_dn_len_l1");

    double ymaxCos = 50, ymaxLen = 70;

    // ── Plot 1: lines only ─────────────────────────────────────────────
    DrawLinesOnly(hDataCos, hMC_cos, h_up_cos_l1, h_dn_cos_l1,
                  "cos(#theta)", "Protons",
                  OUT_DIR + "/lines_cosTheta_truthSplit.png", ymaxCos);
    DrawLinesOnly(hDataLen, hMC_len, h_up_len_l1, h_dn_len_l1,
                  "Track Length [cm]", "Protons",
                  OUT_DIR + "/lines_length_truthSplit.png", ymaxLen);

    // ── Plot 2: band only ─────────────────────────────────────────────
    DrawBandOnly(hDataCos, hMC_cos, gr_cos,
                 "cos(#theta)", "Protons",
                 OUT_DIR + "/band_cosTheta_truthSplit.png", ymaxCos);
    DrawBandOnly(hDataLen, hMC_len, gr_len,
                 "Track Length [cm]", "Protons",
                 OUT_DIR + "/band_length_truthSplit.png", ymaxLen);

    // ── Plot 3: band + lines (h_up/h_dn_line used directly -- these ARE
    //    the band's edges, no separate cloning needed for correctness) ────
    DrawBandPlusLines(hDataCos, hMC_cos, gr_cos, h_up_cos_line, h_dn_cos_line,
                       "cos(#theta)", "Protons",
                       OUT_DIR + "/bandLines_cosTheta_truthSplit.png", ymaxCos);
    DrawBandPlusLines(hDataLen, hMC_len, gr_len, h_up_len_line, h_dn_len_line,
                       "Track Length [cm]", "Protons",
                       OUT_DIR + "/bandLines_length_truthSplit.png", ymaxLen);

    std::cout << "\nDone. All eight plots saved in:\n" << OUT_DIR << "\n";
    fCAF->Close();
    fData->Close();
}

int main()
{
    plot_dataMC_recoProton_truthSplit_all();
    return 0;
}