///////////////////////////////////////////////////////////////////////////////
// plot_dataMC_proton_combined_band.C
//
// Data vs MC comparison with three systematic bands:
//   - GENIE (sum-of-dials, symmetric)
//   - SI / G4RW patched (symmetric, same treatment as GENIE)
//   - Detector (MR6.5 baseline + 23 variations, RMS band)
//   - Total = GENIE (+) SI (+) Detector added in quadrature
//
// Part of the total_systematics workstream, alongside:
//   plot_dataMC_proton_GENIE.C                   (GENIE only, both observables)
//   plot_dataMC_proton_detector_band_MR65sys.C   (detector only, both observables)
//   plot_dataMC_recoProton_truthSplit_all.C      (SI/G4RW only, truth-split)
//
// What changed vs. the original combined script:
//   - CAF_FILE (nominal MC): old 9_May_2026/slides/proton_systematics/output.root
//     -> canonical proton_plots/output.root, same nominal used everywhere
//     else in this workstream.
//   - Detector-systematics band: old MiniRun6.1 production (output_mr6p1.root
//     + 16 detector_systematics/physics_systematics variation files,
//     including an E_field/V_drift systematic that does not exist in the
//     new production) -> new MiniRun6.5 baseline + 23-variation production
//     (anl_2026_samples/), same set used in
//     plot_dataMC_proton_detector_band_MR65sys.C. Since the detector-band
//     nominal is now the SAME file as CAF_FILE, hMC_cos_raw/hMC_len_raw
//     (already loaded from CAF_FILE) are reused directly as the band
//     nominal -- no separate file open needed.
//   - SI/G4RW band: old single "patched" weight file
//     (patched_withRock_genieNu/..._P_cascade_0_4000_5weights.root) applied
//     to EVERY reco entry regardless of true species -> species-split
//     weighting (proton/pi+/pi- via truePDG, everything else w=1) using the
//     current eq3p10 weight files, matching
//     plot_dataMC_recoProton_truthSplit_all.C exactly.
//   - OUT_DIR: updated to total_systematics/, matching the other scripts
//     in this workstream.
//   - GENIE band (SYST_FILE, gr_band_sumDials_p*) and DATA_FILE unchanged.
//
// Two output plots (Total band only, matching the original script's active
// DrawPlot calls -- DrawCombinedPlot, which also overlays the individual
// GENIE/SI bands, is kept below and can be swapped in the same way it
// was commented out originally):
//   band_dataMC_cosTheta_total.png
//   band_dataMC_length_total.png
//
// Also prints a per-bin table (stdout) of each systematic's band half-width
// (GENIE, SI, Detector) and the Total, in both absolute counts and
// fractional %, for both observables.
//
// Compile and run (inside shifter + DUNE env, same as the other scripts):
//   g++ -std=c++17 -O2 -g plot_dataMC_proton_combined_band.C -o plot_dataMC_proton_combined_band `root-config --cflags --libs`
//   ./plot_dataMC_proton_combined_band
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
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <string>
#include <vector>

// =============================================================
// Paths
// =============================================================

// Canonical nominal MC (1E19 POT) -- also the CAF tree source for the
// SI/G4RW species-split reweighting below, and the detector-band nominal.
const std::string CAF_FILE =
    "/global/cfs/cdirs/dune/users/mazam/softwares/patches/2x2_g4rw/validation_proton/codes_for_slides/proton_plots/output.root";

const std::string DATA_FILE =
    "/global/u1/m/mazam/sandbox/v11/sandbox_v11_minervaOff_geomContainRecoOnly_Lgt3cm_300cm_endZneg0to5cmVeto_cosZgtNeg0p9.root";

const std::string SYST_FILE =
    "/global/u1/m/mazam/Geant4ReWeighting/2x2_G4RW/test/on_CAFS/step2_reco/nusyst_step2_reco_output.root";

// SI/G4RW: species-split weight files (current eq3p10 version), matching
// plot_dataMC_recoProton_truthSplit_all.C exactly.
const std::string WFILE_PROTON =
    "/global/cfs/cdirs/dune/users/mazam/G4RW/output_weights/patched_withRock_genieNu_eq3p10/weights_1000_all_LAr_P_cascade_0_9000_5weights.root";

const std::string WFILE_PIPLUS =
    "/global/cfs/cdirs/dune/users/mazam/G4RW/output_weights/patched_withRock_genieNu_eq3p10/weights_1000_all_LAr_PiPlus_cascade_0_9000_5weights.root";

const std::string WFILE_PIMINUS =
    "/global/cfs/cdirs/dune/users/mazam/G4RW/output_weights/patched_withRock_genieNu_eq3p10/weights_1000_all_LAr_PiMinus_cascade_0_9000_5weights.root";

const std::string OUT_DIR =
    "/global/cfs/cdirs/dune/users/mazam/softwares/patches/2x2_g4rw/validation_proton/codes_for_slides/total_systematics";

// Detector systematics: new MiniRun6.5 baseline + 23-variation production.
// Nominal is CAF_FILE itself (reused, see hMC_cos_raw/hMC_len_raw below) --
// no separate nominal file here, unlike the old MiniRun6.1 block.
const std::string SYS_ROOT =
    "/global/u1/m/mazam/MiniRuns/detector_systematics/anl_2026_samples";

const std::vector<std::string> VAR_FILES = {
    SYS_ROOT + "/Birks_A_b/1_down/MR6p5_Birks_A_b_1_down_output.root",
    SYS_ROOT + "/Birks_A_b/1_up/MR6p5_Birks_A_b_1_up_output.root",

    SYS_ROOT + "/Birks_k_b/1_down/MR6p5_Birks_k_b_1_down_output.root",
    SYS_ROOT + "/Birks_k_b/1_up/MR6p5_Birks_k_b_1_up_output.root",

    // SYS_ROOT + "/Box_Model/alpha_down/MR6p5_Box_Model_alpha_down_output.root",
    // SYS_ROOT + "/Box_Model/alpha_up/MR6p5_Box_Model_alpha_up_output.root",
    // SYS_ROOT + "/Box_Model/beta_down/MR6p5_Box_Model_beta_down_output.root",
    // SYS_ROOT + "/Box_Model/beta_up/MR6p5_Box_Model_beta_up_output.root",
    // SYS_ROOT + "/Box_Model/Switch/MR6p5_Box_Model_Switch_output.root",

    // SYS_ROOT + "/e_lifetime/1.25/MR6p5_e_lifetime_1.25_output.root",
    // SYS_ROOT + "/e_lifetime/1.5/MR6p5_e_lifetime_1.5_output.root",
    // SYS_ROOT + "/e_lifetime/1.75/MR6p5_e_lifetime_1.75_output.root",
    // SYS_ROOT + "/e_lifetime/2/MR6p5_e_lifetime_2_output.root",

    SYS_ROOT + "/Long_Diff/1_down/MR6p5_Long_Diff_1_down_output.root",
    SYS_ROOT + "/Long_Diff/1_up/MR6p5_Long_Diff_1_up_output.root",

    SYS_ROOT + "/Temperature/1_down/MR6p5_Temperature_1_down_output.root",
    SYS_ROOT + "/Temperature/1_up/MR6p5_Temperature_1_up_output.root",

    // SYS_ROOT + "/Tran_Diff/10p_down/MR6p5_Tran_Diff_10p_down_output.root",
    // SYS_ROOT + "/Tran_Diff/10p_up/MR6p5_Tran_Diff_10p_up_output.root",
    // SYS_ROOT + "/Tran_Diff/25p_down/MR6p5_Tran_Diff_25p_down_output.root",
    // SYS_ROOT + "/Tran_Diff/25p_up/MR6p5_Tran_Diff_25p_up_output.root",

    SYS_ROOT + "/W_ion/1_down/MR6p5_W_ion_1_down_output.root",
    SYS_ROOT + "/W_ion/1_up/MR6p5_W_ion_1_up_output.root",
};

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
// Weight map -- [0]=w_1p33, [1]=w_0p67. Matches
// plot_dataMC_recoProton_truthSplit_all.C's LoadWeights exactly.
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
    hR->SetDirectory(0);
    return hR;
}

// =============================================================
// Build GENIE band (symmetric) -- unchanged
// =============================================================
TGraphAsymmErrors* BuildGENIEBand(TGraphAsymmErrors* grBand, TH1D* hMC)
{
    TGraphAsymmErrors* gr = new TGraphAsymmErrors();
    for (int i = 0; i < grBand->GetN(); ++i) {
        double x, yNom;
        grBand->GetPoint(i, x, yNom);
        int b = hMC->FindBin(x);
        double yMC  = hMC->GetBinContent(b);
        double sigma = grBand->GetErrorYhigh(i);
        double frac  = (yNom > 0.0) ? sigma / yNom : 0.0;
        double exl   = grBand->GetErrorXlow(i);
        double exh   = grBand->GetErrorXhigh(i);
        gr->SetPoint(i, x, yMC);
        gr->SetPointError(i, exl, exh, frac*yMC, frac*yMC);
    }
    gr->SetFillColorAlpha(kOrange+1, 0.50);
    gr->SetLineColor(kOrange+1);
    return gr;
}

// =============================================================
// Build SI band (symmetric) -- h_nom/h_up/h_dn now come from the
// SPECIES-SPLIT fill loop below (proton/pi+/pi- weight maps via truePDG),
// same treatment as plot_dataMC_recoProton_truthSplit_all.C's
// ComputeFractions + BuildG4RWBand, collapsed into one function since
// this script doesn't need the separate up/down lines.
// sigma = 0.5*(eup + edn), applied symmetrically.
// =============================================================
TGraphAsymmErrors* BuildSIBand(TH1D* h_nom, TH1D* h_up, TH1D* h_dn,
                                 TH1D* hMC)
{
    TGraphAsymmErrors* gr = new TGraphAsymmErrors();
    int nB = hMC->GetNbinsX();
    for (int i = 1; i <= nB; i++) {
        double x    = hMC->GetBinCenter(i);
        double yMC  = hMC->GetBinContent(i);
        double ex   = 0.5 * hMC->GetBinWidth(i);
        double yNom = h_nom->GetBinContent(i);
        double eup  = std::max(0.0, h_up->GetBinContent(i) - yNom);
        double edn  = std::max(0.0, yNom - h_dn->GetBinContent(i));
        double sigma = 0.5 * (eup + edn);
        double frac  = (yNom > 0.0) ? sigma / yNom : 0.0;
        gr->SetPoint(i-1, x, yMC);
        gr->SetPointError(i-1, ex, ex, frac*yMC, frac*yMC);
    }
    gr->SetFillColorAlpha(kBlue+1, 0.35);
    gr->SetLineColor(kBlue+1);
    return gr;
}

// =============================================================
// Build detector band -- h_nomBase is now CAF_FILE's own
// reco_matched_proton_cosL/length (reused via hMC_cos_raw/hMC_len_raw),
// h_vars are the 23 new MiniRun6.5 variations. Formula unchanged.
// =============================================================
TGraphAsymmErrors* BuildDetectorBand(TH1D* h_nomBase, const std::vector<TH1D*>& h_vars,
                                      TH1D* hMC)
{
    int Nvars = h_vars.size();
    TGraphAsymmErrors* gr = new TGraphAsymmErrors();
    int nB = hMC->GetNbinsX();
    for (int i = 1; i <= nB; i++) {
        double x    = hMC->GetBinCenter(i);
        double yDisp = hMC->GetBinContent(i);
        double ex   = 0.5 * hMC->GetBinWidth(i);
        double yNom = h_nomBase->GetBinContent(i);
        double sum = 0.0;
        for (int v = 0; v < Nvars; ++v) {
            double delta = h_vars[v]->GetBinContent(i) - yNom;
            sum += delta * delta;
        }
        double sigmaNom = (Nvars > 0) ? std::sqrt(sum / Nvars) : 0.0;
        double frac      = (yNom > 0.0) ? sigmaNom / yNom : 0.0;
        gr->SetPoint(i-1, x, yDisp);
        gr->SetPointError(i-1, ex, ex, frac*yDisp, frac*yDisp);
    }
    gr->SetFillColorAlpha(kGreen+2, 0.35);
    gr->SetLineColor(kGreen+2);
    return gr;
}

// =============================================================
// Build Total band: GENIE (+) SI (+) Detector in quadrature -- unchanged
// =============================================================
TGraphAsymmErrors* BuildTotalBand(TGraphAsymmErrors* grGENIE,
                                   TGraphAsymmErrors* grSI,
                                   TGraphAsymmErrors* grDet)
{
    TGraphAsymmErrors* gr = new TGraphAsymmErrors();
    int nB = std::min({grGENIE->GetN(), grSI->GetN(), grDet->GetN()});
    for (int i = 0; i < nB; i++) {
        double x, y;
        grGENIE->GetPoint(i, x, y);
        double sig_genie = grGENIE->GetErrorYhigh(i);
        double sig_si   = grSI->GetErrorYhigh(i);
        double sig_det   = grDet->GetErrorYhigh(i);
        double sig_tot   = std::sqrt(sig_genie*sig_genie + sig_si*sig_si + sig_det*sig_det);
        double exl = grGENIE->GetErrorXlow(i);
        double exh = grGENIE->GetErrorXhigh(i);
        gr->SetPoint(i, x, y);
        gr->SetPointError(i, exl, exh, sig_tot, sig_tot);
    }
    gr->SetFillColorAlpha(kOrange+1, 0.50);
    gr->SetLineColor(kOrange+1);
    return gr;
}

// =============================================================
// Print a per-bin table: absolute (+/-) band width and fractional (%)
// width for each systematic type, and the total. "Band width" here is
// the band's half-width (sigma) at each bin -- since every band here is
// symmetric, the full shaded width is just 2x this value.
// =============================================================
void PrintBandTable(const std::string& label,
                     TGraphAsymmErrors* grGENIE,
                     TGraphAsymmErrors* grSI,
                     TGraphAsymmErrors* grDet,
                     TGraphAsymmErrors* grTotal,
                     double xmin, double binWidth)
{
    std::cout << "\n=== Band width table: " << label
              << " (sigma = band half-width; full shaded width = 2 x sigma) ===\n";
    std::cout << std::fixed << std::setprecision(3);
    std::cout << std::setw(4)  << "Bin"
              << std::setw(16) << "Range"
              << std::setw(13) << "GENIE sig"
              << std::setw(10) << "GENIE %"
              << std::setw(13) << "SI sig"
              << std::setw(10) << "SI %"
              << std::setw(13) << "Det sig"
              << std::setw(10) << "Det %"
              << std::setw(13) << "Total sig"
              << std::setw(10) << "Total %"
              << "\n";

    int nB = std::min({grGENIE->GetN(), grSI->GetN(), grDet->GetN(), grTotal->GetN()});
    for (int i = 0; i < nB; i++) {
        double x, yG, yF, yD, yT;
        grGENIE->GetPoint(i, x, yG);
        grSI->GetPoint(i, x, yF);
        grDet->GetPoint(i, x, yD);
        grTotal->GetPoint(i, x, yT);

        double sG = grGENIE->GetErrorYhigh(i);
        double sF = grSI->GetErrorYhigh(i);
        double sD = grDet->GetErrorYhigh(i);
        double sT = grTotal->GetErrorYhigh(i);

        double pG = (yG > 0.0) ? 100.0 * sG / yG : 0.0;
        double pF = (yF > 0.0) ? 100.0 * sF / yF : 0.0;
        double pD = (yD > 0.0) ? 100.0 * sD / yD : 0.0;
        double pT = (yT > 0.0) ? 100.0 * sT / yT : 0.0;

        double lo = xmin + i * binWidth;
        double hi = lo + binWidth;
        std::ostringstream rangeStr;
        rangeStr << std::fixed << std::setprecision(1) << "[" << lo << ", " << hi << ")";

        std::cout << std::setw(4)  << (i + 1)
                  << std::setw(16) << rangeStr.str()
                  << std::setprecision(3)
                  << std::setw(13) << sG
                  << std::setw(9)  << pG << "%"
                  << std::setw(13) << sF
                  << std::setw(9)  << pF << "%"
                  << std::setw(13) << sD
                  << std::setw(9)  << pD << "%"
                  << std::setw(13) << sT
                  << std::setw(9)  << pT << "%"
                  << "\n";
    }
}

// =============================================================
// Draw plot: Data, MC, Total band -- unchanged
// =============================================================
void DrawPlot(TH1D* hData, TH1D* hMC,
              TGraphAsymmErrors* grBand,
              const char* xTitle,
              const char* title,
              const char* bandLabel,
              const std::string& outPng,
              double ymax = -1)
{
    TCanvas* c = new TCanvas("c", "c", 800, 600);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.15);

    hData->SetLineColor(kBlack);  hData->SetLineWidth(3);
    hMC->SetLineColor(kRed);      hMC->SetLineWidth(3);

    hData->SetTitle(title);
    hData->GetXaxis()->SetTitle(xTitle);
    hData->GetYaxis()->SetTitle("Counts");

    double maxY = (ymax > 0) ? ymax : 1.4 * std::max(hData->GetMaximum(), hMC->GetMaximum());
    hData->GetYaxis()->SetRangeUser(0, maxY);

    hData->Draw("HIST E");
    grBand->Draw("2 SAME");
    hMC->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.45, 0.75, 0.82, 0.93);
    leg->SetTextFont(133); leg->SetTextSize(18);
    leg->SetFillStyle(0);  leg->SetBorderSize(0);
    leg->AddEntry(hData,   "Data Reco",               "l");
    leg->AddEntry(hMC,     "MC Reco (norm. to Data)", "l");
    leg->AddEntry(grBand, "Total (GENIE #oplus SI #oplus Det)",    "f");
    leg->Draw();

    TLatex tL; tL.SetNDC();
    tL.DrawLatex(0.20, 0.94, "#bf{DUNE:ND-LAr 2x2}");

    c->SaveAs(outPng.c_str());
    delete c;
}

// =============================================================
// Draw combined plot: Data, MC, GENIE band, SI band, Total band -- unchanged
// =============================================================
void DrawCombinedPlot(TH1D* hData, TH1D* hMC,
                      TGraphAsymmErrors* grGENIE,
                      TGraphAsymmErrors* grSI,
                      TGraphAsymmErrors* grTotal,
                      const char* xTitle,
                      const char* title,
                      const std::string& outPng,
                      double ymax = -1)
{
    TCanvas* c = new TCanvas("c", "c", 800, 600);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.15);

    hData->SetLineColor(kBlack);
    hData->SetLineWidth(3);
    hMC->SetLineColor(kRed);
    hMC->SetLineWidth(3);

    hData->SetTitle(title);
    hData->GetXaxis()->SetTitle(xTitle);
    hData->GetYaxis()->SetTitle("Counts");

    double maxY = (ymax > 0) ? ymax : 1.4 * std::max(hData->GetMaximum(), hMC->GetMaximum());
    hData->GetYaxis()->SetRangeUser(0, maxY);

    // Draw order: total (widest) first, then GENIE, then SI, then MC, then data on top
    hData->Draw("HIST E");
    grTotal->Draw("2 SAME");
    grGENIE->Draw("2 SAME");
    grSI->Draw("2 SAME");
    hMC->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.45, 0.68, 0.82, 0.93);
    leg->SetTextFont(133);
    leg->SetTextSize(18);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hData,   "Data Reco",                   "l");
    leg->AddEntry(hMC,     "MC Reco (norm. to Data)",     "l");
    leg->AddEntry(grGENIE, "GENIE Systematics",           "f");
    leg->AddEntry(grSI,   "SI Systematics",             "f");
    leg->AddEntry(grTotal, "Total (GENIE #oplus SI #oplus Det)",    "f");
    leg->Draw();

    TLatex tL; tL.SetNDC();
    tL.DrawLatex(0.20, 0.94, "#bf{DUNE:ND-LAr 2x2}");

    c->SaveAs(outPng.c_str());
    delete c;
}

// =============================================================
// Main
// =============================================================
void plot_dataMC_proton_combined_band()
{
    gROOT->LoadMacro("protoDUNEStyle.C");
    gROOT->SetStyle("protoDUNEStyle");
    gROOT->ForceStyle();
    gStyle->SetTitleX(0.25);
    gSystem->mkdir(OUT_DIR.c_str(), true);

    // ── Load species-split weight maps ─────────────────────────────────────
    std::map<WKey,W2> wmap_proton, wmap_piplus, wmap_piminus;
    LoadWeights(WFILE_PROTON,  wmap_proton);
    LoadWeights(WFILE_PIPLUS,  wmap_piplus);
    LoadWeights(WFILE_PIMINUS, wmap_piminus);

    // ── Open files ────────────────────────────────────────────────────────
    TFile* fCAF  = TFile::Open(CAF_FILE.c_str());
    TFile* fData = TFile::Open(DATA_FILE.c_str());
    TFile* fSyst = TFile::Open(SYST_FILE.c_str());

    if (!fCAF || fCAF->IsZombie())  { std::cerr << "ERROR: CAF file\n";  return; }
    if (!fData || fData->IsZombie()) { std::cerr << "ERROR: Data file\n"; return; }
    if (!fSyst || fSyst->IsZombie()) { std::cerr << "ERROR: Syst file\n"; return; }

    TTree* tReco = (TTree*)fCAF->Get("syst_g4rw_proton_reco");
    if (!tReco) { std::cerr << "ERROR: Cannot find syst_g4rw_proton_reco\n"; return; }

    // ── Data histograms ───────────────────────────────────────────────────
    TH1D* hDataCos = RebinPreserveErrors((TH1D*)fData->Get("recoProtonCosLAll"),
                                          kNBinsC, kXMinC, kXMaxC, "hDataCos");
    TH1D* hDataLen = RebinPreserveErrors((TH1D*)fData->Get("recoProtonLengthAll"),
                                          kNBinsL, kXMinL, kXMaxL, "hDataLen");

    // ── MC histograms from CAF (for scaling/display AND reused below as
    //    the detector-band nominal, since CAF_FILE is now that nominal) ────
    TH1D* hMC_cos_raw = RebinPreserveErrors((TH1D*)fCAF->Get("reco_matched_proton_cosL"),
                                             kNBinsC, kXMinC, kXMaxC, "hMC_cos_raw");
    TH1D* hMC_len_raw = RebinPreserveErrors((TH1D*)fCAF->Get("reco_matched_proton_length"),
                                             kNBinsL, kXMinL, kXMaxL, "hMC_len_raw");

    // ── Nominal + weighted histograms from tree (for SI band), filled
    //    with SPECIES-SPLIT weights via truePDG -- matches
    //    plot_dataMC_recoProton_truthSplit_all.C exactly. ────────────────
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

    for (long i = 0; i < N; i++) {
        tReco->GetEntry(i);
        WKey key = {cafFileIdx, (int)nu_id, G4ID};

        W2 w = {1.0, 1.0};
        if (truePDG == 2212) {
            if (wmap_proton.count(key)) w = wmap_proton[key];
            nProton++;
        } else if (truePDG == 211) {
            if (wmap_piplus.count(key)) w = wmap_piplus[key];
            nPiPlus++;
        } else if (truePDG == -211) {
            if (wmap_piminus.count(key)) w = wmap_piminus[key];
            nPiMinus++;
        } else {
            nOther++;
        }

        h_nom_cos->Fill(cosT);
        h_up_cos->Fill(cosT, w[0]); h_dn_cos->Fill(cosT, w[1]);

        h_nom_len->Fill(length);
        h_up_len->Fill(length, w[0]); h_dn_len->Fill(length, w[1]);
    }

    std::cout << "\n--- Truth-species breakdown ---\n";
    std::cout << "Proton (2212): " << nProton  << " / " << N << "\n";
    std::cout << "pi+     (211): " << nPiPlus  << " / " << N << "\n";
    std::cout << "pi-    (-211): " << nPiMinus << " / " << N << "\n";
    std::cout << "Muon/Other (w=1, no G4RW weight applied): " << nOther << " / " << N << "\n";

    // ── Scale MC to data ──────────────────────────────────────────────────
    double scaleCos = (hMC_cos_raw->Integral("width") > 0.0)
                    ? hDataCos->Integral("width") / hMC_cos_raw->Integral("width") : 1.0;
    TH1D* hMC_cos = (TH1D*)hMC_cos_raw->Clone("hMC_cos");
    hMC_cos->Scale(scaleCos);

    double scaleLen = (hMC_len_raw->Integral("width") > 0.0)
                    ? hDataLen->Integral("width") / hMC_len_raw->Integral("width") : 1.0;
    TH1D* hMC_len = (TH1D*)hMC_len_raw->Clone("hMC_len");
    hMC_len->Scale(scaleLen);

    // ── GENIE bands ───────────────────────────────────────────────────────
    TGraphAsymmErrors* grGENIE_cos_raw =
        (TGraphAsymmErrors*)fSyst->Get("gr_band_sumDials_pCosL");
    if (!grGENIE_cos_raw) { std::cerr << "ERROR: GENIE cosL band not found\n"; return; }
    grGENIE_cos_raw = (TGraphAsymmErrors*)grGENIE_cos_raw->Clone("grGENIE_cos_raw");

    TGraphAsymmErrors* grGENIE_len_raw =
        (TGraphAsymmErrors*)fSyst->Get("gr_band_sumDials_pLength");
    if (!grGENIE_len_raw) { std::cerr << "ERROR: GENIE length band not found\n"; return; }
    grGENIE_len_raw = (TGraphAsymmErrors*)grGENIE_len_raw->Clone("grGENIE_len_raw");

    auto* grGENIE_cos = BuildGENIEBand(grGENIE_cos_raw, hMC_cos);
    auto* grGENIE_len = BuildGENIEBand(grGENIE_len_raw, hMC_len);

    // ── SI bands (species-split) ────────────────────────────────────────
    auto* grSI_cos = BuildSIBand(h_nom_cos, h_up_cos, h_dn_cos, hMC_cos);
    auto* grSI_len = BuildSIBand(h_nom_len, h_up_len, h_dn_len, hMC_len);

    // ── Detector bands -- nominal reused from CAF_FILE (hMC_cos_raw/
    //    hMC_len_raw), no separate MR6.1 file. ──────────────────────────
    std::vector<TH1D*> h_vars_cos, h_vars_len;
    int Nvars = (int)VAR_FILES.size();
    for (int v = 0; v < Nvars; v++) {
        TFile* fVar = TFile::Open(VAR_FILES[v].c_str());
        if (!fVar || fVar->IsZombie()) { std::cerr << "ERROR: var file " << v << ": " << VAR_FILES[v] << "\n"; return; }
        TH1D* hc = RebinPreserveErrors((TH1D*)fVar->Get("reco_matched_proton_cosL"),
                                        kNBinsC, kXMinC, kXMaxC, Form("h_var%d_cos", v));
        TH1D* hl = RebinPreserveErrors((TH1D*)fVar->Get("reco_matched_proton_length"),
                                        kNBinsL, kXMinL, kXMaxL, Form("h_var%d_len", v));
        h_vars_cos.push_back(hc);
        h_vars_len.push_back(hl);
        fVar->Close();
    }
    // Normalize variations to the CAF_FILE nominal (POT equalisation)
    double norm_cos = hMC_cos_raw->Integral();
    double norm_len = hMC_len_raw->Integral();
    for (int v = 0; v < Nvars; v++) {
        if (h_vars_cos[v]->Integral() > 0.0) h_vars_cos[v]->Scale(norm_cos / h_vars_cos[v]->Integral());
        if (h_vars_len[v]->Integral() > 0.0) h_vars_len[v]->Scale(norm_len / h_vars_len[v]->Integral());
    }

    auto* grDet_cos = BuildDetectorBand(hMC_cos_raw, h_vars_cos, hMC_cos);
    auto* grDet_len = BuildDetectorBand(hMC_len_raw, h_vars_len, hMC_len);

    // ── Total bands ───────────────────────────────────────────────────────
    auto* grTotal_cos = BuildTotalBand(grGENIE_cos, grSI_cos, grDet_cos);
    auto* grTotal_len = BuildTotalBand(grGENIE_len, grSI_len, grDet_len);

    // ── Band width tables ────────────────────────────────────────────────
    PrintBandTable("cos(theta)", grGENIE_cos, grSI_cos, grDet_cos, grTotal_cos,
                    kXMinC, (kXMaxC - kXMinC) / kNBinsC);
    PrintBandTable("Track Length [cm]", grGENIE_len, grSI_len, grDet_len, grTotal_len,
                    kXMinL, (kXMaxL - kXMinL) / kNBinsL);

    // ── Draw ──────────────────────────────────────────────────────────────
    // DrawCombinedPlot(hDataCos, hMC_cos, grGENIE_cos, grSI_cos, grTotal_cos,
    //                  "cos(#theta)", "Protons",
    //                  OUT_DIR + "/band_dataMC_cosTheta_combined.png", 50);

    // DrawCombinedPlot(hDataLen, hMC_len, grGENIE_len, grSI_len, grTotal_len,
    //                  "Track Length [cm]", "Protons",
    //                  OUT_DIR + "/band_dataMC_length_combined.png", 70);

    DrawPlot(hDataCos, hMC_cos, grTotal_cos,
             "cos(#theta)", "Protons", "Total Systematics",
             OUT_DIR + "/band_dataMC_cosTheta_total.png", 50);

    DrawPlot(hDataLen, hMC_len, grTotal_len,
             "Track Length [cm]", "Protons", "Total Systematics",
             OUT_DIR + "/band_dataMC_length_total.png", 70);

    std::cout << "\nDone. Plots saved in:\n" << OUT_DIR << "\n";

    fCAF->Close();
    fData->Close();
    fSyst->Close();
}

int main()
{
    plot_dataMC_proton_combined_band();
    return 0;
}