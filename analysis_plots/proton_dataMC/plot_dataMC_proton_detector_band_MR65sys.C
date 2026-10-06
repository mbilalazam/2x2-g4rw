///////////////////////////////////////////////////////////////////////////////
// plot_dataMC_proton_detector_band_MR65sys.C
//
// Data vs MC comparison with detector systematic band for reco-matched protons.
// Observables: cosTheta and track length.
//
// Updated to use the new MiniRun6.5_1E18_Systematics production
// (23 detector/physics variations run via purity_dlp_multiplicity.cc,
// see anl_2026_samples/, Sep 2026).
//
// Machinery: SAME as plotSPINE_systematics_MR65_band_v4_dataNorm.C /
//            plot_dataMC_proton_detector_band.C (unchanged)
//   - Band nominal + N variations
//   - Per-bin RMS band: sigma_i = sqrt(1/N * sum_v (var_v - nom)^2)
//   - Fractional uncertainty transferred to the MR6.5 display nominal
//   - MR6.5 display nominal normalized to data
//
// What changed vs. the original script:
//   - Band-RMS nominal (was MR6.1 output_mr6p1.root) -> new MiniRun6.5
//     baseline output (MR6p5_baseline_output.root), produced with the
//     same code/selection/binning as the new variations.
//   - 16 MiniRun6.1 variation files -> 23 new MiniRun6.5 variation files
//     (Birks_A_b, Birks_k_b, Box_Model x5, e_lifetime x4, Long_Diff,
//      Temperature, Tran_Diff x4, W_ion).
//   - FILE_MR65 (display nominal, scaled to data) and FILE_DATA unchanged.
//   - OUT_DIR updated per request.
//
// Formatting: same as GENIE/G4RW scripts
//   - protoDUNEStyle, orange band, fixed ymax
//
// Two output plots:
//   band_dataMC_cosTheta_detector.png
//   band_dataMC_length_detector.png
//
// Compile and run (inside shifter + DUNE env, same as purity_dlp_multiplicity.cc):
//   g++ -std=c++17 -O2 -g plot_dataMC_proton_detector_band_MR65sys.C -o plot_dataMC_proton_detector_band_MR65sys `root-config --cflags --libs`
//   ./plot_dataMC_proton_detector_band_MR65sys
///////////////////////////////////////////////////////////////////////////////

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TGraphAsymmErrors.h"

#include <vector>
#include <string>
#include <cmath>
#include <iostream>

// =============================================================
// Paths
// =============================================================

// New MiniRun6.5_1E18_Systematics production (24 folders: baseline + 23 vars)
const std::string SYS_ROOT =
    "/global/u1/m/mazam/MiniRuns/detector_systematics/anl_2026_samples";

// Band-RMS nominal: new baseline, same code/binning as the new variations.
// const std::string FILE_NOM = SYS_ROOT + "/baseline/MR6p5_baseline_output.root";
const std::string FILE_NOM = "/global/cfs/cdirs/dune/users/mazam/softwares/patches/2x2_g4rw/validation_proton/codes_for_slides/proton_plots/output.root";

// Display nominal (scaled to data) -- as specified, unchanged from before.
const std::string FILE_MR65 =
    "/global/cfs/cdirs/dune/users/mazam/softwares/patches/2x2_g4rw/validation_proton/codes_for_slides/proton_plots/output.root";

// Data -- as specified, unchanged from before.
const std::string FILE_DATA =
    "/global/u1/m/mazam/sandbox/v11/sandbox_v11_minervaOff_geomContainRecoOnly_Lgt3cm_300cm_endZneg0to5cmVeto_cosZgtNeg0p9.root";

const std::string OUT_DIR =
    "/global/cfs/cdirs/dune/users/mazam/softwares/patches/2x2_g4rw/validation_proton/codes_for_slides/total_systematics";

// 23 new variation files (baseline excluded -- it's FILE_NOM above)
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
// Build detector band:
//   1. Compute fractional RMS from the band-nominal + variations
//   2. Apply fractional uncertainty to the MR6.5 (data-scaled) display MC
// Same machinery as plotSPINE_systematics_MR65_band_v4_dataNorm.C /
// plot_dataMC_proton_detector_band.C -- UNCHANGED.
// =============================================================
TGraphAsymmErrors* BuildDetectorBand(TH1D* h_nomBase,
                                      const std::vector<TH1D*>& h_vars,
                                      TH1D* h_mc65_scaled)
{
    int Nvars = h_vars.size();
    TGraphAsymmErrors* gr = new TGraphAsymmErrors();
    int nB = h_mc65_scaled->GetNbinsX();

    for (int i = 1; i <= nB; i++) {
        double x    = h_mc65_scaled->GetBinCenter(i);
        double y65  = h_mc65_scaled->GetBinContent(i);
        double ex   = 0.5 * h_mc65_scaled->GetBinWidth(i);
        double yNom = h_nomBase->GetBinContent(i);

        // RMS of variations around the band-nominal
        double sum = 0.0;
        for (int v = 0; v < Nvars; ++v) {
            double delta = h_vars[v]->GetBinContent(i) - yNom;
            sum += delta * delta;
        }
        double sigmaNom = (Nvars > 0) ? std::sqrt(sum / Nvars) : 0.0;

        // Transfer fractional uncertainty to the MR6.5 display nominal
        double frac    = (yNom > 0.0) ? sigmaNom / yNom : 0.0;
        double sigma65 = frac * y65;

        gr->SetPoint(i-1, x, y65);
        gr->SetPointError(i-1, ex, ex, sigma65, sigma65);
    }
    gr->SetFillColorAlpha(kOrange+1, 0.50);
    gr->SetLineColor(kOrange+1);
    gr->SetMarkerSize(0);
    return gr;
}

// =============================================================
// Draw plot: same formatting as GENIE/G4RW scripts
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

    hData->SetLineColor(kBlack);
    hData->SetLineWidth(3);
    hMC->SetLineColor(kRed);
    hMC->SetLineWidth(3);

    hData->SetTitle(title);
    hData->GetXaxis()->SetTitle(xTitle);
    hData->GetYaxis()->SetTitle("Counts");

    double maxY = (ymax > 0) ? ymax : 1.4 * std::max(hData->GetMaximum(), hMC->GetMaximum());
    hData->GetYaxis()->SetRangeUser(0, maxY);

    hData->Draw("HIST E");
    grBand->Draw("2 SAME");
    hMC->Draw("HIST SAME");

	TLegend* leg = new TLegend(0.50, 0.75, 0.80, 0.93);
    leg->SetTextFont(133);
    leg->SetTextSize(18);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hData,   "Data Reco",               "l");
    leg->AddEntry(hMC,     "MC Reco (norm. to Data)", "l");
    leg->AddEntry(grBand,  bandLabel,                 "f");
    leg->Draw();

    TLatex tL; tL.SetNDC();
    tL.DrawLatex(0.20, 0.94, "#bf{DUNE:ND-LAr 2x2}");

    c->SaveAs(outPng.c_str());
    delete c;
}

// =============================================================
// Main
// =============================================================
void plot_dataMC_proton_detector_band_MR65sys()
{
    gROOT->LoadMacro("protoDUNEStyle.C");
    gROOT->SetStyle("protoDUNEStyle");
    gROOT->ForceStyle();
    gStyle->SetTitleX(0.25);
    gSystem->mkdir(OUT_DIR.c_str(), true);

    // ── Open files ────────────────────────────────────────────────────────
    TFile* fNom  = TFile::Open(FILE_NOM.c_str());
    TFile* fMR65 = TFile::Open(FILE_MR65.c_str());
    TFile* fData = TFile::Open(FILE_DATA.c_str());

    if (!fNom  || fNom->IsZombie())  { std::cerr << "ERROR: band-nominal (new baseline) file: " << FILE_NOM << "\n"; return; }
    if (!fMR65 || fMR65->IsZombie()) { std::cerr << "ERROR: MR6.5 display file: " << FILE_MR65 << "\n";         return; }
    if (!fData || fData->IsZombie()) { std::cerr << "ERROR: Data file: " << FILE_DATA << "\n";           return; }

    // ── Load band-nominal histograms (new baseline) ────────────────────────
    TH1D* h_nomBase_cos_raw = (TH1D*)fNom->Get("reco_matched_proton_cosL");
    TH1D* h_nomBase_len_raw = (TH1D*)fNom->Get("reco_matched_proton_length");
    if (!h_nomBase_cos_raw || !h_nomBase_len_raw) {
        std::cerr << "ERROR: Missing band-nominal histograms in " << FILE_NOM << "\n"; return;
    }
    TH1D* h_nomBase_cos = RebinPreserveErrors(h_nomBase_cos_raw, kNBinsC, kXMinC, kXMaxC, "h_nomBase_cos");
    TH1D* h_nomBase_len = RebinPreserveErrors(h_nomBase_len_raw, kNBinsL, kXMinL, kXMaxL, "h_nomBase_len");

    // ── Load MR6.5 display histograms ───────────────────────────────────────
    TH1D* h_mr65_cos_raw = (TH1D*)fMR65->Get("reco_matched_proton_cosL");
    TH1D* h_mr65_len_raw = (TH1D*)fMR65->Get("reco_matched_proton_length");
    if (!h_mr65_cos_raw || !h_mr65_len_raw) {
        std::cerr << "ERROR: Missing MR6.5 histograms in " << FILE_MR65 << "\n"; return;
    }
    TH1D* h_mr65_cos = RebinPreserveErrors(h_mr65_cos_raw, kNBinsC, kXMinC, kXMaxC, "h_mr65_cos");
    TH1D* h_mr65_len = RebinPreserveErrors(h_mr65_len_raw, kNBinsL, kXMinL, kXMaxL, "h_mr65_len");

    // ── Load data histograms ──────────────────────────────────────────────
    TH1D* hDataCos_raw = (TH1D*)fData->Get("recoProtonCosLAll");
    TH1D* hDataLen_raw = (TH1D*)fData->Get("recoProtonLengthAll");
    if (!hDataCos_raw || !hDataLen_raw) {
        std::cerr << "ERROR: Missing data histograms in " << FILE_DATA << "\n"; return;
    }
    TH1D* hDataCos = RebinPreserveErrors(hDataCos_raw, kNBinsC, kXMinC, kXMaxC, "hDataCos");
    TH1D* hDataLen = RebinPreserveErrors(hDataLen_raw, kNBinsL, kXMinL, kXMaxL, "hDataLen");

    // ── Load 23 variation histograms ────────────────────────────────────────
    std::vector<TH1D*> h_vars_cos, h_vars_len;
    int Nvars = (int)VAR_FILES.size();

    for (int v = 0; v < Nvars; v++) {
        TFile* fVar = TFile::Open(VAR_FILES[v].c_str());
        if (!fVar || fVar->IsZombie()) {
            std::cerr << "ERROR: Cannot open variation file: " << VAR_FILES[v] << "\n";
            return;
        }
        TH1D* hc_raw = (TH1D*)fVar->Get("reco_matched_proton_cosL");
        TH1D* hl_raw = (TH1D*)fVar->Get("reco_matched_proton_length");
        if (!hc_raw || !hl_raw) {
            std::cerr << "ERROR: Missing histograms in: " << VAR_FILES[v] << "\n";
            return;
        }
        // Clone before closing file
        TH1D* hc = RebinPreserveErrors(hc_raw, kNBinsC, kXMinC, kXMaxC,
                                        Form("h_var%d_cos", v));
        TH1D* hl = RebinPreserveErrors(hl_raw, kNBinsL, kXMinL, kXMaxL,
                                        Form("h_var%d_len", v));
        h_vars_cos.push_back(hc);
        h_vars_len.push_back(hl);
        fVar->Close();
        std::cout << "Loaded variation " << v << ": " << VAR_FILES[v] << "\n";
    }

    // ── Scale MR6.5 display MC to data ─────────────────────────────────────
    double scaleCos = (h_mr65_cos->Integral("width") > 0.0)
                    ? hDataCos->Integral("width") / h_mr65_cos->Integral("width") : 1.0;
    TH1D* hMC_cos = (TH1D*)h_mr65_cos->Clone("hMC_cos");
    hMC_cos->Scale(scaleCos);

    double scaleLen = (h_mr65_len->Integral("width") > 0.0)
                    ? hDataLen->Integral("width") / h_mr65_len->Integral("width") : 1.0;
    TH1D* hMC_len = (TH1D*)h_mr65_len->Clone("hMC_len");
    hMC_len->Scale(scaleLen);

    // ── Normalize variations to the band-nominal (POT equalisation) ───────
    double norm_cos = h_nomBase_cos->Integral();
    double norm_len = h_nomBase_len->Integral();
    for (int v = 0; v < Nvars; v++) {
        double vi_cos = h_vars_cos[v]->Integral();
        double vi_len = h_vars_len[v]->Integral();
        if (vi_cos > 0.0) h_vars_cos[v]->Scale(norm_cos / vi_cos);
        if (vi_len > 0.0) h_vars_len[v]->Scale(norm_len / vi_len);
    }

    // ── Build detector bands ──────────────────────────────────────────────
    auto* gr_cos = BuildDetectorBand(h_nomBase_cos, h_vars_cos, hMC_cos);
    auto* gr_len = BuildDetectorBand(h_nomBase_len, h_vars_len, hMC_len);

    // ── Draw plots ────────────────────────────────────────────────────────
    DrawPlot(hDataCos, hMC_cos, gr_cos,
             "cos(#theta)", "Protons",
             "Detector Systematics",
             OUT_DIR + "/band_dataMC_cosTheta_detector.png", 50);

    DrawPlot(hDataLen, hMC_len, gr_len,
             "Track Length [cm]", "Protons",
             "Detector Systematics",
             OUT_DIR + "/band_dataMC_length_detector.png", 70);

    std::cout << "\nDone. Plots saved in:\n" << OUT_DIR << "\n";
    fNom->Close();
    fMR65->Close();
    fData->Close();
}

int main()
{
    plot_dataMC_proton_detector_band_MR65sys();
    return 0;
}