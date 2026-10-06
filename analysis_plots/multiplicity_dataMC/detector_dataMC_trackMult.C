///////////////////////////////////////////////////////////////////////////////
// detector_dataMC_trackMult.C
//
// Data vs MC comparison for reco track multiplicity with a detector
// systematics band (MiniRun6.5 baseline + 23 variations).
//
// Adapted from the old plotSPINE_systematics_MR65_band_v4_dataNorm.C.
// What changed vs. that version:
//   - Two nominals (MR6.1 for the RMS band, MR6.5 for display) collapsed
//     into ONE: nominal_mc.root's track_multGood, used for both roles.
//     This matches the final approach settled on for the proton study --
//     it's the correct choice here specifically because nominal_mc.root
//     and all 23 variation files below were produced by the exact same
//     code (nominal_mc.cc), so there's no cross-production mismatch to
//     account for.
//   - 16 MiniRun6.1 detector/physics variation files -> 23 new MiniRun6.5
//     variation files (Birks_A_b, Birks_k_b, Box_Model x5, e_lifetime x4,
//     Long_Diff, Temperature, Tran_Diff x4, W_ion), from
//     detector_systematic_output_root/. All 23 are included for now;
//     comment out individual entries in VAR_FILES to trim later.
//   - Dropped: the full covariance-matrix/correlation-printing block and
//     the three extra plots (all-variations-lines, MR6.5-only,
//     MR6.5-only-splitpad) from the old script -- none of those fed the
//     final Data+MC+band plots this script produces. Kept: the RMS band
//     logic, the mult=0 bin suppression (presentation choice), and the
//     "total yield systematic" summary print.
//   - DATA_FILE/NOM_FILE updated to nominal_data.root/nominal_mc.root.
//
// Produces 2 plots:
//   detector_dataMC_trackMult.png            (no fractional panel)
//   detector_dataMC_trackMult_splitPad.png   (with sigma/N panel below)
//
// Compile:
//   g++ -std=c++17 -O2 -g detector_dataMC_trackMult.C \
//       -o detector_dataMC_trackMult `root-config --cflags --libs`
//
// Run:
//   ./detector_dataMC_trackMult
///////////////////////////////////////////////////////////////////////////////

#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TGraphAsymmErrors.h"
#include "TPad.h"

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <algorithm>

// =============================================================
// Paths
// =============================================================
const std::string NOM_FILE =
    "/global/cfs/cdirs/dune/users/mazam/talks/workshop_Sep2026/nominal_mc.root";

const std::string DATA_FILE =
    "/global/cfs/cdirs/dune/users/mazam/talks/workshop_Sep2026/nominal_data.root";

const std::string SYS_ROOT =
    "/global/cfs/cdirs/dune/users/mazam/talks/workshop_Sep2026/multiplicity_systematics/detector_systematic_output_root";

const std::string OUT_DIR =
    "/global/cfs/cdirs/dune/users/mazam/talks/workshop_Sep2026/multiplicity_systematics";

// Locked legend position -- use these SAME four numbers in the GENIE and
// SI scripts too, for a consistent legend box across all three plot types.
const double LEG_X1 = 0.55, LEG_Y1 = 0.60, LEG_X2 = 0.88, LEG_Y2 = 0.88;

// Fixed y-axis maximum (was dynamic 1.35*max before) -- locked at 80 so
// all three systematic-source plots share the same y-scale.
const double Y_MAX = 80.0;

// All 23 variation files (baseline excluded -- NOM_FILE serves as nominal
// for both the RMS band and the display curve). Comment out entries here
// to trim later, same pattern as the proton detector script.
const std::vector<std::string> VAR_FILES = {
    SYS_ROOT + "/nominal_mc_MR6p5_Birks_A_b_1_down_output.root",
    SYS_ROOT + "/nominal_mc_MR6p5_Birks_A_b_1_up_output.root",

    SYS_ROOT + "/nominal_mc_MR6p5_Birks_k_b_1_down_output.root",
    SYS_ROOT + "/nominal_mc_MR6p5_Birks_k_b_1_up_output.root",

    SYS_ROOT + "/nominal_mc_MR6p5_Box_Model_alpha_down_output.root",
    SYS_ROOT + "/nominal_mc_MR6p5_Box_Model_alpha_up_output.root",
    SYS_ROOT + "/nominal_mc_MR6p5_Box_Model_beta_down_output.root",
    SYS_ROOT + "/nominal_mc_MR6p5_Box_Model_beta_up_output.root",
    SYS_ROOT + "/nominal_mc_MR6p5_Box_Model_Switch_output.root",

    SYS_ROOT + "/nominal_mc_MR6p5_e_lifetime_1.25_output.root",
    SYS_ROOT + "/nominal_mc_MR6p5_e_lifetime_1.5_output.root",
    SYS_ROOT + "/nominal_mc_MR6p5_e_lifetime_1.75_output.root",
    SYS_ROOT + "/nominal_mc_MR6p5_e_lifetime_2_output.root",

    SYS_ROOT + "/nominal_mc_MR6p5_Long_Diff_1_down_output.root",
    SYS_ROOT + "/nominal_mc_MR6p5_Long_Diff_1_up_output.root",

    SYS_ROOT + "/nominal_mc_MR6p5_Temperature_1_down_output.root",
    SYS_ROOT + "/nominal_mc_MR6p5_Temperature_1_up_output.root",

    SYS_ROOT + "/nominal_mc_MR6p5_Tran_Diff_10p_down_output.root",
    SYS_ROOT + "/nominal_mc_MR6p5_Tran_Diff_10p_up_output.root",
    SYS_ROOT + "/nominal_mc_MR6p5_Tran_Diff_25p_down_output.root",
    SYS_ROOT + "/nominal_mc_MR6p5_Tran_Diff_25p_up_output.root",

    SYS_ROOT + "/nominal_mc_MR6p5_W_ion_1_down_output.root",
    SYS_ROOT + "/nominal_mc_MR6p5_W_ion_1_up_output.root",
};

// =============================================================
// Main
// =============================================================
void detector_dataMC_trackMult()
{
    gROOT->LoadMacro("protoDUNEStyle.C");
    gROOT->SetStyle("protoDUNEStyle");
    gROOT->ForceStyle();
    gStyle->SetOptStat(0);
    gStyle->SetTitleX(0.35);
    gSystem->mkdir(OUT_DIR.c_str(), true);

    // ── Open files ────────────────────────────────────────────────────────
    TFile* fNom  = TFile::Open(NOM_FILE.c_str());
    TFile* fData = TFile::Open(DATA_FILE.c_str());

    if (!fNom  || fNom->IsZombie())  { std::cerr << "ERROR: nominal file: " << NOM_FILE << "\n"; return; }
    if (!fData || fData->IsZombie()) { std::cerr << "ERROR: data file: " << DATA_FILE << "\n"; return; }

    TH1D* h_nom = (TH1D*)fNom->Get("track_multGood");
    if (!h_nom) { std::cerr << "ERROR: track_multGood not found in " << NOM_FILE << "\n"; return; }
    h_nom->SetDirectory(0);

    TH1D* h_data = (TH1D*)fData->Get("track_mult");
    if (!h_data) { std::cerr << "ERROR: track_mult not found in " << DATA_FILE << "\n"; return; }
    h_data->SetDirectory(0);

    // ── Load variation histograms ────────────────────────────────────────
    std::vector<TH1D*> h_vars;
    for (const auto& fname : VAR_FILES) {
        TFile* f = TFile::Open(fname.c_str());
        if (!f || f->IsZombie()) {
            std::cerr << "ERROR: cannot open variation file: " << fname << "\n";
            continue;
        }
        TH1D* h = (TH1D*)f->Get("track_multGood");
        if (h) {
            h->SetDirectory(0);
            h_vars.push_back(h);
        } else {
            std::cerr << "ERROR: track_multGood not found in " << fname << "\n";
        }
        f->Close();
    }
    std::cout << "Loaded " << h_vars.size() << " / " << VAR_FILES.size() << " variation histograms.\n";

    const int nbins = h_nom->GetNbinsX();
    const int Nvars = (int)h_vars.size();

    // ── POT-normalize variations to the nominal's integral ───────────────
    double norm_int = h_nom->Integral();
    for (auto* h : h_vars) {
        double var_int = h->Integral();
        if (var_int > 0.0) h->Scale(norm_int / var_int);
    }

    // ── Suppress multiplicity = 0 bin (presentation choice, as before) ──
    h_nom->SetBinContent(1, 0.0);
    h_nom->SetBinError(1, 0.0);
    h_data->SetBinContent(1, 0.0);
    h_data->SetBinError(1, 0.0);
    for (auto* h : h_vars) {
        h->SetBinContent(1, 0.0);
        h->SetBinError(1, 0.0);
    }

    // ── Per-bin fractional RMS uncertainty relative to the nominal ──────
    std::vector<double> frac_syst(nbins + 1, 0.0);
    for (int i = 2; i <= nbins; ++i) {
        double y = h_nom->GetBinContent(i);
        double sum = 0.0;
        for (int v = 0; v < Nvars; ++v) {
            double delta = h_vars[v]->GetBinContent(i) - y;
            sum += delta * delta;
        }
        double sigma = (Nvars > 0) ? std::sqrt(sum / Nvars) : 0.0;
        frac_syst[i] = (y > 0.0) ? sigma / y : 0.0;

        if (i <= 11) {  // print up to multiplicity = 10
            double rel = 100.0 * frac_syst[i];
            std::cout << "Bin " << i << " (mult = " << (i - 1) << "): "
                      << "Nominal = " << y << ", Syst = ±" << sigma
                      << " (" << rel << "%)\n";
        }
    }

    // ── Total yield systematic (mult 1-10) ───────────────────────────────
    double total_nominal = 0.0;
    for (int i = 2; i <= 11; ++i) total_nominal += h_nom->GetBinContent(i);

    double sum_sq = 0.0;
    for (int v = 0; v < Nvars; ++v) {
        double delta_sum = 0.0;
        for (int i = 2; i <= 11; ++i)
            delta_sum += h_vars[v]->GetBinContent(i) - h_nom->GetBinContent(i);
        sum_sq += delta_sum * delta_sum;
    }
    double sigma_total = (Nvars > 0) ? std::sqrt(sum_sq / Nvars) : 0.0;
    double rel_total   = (total_nominal > 0) ? 100.0 * sigma_total / total_nominal : 0.0;

    std::cout << "\n=== Total yield systematic (mult 1-10) ===\n";
    std::cout << "Nominal total = " << total_nominal << "\n";
    std::cout << "Systematic    = ±" << sigma_total << " (" << rel_total << "%)\n";

    // ── Scale MC (nominal) to data ────────────────────────────────────────
    TH1D* hMC = (TH1D*)h_nom->Clone("hMC_dataNorm_trackMult");
    double dataInt = h_data->Integral();
    double mcInt   = hMC->Integral();
    double scale   = (mcInt > 0.0) ? dataInt / mcInt : 1.0;
    hMC->Scale(scale);

    std::cout << "\nData integral: " << dataInt << "\n";
    std::cout << "MC integral (pre-scale): " << mcInt << "\n";
    std::cout << "Scale factor: " << scale << "\n";

    // ── Build the detector band on the data-normalized MC ────────────────
    TGraphAsymmErrors* gr_band = new TGraphAsymmErrors();
    int ip = 0;
    for (int i = 2; i <= nbins; ++i) {
        double x     = hMC->GetBinCenter(i);
        double y     = hMC->GetBinContent(i);
        double ex    = 0.5 * hMC->GetBinWidth(i);
        double sigma = frac_syst[i] * y;
        gr_band->SetPoint(ip, x, y);
        gr_band->SetPointError(ip, ex, ex, sigma, sigma);
        ip++;
    }
    gr_band->SetFillColorAlpha(kOrange + 1, 0.30);
    gr_band->SetLineColor(kOrange + 1);
    gr_band->SetMarkerSize(0);

    // ── Style ─────────────────────────────────────────────────────────────
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(1.2);
    h_data->SetLineColor(kBlack);
    h_data->SetMarkerColor(kBlack);
    h_data->SetTitle("Track Multiplicity");
    h_data->GetXaxis()->SetTitle("Number of Reconstructed Tracks");
    h_data->GetYaxis()->SetTitle("Number of Interactions");
    h_data->GetXaxis()->SetRangeUser(0, 10);
    h_data->GetYaxis()->SetTitleOffset(1.1);

    hMC->SetLineColor(kRed);
    hMC->SetLineWidth(3);

    double yMax = Y_MAX;

    TLatex tL;
    tL.SetNDC();

    // =========================================================================
    // Plot 1: Data + MC(norm. to data) + Detector band -- no fractional panel
    // =========================================================================
    {
        TCanvas* c = new TCanvas("c_detector_trackMult", "", 800, 600);
        c->SetLeftMargin(0.14);

        h_data->GetYaxis()->SetRangeUser(0, yMax);
        h_data->Draw("E1");
        gr_band->Draw("2 SAME");
        hMC->Draw("HIST SAME");
        h_data->Draw("E1 SAME");

        TLegend* leg = new TLegend(LEG_X1, LEG_Y1, LEG_X2, LEG_Y2);
        leg->SetTextFont(43); leg->SetTextSize(18);
        leg->SetFillStyle(0); leg->SetBorderSize(0);
        leg->AddEntry(h_data, "Data",                    "lep");
        leg->AddEntry(hMC,    "MC Reco (norm. to Data)",  "l");
        leg->AddEntry(gr_band,"Detector Systematics",     "f");
        leg->Draw();

        tL.DrawLatex(0.18, 0.94, "#bf{DUNE:ND-LAr 2x2}");

        c->SaveAs((OUT_DIR + "/detector_dataMC_trackMult.png").c_str());
        delete c;
    }

    // =========================================================================
    // Plot 2: same, with a fractional (sigma/N) panel below
    // =========================================================================
    {
        TH1D* hFrac = (TH1D*)hMC->Clone("hFrac_detector_trackMult");
        hFrac->Reset();
        for (int i = 2; i <= nbins; ++i) hFrac->SetBinContent(i, frac_syst[i]);

        hFrac->SetLineColor(kOrange + 1);
        hFrac->SetLineWidth(2);
        hFrac->SetFillColorAlpha(kOrange + 1, 0.30);
        hFrac->SetTitle("");
        hFrac->GetXaxis()->SetTitle("Number of Reconstructed Tracks");
        hFrac->GetYaxis()->SetTitle("#sigma/N");
        hFrac->GetXaxis()->SetRangeUser(0, 10);
        hFrac->GetYaxis()->SetRangeUser(0, 0.3);
        hFrac->GetXaxis()->SetLabelSize(0.12);
        hFrac->GetXaxis()->SetTitleSize(0.14);
        hFrac->GetXaxis()->SetTitleOffset(1.0);
        hFrac->GetYaxis()->SetLabelSize(0.10);
        hFrac->GetYaxis()->SetTitleSize(0.12);
        hFrac->GetYaxis()->SetTitleOffset(0.55);
        hFrac->GetYaxis()->SetNdivisions(505);

        TCanvas* c = new TCanvas("c_detector_trackMult_split", "", 800, 750);

        TPad* p1 = new TPad("p1_det", "", 0, 0.3, 1, 1.0);
        p1->SetLeftMargin(0.14); p1->SetRightMargin(0.05);
        p1->SetBottomMargin(0.02); p1->SetTopMargin(0.08);
        p1->Draw(); p1->cd();

        h_data->GetXaxis()->SetLabelSize(0);
        h_data->GetXaxis()->SetTitleSize(0);
        h_data->GetYaxis()->SetRangeUser(0, yMax);
        h_data->Draw("E1");
        gr_band->Draw("2 SAME");
        hMC->Draw("HIST SAME");
        h_data->Draw("E1 SAME");

        TLegend* leg = new TLegend(LEG_X1, LEG_Y1, LEG_X2, LEG_Y2);
        leg->SetTextFont(43); leg->SetTextSize(18);
        leg->SetFillStyle(0); leg->SetBorderSize(0);
        leg->AddEntry(h_data, "Data",                    "lep");
        leg->AddEntry(hMC,    "MC Reco (norm. to Data)",  "l");
        leg->AddEntry(gr_band,"Detector Systematics",     "f");
        leg->Draw();

        tL.DrawLatex(0.15, 0.94, "#bf{DUNE:ND-LAr 2x2}");

        c->cd();
        TPad* p2 = new TPad("p2_det", "", 0, 0.0, 1, 0.3);
        p2->SetLeftMargin(0.14); p2->SetRightMargin(0.05);
        p2->SetTopMargin(0.02); p2->SetBottomMargin(0.35);
        p2->Draw(); p2->cd();
        hFrac->Draw("HIST");

        c->SaveAs((OUT_DIR + "/detector_dataMC_trackMult_splitPad.png").c_str());
        delete hFrac;
        delete c;
    }

    std::cout << "\nDone. Plots saved in:\n" << OUT_DIR << "\n";
}

int main()
{
    detector_dataMC_trackMult();
    return 0;
}