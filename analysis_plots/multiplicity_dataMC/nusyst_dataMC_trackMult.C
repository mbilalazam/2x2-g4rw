///////////////////////////////////////////////////////////////////////////////
// nusyst_dataMC_trackMult.C
//
// Data vs MC comparison for reco track multiplicity with GENIE nusystematics
// (sum-of-dials band only — matches the treatment used for the proton study).
//
// MC is normalized to data (scale = data integral / MC nominal integral).
// The same scale factor is applied to all universe histograms so that
// systematic shape variations are preserved correctly.
//
// Reads:
//   1) Data ROOT file  -> track_mult histogram
//   2) MC ROOT file    -> track_multGood histogram (nominal)
//                      -> syst_reco_events TTree
//                         branches: genieIdx, trackMult, cafFileIdx
//   3) nusyst_output_XXXXXXX.root files (SystWeights tree)
//      Matched via cafFileIdx -> nusyst file index
//      Matched via genieIdx   -> row in SystWeights tree
//
// Produces 1 plot:
//   - sumDials_withFrac : data points + MC nominal + sum-of-dials band,
//                         with a fractional-uncertainty (sigma/N) panel below
//                         -- title "Track Multiplicity", axes matching the
//                         proton-study style.
//
// The totWeight/multisim band (gr_band_tot) is still computed and written to
// the output ROOT file for reference, but is no longer plotted -- this study
// uses the sum-of-dials approach only, matching the proton study.
//
// NOTE: ROOT only — no GENIE libraries required.
//
// Compile:
//   g++ -std=c++17 -O2 -g nusyst_dataMC_trackMult.C -o nusyst_dataMC_trackMult `root-config --cflags --libs`
//
// Run:
//   ./nusyst_dataMC_trackMult
///////////////////////////////////////////////////////////////////////////////

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TObjArray.h"
#include "TROOT.h"
#include "TGraphAsymmErrors.h"

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <cmath>
#include <algorithm>

// =============================================================================
// Helpers
// =============================================================================
int ParseArrayLenFromBranchTitle(TBranch* br)
{
    if (!br) return 0;
    TString title = br->GetTitle();
    int i1 = title.Index("[");
    int i2 = title.Index("]");
    if (i1 < 0 || i2 < 0 || i2 <= i1) return 0;
    return TString(title(i1 + 1, i2 - i1 - 1)).Atoi();
}

TGraphAsymmErrors* BuildBand(const std::vector<TH1D*>& h_univ,
                              const TH1D* h_nom,
                              int nBins, int nUniverses,
                              const char* name,
                              int fillColor, float alpha)
{
    // Build covariance diagonal then make asymmetric error band
    std::vector<double> sigma(nBins + 1, 0.0);
    double norm = (nUniverses > 1) ? 1.0 / nUniverses : 0.0;
    for (int i = 1; i <= nBins; ++i) {
        double c = 0.0;
        for (int u = 0; u < nUniverses; ++u) {
            double di = h_univ[u]->GetBinContent(i) - h_nom->GetBinContent(i);
            c += di * di;
        }
        sigma[i] = std::sqrt(c * norm);
    }

    TGraphAsymmErrors* gr = new TGraphAsymmErrors();
    gr->SetName(name);
    for (int i = 1; i <= nBins; ++i) {
        double x  = h_nom->GetBinCenter(i);
        double y  = h_nom->GetBinContent(i);
        double ex = 0.5 * h_nom->GetBinWidth(i);
        gr->SetPoint(i-1, x, y);
        gr->SetPointError(i-1, ex, ex, sigma[i], sigma[i]);
    }
    gr->SetFillColorAlpha(fillColor, alpha);
    gr->SetLineColor(fillColor);
    return gr;
}

TGraphAsymmErrors* BuildSumDialsBand(
    const std::map<TString, std::vector<TH1D*>>& h_dial,
    const TH1D* h_nom,
    int nBins, int nUniverses,
    const std::vector<TString>& dialNames,
    const char* name,
    int fillColor, float alpha)
{
    double norm = (nUniverses > 1) ? 1.0 / nUniverses : 0.0;
    std::vector<double> sigma2(nBins + 1, 0.0);

    for (const auto& dial : dialNames) {
        const auto& vec = h_dial.at(dial);
        for (int i = 1; i <= nBins; ++i) {
            double c = 0.0;
            for (int u = 0; u < nUniverses; ++u) {
                double di = vec[u]->GetBinContent(i) - h_nom->GetBinContent(i);
                c += di * di;
            }
            sigma2[i] += c * norm;
        }
    }

    TGraphAsymmErrors* gr = new TGraphAsymmErrors();
    gr->SetName(name);
    for (int i = 1; i <= nBins; ++i) {
        double x  = h_nom->GetBinCenter(i);
        double y  = h_nom->GetBinContent(i);
        double ex = 0.5 * h_nom->GetBinWidth(i);
        double ey = std::sqrt(sigma2[i]);
        gr->SetPoint(i-1, x, y);
        gr->SetPointError(i-1, ex, ex, ey, ey);
    }
    gr->SetFillColorAlpha(fillColor, alpha);
    gr->SetLineColor(fillColor);
    return gr;
}

// =============================================================================
// Main
// =============================================================================
void nusyst_dataMC_trackMult()
{
    gROOT->LoadMacro("protoDUNEStyle.C");
    gROOT->SetStyle("protoDUNEStyle");
    gROOT->ForceStyle();
    gStyle->SetTitleX(0.35);
    gStyle->SetOptFit(111);

    // =========================================================================
    // Configuration
    // =========================================================================
    const std::string dataFile   = "/global/cfs/cdirs/dune/users/mazam/talks/workshop_Sep2026/nominal_data.root";
    const std::string mcFile     = "/global/cfs/cdirs/dune/users/mazam/talks/workshop_Sep2026/nominal_mc.root";

    const std::string nusystDir     = "/global/cfs/cdirs/dune/users/zdjurcic/nd_nusyst_MR6p5/test/ndnusyst_root_test0/";
    const std::string outputDir     = "/global/cfs/cdirs/dune/users/mazam/talks/workshop_Sep2026/multiplicity_systematics";
    const std::string outputFile    = outputDir + "/nusyst_dataMC_trackMult_output.root";

    const int minCafFileIdx = 0;
    const int maxCafFileIdx = 999;

    // Binning — match track_multGood binning in MC CAF output
    const int    nBinsM = 15;
    const double xMinM  = 0.0, xMaxM = 15.0;

    // Locked cosmetics -- same values used in detector_dataMC_trackMult.C
    // and g4rw_dataMC_trackMult_patched_splitpad_pm1sigma.C.
    const double Y_MAX = 80.0;
    const double LEG_X1 = 0.55, LEG_Y1 = 0.60, LEG_X2 = 0.88, LEG_Y2 = 0.88;

    // =========================================================================
    // Step 1: Discover nUniverses and dial names
    // =========================================================================
    int nUniverses = 0;
    std::vector<TString> dialNames;
    {
        TFile* f0 = TFile::Open(Form("%s/nusyst_output_%07d.root",
                                     nusystDir.c_str(), minCafFileIdx));
        if (!f0 || f0->IsZombie()) {
            std::cerr << "ERROR: Cannot open first nusyst file.\n"; return;
        }
        TTree* t0 = (TTree*)f0->Get("SystWeights");
        if (!t0) {
            std::cerr << "ERROR: No SystWeights tree.\n"; f0->Close(); return;
        }
        TBranch* brTot = t0->GetBranch("totWeight");
        nUniverses = ParseArrayLenFromBranchTitle(brTot);
        if (nUniverses <= 0) {
            std::cerr << "ERROR: Could not parse nUniverses.\n"; f0->Close(); return;
        }
        TObjArray* branches = t0->GetListOfBranches();
        for (int i = 0; i < branches->GetEntries(); ++i) {
            TBranch* br = (TBranch*)branches->At(i);
            TString name = br->GetName();
            if (name == "totWeight" || name == "EventID" || name == "genieIdx") continue;
            if (ParseArrayLenFromBranchTitle(br) == nUniverses)
                dialNames.push_back(name);
        }
        f0->Close();
    }
    std::cout << "nUniverses = " << nUniverses << "\n";
    std::cout << "nDials     = " << dialNames.size() << "\n";
    for (size_t i = 0; i < dialNames.size(); ++i)
        std::cout << "  dial[" << i << "] = " << dialNames[i] << "\n";


    // =========================================================================
    // Step 2: Load data and MC nominal histograms
    // =========================================================================
    TFile* fData = TFile::Open(dataFile.c_str());
    TFile* fMC   = TFile::Open(mcFile.c_str());

    if (!fData || fData->IsZombie()) {
        std::cerr << "ERROR: Cannot open data file.\n"; return;
    }
    if (!fMC || fMC->IsZombie()) {
        std::cerr << "ERROR: Cannot open MC file.\n"; return;
    }

    TH1D* hData = (TH1D*)fData->Get("track_mult");
    TH1D* hMC   = (TH1D*)fMC->Get("track_multGood");

    if (!hData || !hMC) {
        std::cerr << "ERROR: Could not retrieve histograms.\n"; return;
    }

    hData->SetDirectory(0);
    hMC->SetDirectory(0);
    fData->Close();

    std::cout << "\nData integral: " << hData->Integral() << "\n";
    std::cout << "MC integral:   " << hMC->Integral()   << "\n";

    // =========================================================================
    // Step 3: Compute normalization scale factor
    //   Applied once to MC nominal and all universes
    // =========================================================================
    double scale = 1.0;
    if (hMC->Integral() > 0)
        scale = hData->Integral() / hMC->Integral();
    std::cout << "Normalization scale factor: " << scale << "\n";

    // Scale the nominal
    hMC->Scale(scale);

    // =========================================================================
    // Step 4: Load syst_reco_events from MC file
    // =========================================================================
    TTree* tReco = (TTree*)fMC->Get("syst_reco_events");
    if (!tReco) {
        std::cerr << "ERROR: No syst_reco_events tree.\n";
        fMC->Close(); return;
    }

    int r_genieIdx, r_trackMult, r_cafFileIdx;
    tReco->SetBranchAddress("genieIdx",   &r_genieIdx);
    tReco->SetBranchAddress("trackMult",  &r_trackMult);
    tReco->SetBranchAddress("cafFileIdx", &r_cafFileIdx);

    long nRecoEntries = tReco->GetEntries();
    std::cout << "\nTotal syst_reco_events entries: " << nRecoEntries << "\n";

    std::map<int, std::vector<long>> eventsPerFile;
    for (long i = 0; i < nRecoEntries; ++i) {
        tReco->GetEntry(i);
        if (r_cafFileIdx >= minCafFileIdx && r_cafFileIdx < maxCafFileIdx)
            eventsPerFile[r_cafFileIdx].push_back(i);
    }
    long nSelected = 0;
    for (const auto& kv : eventsPerFile) nSelected += kv.second.size();
    std::cout << "  Selected: " << nSelected << "\n";

    // =========================================================================
    // Step 5: Create universe histograms
    // =========================================================================
    std::vector<TH1D*> h_univ_tot(nUniverses);
    for (int u = 0; u < nUniverses; ++u) {
        h_univ_tot[u] = new TH1D(Form("h_univ_totWeight_trackMult_%d", u),
                                  "", nBinsM, xMinM, xMaxM);
        h_univ_tot[u]->Sumw2();
    }

    std::map<TString, std::vector<TH1D*>> h_univ_dial;
    for (const auto& dial : dialNames) {
        std::vector<TH1D*> vec(nUniverses);
        for (int u = 0; u < nUniverses; ++u) {
            vec[u] = new TH1D(Form("h_univ_%s_trackMult_%d", dial.Data(), u),
                              "", nBinsM, xMinM, xMaxM);
            vec[u]->Sumw2();
        }
        h_univ_dial[dial] = vec;
    }

    // =========================================================================
    // Step 6: Fill universe histograms
    // =========================================================================
    long totalFilled = 0;
    int  nFP         = 0;

    for (const auto& kv : eventsPerFile) {
        int cafIdx          = kv.first;
        const auto& entries = kv.second;

        TString nusystFileName = Form("%s/nusyst_output_%07d.root",
                                      nusystDir.c_str(), cafIdx);
        TFile* fNusyst = TFile::Open(nusystFileName);
        if (!fNusyst || fNusyst->IsZombie()) { continue; }

        TTree* tWeights = (TTree*)fNusyst->Get("SystWeights");
        if (!tWeights) { fNusyst->Close(); continue; }

        Double_t* totWeight = new Double_t[nUniverses];
        int w_genieIdx = -1;
        tWeights->SetBranchAddress("totWeight", totWeight);
        tWeights->SetBranchAddress("genieIdx",  &w_genieIdx);

        std::map<TString, std::vector<double>> dialBuffer;
        for (const auto& dial : dialNames) {
            dialBuffer[dial] = std::vector<double>(nUniverses);
            tWeights->SetBranchAddress(dial, dialBuffer[dial].data());
        }

        std::map<int, std::vector<double>>                    wMap_tot;
        std::map<int, std::map<TString, std::vector<double>>> wMap_dial;
        for (long iW = 0; iW < tWeights->GetEntries(); ++iW) {
            tWeights->GetEntry(iW);
            std::vector<double> wt(nUniverses);
            for (int u = 0; u < nUniverses; ++u) wt[u] = totWeight[u];
            wMap_tot[w_genieIdx] = wt;
            for (const auto& dial : dialNames)
                wMap_dial[w_genieIdx][dial] = dialBuffer[dial];
        }

        for (long entryIdx : entries) {
            tReco->GetEntry(entryIdx);
            double mult = (double)r_trackMult;

            auto itTot = wMap_tot.find(r_genieIdx);
            if (itTot != wMap_tot.end()) {
                for (int u = 0; u < nUniverses; ++u)
                    h_univ_tot[u]->Fill(mult, itTot->second[u]);
            } else {
                for (int u = 0; u < nUniverses; ++u)
                    h_univ_tot[u]->Fill(mult, 1.0);
            }

            auto itDial = wMap_dial.find(r_genieIdx);
            if (itDial != wMap_dial.end()) {
                for (const auto& dial : dialNames)
                    for (int u = 0; u < nUniverses; ++u)
                        h_univ_dial[dial][u]->Fill(mult, itDial->second[dial][u]);
            } else {
                for (const auto& dial : dialNames)
                    for (int u = 0; u < nUniverses; ++u)
                        h_univ_dial[dial][u]->Fill(mult, 1.0);
            }
            totalFilled++;
        }

        delete[] totWeight;
        fNusyst->Close();
        ++nFP;
        if (nFP % 100 == 0)
            std::cout << "  trackMult: " << nFP << " files\n";
    }
    std::cout << "trackMult events filled: " << totalFilled << "\n";

    // =========================================================================
    // Step 7: Apply scale factor to all universe histograms
    // =========================================================================
    for (int u = 0; u < nUniverses; ++u)
        h_univ_tot[u]->Scale(scale);

    for (const auto& dial : dialNames)
        for (int u = 0; u < nUniverses; ++u)
            h_univ_dial[dial][u]->Scale(scale);

    // =========================================================================
    // Step 8: Build uncertainty bands
    //   gr_band_tot  : totWeight/multisim -- kept for reference (written to
    //                  the output ROOT file below), not plotted.
    //   gr_band_sum  : sum-of-dials -- the band actually used/plotted,
    //                  matching the proton study.
    // =========================================================================
    TGraphAsymmErrors* gr_band_tot = BuildBand(
        h_univ_tot, hMC, nBinsM, nUniverses,
        "gr_band_totWeight_trackMult", kBlue+1, 0.3);

    TGraphAsymmErrors* gr_band_sum = BuildSumDialsBand(
        h_univ_dial, hMC, nBinsM, nUniverses, dialNames,
        "gr_band_sumDials_trackMult", kOrange+1, 0.3);

    // =========================================================================
    // Step 9: Save output ROOT file
    // =========================================================================
    gSystem->mkdir(outputDir.c_str(), true);
    TFile* fOut = TFile::Open(outputFile.c_str(), "RECREATE");
    hData->Write("hData_trackMult");
    hMC->Write("hMC_trackMult_normToData");
    gr_band_tot->Write();
    gr_band_sum->Write();
    for (int u = 0; u < nUniverses; ++u) h_univ_tot[u]->Write();
    fOut->Close();
    std::cout << "\nOutput saved to: " << outputFile << "\n";

    // =========================================================================
    // Step 10: Style
    // =========================================================================
    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(1.2);
    hData->SetLineColor(kBlack);
    hData->SetMarkerColor(kBlack);
    hData->GetXaxis()->SetTitle("Number of Reconstructed Tracks");
    hData->GetXaxis()->SetRangeUser(0, 10);
    hData->GetYaxis()->SetTitle("Number of Interactions");
    hData->GetYaxis()->SetTitleOffset(1.1);

    hMC->SetLineColor(kRed);
    hMC->SetLineWidth(3);

    double yMax = std::max(hData->GetMaximum(), hMC->GetMaximum()) * 1.5;

    TLatex tL;
    tL.SetNDC();

    // =========================================================================
    // Step 11: Plot — sum-of-dials band, single canvas (no fractional panel)
    // =========================================================================
    {
        TCanvas* c = new TCanvas("c_dataMC_sumDials", "", 900, 650);
        c->SetLeftMargin(0.14);

        hData->SetTitle("Track Multiplicity");
        hData->GetYaxis()->SetRangeUser(0, yMax);
        hData->Draw("E1");
        gr_band_sum->Draw("2 SAME");
        hMC->Draw("HIST SAME");
        hData->Draw("E1 SAME");

        TLegend* leg = new TLegend(0.55, 0.60, 0.88, 0.88);
        leg->SetTextFont(43); leg->SetTextSize(18);
        leg->SetFillStyle(0); leg->SetBorderSize(0);
        leg->AddEntry(hData,       "Data",                    "lep");
        leg->AddEntry(hMC,         "MC Reco (norm. to Data)", "l");
        leg->AddEntry(gr_band_sum, "GENIE Systematics",       "f");
        leg->Draw();
        tL.DrawLatex(0.20, 0.94, "#bf{DUNE:ND-LAr 2x2}");

        c->SaveAs((outputDir + "/nusyst_dataMC_trackMult_sumDials.png").c_str());
        delete c;
    }

    // =========================================================================
    // Step 12: Plot — sum-of-dials band with fractional uncertainty panel
    // =========================================================================
    {
        // Build fractional uncertainty histogram from gr_band_sum
        TH1D* hFrac = (TH1D*)hMC->Clone("hFrac_sumDials");
        hFrac->Reset();
        for (int i = 0; i < gr_band_sum->GetN(); ++i) {
            double x, y;
            gr_band_sum->GetPoint(i, x, y);
            int b    = hFrac->FindBin(x);
            double f = (y > 0) ? gr_band_sum->GetErrorYhigh(i) / y : 0.0;
            hFrac->SetBinContent(b, f);
        }
        hFrac->SetLineColor(kOrange + 1);
        hFrac->SetLineWidth(2);
        hFrac->SetFillColorAlpha(kOrange + 1, 0.3);
        hFrac->SetTitle("");
        hFrac->GetXaxis()->SetTitle("Number of Reconstructed Tracks");
        hFrac->GetYaxis()->SetTitle("#sigma/N");
        hFrac->GetXaxis()->SetRangeUser(0, 10);
        hFrac->GetYaxis()->SetRangeUser(0, 0.30);
        hFrac->GetXaxis()->SetLabelSize(0.12);
        hFrac->GetXaxis()->SetTitleSize(0.14);
        hFrac->GetXaxis()->SetTitleOffset(1.0);
        hFrac->GetYaxis()->SetLabelSize(0.10);
        hFrac->GetYaxis()->SetTitleSize(0.12);
        hFrac->GetYaxis()->SetTitleOffset(0.55);
        hFrac->GetYaxis()->SetNdivisions(505);

        TCanvas* c = new TCanvas("c_dataMC_sumDials_withFrac", "", 800, 750);

        TPad* p1 = new TPad("p1_sumDials", "", 0, 0.3, 1, 1.0);
        p1->SetLeftMargin(0.14);
        p1->SetRightMargin(0.05);
        p1->SetBottomMargin(0.02);
        p1->SetTopMargin(0.08);
        p1->Draw();
        p1->cd();

        hData->SetTitle("Track Multiplicity");
        hData->GetYaxis()->SetRangeUser(0, yMax);
        hData->GetXaxis()->SetLabelSize(0);
        hData->GetXaxis()->SetTitleSize(0);
        hData->Draw("E1");
        gr_band_sum->Draw("2 SAME");
        hMC->Draw("HIST SAME");
        hData->Draw("E1 SAME");

        TLegend* leg = new TLegend(0.55, 0.60, 0.88, 0.88);
        leg->SetTextFont(43); leg->SetTextSize(18);
        leg->SetFillStyle(0); leg->SetBorderSize(0);
        leg->AddEntry(hData,        "Data",                     "lep");
        leg->AddEntry(hMC,          "MC Reco (norm. to Data)",  "l");
        leg->AddEntry(gr_band_sum,  "GENIE Systematics",        "f");
        leg->Draw();
        tL.DrawLatex(0.15, 0.94, "#bf{DUNE:ND-LAr 2x2}");

        c->cd();
        TPad* p2 = new TPad("p2_sumDials", "", 0, 0.0, 1, 0.3);
        p2->SetLeftMargin(0.14);
        p2->SetRightMargin(0.05);
        p2->SetTopMargin(0.02);
        p2->SetBottomMargin(0.35);
        p2->Draw();
        p2->cd();
        hFrac->Draw("HIST");

        c->SaveAs((outputDir + "/nusyst_dataMC_trackMult_sumDials_splitPad.png").c_str());
        delete hFrac;
        delete c;
    }

    // =========================================================================
    // Done
    // =========================================================================
    fMC->Close();
    std::cout << "\n========================================================\n";
    std::cout << "nusyst_dataMC_trackMult complete.\n";
    std::cout << "  Output ROOT: " << outputFile << "\n";
    std::cout << "  Plots:       " << outputDir  << "/\n";
    std::cout << "========================================================\n";
}

int main()
{
    nusyst_dataMC_trackMult();
    return 0;
}