///////////////////////////////////////////////////////////////////////////////
// plot_dataMC_proton_GENIE.C
//
// Data vs MC comparison with GENIE (nusystematics sum-of-dials) systematic
// band for reco-matched protons. Combines the two former separate scripts:
//   plot_dataMC_proton_length_GENIE.C
//   plot_dataMC_proton_cosTheta_GENIE.C
// into one, sharing a single MakeGeniePlot() helper parameterized by
// observable, so length and cosTheta can't drift out of sync (they had:
// the old cosTheta script pointed fMC at the stale May 2026 nominal
// (9_May_2026/slides/proton_systematics/output.root) while the length
// script had already been updated to the canonical nominal below).
//
// Part of the total_systematics workstream, alongside:
//   plot_dataMC_proton_detector_band_MR65sys.C  (detector systematics)
//   G4RW covariance/band scripts                (hadronic/G4RW systematics)
// Eventually combined into a total (GENIE + G4RW + detector) uncertainty.
//
// Machinery (unchanged from the original scripts):
//   - GENIE bands (gr_band_sumDials_pLength / gr_band_sumDials_pCosL) are
//     pre-computed elsewhere (nusyst_step2_reco_output.root) as ABSOLUTE
//     sigma vs. their own nominal.
//   - Here, each absolute band is converted to a FRACTIONAL uncertainty
//     (sigma / yNom, using the band's own internal nominal points -- NOT
//     hMC), then re-applied onto hMC after hMC is normalized to data.
//   - i.e. each band's fractional SIZE is independent of which nominal MC
//     file hMC points to; hMC only sets the central curve compared to data
//     and what the band gets re-centered on.
//
// Inputs (shared across both observables):
//   FILE_MC   : nominal MC -- canonical MR6.5 nominal (1E19 POT), same file
//               used as FILE_NOM/FILE_MR65 in
//               plot_dataMC_proton_detector_band_MR65sys.C.
//   FILE_DATA : data (sandbox_v11, same file used everywhere else in this
//               workstream)
//   FILE_SYST : pre-computed GENIE bands (both observables, one file)
//
// Output:
//   proton_length_data_reco_GENIE_band.png
//   proton_cosTheta_data_reco_GENIE_band.png
//   (relative paths -- no OUT_DIR, matching the original scripts; run this
//   from inside total_systematics/ so the PNGs land with the rest.)
//
// Compile and run (inside shifter + DUNE env, same as the detector script):
//   g++ -std=c++17 -O2 -g plot_dataMC_proton_GENIE.C -o plot_dataMC_proton_GENIE `root-config --cflags --libs`
//   ./plot_dataMC_proton_GENIE
///////////////////////////////////////////////////////////////////////////////

#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGraphAsymmErrors.h"
#include <iostream>
#include <cmath>
#include <string>

// ------------------------------------------------------------
// Paths (shared by both observables)
// ------------------------------------------------------------

// Canonical nominal MC -- same file used as FILE_NOM/FILE_MR65 in
// plot_dataMC_proton_detector_band_MR65sys.C (1E19 POT).
const std::string FILE_MC =
    "/global/cfs/cdirs/dune/users/mazam/softwares/patches/2x2_g4rw/validation_proton/codes_for_slides/proton_plots/output.root";

const std::string FILE_DATA =
    "/global/u1/m/mazam/sandbox/v11/sandbox_v11_minervaOff_geomContainRecoOnly_Lgt3cm_300cm_endZneg0to5cmVeto_cosZgtNeg0p9.root";

const std::string FILE_SYST =
    "/global/u1/m/mazam/Geant4ReWeighting/2x2_G4RW/test/on_CAFS/step2_reco/nusyst_step2_reco_output.root";

// ------------------------------------------------------------
// Helper: Rebin preserving errors
// ------------------------------------------------------------
TH1D* RebinPreserveErrors(TH1D* h, int nbins, double xmin, double xmax, const char* name)
{
	TH1D* hR = new TH1D(name, "", nbins, xmin, xmax);
	hR->Sumw2();

	for (int i = 1; i <= h->GetNbinsX(); ++i)
	{
		double x = h->GetBinCenter(i);
		if (x < xmin || x >= xmax) continue;

		int b = hR->FindBin(x);

		double c_old = hR->GetBinContent(b);
		double e_old = hR->GetBinError(b);

		double c = h->GetBinContent(i);
		double e = h->GetBinError(i);

		hR->SetBinContent(b, c_old + c);
		hR->SetBinError(b, std::sqrt(e_old*e_old + e*e));
	}
	hR->SetDirectory(0);
	return hR;
}

// ------------------------------------------------------------
// One observable's worth of work: load, rebin, normalize to data,
// convert the GENIE band to fractional and re-center it, draw, save.
// ------------------------------------------------------------
void MakeGeniePlot(TFile& fMC, TFile& fData, TFile& fSyst,
                    const char* dataHistName,
                    const char* mcHistName,
                    const char* bandGraphName,
                    int nbins, double xmin, double xmax,
                    const char* xTitle,
                    double yMaxUser,
                    const std::string& outPng)
{
	TH1D* hData_raw = (TH1D*)fData.Get(dataHistName);
	TH1D* hMC_raw   = (TH1D*)fMC.Get(mcHistName);

	if (!hData_raw || !hMC_raw) {
		std::cerr << "Missing Data or MC hist for " << outPng << "\n";
		return;
	}

	TH1D* hData = RebinPreserveErrors(hData_raw, nbins, xmin, xmax, "hData_local");
	TH1D* hMC   = RebinPreserveErrors(hMC_raw,   nbins, xmin, xmax, "hMC_local");

	// --------------------------------------------------------
	// Load GENIE sum-of-dials band
	// --------------------------------------------------------
	TGraphAsymmErrors* grBand = (TGraphAsymmErrors*)fSyst.Get(bandGraphName);
	if (!grBand) {
		std::cerr << "GENIE band '" << bandGraphName << "' not found\n";
		return;
	}
	grBand = (TGraphAsymmErrors*)grBand->Clone("grBand_local");

	// --------------------------------------------------------
	// Normalize MC to Data
	// --------------------------------------------------------
	double intData = hData->Integral("width");
	double intMC   = hMC->Integral("width");
	double scale   = (intMC > 0.0) ? intData / intMC : 1.0;

	hMC->Scale(scale);

	// --------------------------------------------------------
	// Convert absolute band -> scaled (fractional, re-centered) band
	// --------------------------------------------------------
	TGraphAsymmErrors* grScaled = new TGraphAsymmErrors();

	for (int i = 0; i < grBand->GetN(); ++i)
	{
		double x, yNom;
		grBand->GetPoint(i, x, yNom);  // nominal absolute (band's own nominal)

		int b = hMC->FindBin(x);
		double yMC = hMC->GetBinContent(b);

		double sigma = grBand->GetErrorYhigh(i);

		// convert absolute -> fractional
		double frac = (yNom > 0.0) ? sigma / yNom : 0.0;

		double exl = grBand->GetErrorXlow(i);
		double exh = grBand->GetErrorXhigh(i);

		grScaled->SetPoint(i, x, yMC);
		grScaled->SetPointError(i, exl, exh, frac*yMC, frac*yMC);
	}

	grScaled->SetFillColorAlpha(kOrange + 1, 0.50);
	grScaled->SetLineColor(kOrange + 1);

	// --------------------------------------------------------
	// Style
	// --------------------------------------------------------
	hData->SetLineColor(kBlack);
	hData->SetLineWidth(3);

	hMC->SetLineColor(kRed);
	hMC->SetLineWidth(3);

	hData->SetTitle("Protons");
	hData->GetXaxis()->SetTitle(xTitle);
	hData->GetYaxis()->SetTitle("Counts");
	hData->GetYaxis()->SetRangeUser(0, yMaxUser);

	// --------------------------------------------------------
	// Draw
	// --------------------------------------------------------
	TCanvas* c = new TCanvas("c", "c", 800, 600);
	gPad->SetLeftMargin(0.18);
	gPad->SetRightMargin(0.15);

	hData->Draw("HIST E");
	grScaled->Draw("2 SAME");
	hMC->Draw("HIST SAME");

	TLegend* leg = new TLegend(0.50, 0.75, 0.80, 0.93);
	leg->SetTextFont(133);
	leg->SetTextSize(18);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->AddEntry(hData, "Data Reco", "l");
	leg->AddEntry(hMC, "MC Reco (norm. to Data)", "l");
	leg->AddEntry(grScaled, "GENIE Systematics", "f");
	leg->Draw();

	TLatex tL;
	tL.SetNDC();
	tL.DrawLatex(0.20, 0.94, "#bf{DUNE:ND-LAr 2x2}");

	c->SaveAs(outPng.c_str());

	delete c;
}

// ------------------------------------------------------------
// MAIN
// ------------------------------------------------------------
void plot_dataMC_proton_GENIE()
{
	gROOT->LoadMacro("protoDUNEStyle.C");
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetTitleX(0.25);

	TFile fMC(FILE_MC.c_str());
	TFile fData(FILE_DATA.c_str());
	TFile fSyst(FILE_SYST.c_str());

	if (fMC.IsZombie())   { std::cerr << "Cannot open MC file: "   << FILE_MC   << "\n"; return; }
	if (fData.IsZombie()) { std::cerr << "Cannot open data file: " << FILE_DATA << "\n"; return; }
	if (fSyst.IsZombie()) { std::cerr << "Cannot open nusyst file: " << FILE_SYST << "\n"; return; }

	// Track length
	MakeGeniePlot(fMC, fData, fSyst,
	              "recoProtonLengthAll",
	              "reco_matched_proton_length",
	              "gr_band_sumDials_pLength",
	              15, 3.0, 83.0,
	              "Track Length [cm]",
	              70,
	              "proton_length_data_reco_GENIE_band.png");

	// cos(theta)
	MakeGeniePlot(fMC, fData, fSyst,
	              "recoProtonCosLAll",
	              "reco_matched_proton_cosL",
	              "gr_band_sumDials_pCosL",
	              15, -1.0, 1.0,
	              "cos(#theta)",
	              50,
	              "proton_cosTheta_data_reco_GENIE_band.png");
}

int main()
{
	plot_dataMC_proton_GENIE();
	return 0;
}