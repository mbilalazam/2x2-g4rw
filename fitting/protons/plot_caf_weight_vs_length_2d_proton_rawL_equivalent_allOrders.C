///////////////////////////////////////////////////////////////////////////////
// plot_caf_weight_vs_length_2d_proton_rawL_equivalent_allOrders.C
//
// Shows that a RAW-L polynomial CAN describe the data -- pol3-pol7 in raw L,
// derived by algebraic expansion of the numerically-stable scaled-u fit
// (NOT a direct raw-L fit, which fails due to matrix ill-conditioning, as
// shown previously). See the original version of this script for the full
// derivation; this revision only changes formatting and adds properly-
// propagated errors -- the underlying method is unchanged.
//
// REVISIONS in this version:
//   1. Legend repositioned to match g4rw_dataMC_trackMult.C's convention
//      (LEG_X1=0.55, LEG_Y1, LEG_X2=0.90, LEG_Y2=0.88).
//   2. TWO separate plots again (matching the original per-band slide
//      layout), each with the SAME full background scatter (every pdg==2212
//      row, both bands, exactly like the original pre-band-split plots
//      always showed -- classification only ever gated the FIT, never the
//      histogram fill). Only that plot's own band gets fit curves drawn;
//      the other band's points are visible but left unfit.
//   4. Errors are now printed everywhere they were previously missing:
//        - Raw-L equation coefficients now include their own propagated
//          error, computed correctly via Cov(raw) = A * Cov(scaled) * A^T,
//          where A is the SAME binomial-expansion matrix used to compute
//          the coefficients themselves (not a separate/approximate
//          calculation) -- this is an exact linear-covariance propagation,
//          valid because the raw<->scaled coefficient map is linear.
//        - The w(L) value table now prints "value +- error" per cell,
//          using the existing (already-correct) full-covariance
//          EvalWithError machinery in the scaled basis.
//
// Output: TWO plots (both showing the full scatter; each with only its
// own band's curves), plus per-band printed equations/tables/cross-checks
// (now with errors):
//   plots/caf_weight_vs_length_upper_2d_1p33_rawL_equivalent.png
//   plots/caf_weight_vs_length_lower_2d_1p33_rawL_equivalent.png
//
// Compile and run:
//   g++ -std=c++17 -O2 -g plot_caf_weight_vs_length_2d_proton_rawL_equivalent_allOrders.C -o plot_caf_weight_vs_length_2d_proton_rawL_equivalent_allOrders `root-config --cflags --libs`
//   ./plot_caf_weight_vs_length_2d_proton_rawL_equivalent_allOrders > fit_output_rawL_equivalent.txt 2>&1
///////////////////////////////////////////////////////////////////////////////

#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TMatrixDSym.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"

#include <string>
#include <sstream>
#include <vector>
#include <array>
#include <utility>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cmath>

const std::string kInFile  = "/global/cfs/cdirs/dune/users/mazam/softwares/patches/2x2_g4rw/validation_proton/caf_proton_list_eq3p10_hadInt.root";
const std::string kOutDir  = "/global/cfs/cdirs/dune/users/mazam/softwares/patches/2x2_g4rw/validation_proton/fitting";

const std::string kBranchName = "w_1p33";
const std::string kLabel      = "f=1.33";
const std::string kFileTag    = "1p33";

const int nBinsX = 80;
const double plotXMin = 0.0, plotXMax = 400.0;
const int nBinsY = 100;
const double yMin = 0.0, yMax = 2.0;

const double kFitXMin = 5.0, kFitXMax = 290.0;
const double kScaledX0    = 0.5 * (kFitXMin + kFitXMax);
const double kScaledXScale = 0.5 * (kFitXMax - kFitXMin);
const int kMedianSplitBins = 40;

const std::array<int,5> kFitOrders = {3, 4, 5, 6, 7};

// Same color palettes for the upper (hadronic-interaction) and the lower (not-4/121) bands
const std::array<int,5> kUpperColors = {kRed, kRed+2, kOrange+7, kOrange-3, kYellow+2};
const std::array<int,5> kLowerColors = {kRed, kRed+2, kOrange+7, kOrange-3, kYellow+2};

// Same legend position/size convention as g4rw_dataMC_trackMult.C.
const double LEG_X1 = 0.55, LEG_Y1 = 0.45, LEG_X2 = 0.90, LEG_Y2 = 0.88;

const std::array<double,4> kTableLengths = {50.0, 150.0, 200.0, 250.0};

enum class Band { kUpper, kLower };
std::string BandName(Band band) { return (band == Band::kUpper) ? "upper" : "lower"; }

void DrawDuneLabel() {
    TLatex tL; tL.SetNDC();
    tL.DrawLatex(0.13, 0.94, "#bf{DUNE:ND-LAr 2x2}");
}

long double BinomCoeff(int k, int j) {
    long double result = 1.0L;
    for (int i = 0; i < j; ++i) result = result * (k - i) / (i + 1);
    return result;
}

// A[j][k] such that raw_c[j] = sum_k A[j][k] * scaled_d[k] -- the SAME
// matrix is reused both to compute the raw coefficients themselves AND
// (via A * Cov(scaled) * A^T) to propagate their errors, so the two are
// guaranteed consistent with each other.
std::vector<std::vector<long double>> BuildExpansionMatrix(double x0, double xScale, int order) {
    std::vector<std::vector<long double>> A(order + 1, std::vector<long double>(order + 1, 0.0L));
    for (int k = 0; k <= order; ++k) {
        long double xScaleK = powl((long double)xScale, k);
        for (int j = 0; j <= k; ++j) {
            long double binom = BinomCoeff(k, j);
            long double x0pow = powl((long double)(-x0), k - j);
            A[j][k] = binom * x0pow / xScaleK;
        }
    }
    return A;
}

std::vector<long double> ExpandCoeffs(const std::vector<std::vector<long double>>& A, const std::vector<double>& d) {
    int order = (int)d.size() - 1;
    std::vector<long double> raw(order + 1, 0.0L);
    for (int j = 0; j <= order; ++j)
        for (int k = j; k <= order; ++k)
            raw[j] += A[j][k] * (long double)d[k];
    return raw;
}

// Cov(raw)[j][j] = sum_{k1,k2>=j} A[j][k1]*A[j][k2]*Cov(scaled)[k1][k2] --
// exact linear-covariance propagation (valid since raw = A*scaled is a
// linear map), using the SAME A matrix as ExpandCoeffs above.
std::vector<double> RawCoeffErrors(const std::vector<std::vector<long double>>& A, const TMatrixDSym& covScaled, int order) {
    std::vector<double> errs(order + 1, 0.0);
    for (int j = 0; j <= order; ++j) {
        long double var = 0.0L;
        for (int k1 = j; k1 <= order; ++k1)
            for (int k2 = j; k2 <= order; ++k2)
                var += A[j][k1] * A[j][k2] * (long double)covScaled(k1, k2);
        double v = (double)var;
        errs[j] = (v > 0.0) ? std::sqrt(v) : 0.0;
    }
    return errs;
}

// Full-covariance evaluation in the (numerically safe) scaled basis --
// same exact math used throughout this whole project's EvalWithError.
std::pair<double,double> EvalWithError(const std::vector<double>& d, const TMatrixDSym& cov, int order, double L,
                                        double x0, double xScale) {
    double u = (L - x0) / xScale;
    std::vector<double> basis(order + 1);
    double up = 1.0;
    for (int k = 0; k <= order; ++k) { basis[k] = up; up *= u; }

    double value = 0.0;
    for (int k = 0; k <= order; ++k) value += basis[k] * d[k];

    double variance = 0.0;
    for (int i = 0; i <= order; ++i)
        for (int j = 0; j <= order; ++j)
            variance += basis[i] * basis[j] * cov(i, j);
    double error = (variance > 0.0) ? std::sqrt(variance) : 0.0;

    return {value, error};
}

double EvalRawExpanded(const std::vector<long double>& rawCoeffs, double L) {
    long double val = 0.0L, Lp = 1.0L;
    for (long double c : rawCoeffs) { val += c * Lp; Lp *= (long double)L; }
    return (double)val;
}

// ---------------------------------------------------------------------
// Fits + prints everything for one band, and returns the drawable raw-L
// TF1* curves (5, one per order) plus the raw (length,weight) points --
// the caller combines both bands' returns onto one shared plot.
// ---------------------------------------------------------------------
void FitAndPrintBand(TTree* t, Band band, const std::array<int,5>& colors,
                      std::vector<TF1*>& outFits, std::vector<double>& outLengths, std::vector<double>& outWeights) {

    int pdg, hadInteraction = -1;
    double length = -999;
    double w = -999;
    t->SetBranchAddress("pdg",            &pdg);
    t->SetBranchAddress("length",         &length);
    t->SetBranchAddress(kBranchName.c_str(), &w);
    t->SetBranchAddress("hadInteraction", &hadInteraction);

    int wantHadInt = (band == Band::kUpper) ? 1 : 0;

    Long64_t N = t->GetEntries();
    for (Long64_t i = 0; i < N; ++i) {
        t->GetEntry(i);
        if (pdg != 2212) continue;
        if (hadInteraction != wantHadInt) continue;
        outLengths.push_back(length);
        outWeights.push_back(w);
    }
    std::cout << "\n\n################ " << BandName(band) << " band: " << outLengths.size() << " points ################\n";

    // ---- Fit range + median/MAD-robust profile (unchanged from before) ----
    std::vector<double> rangeL, rangeW;
    for (size_t i = 0; i < outLengths.size(); ++i)
        if (outLengths[i] >= kFitXMin && outLengths[i] <= kFitXMax) { rangeL.push_back(outLengths[i]); rangeW.push_back(outWeights[i]); }

    double fitBinWidth = (kFitXMax - kFitXMin) / kMedianSplitBins;
    std::vector<std::vector<double>> binX(kMedianSplitBins), binY(kMedianSplitBins);
    for (size_t i = 0; i < rangeL.size(); ++i) {
        int b = (int)((rangeL[i] - kFitXMin) / fitBinWidth);
        if (b < 0 || b >= kMedianSplitBins) continue;
        binX[b].push_back(rangeL[i]); binY[b].push_back(rangeW[i]);
    }
    std::vector<double> profX, profY, profXErr, profYErr;
    for (int b = 0; b < kMedianSplitBins; ++b) {
        size_t n = binX[b].size();
        if (n < 2) continue;
        double meanX = 0; for (size_t i = 0; i < n; ++i) meanX += binX[b][i]; meanX /= n;
        std::vector<double> ySorted = binY[b]; std::sort(ySorted.begin(), ySorted.end());
        size_t mid = n / 2;
        double medianY = (n % 2 == 0) ? 0.5*(ySorted[mid-1]+ySorted[mid]) : ySorted[mid];
        std::vector<double> absDev(n);
        for (size_t i = 0; i < n; ++i) absDev[i] = std::fabs(binY[b][i]-medianY);
        std::sort(absDev.begin(), absDev.end());
        double mad = (n % 2 == 0) ? 0.5*(absDev[mid-1]+absDev[mid]) : absDev[mid];
        double robustErrY = 1.2533 * (1.4826*mad) / std::sqrt((double)n);
        double varX = 0; for (size_t i = 0; i < n; ++i) varX += (binX[b][i]-meanX)*(binX[b][i]-meanX);
        varX /= (n-1);
        profX.push_back(meanX); profY.push_back(medianY);
        profXErr.push_back(std::sqrt(varX/n)); profYErr.push_back(robustErrY);
    }
    TGraphErrors* gProfile = new TGraphErrors((int)profX.size(), profX.data(), profY.data(), profXErr.data(), profYErr.data());

    // ---- Fit scaled-u for every order (WITH covariance this time -- "S") ----
    std::vector<std::vector<double>> scaledCoeffsByOrder;
    std::vector<TMatrixDSym> covByOrder;
    std::vector<std::vector<long double>> rawCoeffsByOrder;
    std::vector<std::vector<double>> rawErrsByOrder;

    for (size_t oi = 0; oi < kFitOrders.size(); ++oi) {
        int order = kFitOrders[oi];
        std::string u = "((x-" + std::to_string(kScaledX0) + ")/" + std::to_string(kScaledXScale) + ")";
        std::string formula = "[0]";
        for (int p = 1; p <= order; ++p) formula += "+[" + std::to_string(p) + "]*" + u + "^" + std::to_string(p);
        TF1* fScaled = new TF1(("fScaled" + std::to_string(order) + "_" + BandName(band)).c_str(),
                                formula.c_str(), kFitXMin, kFitXMax);
        TFitResultPtr fr = gProfile->Fit(fScaled, "QNS");   // S: keep covariance matrix

        std::vector<double> d(order + 1);
        for (int p = 0; p <= order; ++p) d[p] = fScaled->GetParameter(p);
        scaledCoeffsByOrder.push_back(d);

        TMatrixDSym cov(order + 1);
        if (fr.Get() && fr->IsValid()) cov = fr->GetCovarianceMatrix();
        else { cov.Zero(); std::cout << "  [WARNING] pol" << order << " (" << BandName(band) << "): invalid covariance -- errors will read 0.\n"; }
        covByOrder.push_back(cov);

        auto A = BuildExpansionMatrix(kScaledX0, kScaledXScale, order);
        rawCoeffsByOrder.push_back(ExpandCoeffs(A, d));
        rawErrsByOrder.push_back(RawCoeffErrors(A, cov, order));

        // Drawable raw-L TF1, built from the expanded coefficients.
        std::string rawFormula = "pol" + std::to_string(order);
        TF1* fRawEq = new TF1(("fRawEq" + std::to_string(order) + "_" + BandName(band)).c_str(),
                               rawFormula.c_str(), kFitXMin, kFitXMax);
        for (int j = 0; j <= order; ++j) fRawEq->SetParameter(j, (double)rawCoeffsByOrder.back()[j]);
        fRawEq->SetLineColor(colors[oi]);
        fRawEq->SetLineWidth(3);
        outFits.push_back(fRawEq);
    }

    // ---- Print equations WITH errors ----
    std::cout << "\n=== Raw-L Fit Equations, with errors (" << BandName(band) << " band) ===\n"
               << "(Coefficients from algebraic expansion of the scaled-u fit; errors from\n"
               << " Cov(raw) = A*Cov(scaled)*A^T, exact linear-covariance propagation.)\n";
    for (size_t oi = 0; oi < kFitOrders.size(); ++oi) {
        int order = kFitOrders[oi];
        auto& raw = rawCoeffsByOrder[oi];
        auto& errs = rawErrsByOrder[oi];
        std::cout << "\npol" << order << ":  w(L) = " << (double)raw[0];
        for (int j = 1; j <= order; ++j) {
            double cj = (double)raw[j];
            std::cout << (cj >= 0 ? " + " : " - ") << std::fabs(cj) << "*L^" << j;
        }
        std::cout << "\n  coefficients:";
        for (int j = 0; j <= order; ++j) std::cout << "  c" << j << "=" << (double)raw[j] << " (+-" << errs[j] << ")";
        std::cout << "\n";
    }

    // ---- Value table WITH errors: rows = L, columns = pol3..pol7 ----
    std::cout << "\n=== w(L) +- error, by order (" << BandName(band) << " band) ===\n";
    std::cout << std::setw(8) << "L (cm)";
    for (int order : kFitOrders) std::cout << std::setw(24) << ("pol" + std::to_string(order));
    std::cout << "\n";
    for (double L : kTableLengths) {
        std::cout << std::setw(8) << L;
        for (size_t oi = 0; oi < kFitOrders.size(); ++oi) {
            auto ve = EvalWithError(scaledCoeffsByOrder[oi], covByOrder[oi], kFitOrders[oi], L, kScaledX0, kScaledXScale);
            std::ostringstream cell; cell << std::setprecision(6) << ve.first << " +- " << std::setprecision(3) << ve.second;
            std::cout << std::setw(24) << cell.str();
        }
        std::cout << "\n";
    }

    // ---- Cross-check: expanded-raw-L eval vs. safe u-based eval (unchanged) ----
    std::cout << "\n=== Cross-check: |raw-L-expanded eval - safe u-based eval| (" << BandName(band) << " band) ===\n";
    std::cout << std::setw(8) << "L (cm)";
    for (int order : kFitOrders) std::cout << std::setw(14) << ("pol" + std::to_string(order));
    std::cout << "\n";
    for (double L : kTableLengths) {
        std::cout << std::setw(8) << L;
        for (size_t oi = 0; oi < kFitOrders.size(); ++oi) {
            auto ve = EvalWithError(scaledCoeffsByOrder[oi], covByOrder[oi], kFitOrders[oi], L, kScaledX0, kScaledXScale);
            double vRaw = EvalRawExpanded(rawCoeffsByOrder[oi], L);
            std::cout << std::setw(14) << std::fabs(vRaw - ve.first);
        }
        std::cout << "\n";
    }
}

// Every pdg==2212 row, regardless of hadInteraction (including -1/unknown)
// -- the SAME full background scatter shown on BOTH plots below, matching
// how the original (pre-band-split) proton plots always looked: the
// classification only ever decided which points feed a FIT, never which
// points get drawn in the histogram itself.
void BuildFullScatter(TTree* t, std::vector<double>& lengths, std::vector<double>& weights) {
    int pdg; double length = -999, w = -999;
    t->SetBranchAddress("pdg",    &pdg);
    t->SetBranchAddress("length", &length);
    t->SetBranchAddress(kBranchName.c_str(), &w);
    Long64_t N = t->GetEntries();
    for (Long64_t i = 0; i < N; ++i) {
        t->GetEntry(i);
        if (pdg != 2212) continue;
        lengths.push_back(length);
        weights.push_back(w);
    }
}

void MakeOneBandPlot(const std::vector<double>& fullL, const std::vector<double>& fullW,
                      const std::vector<TF1*>& fits, Band band) {

    TH2D* hist = new TH2D(("h_rawLequiv_" + BandName(band)).c_str(), ";Track Length [cm];G4RW Weight",
                           nBinsX, plotXMin, plotXMax, nBinsY, yMin, yMax);
    for (size_t i = 0; i < fullL.size(); ++i) hist->Fill(fullL[i], fullW[i]);   // full scatter, both bands, on every plot

    TCanvas* c = new TCanvas(("c_rawLequiv_" + BandName(band)).c_str(), "", 800, 650);
    c->SetLeftMargin(0.14);
    c->SetRightMargin(0.14);
    hist->GetYaxis()->SetTitleOffset(1.3);
    hist->Draw("COLZ");
    for (int k = (int)fits.size() - 1; k >= 0; --k) fits[k]->Draw("L SAME");   // ONLY this band's curves

    TLegend* leg = new TLegend(LEG_X1, LEG_Y1, LEG_X2, LEG_Y2);
    leg->SetTextFont(43); leg->SetTextSize(15);
    leg->SetFillStyle(0); leg->SetBorderSize(0);
    for (size_t oi = 0; oi < kFitOrders.size(); ++oi)
        leg->AddEntry(fits[oi], ("pol" + std::to_string(kFitOrders[oi])).c_str(), "l");
    leg->Draw();

    DrawDuneLabel();
    TLatex tR; tR.SetNDC(); tR.SetTextSize(0.035);
    tR.DrawLatex(0.50, 0.94, ("Truth Protons (CAF), " + kLabel).c_str());

    std::string outPng = kOutDir + "/caf_weight_vs_length_" + BandName(band) + "_2d_" + kFileTag + "_rawL_equivalent.png";
    c->SaveAs(outPng.c_str());
    delete c;
}

void MakePerBandPlots(TTree* t) {

    std::vector<double> fullL, fullW;
    BuildFullScatter(t, fullL, fullW);

    std::vector<TF1*> upperFits, lowerFits;
    std::vector<double> upperL, upperW, lowerL, lowerW;   // unused for plotting now, kept for FitAndPrintBand's own fit-range logic

    FitAndPrintBand(t, Band::kUpper, kUpperColors, upperFits, upperL, upperW);
    FitAndPrintBand(t, Band::kLower, kLowerColors, lowerFits, lowerL, lowerW);

    MakeOneBandPlot(fullL, fullW, upperFits, Band::kUpper);   // upper curves fit + drawn; lower band's points visible, unfit
    MakeOneBandPlot(fullL, fullW, lowerFits, Band::kLower);   // lower curves fit + drawn; upper band's points visible, unfit
}

void plot_caf_weight_vs_length_2d_proton_rawL_equivalent_allOrders() {

    gROOT->LoadMacro("protoDUNEStyle.C");
    gROOT->SetStyle("protoDUNEStyle");
    gROOT->ForceStyle();
    gSystem->mkdir(kOutDir.c_str(), true);

    TFile* f = TFile::Open(kInFile.c_str());
    if (!f || f->IsZombie()) { std::cerr << "ERROR: cannot open " << kInFile << "\n"; return; }
    TTree* t = (TTree*)f->Get("caf_proton_list");
    if (!t) { std::cerr << "ERROR: no caf_proton_list tree in " << kInFile << "\n"; return; }

    MakePerBandPlots(t);

    std::cout << "\nDone. Wrote both raw-L-equivalent plots to " << kOutDir << "/\n";
    f->Close();
}

int main() {
    plot_caf_weight_vs_length_2d_proton_rawL_equivalent_allOrders();
    return 0;
}