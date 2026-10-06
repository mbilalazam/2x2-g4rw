///////////////////////////////////////////////////////////////////////////////
// find_example_protons_from_caf_list.C
//
// Reads caf_proton_list.root (produced by dump_caf_proton_list.C) and
// classifies each pdg==2212 entry into one of:
//   UPPER  : self-consistent with w = f * exp((1-f)*S)      (correctly tagged)
//   LOWER  : self-consistent with w = exp((1-f)*S), no f factor (missing tag)
//   FLAT   : all five weights ~equal (spread < 0.01) -- never went through
//            the reweighting formula at all
//   (ambiguous entries -- neither hypothesis fits w_2.00 within 0.05 -- are
//    silently excluded from upper/lower/flat, not force-classified)
//
// Classification method: invert S from w_1.33 under BOTH hypotheses, predict
// w_2.00 from each, and assign to whichever hypothesis's prediction is
// closer to the actual saved w_2.00. This is the SAME self-consistency test
// used for pions, just applied per-proton here -- it does NOT assume protons
// share the pion mechanism, it only sorts by which simple formula the saved
// weights happen to match.
//
// For each category, prints up to kMaxPerBin candidates per length bin
// (2-100, 100-200, 200-400 cm), plus a separate "long tail" listing
// (length > 200 cm) drawn from the already-classified upper/lower sets.
//
// Every candidate printed is, by construction, already in caf_proton_list.root
// and therefore already CAF-traceable -- safe to use directly as targets for
// a proton EDepSim trajectory check (e.g. a validate_proton_example script).
//
// Run with:
//   root -l -b -q find_example_protons_from_caf_list.C
///////////////////////////////////////////////////////////////////////////////

#include "TFile.h"
#include "TTree.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

const std::string kCafProtonListFile = "caf_proton_list.root";

const double kMinLen = 2.0;      // exclude trivial near-zero-length tracks
const double kMaxLen = 400.0;
const double kFlatSpread = 0.01; // max(w)-min(w) below this => "flat" (never reweighted)
const double kMaxFitErr  = 0.05; // worse-than-this fit to w_2.00 under BOTH hypotheses => ambiguous, excluded
const int    kMaxPerBin  = 2;    // candidates printed per length bin per category
const double kLongTailMinLen = 200.0;

struct Candidate {
    int cafFileIdx, RunId, EventId, TrackId;
    double len;
    double w034, w067, w133, w166, w200;
};

enum class Category { UPPER, LOWER, FLAT, AMBIGUOUS };

Category ClassifyProton(double w034, double w067, double w133, double w166, double w200) {
    double vals[5] = {w034, w067, w133, w166, w200};
    double wmax = *std::max_element(vals, vals+5);
    double wmin = *std::min_element(vals, vals+5);
    if (wmax - wmin < kFlatSpread) return Category::FLAT;

    if (w133 <= 0.02 || w133 >= 1.32) return Category::AMBIGUOUS;   // too close to formula's edge to invert safely

    double S_upper = std::log(w133 / 1.33) / (1.0 - 1.33);
    double S_lower = std::log(w133)        / (1.0 - 1.33);

    double predUpper200 = 2.00 * std::exp((1.0 - 2.00) * S_upper);
    double predLower200 = 1.00 * std::exp((1.0 - 2.00) * S_lower);

    double errUpper = std::fabs(predUpper200 - w200);
    double errLower = std::fabs(predLower200 - w200);

    if (std::min(errUpper, errLower) > kMaxFitErr) return Category::AMBIGUOUS;
    return (errUpper < errLower) ? Category::UPPER : Category::LOWER;
}

void PrintCategory(const std::string& name, std::vector<Candidate>& cands, bool binByLength = true) {
    std::cout << "\n================================================================\n";
    std::cout << " " << name << "  --  " << cands.size() << " candidate(s) total\n";
    std::cout << "================================================================\n";

    if (cands.empty()) {
        std::cout << "  (none found)\n";
        return;
    }

    std::cout << std::setw(10) << "cafFileIdx" << std::setw(10) << "RunId" << std::setw(12) << "EventId"
               << std::setw(10) << "TrackId" << std::setw(10) << "len[cm]"
               << std::setw(10) << "w_0.34" << std::setw(10) << "w_0.67" << std::setw(10) << "w_1.33"
               << std::setw(10) << "w_1.66" << std::setw(10) << "w_2.00" << "\n";

    std::vector<std::pair<double,double>> bins;
    if (binByLength) bins = {{kMinLen,100.0}, {100.0,200.0}, {200.0,kMaxLen}};
    else             bins = {{kMinLen,kMaxLen}};

    for (auto& bin : bins) {
        std::vector<Candidate*> inBin;
        for (auto& c : cands) if (c.len >= bin.first && c.len < bin.second) inBin.push_back(&c);
        std::sort(inBin.begin(), inBin.end(), [](Candidate* a, Candidate* b){ return a->len < b->len; });

        size_t nPrint = std::min((size_t)kMaxPerBin, inBin.size());
        for (size_t i = 0; i < nPrint; ++i) {
            Candidate* c = inBin[i];
            std::cout << std::setw(10) << c->cafFileIdx << std::setw(10) << c->RunId << std::setw(12) << c->EventId
                       << std::setw(10) << c->TrackId << std::setw(10) << std::fixed << std::setprecision(1) << c->len
                       << std::setw(10) << std::setprecision(4) << c->w034 << std::setw(10) << c->w067
                       << std::setw(10) << c->w133 << std::setw(10) << c->w166 << std::setw(10) << c->w200 << "\n";
        }
    }
}

void find_example_protons_from_caf_list() {

    TFile* f = TFile::Open(kCafProtonListFile.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: cannot open " << kCafProtonListFile << " -- run dump_caf_proton_list.C first.\n";
        return;
    }
    TTree* t = (TTree*)f->Get("caf_proton_list");
    if (!t) {
        std::cerr << "ERROR: no caf_proton_list tree in " << kCafProtonListFile << "\n";
        return;
    }

    int cafFileIdx, RunId, EventId, TrackId, pdg;
    double length, w034, w067, w133, w166, w200;
    t->SetBranchAddress("cafFileIdx", &cafFileIdx);
    t->SetBranchAddress("RunId",      &RunId);
    t->SetBranchAddress("EventId",    &EventId);
    t->SetBranchAddress("TrackId",    &TrackId);
    t->SetBranchAddress("pdg",        &pdg);
    t->SetBranchAddress("length",     &length);
    t->SetBranchAddress("w_0p34",     &w034);
    t->SetBranchAddress("w_0p67",     &w067);
    t->SetBranchAddress("w_1p33",     &w133);
    t->SetBranchAddress("w_1p66",     &w166);
    t->SetBranchAddress("w_2p00",     &w200);

    std::vector<Candidate> upper, lower, flat;
    long nTotal = 0, nSkippedPdg = 0, nSkippedLen = 0, nAmbiguous = 0;

    Long64_t N = t->GetEntries();
    for (Long64_t i = 0; i < N; ++i) {
        t->GetEntry(i);
        nTotal++;
        if (pdg != 2212) { nSkippedPdg++; continue; }
        if (length < kMinLen || length > kMaxLen) { nSkippedLen++; continue; }

        Category cat = ClassifyProton(w034, w067, w133, w166, w200);
        Candidate c{cafFileIdx, RunId, EventId, TrackId, length, w034, w067, w133, w166, w200};

        if (cat == Category::UPPER) upper.push_back(c);
        else if (cat == Category::LOWER) lower.push_back(c);
        else if (cat == Category::FLAT) flat.push_back(c);
        else nAmbiguous++;
    }

    std::cout << "Scanned " << nTotal << " rows (skipped " << nSkippedPdg << " non-proton, "
               << nSkippedLen << " outside length range).\n";
    std::cout << "Classified: " << upper.size() << " upper, " << lower.size() << " lower, "
               << flat.size() << " flat, " << nAmbiguous << " ambiguous (excluded).\n";

    PrintCategory("UPPER TREND (w = f*exp((1-f)S), correctly tagged)", upper);
    PrintCategory("LOWER TREND (w = exp((1-f)S), missing tag)", lower);
    PrintCategory("FLAT w=1.00 (never went through reweighting formula)", flat, /*binByLength=*/false);

    // Long-tail: pull from already-classified upper+lower, length > kLongTailMinLen
    std::vector<Candidate> longTail;
    for (auto& c : upper) if (c.len > kLongTailMinLen) longTail.push_back(c);
    for (auto& c : lower) if (c.len > kLongTailMinLen) longTail.push_back(c);
    std::sort(longTail.begin(), longTail.end(), [](const Candidate& a, const Candidate& b){ return a.len < b.len; });
    if (longTail.size() > 5) longTail.resize(5);
    PrintCategory("LONG-TRACK TAIL (length > " + std::to_string((int)kLongTailMinLen) + " cm, from upper+lower)",
                  longTail, /*binByLength=*/false);

    std::cout << "\nDone. Every candidate above is already confirmed present in the CAF --\n"
               << "use (cafFileIdx, RunId, EventId, TrackId) directly as targets for a\n"
               << "proton EDepSim trajectory check.\n";

    f->Close();
}

int main() {
    find_example_protons_from_caf_list();
    return 0;
}