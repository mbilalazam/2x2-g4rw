///////////////////////////////////////////////////////////////////////////////
// find_example_pions_from_caf_list.C
//
// Same purpose as find_example_pions_by_weight.c++ (find candidate pi+
// tracks in three w_1.33 bands: low tail, mid-range ~0.5, near w=1), but
// scanning caf_pion_list.root (produced by dump_caf_pion_list.C) instead
// of the raw g4rwEvt weight tree.
//
// WHY: the raw weight file has no guarantee a given entry is traceable
// back through the CAF -- confirmed directly (file_idx=488/RunId=589/
// EventId=32718/TrackId=3 exists in the raw weight file but NOT in the
// CAF's syst_g4rw_pion_truth tree). Every row in caf_pion_list.root is,
// by construction, already a CAF-matched pion, so any candidate picked
// from here is guaranteed traceable to EDepSim and the weight files via
// the direct (cafFileIdx, RunId, EventId, TrackId) keys -- no offset,
// no guessing.
//
// Same filters as before, reapplied here:
//   - kMinLen / kMaxLen: excludes trivial near-zero-length tracks
//   - kMinSpread: excludes tracks where all five weights are ~equal
//     (never actually went through the reweighting formula)
//
// Run with:
//   root -l -b -q find_example_pions_from_caf_list.C

#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

const std::string kCafPionListFile = "caf_pion_list.root";   // output of dump_caf_pion_list.C

const int kMaxPerBand = 15;

const double kMinLen = 2.0;     // cm -- excludes trivial/degenerate tracks
const double kMaxLen = 300.0;   // cm -- generous sanity bound

const double kMinSpread = 0.02; // excludes flat weight=1-for-all-f tracks (never actually reweighted)

struct Band { std::string name; double lo; double hi; double center; };
const std::vector<Band> kBands = {
    {"BAND A (w_1.33 ~ 0.05-0.30, low tail)",   0.05, 0.30, 0.175},
    {"BAND B (w_1.33 ~ 0.45-0.55, mid-range)",   0.45, 0.55, 0.50},
    {"BAND C (w_1.33 ~ 0.95-1.05, near w=1)",    0.95, 1.05, 1.00}
};

struct Candidate {
    int cafFileIdx, RunId, EventId, TrackId, pdg;
    double len, cosTheta;
    double w034, w067, w133, w166, w200;
    double distToCenter, spread;
};

void find_example_pions_from_caf_list() {

    TFile* f = TFile::Open(kCafPionListFile.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: could not open " << kCafPionListFile
                   << " -- run dump_caf_pion_list.C first.\n";
        return;
    }
    TTree* t = (TTree*)f->Get("caf_pion_list");
    if (!t) {
        std::cerr << "ERROR: no caf_pion_list tree in " << kCafPionListFile << "\n";
        return;
    }

    int cafFileIdx, RunId, EventId, TrackId, pdg;
    double cosTheta, length, w034, w067, w133, w166, w200;
    t->SetBranchAddress("cafFileIdx", &cafFileIdx);
    t->SetBranchAddress("RunId",      &RunId);
    t->SetBranchAddress("EventId",    &EventId);
    t->SetBranchAddress("TrackId",    &TrackId);
    t->SetBranchAddress("pdg",        &pdg);
    t->SetBranchAddress("cosTheta",   &cosTheta);
    t->SetBranchAddress("length",     &length);
    t->SetBranchAddress("w_0p34",     &w034);
    t->SetBranchAddress("w_0p67",     &w067);
    t->SetBranchAddress("w_1p33",     &w133);
    t->SetBranchAddress("w_1p66",     &w166);
    t->SetBranchAddress("w_2p00",     &w200);

    std::vector<std::vector<Candidate>> bandCandidates(kBands.size());

    Long64_t nEntries = t->GetEntries();
    std::cout << "Scanning " << nEntries << " CAF-matched entries in " << kCafPionListFile << " ...\n";

    int nSkippedPdg = 0, nSkippedLen = 0, nSkippedSpread = 0;

    for (Long64_t i = 0; i < nEntries; ++i) {
        t->GetEntry(i);
        if (pdg != 211) { nSkippedPdg++; continue; }
        if (length < kMinLen || length > kMaxLen) { nSkippedLen++; continue; }

        double wmax = std::max({w034, w067, w133, w166, w200});
        double wmin = std::min({w034, w067, w133, w166, w200});
        double spread = wmax - wmin;
        if (spread < kMinSpread) { nSkippedSpread++; continue; }

        for (size_t b = 0; b < kBands.size(); ++b) {
            const auto& band = kBands[b];
            if (w133 < band.lo || w133 > band.hi) continue;

            Candidate c;
            c.cafFileIdx = cafFileIdx; c.RunId = RunId; c.EventId = EventId; c.TrackId = TrackId;
            c.pdg = pdg; c.len = length; c.cosTheta = cosTheta;
            c.w034 = w034; c.w067 = w067; c.w133 = w133; c.w166 = w166; c.w200 = w200;
            c.distToCenter = std::fabs(w133 - band.center);
            c.spread = spread;

            bandCandidates[b].push_back(c);
        }
    }

    std::cout << "Skipped: " << nSkippedPdg << " (pdg!=211), " << nSkippedLen
               << " (length out of [" << kMinLen << "," << kMaxLen << "]), "
               << nSkippedSpread << " (spread < " << kMinSpread << ")\n";

    for (size_t b = 0; b < kBands.size(); ++b) {
        auto& cands = bandCandidates[b];
        std::cout << "\n================================================================\n";
        std::cout << " " << kBands[b].name << "  --  " << cands.size() << " candidate(s) found total\n";
        std::cout << "================================================================\n";

        if (cands.empty()) {
            std::cout << "  (none found -- consider widening the band)\n";
            continue;
        }

        std::sort(cands.begin(), cands.end(),
                  [](const Candidate& a, const Candidate& c) { return a.distToCenter < c.distToCenter; });

        size_t nPrint = std::min((size_t)kMaxPerBand, cands.size());

        std::cout << std::setw(10) << "cafFileIdx" << std::setw(10) << "RunId" << std::setw(12) << "EventId"
                   << std::setw(10) << "TrackId" << std::setw(10) << "len[cm]"
                   << std::setw(10) << "w_0.34" << std::setw(10) << "w_0.67" << std::setw(10) << "w_1.33"
                   << std::setw(10) << "w_1.66" << std::setw(10) << "w_2.00" << "\n";

        for (size_t k = 0; k < nPrint; ++k) {
            const auto& c = cands[k];
            std::cout << std::setw(10) << c.cafFileIdx << std::setw(10) << c.RunId << std::setw(12) << c.EventId
                       << std::setw(10) << c.TrackId << std::setw(10) << std::fixed << std::setprecision(3) << c.len
                       << std::setw(10) << c.w034 << std::setw(10) << c.w067 << std::setw(10) << c.w133
                       << std::setw(10) << c.w166 << std::setw(10) << c.w200 << "\n";
        }
    }

    std::cout << "\nDone. Every candidate above is already confirmed present in the CAF\n"
               << "(that's what caf_pion_list.root is built from) -- use (cafFileIdx, RunId,\n"
               << "EventId, TrackId) directly as the target constants in a copy of\n"
               << "validate_pion_example_2.c++ / caf_edepsim_weights_matcher.C.\n";

    f->Close();
}

int main() {
    find_example_pions_from_caf_list();
    return 0;
}