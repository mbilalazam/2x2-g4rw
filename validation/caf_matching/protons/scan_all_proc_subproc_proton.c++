// scan_all_proc_subproc_proton.c++
//
// Proton version of scan_all_proc_subproc.c++ (pion). Same structure:
// for EVERY proton in caf_proton_list.root, opens its patched EDepSim
// trajectory and records the LAST point's (proc, subproc) pair -- the
// exact quantity do_weighting_v2_protons.c++'s last_is_inel check reads.
//
// EXTENDED beyond the pion version: also reads the 5 saved weights
// (already present in caf_proton_list.root, no separate weight-file
// join needed) and computes their spread (max-min). This catches a
// distinct anomaly found in manual proton spot-checks (examples 3, 5, 6)
// that a tag-only scan CANNOT see: tracks whose last point IS correctly
// tagged 4/121, yet whose saved weight is suspiciously flat (near-equal)
// across all five scale factors -- inconsistent with the real,
// non-negligible Sigma_I(p) over these tracks' lengths. Cross-tabulates
// tag vs flat/varies so we can see how common this is per tag category,
// at full-sample scale.
//
// RESUMABLE: output file/trees are opened at the START (UPDATE mode if
// they already exist from a prior partial run), already-processed
// tracks are skipped, and progress is checkpointed to disk periodically
// -- so a crash/restart loses only the last few files' work, not
// everything.
//
// Efficiency note: same as pion version -- groups targets by cafFileIdx
// and opens each EDepSim file exactly once.
//
// Output:
//   scan_all_proc_subproc_proton.root   -- TH1 (tag counts) + trees
//   2d_plots/proc_subproc_hist_proton.png
//
// Compile:
//   g++ -o scan_all_proc_subproc_proton scan_all_proc_subproc_proton.c++ \
//     -I ${EDEPSIM_INC} -L ${EDEPSIM_LIB} -ledepsim_io -I ${ROOT_INC} \
//     -L ${ROOTSYS}/lib \
//     -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree \
//     -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread \
//     -lMultiProc -lROOTDataFrame -pthread -lm -ldl -rdynamic \
//     --std=c++17
//
// Run (this will take a while -- up to ~1000 file opens):
//   ./scan_all_proc_subproc_proton > scan_all_proc_subproc_proton_output.txt 2>&1
// Just re-run the same command to resume after an interruption.

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TLatex.h"
#include "EDepSim/TG4Event.h"
#include "EDepSim/TG4Trajectory.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <tuple>
#include <cmath>
#include <algorithm>

const std::string kCafProtonListFile = "caf_proton_list.root";
const std::string kOutRoot = "scan_all_proc_subproc_proton.root";
const std::string kOutDir  = "2d_plots";
const double kFlatSpread = 0.01;   // same threshold used in find_example_protons_from_caf_list.C
const int kMaxFilesToProcess = -1;  // -1 = no limit (full run). Set small for a quick test first.
const int kCheckpointEvery = 10;    // flush output to disk every N files processed

std::string Pad7(int idx) {
    char buf[8];
    snprintf(buf, sizeof(buf), "%07d", idx);
    return std::string(buf);
}
std::string PatchedNewPath(int idx) {
    return "/global/cfs/cdirs/dune/users/mazam/EdepSim_G4RW_nu/"
           "patched_withRock_genieNu_OLD_pending_delete/"
           "MiniRun5_1E19_RHC.nu." + Pad7(idx) + ".G4RW.EDEPSIM.root";
}

struct TargetInfo {
    int RunId, EventId, TrackId;
    double w034, w067, w133, w166, w200;
    bool found = false;
};

int main() {

    TFile* fIn = TFile::Open(kCafProtonListFile.c_str());
    if (!fIn || fIn->IsZombie()) {
        std::cerr << "ERROR: cannot open " << kCafProtonListFile << "\n";
        return 1;
    }
    TTree* t = (TTree*)fIn->Get("caf_proton_list");
    if (!t) {
        std::cerr << "ERROR: no caf_proton_list tree.\n";
        return 1;
    }

    int cafFileIdx, RunId, EventId, TrackId, pdg;
    double w034, w067, w133, w166, w200;
    t->SetBranchAddress("cafFileIdx", &cafFileIdx);
    t->SetBranchAddress("RunId",      &RunId);
    t->SetBranchAddress("EventId",    &EventId);
    t->SetBranchAddress("TrackId",    &TrackId);
    t->SetBranchAddress("pdg",        &pdg);
    t->SetBranchAddress("w_0p34",     &w034);
    t->SetBranchAddress("w_0p67",     &w067);
    t->SetBranchAddress("w_1p33",     &w133);
    t->SetBranchAddress("w_1p66",     &w166);
    t->SetBranchAddress("w_2p00",     &w200);

    // Group targets by cafFileIdx, carrying the saved weights along.
    std::map<int, std::vector<TargetInfo>> targetsByFile;
    Long64_t N = t->GetEntries();
    long nProtons = 0;
    for (Long64_t i = 0; i < N; ++i) {
        t->GetEntry(i);
        if (pdg != 2212) continue;
        TargetInfo ti;
        ti.RunId = RunId; ti.EventId = EventId; ti.TrackId = TrackId;
        ti.w034 = w034; ti.w067 = w067; ti.w133 = w133; ti.w166 = w166; ti.w200 = w200;
        targetsByFile[cafFileIdx].push_back(ti);
        nProtons++;
    }
    fIn->Close();

    std::cout << "Grouped " << nProtons << " proton tracks across " << targetsByFile.size()
               << " distinct EDepSim files.\n\n";

    // --- Open output file/trees NOW (not at the end) so partial progress
    // survives a crash. If kOutRoot already exists from a prior run,
    // open it in UPDATE mode and read back what's already there.
    gROOT->LoadMacro("protoDUNEStyle.C");
    gROOT->SetStyle("protoDUNEStyle");
    gROOT->ForceStyle();
    gSystem->mkdir(kOutDir.c_str(), true);

    bool resuming = (gSystem->AccessPathName(kOutRoot.c_str()) == 0);   // 0 = file exists
    TFile* fOut = new TFile(kOutRoot.c_str(), resuming ? "UPDATE" : "RECREATE");

    // (proc, subproc) -> total count, and -> count of those ALSO flat
    std::map<std::pair<int,int>, long> tagCounts;
    std::map<std::pair<int,int>, long> tagFlatCounts;
    long nFound = 0, nNotFound = 0, nFileOpenFail = 0, nEmptyPoints = 0;

    // Set of already-done (cafFileIdx,RunId,EventId,TrackId), from a prior run.
    std::set<std::tuple<int,int,int,int>> alreadyDone;

    TTree* tOut = nullptr;
    TTree* tUnmatched = nullptr;
    int out_cafFileIdx, out_RunId, out_EventId, out_TrackId, out_lastProc, out_lastSubproc, out_isFlat;
    double out_spread;
    int u_cafFileIdx, u_RunId, u_EventId, u_TrackId;
    std::string u_reason;

    if (resuming) {
        tOut = (TTree*)fOut->Get("proc_subproc_by_track");
        tUnmatched = (TTree*)fOut->Get("unmatched_tracks");
        if (tOut && tUnmatched) {
            std::cout << "Resuming from existing " << kOutRoot << " ("
                       << tOut->GetEntries() << " tracks + " << tUnmatched->GetEntries()
                       << " unmatched already recorded).\n";

            // Read back tOut's existing entries to rebuild tagCounts/tagFlatCounts
            // and the skip-set, using LOCAL read variables (not the out_* ones
            // we'll reuse for writing new entries below).
            int r_cafFileIdx, r_RunId, r_EventId, r_TrackId, r_lastProc, r_lastSubproc, r_isFlat;
            double r_spread;
            tOut->SetBranchAddress("cafFileIdx", &r_cafFileIdx);
            tOut->SetBranchAddress("RunId",      &r_RunId);
            tOut->SetBranchAddress("EventId",    &r_EventId);
            tOut->SetBranchAddress("TrackId",    &r_TrackId);
            tOut->SetBranchAddress("lastProc",   &r_lastProc);
            tOut->SetBranchAddress("lastSubproc",&r_lastSubproc);
            tOut->SetBranchAddress("spread",     &r_spread);
            tOut->SetBranchAddress("isFlat",     &r_isFlat);
            for (Long64_t i = 0; i < tOut->GetEntries(); ++i) {
                tOut->GetEntry(i);
                tagCounts[{r_lastProc, r_lastSubproc}]++;
                if (r_isFlat) tagFlatCounts[{r_lastProc, r_lastSubproc}]++;
                alreadyDone.insert({r_cafFileIdx, r_RunId, r_EventId, r_TrackId});
                nFound++;
            }

            int ur_cafFileIdx, ur_RunId, ur_EventId, ur_TrackId;
            std::string* ur_reason = nullptr;
            tUnmatched->SetBranchAddress("cafFileIdx", &ur_cafFileIdx);
            tUnmatched->SetBranchAddress("RunId",      &ur_RunId);
            tUnmatched->SetBranchAddress("EventId",    &ur_EventId);
            tUnmatched->SetBranchAddress("TrackId",    &ur_TrackId);
            tUnmatched->SetBranchAddress("reason",     &ur_reason);
            for (Long64_t i = 0; i < tUnmatched->GetEntries(); ++i) {
                tUnmatched->GetEntry(i);
                alreadyDone.insert({ur_cafFileIdx, ur_RunId, ur_EventId, ur_TrackId});
                if (*ur_reason == "never_matched") nNotFound++;
                else if (*ur_reason == "empty_points") nEmptyPoints++;
                else nFileOpenFail++;
            }

            std::cout << "Rebuilt " << alreadyDone.size() << " already-done tracks to skip.\n\n";
        } else {
            std::cout << kOutRoot << " exists but doesn't look like a valid partial run "
                       << "(missing expected trees) -- starting fresh.\n";
            tOut = nullptr; tUnmatched = nullptr;
            resuming = false;
        }
    }

    if (!tOut) {
        tOut = new TTree("proc_subproc_by_track",
            "Last-point proc/subproc + weight spread per proton track, joinable to caf_proton_list via cafFileIdx/RunId/EventId/TrackId");
        tOut->Branch("cafFileIdx", &out_cafFileIdx, "cafFileIdx/I");
        tOut->Branch("RunId",      &out_RunId,      "RunId/I");
        tOut->Branch("EventId",    &out_EventId,    "EventId/I");
        tOut->Branch("TrackId",    &out_TrackId,    "TrackId/I");
        tOut->Branch("lastProc",   &out_lastProc,   "lastProc/I");
        tOut->Branch("lastSubproc",&out_lastSubproc,"lastSubproc/I");
        tOut->Branch("spread",     &out_spread,     "spread/D");
        tOut->Branch("isFlat",     &out_isFlat,     "isFlat/I");
    } else {
        // Re-point branch addresses at the out_* variables for future Fill() calls.
        tOut->SetBranchAddress("cafFileIdx", &out_cafFileIdx);
        tOut->SetBranchAddress("RunId",      &out_RunId);
        tOut->SetBranchAddress("EventId",    &out_EventId);
        tOut->SetBranchAddress("TrackId",    &out_TrackId);
        tOut->SetBranchAddress("lastProc",   &out_lastProc);
        tOut->SetBranchAddress("lastSubproc",&out_lastSubproc);
        tOut->SetBranchAddress("spread",     &out_spread);
        tOut->SetBranchAddress("isFlat",     &out_isFlat);
    }

    if (!tUnmatched) {
        tUnmatched = new TTree("unmatched_tracks",
            "Tracks that failed the EDepSim lookup (never_matched/file_open_failed/no_edepsim_tree) or had zero trajectory points (empty_points)");
        tUnmatched->Branch("cafFileIdx", &u_cafFileIdx, "cafFileIdx/I");
        tUnmatched->Branch("RunId",      &u_RunId,      "RunId/I");
        tUnmatched->Branch("EventId",    &u_EventId,    "EventId/I");
        tUnmatched->Branch("TrackId",    &u_TrackId,    "TrackId/I");
        tUnmatched->Branch("reason",     &u_reason);
    } else {
        tUnmatched->SetBranchAddress("cafFileIdx", &u_cafFileIdx);
        tUnmatched->SetBranchAddress("RunId",      &u_RunId);
        tUnmatched->SetBranchAddress("EventId",    &u_EventId);
        tUnmatched->SetBranchAddress("TrackId",    &u_TrackId);
        tUnmatched->SetBranchAddress("reason",     &u_reason);
    }

    long nFilesProcessed = 0;
    long nFilesSkippedEntirely = 0;
    for (const auto& kv : targetsByFile) {
        if (kMaxFilesToProcess > 0 && nFilesProcessed >= kMaxFilesToProcess) break;
        int fileIdx = kv.first;
        const auto& targets = kv.second;

        // Lookup keyed by (RunId,EventId,TrackId), carrying weights + found flag.
        // Skip anything already recorded from a prior run.
        std::map<std::tuple<int,int,int>, TargetInfo> needed;
        for (const auto& ti : targets) {
            if (alreadyDone.count({fileIdx, ti.RunId, ti.EventId, ti.TrackId})) continue;
            needed[{ti.RunId, ti.EventId, ti.TrackId}] = ti;
        }
        if (needed.empty()) {
            // Every target in this file was already handled in a prior run --
            // don't even bother opening its EDepSim file.
            nFilesSkippedEntirely++;
            nFilesProcessed++;
            continue;
        }

        std::string path = PatchedNewPath(fileIdx);
        TFile* f = TFile::Open(path.c_str());
        if (!f || f->IsZombie()) {
            nFileOpenFail++;
            nNotFound += targets.size();
            for (const auto& ti : targets) {
                u_cafFileIdx = fileIdx; u_RunId = ti.RunId; u_EventId = ti.EventId; u_TrackId = ti.TrackId;
                u_reason = "file_open_failed";
                tUnmatched->Fill();
            }
            continue;
        }
        TTree* tree = (TTree*)f->Get("EDepSimEvents");
        if (!tree) {
            nFileOpenFail++;
            nNotFound += targets.size();
            for (const auto& ti : targets) {
                u_cafFileIdx = fileIdx; u_RunId = ti.RunId; u_EventId = ti.EventId; u_TrackId = ti.TrackId;
                u_reason = "no_edepsim_tree";
                tUnmatched->Fill();
            }
            f->Close();
            continue;
        }
        TG4Event* event = nullptr;
        tree->SetBranchAddress("Event", &event);

        for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
            tree->GetEntry(i);
            for (const auto& traj : event->Trajectories) {
                if (traj.GetPDGCode() != 2212) continue;
                std::tuple<int,int,int> key = {event->RunId, event->EventId, traj.GetTrackId()};
                auto it = needed.find(key);
                if (it == needed.end() || it->second.found) continue;

                if (traj.Points.empty()) {
                    nEmptyPoints++;
                    u_cafFileIdx = fileIdx; u_RunId = event->RunId; u_EventId = event->EventId;
                    u_TrackId = traj.GetTrackId(); u_reason = "empty_points";
                    tUnmatched->Fill();
                    it->second.found = true;
                    continue;
                }
                const auto& lastPt = traj.Points.back();
                int proc = lastPt.GetProcess();
                int sub  = lastPt.GetSubprocess();
                tagCounts[{proc, sub}]++;

                const TargetInfo& ti = it->second;
                double wvals[5] = {ti.w034, ti.w067, ti.w133, ti.w166, ti.w200};
                double wmax = *std::max_element(wvals, wvals+5);
                double wmin = *std::min_element(wvals, wvals+5);
                double spread = wmax - wmin;
                bool isFlat = (spread < kFlatSpread);
                if (isFlat) tagFlatCounts[{proc, sub}]++;

                out_cafFileIdx  = fileIdx;
                out_RunId       = event->RunId;
                out_EventId     = event->EventId;
                out_TrackId     = traj.GetTrackId();
                out_lastProc    = proc;
                out_lastSubproc = sub;
                out_spread      = spread;
                out_isFlat      = isFlat ? 1 : 0;
                tOut->Fill();

                it->second.found = true;
                nFound++;
            }
        }
        f->Close();

        for (const auto& kv2 : needed) {
            if (kv2.second.found) continue;
            nNotFound++;
            u_cafFileIdx = fileIdx; u_RunId = std::get<0>(kv2.first);
            u_EventId = std::get<1>(kv2.first); u_TrackId = std::get<2>(kv2.first);
            u_reason = "never_matched";
            tUnmatched->Fill();
        }

        nFilesProcessed++;
        if (nFilesProcessed % kCheckpointEvery == 0) {
            tOut->AutoSave();
            tUnmatched->AutoSave();
            std::cout << "  ...checkpointed at " << nFilesProcessed << " files processed ("
                       << nFilesSkippedEntirely << " already done, skipped)\n";
        }
        if (nFilesProcessed % 100 == 0) {
            std::cout << "  ...processed " << nFilesProcessed << " / " << targetsByFile.size() << " files\n";
        }
    }

    std::cout << "\n================================================================\n";
    std::cout << " RESULTS: last-point (proc, subproc) counts across " << nFound
               << " found tracks\n";
    std::cout << " (" << nNotFound << " tracks never matched in their file, "
               << nEmptyPoints << " matched but with zero trajectory points, "
               << nFileOpenFail << " file open failures)\n";
    if (resuming) {
        std::cout << " (" << nFilesSkippedEntirely << " files skipped entirely -- already fully done in a prior run)\n";
    }
    std::cout << "================================================================\n";

    std::vector<std::pair<std::pair<int,int>, long>> sortedCounts(tagCounts.begin(), tagCounts.end());
    std::sort(sortedCounts.begin(), sortedCounts.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });

    std::cout << std::setw(10) << "proc" << std::setw(10) << "subproc" << std::setw(12) << "count"
               << std::setw(11) << "fraction" << std::setw(14) << "alsoFlat" << std::setw(12) << "flatFrac\n";
    for (const auto& kv2 : sortedCounts) {
        long flatN = tagFlatCounts.count(kv2.first) ? tagFlatCounts[kv2.first] : 0;
        std::cout << std::setw(10) << kv2.first.first << std::setw(10) << kv2.first.second
                   << std::setw(12) << kv2.second
                   << std::setw(10) << std::fixed << std::setprecision(1)
                   << (100.0 * kv2.second / nFound) << "%"
                   << std::setw(14) << flatN
                   << std::setw(11) << std::setprecision(1)
                   << (100.0 * flatN / kv2.second) << "%\n";
    }
    std::cout << "\n(\"alsoFlat\"/\"flatFrac\": of the tracks with THIS tag, how many have\n"
               << " saved-weight spread < " << kFlatSpread << " across all 5 scale factors --\n"
               << " this is the tagged-but-flat anomaly found in manual spot checks.)\n";

    // --- Histogram with string bin labels "proc/subproc" ---
    TH1D* h = new TH1D("h_proc_subproc", "Last Trajectory Point: proc/subproc;;Number of Protons",
                        (int)sortedCounts.size(), 0, (int)sortedCounts.size());
    for (size_t i = 0; i < sortedCounts.size(); ++i) {
        std::string label = std::to_string(sortedCounts[i].first.first) + "/" +
                             std::to_string(sortedCounts[i].first.second);
        h->GetXaxis()->SetBinLabel(i+1, label.c_str());
        h->SetBinContent(i+1, sortedCounts[i].second);
    }
    h->GetXaxis()->SetLabelSize(0.05);
    h->SetFillColor(kAzure+1);
    h->SetTitle("");
    h->Write();

    TCanvas* c = new TCanvas("c_proc_subproc", "", 900, 650);
    c->SetLeftMargin(0.13);
    c->SetBottomMargin(0.15);
    h->Draw("HIST");
    TLatex tL; tL.SetNDC();
    tL.DrawLatex(0.13, 0.94, "#bf{DUNE:ND-LAr 2x2}");
    TLatex tR; tR.SetNDC(); tR.SetTextSize(0.035);
    tR.DrawLatex(0.50, 0.94, "Truth protons (CAF-matched), last-point tag");
    c->SaveAs((kOutDir + "/proc_subproc_hist_proton.png").c_str());

    fOut->cd();
    tOut->Write("", TObject::kOverwrite);
    tUnmatched->Write("", TObject::kOverwrite);
    fOut->Close();

    std::cout << "\nDone. Wrote " << kOutRoot << " and " << kOutDir << "/proc_subproc_hist_proton.png\n";
    return 0;
}
