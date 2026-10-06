// scan_all_proc_subproc.c++
//
// For EVERY pi+ in caf_pion_list.root, opens its patched (NEW) EDepSim
// trajectory and records the LAST point's (proc, subproc) pair -- the
// exact quantity do_weighting_v2_PiPlus.c++'s last_is_inel check reads.
// Fills a labeled bar-chart histogram of how often each (proc,subproc)
// combination appears as the terminal tag, across the whole sample.
//
// Efficiency note: with ~6400 pi+ tracks spread across up to 1000
// EDepSim files, this groups targets by cafFileIdx FIRST and opens each
// file exactly once, scanning its tree in a single pass for every track
// that belongs to it -- not one file-open per track.
//
// Output:
//   scan_all_proc_subproc_output.txt   -- printed counts (redirect manually)
//   scan_all_proc_subproc.root         -- TH1 with bin labels "proc/subproc"
//   2d_plots/proc_subproc_hist.png        -- same histogram, drawn
//
// Compile:
//   g++ -o scan_all_proc_subproc scan_all_proc_subproc.c++ \
//     -I ${EDEPSIM_INC} -L ${EDEPSIM_LIB} -ledepsim_io -I ${ROOT_INC} \
//     -L ${ROOTSYS}/lib \
//     -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree \
//     -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread \
//     -lMultiProc -lROOTDataFrame -pthread -lm -ldl -rdynamic \
//     --std=c++17
//
// Run (this will take a while -- up to ~1000 file opens):
//   ./scan_all_proc_subproc > scan_all_proc_subproc_output.txt 2>&1

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
#include <vector>
#include <tuple>

const std::string kCafPionListFile = "caf_pion_list.root";
const std::string kOutRoot = "scan_all_proc_subproc.root";
const std::string kOutDir  = "2d_plots";

std::string Pad7(int idx) {
    char buf[8];
    snprintf(buf, sizeof(buf), "%07d", idx);
    return std::string(buf);
}
std::string PatchedNewPath(int idx) {
    return "/global/cfs/cdirs/dune/users/mazam/EdepSim_G4RW_nu/"
           "patched_withRock_genieNu/"
           "MiniRun5_1E19_RHC.nu." + Pad7(idx) + ".G4RW.EDEPSIM.root";
}

int main() {

    TFile* fIn = TFile::Open(kCafPionListFile.c_str());
    if (!fIn || fIn->IsZombie()) {
        std::cerr << "ERROR: cannot open " << kCafPionListFile << "\n";
        return 1;
    }
    TTree* t = (TTree*)fIn->Get("caf_pion_list");
    if (!t) {
        std::cerr << "ERROR: no caf_pion_list tree.\n";
        return 1;
    }

    int cafFileIdx, RunId, EventId, TrackId, pdg;
    t->SetBranchAddress("cafFileIdx", &cafFileIdx);
    t->SetBranchAddress("RunId",      &RunId);
    t->SetBranchAddress("EventId",    &EventId);
    t->SetBranchAddress("TrackId",    &TrackId);
    t->SetBranchAddress("pdg",        &pdg);

    // Group targets by cafFileIdx.
    std::map<int, std::vector<std::tuple<int,int,int>>> targetsByFile;
    Long64_t N = t->GetEntries();
    long nPiPlus = 0;
    for (Long64_t i = 0; i < N; ++i) {
        t->GetEntry(i);
        if (pdg != 211) continue;
        targetsByFile[cafFileIdx].push_back({RunId, EventId, TrackId});
        nPiPlus++;
    }
    fIn->Close();

    std::cout << "Grouped " << nPiPlus << " pi+ tracks across " << targetsByFile.size()
               << " distinct EDepSim files.\n\n";

    // (proc, subproc) -> count
    std::map<std::pair<int,int>, long> tagCounts;
    long nFound = 0, nNotFound = 0, nFileOpenFail = 0, nEmptyPoints = 0;
    
    TTree* tOut = new TTree("proc_subproc_by_track",
        "Last-point proc/subproc per pi+ track, joinable to caf_pion_list via cafFileIdx/RunId/EventId/TrackId");
    int out_cafFileIdx, out_RunId, out_EventId, out_TrackId, out_lastProc, out_lastSubproc;
    tOut->Branch("cafFileIdx", &out_cafFileIdx, "cafFileIdx/I");
    tOut->Branch("RunId",      &out_RunId,      "RunId/I");
    tOut->Branch("EventId",    &out_EventId,    "EventId/I");
    tOut->Branch("TrackId",    &out_TrackId,    "TrackId/I");
    tOut->Branch("lastProc",   &out_lastProc,   "lastProc/I");
    tOut->Branch("lastSubproc",&out_lastSubproc,"lastSubproc/I");

    
    TTree* tUnmatched = new TTree("unmatched_tracks",
        "Tracks that failed the EDepSim lookup (never_matched/file_open_failed/no_edepsim_tree) or had zero trajectory points (empty_points)");
    int u_cafFileIdx, u_RunId, u_EventId, u_TrackId;
    std::string u_reason;
    tUnmatched->Branch("cafFileIdx", &u_cafFileIdx, "cafFileIdx/I");
    tUnmatched->Branch("RunId",      &u_RunId,      "RunId/I");
    tUnmatched->Branch("EventId",    &u_EventId,    "EventId/I");
    tUnmatched->Branch("TrackId",    &u_TrackId,    "TrackId/I");
    tUnmatched->Branch("reason",     &u_reason);

    long nFilesProcessed = 0;
    for (const auto& kv : targetsByFile) {
        int fileIdx = kv.first;
        const auto& targets = kv.second;

        // Build a lookup set for this file's needed (RunId,EventId,TrackId).
        std::map<std::tuple<int,int,int>, bool> needed;
        for (const auto& tup : targets) needed[tup] = false;

        std::string path = PatchedNewPath(fileIdx);
        TFile* f = TFile::Open(path.c_str());
        if (!f || f->IsZombie()) {
            nFileOpenFail++;
            nNotFound += targets.size();
            for (const auto& tup : targets) {
                u_cafFileIdx = fileIdx; u_RunId = std::get<0>(tup);
                u_EventId = std::get<1>(tup); u_TrackId = std::get<2>(tup);
                u_reason = "file_open_failed";
                tUnmatched->Fill();
            }
            continue;
        }
        TTree* tree = (TTree*)f->Get("EDepSimEvents");
        if (!tree) {
            nFileOpenFail++;
            nNotFound += targets.size();
            for (const auto& tup : targets) {
                u_cafFileIdx = fileIdx; u_RunId = std::get<0>(tup);
                u_EventId = std::get<1>(tup); u_TrackId = std::get<2>(tup);
                u_reason = "no_edepsim_tree";
                tUnmatched->Fill();
            }
            f->Close();
            continue;
        }
        TG4Event* event = nullptr;
        tree->SetBranchAddress("Event", &event);

        // Single pass through this file's tree; check every event against
        // every (RunId,EventId) still needed from this file.
        for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
            tree->GetEntry(i);
            for (const auto& traj : event->Trajectories) {
                if (traj.GetPDGCode() != 211) continue;
                std::tuple<int,int,int> key = {event->RunId, event->EventId, traj.GetTrackId()};
                auto it = needed.find(key);
                if (it == needed.end() || it->second) continue;   // not needed, or already found

                if (traj.Points.empty()) {
                    nEmptyPoints++;
                    u_cafFileIdx = fileIdx; u_RunId = event->RunId; u_EventId = event->EventId;
                    u_TrackId = traj.GetTrackId(); u_reason = "empty_points";
                    tUnmatched->Fill();
                    it->second = true;
                    continue;
                }
                const auto& lastPt = traj.Points.back();
                int proc = lastPt.GetProcess();
                int sub  = lastPt.GetSubprocess();
                tagCounts[{proc, sub}]++;

                out_cafFileIdx  = fileIdx;
                out_RunId       = event->RunId;
                out_EventId     = event->EventId;
                out_TrackId     = traj.GetTrackId();
                out_lastProc    = proc;
                out_lastSubproc = sub;
                tOut->Fill();

                it->second = true;
                nFound++;
            }
        }
        f->Close();

        for (const auto& kv2 : needed) {
            if (kv2.second) continue;
            nNotFound++;
            u_cafFileIdx = fileIdx; u_RunId = std::get<0>(kv2.first);
            u_EventId = std::get<1>(kv2.first); u_TrackId = std::get<2>(kv2.first);
            u_reason = "never_matched";
            tUnmatched->Fill();
        }

        nFilesProcessed++;
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
    std::cout << "================================================================\n";

    // Sort by count descending for the printout.
    std::vector<std::pair<std::pair<int,int>, long>> sortedCounts(tagCounts.begin(), tagCounts.end());
    std::sort(sortedCounts.begin(), sortedCounts.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });

    std::cout << std::setw(10) << "proc" << std::setw(10) << "subproc" << std::setw(12) << "count"
               << std::setw(12) << "fraction\n";
    for (const auto& kv2 : sortedCounts) {
        std::cout << std::setw(10) << kv2.first.first << std::setw(10) << kv2.first.second
                   << std::setw(12) << kv2.second
                   << std::setw(11) << std::fixed << std::setprecision(1)
                   << (100.0 * kv2.second / nFound) << "%\n";
    }

    // --- Histogram with string bin labels "proc/subproc" ---
    gROOT->LoadMacro("protoDUNEStyle.C");
    gROOT->SetStyle("protoDUNEStyle");
    gROOT->ForceStyle();
    gSystem->mkdir(kOutDir.c_str(), true);

    TFile* fOut = new TFile(kOutRoot.c_str(), "RECREATE");
    TH1D* h = new TH1D("h_proc_subproc", "Last Trajectory Point: proc/subproc;;Number of Pions",
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
    tR.DrawLatex(0.55, 0.94, "Truth #pi^{+} (CAF-matched), last-point tag");
    c->SaveAs((kOutDir + "/proc_subproc_hist.png").c_str());

    fOut->cd();
    tOut->Write();
    tUnmatched->Write();
    fOut->Close();

    std::cout << "\nDone. Wrote " << kOutRoot << " and " << kOutDir << "/proc_subproc_hist.png\n";
    return 0;
}
