// check_patched_traj_nkaon_all1000.C
//
// Loops over all 1000 files in the extended-patch (pi/p/n/K) EDepSim sample
// and accumulates trajectory point-count statistics per particle species.
// Prints progress every N files since this processes ~880k events total
// and can take a while.
//
// Run with:
//   root -l -b -q check_patched_traj_nkaon_all1000.C

#ifdef __CLING__
R__LOAD_LIBRARY(libedepsim_io)
#endif

#include "TFile.h"
#include "TTree.h"
#include "TG4Event.h"
#include <iostream>
#include <algorithm>
#include <map>
#include <string>

// ── Configure here ────────────────────────────────────────────────────────
const std::string DIR    = "/global/cfs/cdirs/dune/users/mazam/EdepSim_G4RW_nu/patched_withRock_genieNu_v2/";
const std::string PREFIX = "MiniRun5_1E19_RHC.nu.";
const std::string SUFFIX = ".G4RW.EDEPSIM.root";
const int N_FILES = 1000;
const int PROGRESS_EVERY = 50;   // print a progress line every N files
// ─────────────────────────────────────────────────────────────────────────

struct ChannelStats {
    long long n     = 0;
    double    sum   = 0;
    int       min   = 999999;
    int       max   = 0;
};

void check_patched_traj_nkaon_all1000()
{
    std::map<std::string, ChannelStats> stats;
    stats["PiPlus (211)"]  = ChannelStats();
    stats["PiMinus(-211)"] = ChannelStats();
    stats["proton (2212)"] = ChannelStats();
    stats["neutron(2112)"] = ChannelStats();
    stats["K+     (321)"]  = ChannelStats();
    stats["K-     (-321)"] = ChannelStats();

    int filesOk = 0, filesMissing = 0;
    long long totalEvents = 0;

    for (int g = 0; g < N_FILES; g++) {
        char gStr[8];
        snprintf(gStr, sizeof(gStr), "%07d", g);
        std::string fname = DIR + PREFIX + gStr + SUFFIX;

        TFile f(fname.c_str());
        if (f.IsZombie()) {
            std::cerr << "WARNING: Cannot open " << fname << " -- skipping." << std::endl;
            filesMissing++;
            continue;
        }

        TTree* tree = (TTree*)f.Get("EDepSimEvents");
        if (!tree) {
            std::cerr << "WARNING: No EDepSimEvents in " << fname << " -- skipping." << std::endl;
            filesMissing++;
            continue;
        }

        TG4Event* event = nullptr;
        tree->SetBranchAddress("Event", &event);

        int nEvents = tree->GetEntries();
        totalEvents += nEvents;
        filesOk++;

        for (int i = 0; i < nEvents; i++) {
            tree->GetEntry(i);
            for (const auto& traj : event->Trajectories) {
                int pdg  = traj.GetPDGCode();
                int npts = (int)traj.Points.size();

                std::string key = "";
                if (pdg == 211)                key = "PiPlus (211)";
                else if (pdg == -211)           key = "PiMinus(-211)";
                else if (pdg == 2212)          key = "proton (2212)";
                else if (pdg == 2112)          key = "neutron(2112)";
                else if (pdg == 321)           key = "K+     (321)";
                else if (pdg == -321)          key = "K-     (-321)";
                else continue;

                auto& s = stats[key];
                s.n++;
                s.sum += npts;
                s.min  = std::min(s.min, npts);
                s.max  = std::max(s.max, npts);
            }
        }

        f.Close();

        if ((g + 1) % PROGRESS_EVERY == 0) {
            std::cout << "  ...processed " << (g + 1) << "/" << N_FILES
                       << " files (" << totalEvents << " events so far)" << std::endl;
        }
    }

    std::cout << "\n=== Summary over " << filesOk << "/" << N_FILES
               << " files (" << filesMissing << " missing/unreadable) ===" << std::endl;
    std::cout << "Total events: " << totalEvents << std::endl;

    std::cout << "\n" << std::left
              << std::setw(16) << "Particle"
              << std::setw(12) << "Count"
              << std::setw(10) << "Avg pts"
              << std::setw(8)  << "Min"
              << std::setw(8)  << "Max" << std::endl;
    std::cout << std::string(54, '-') << std::endl;

    for (auto& [key, s] : stats) {
        if (s.n > 0) {
            std::cout << std::left << std::setw(16) << key
                       << std::setw(12) << s.n
                       << std::setw(10) << (s.sum / s.n)
                       << std::setw(8)  << s.min
                       << std::setw(8)  << s.max << std::endl;
        } else {
            std::cout << std::left << std::setw(16) << key
                       << "none found" << std::endl;
        }
    }
}