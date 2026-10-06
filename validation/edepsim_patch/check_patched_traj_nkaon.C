// check_patched_traj_nkaon.C
//
// Checks average, min, max nPoints for pions, protons, neutrons, and
// kaons (K+/K-) in a patched edep-sim output file, to verify the extended
// bypass condition (std::abs(pdg)==211 || pdg==2212 || pdg==2112 ||
// std::abs(pdg)==321) is actually taking effect for the newly-added PDG
// codes (2112, 321, -321), not just the originally-patched ones (211, 2212).
//
// Run with:
//   root -l -b -q 'check_patched_traj_nkaon.C("/path/to/file.root")'

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

struct ChannelStats {
    int    n     = 0;
    double sum   = 0;
    int    min   = 999999;
    int    max   = 0;
};

void check_patched_traj_nkaon(const std::string& fname)
{
    TFile f(fname.c_str());
    if (f.IsZombie()) {
        std::cerr << "ERROR: Cannot open " << fname << std::endl;
        return;
    }

    TTree* tree = (TTree*)f.Get("EDepSimEvents");
    if (!tree) {
        std::cerr << "ERROR: No EDepSimEvents in " << fname << std::endl;
        return;
    }

    TG4Event* event = nullptr;
    tree->SetBranchAddress("Event", &event);

    int nEvents = tree->GetEntries();
    std::cout << "File: " << fname << "  (" << nEvents << " events)" << std::endl;

    std::map<std::string, ChannelStats> stats;
    stats["pi+/-  (211)"]  = ChannelStats();
    stats["proton (2212)"] = ChannelStats();
    stats["neutron(2112)"] = ChannelStats();
    stats["K+     (321)"]  = ChannelStats();
    stats["K-     (-321)"] = ChannelStats();

    for (int i = 0; i < nEvents; i++) {
        tree->GetEntry(i);
        for (const auto& traj : event->Trajectories) {
            int pdg  = traj.GetPDGCode();
            int npts = (int)traj.Points.size();

            std::string key = "";
            if (std::abs(pdg) == 211)      key = "pi+/-  (211)";
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

    std::cout << "\n" << std::left
              << std::setw(16) << "Particle"
              << std::setw(10) << "Count"
              << std::setw(10) << "Avg pts"
              << std::setw(8)  << "Min"
              << std::setw(8)  << "Max" << std::endl;
    std::cout << std::string(52, '-') << std::endl;

    for (auto& [key, s] : stats) {
        if (s.n > 0) {
            std::cout << std::left << std::setw(16) << key
                       << std::setw(10) << s.n
                       << std::setw(10) << (s.sum / s.n)
                       << std::setw(8)  << s.min
                       << std::setw(8)  << s.max << std::endl;
        } else {
            std::cout << std::left << std::setw(16) << key
                       << "none found in this file" << std::endl;
        }
    }

    std::cout << "\nExpected: pi+/- ~45 pts, proton ~9 pts (already-patched,"
                 " should be unaffected by this change)." << std::endl;
    std::cout << "If the extension worked: neutron/K+/K- averages should be"
                 " similarly large (tens of points), NOT ~2 pts (which would"
                 " indicate the bypass condition did not take effect for"
                 " these PDG codes)." << std::endl;
}