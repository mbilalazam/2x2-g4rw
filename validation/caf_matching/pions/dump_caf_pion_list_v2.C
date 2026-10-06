// dump_caf_pion_list_v2.C
//
// FIX vs dump_caf_pion_list.C: that script decoded RunId/EventId from nu_id
// via nu_id/1000000 and nu_id%1000000. That decode does NOT match the real
// matching convention used by purity_dlp_multiplicity.cc itself when it
// filled w_patched:
//
//   WKey5 g4rwKey = {g4rw_fileIdx, (int)g4rw_nu_id, g4rw_G4ID};
//   // matched directly against g4rwEvt's (file_idx, event_id, track_id)
//
// i.e. nu_id IS event_id, undecoded, and cafFileIdx IS file_idx. The old
// decode was an invented convention with no basis in the actual weight
// lookup logic, and was lossy for any nu_id >= 1e6 (only nu_id % 1e6
// survived under the wrong name "EventId"; the true value needed for
// matching was discarded).
//
// This version stores file_idx and event_id directly, matching the real
// convention, so a join against any g4rwEvt file (old formula, Eq. 3.10,
// or otherwise) uses the exact same key purity_dlp_multiplicity.cc used.
//
// Output:
//   caf_pion_list_v2.root   -- TTree "caf_pion_list" with corrected keys
//   caf_pion_list_v2.txt    -- same content, plain text
//
// Run with:
//   root -l -b -q dump_caf_pion_list_v2.C

#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <map>
#include <tuple>

const std::string kCafFile =
    "/global/u1/m/mazam/Geant4ReWeighting/2x2_G4RW/test/on_CAFS/31Aug2026/output.root";

const std::string kOutRoot = "caf_pion_list_v2.root";
const std::string kOutTxt  = "caf_pion_list_v2.txt";

bool HasBranch(TTree* t, const std::string& name) {
    return t->GetBranch(name.c_str()) != nullptr;
}

void dump_caf_pion_list_v2() {

    TFile* fIn = TFile::Open(kCafFile.c_str());
    if (!fIn || fIn->IsZombie()) {
        std::cerr << "ERROR: cannot open CAF file: " << kCafFile << "\n";
        return;
    }
    TTree* tIn = (TTree*)fIn->Get("syst_g4rw_pion_truth");
    if (!tIn) {
        std::cerr << "ERROR: no syst_g4rw_pion_truth tree in CAF file.\n";
        fIn->Close();
        return;
    }

    if (!HasBranch(tIn, "cafFileIdx") || !HasBranch(tIn, "nu_id") || !HasBranch(tIn, "G4ID")) {
        std::cerr << "ERROR: one of the required branches (cafFileIdx, nu_id, G4ID) is missing.\n";
        fIn->Close();
        return;
    }

    int cafFileIdx, G4ID;
    long long nu_id;
    tIn->SetBranchAddress("cafFileIdx", &cafFileIdx);
    tIn->SetBranchAddress("nu_id",      &nu_id);
    tIn->SetBranchAddress("G4ID",       &G4ID);

    bool hasPdg      = HasBranch(tIn, "pdg");
    bool hasCosTheta = HasBranch(tIn, "cosTheta");
    bool hasLength   = HasBranch(tIn, "length");
    bool hasWPatched = HasBranch(tIn, "w_patched");

    int pdg = -999;
    double cosTheta = -999, length = -999;
    double w_patched[5] = {-999,-999,-999,-999,-999};

    if (hasPdg)      tIn->SetBranchAddress("pdg", &pdg);
    if (hasCosTheta) tIn->SetBranchAddress("cosTheta", &cosTheta);
    if (hasLength)   tIn->SetBranchAddress("length", &length);
    if (hasWPatched) tIn->SetBranchAddress("w_patched", w_patched);

    bool hasEndXYZ = HasBranch(tIn, "endX") && HasBranch(tIn, "endY") && HasBranch(tIn, "endZ");
    double endX = -9999, endY = -9999, endZ = -9999;
    if (hasEndXYZ) {
        tIn->SetBranchAddress("endX", &endX);
        tIn->SetBranchAddress("endY", &endY);
        tIn->SetBranchAddress("endZ", &endZ);
    }

    // --- Output tree ---
    TFile* fOut = new TFile(kOutRoot.c_str(), "RECREATE");
    TTree* tOut = new TTree("caf_pion_list", "Flat list of CAF-matched truth pions (corrected keys)");

    // NOTE: file_idx / event_id now match g4rwEvt's naming and semantics
    // directly -- this IS the real matching convention, not a derived one.
    int out_file_idx, out_event_id, out_track_id, out_pdg;
    double out_cosTheta, out_length;
    double out_w034, out_w067, out_w133, out_w166, out_w200;

    tOut->Branch("file_idx",  &out_file_idx,  "file_idx/I");
    tOut->Branch("event_id",  &out_event_id,  "event_id/I");
    tOut->Branch("track_id",  &out_track_id,  "track_id/I");
    tOut->Branch("pdg",       &out_pdg,       "pdg/I");
    tOut->Branch("cosTheta",  &out_cosTheta,  "cosTheta/D");
    tOut->Branch("length",    &out_length,    "length/D");
    tOut->Branch("w_0p34",    &out_w034,      "w_0p34/D");
    tOut->Branch("w_0p67",    &out_w067,      "w_0p67/D");
    tOut->Branch("w_1p33",    &out_w133,      "w_1p33/D");
    tOut->Branch("w_1p66",    &out_w166,      "w_1p66/D");
    tOut->Branch("w_2p00",    &out_w200,      "w_2p00/D");
    double out_endX, out_endY, out_endZ;
    tOut->Branch("endX", &out_endX, "endX/D");
    tOut->Branch("endY", &out_endY, "endY/D");
    tOut->Branch("endZ", &out_endZ, "endZ/D");

    std::ofstream txt(kOutTxt);
    txt << std::setw(10) << "file_idx" << std::setw(12) << "event_id"
        << std::setw(10) << "track_id" << std::setw(8) << "pdg"
        << std::setw(12) << "cosTheta" << std::setw(12) << "length"
        << std::setw(10) << "w_0.34" << std::setw(10) << "w_0.67" << std::setw(10) << "w_1.33"
        << std::setw(10) << "w_1.66" << std::setw(10) << "w_2.00" << "\n";

    Long64_t N = tIn->GetEntries();
    std::cout << "Dumping " << N << " entries from syst_g4rw_pion_truth ...\n";

    long n_overflow = 0;   // nu_id that doesn't fit in an int -- worth knowing about

    for (Long64_t i = 0; i < N; ++i) {
        tIn->GetEntry(i);

        out_file_idx = cafFileIdx;               // matches g4rwEvt's file_idx directly
        if (nu_id > INT_MAX || nu_id < INT_MIN) n_overflow++;
        out_event_id = (int)nu_id;               // matches g4rwEvt's event_id directly, UNDECODED
        out_track_id = G4ID;                     // matches g4rwEvt's track_id directly
        out_pdg      = pdg;
        out_cosTheta = cosTheta;
        out_length   = length;
        out_w034 = w_patched[0]; out_w067 = w_patched[1]; out_w133 = w_patched[2];
        out_w166 = w_patched[3]; out_w200 = w_patched[4];
        out_endX = endX; out_endY = endY; out_endZ = endZ;

        tOut->Fill();

        txt << std::setw(10) << out_file_idx << std::setw(12) << out_event_id
            << std::setw(10) << out_track_id << std::setw(8) << out_pdg
            << std::setw(12) << std::fixed << std::setprecision(4) << out_cosTheta
            << std::setw(12) << std::setprecision(3) << out_length
            << std::setw(10) << std::setprecision(4) << out_w034 << std::setw(10) << out_w067
            << std::setw(10) << out_w133 << std::setw(10) << out_w166 << std::setw(10) << out_w200
            << "\n";
    }

    tOut->Write();
    fOut->Close();
    txt.close();
    fIn->Close();

    std::cout << "\nDone. Wrote " << N << " rows to " << kOutRoot << " and " << kOutTxt << "\n";
    if (n_overflow > 0) {
        std::cout << "WARNING: " << n_overflow << " entries had nu_id outside int range "
                  << "-- these may need long long handling if this fires.\n";
    }
    std::cout << "\nKeys now match purity_dlp_multiplicity.cc's real lookup convention:\n"
              << "  file_idx = cafFileIdx, event_id = nu_id (undecoded), track_id = G4ID\n"
              << "-- the same (file_idx, event_id, track_id) semantics used by every\n"
              << "g4rwEvt weight file (old formula, Eq. 3.10, etc).\n";
}

int main() {
    dump_caf_pion_list_v2();
    return 0;
}