// dump_caf_proton_list_eq3p10.C
//
// Dumps EVERY entry in the CAF's truth-proton tree (syst_g4rw_proton_truth)
// to a flat output ROOT tree + a printed table, with every field needed
// to trace a proton back through EDepSim and the G4RW weight files:
//   cafFileIdx, RunId, EventId (decoded from nu_id), TrackId (G4ID),
//   pdg, length, cosTheta, and w_patched[5] if present.
//
// Output:
//   caf_proton_list_eq3p10.root   -- TTree "caf_proton_list" with one row per proton
//   caf_proton_list_eq3p10.txt    -- same content, plain text (easy to grep/sort)
//
// Run with:
//   root -l -b -q dump_caf_proton_list_eq3p10.C

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

const std::string kOutRoot = "caf_proton_list_eq3p10.root";
const std::string kOutTxt  = "caf_proton_list_eq3p10.txt";

bool HasBranch(TTree* t, const std::string& name) {
    return t->GetBranch(name.c_str()) != nullptr;
}

void dump_caf_proton_list_eq3p10() {

    TFile* fIn = TFile::Open(kCafFile.c_str());
    if (!fIn || fIn->IsZombie()) {
        std::cerr << "ERROR: cannot open CAF file: " << kCafFile << "\n";
        return;
    }
    TTree* tIn = (TTree*)fIn->Get("syst_g4rw_proton_truth");
    if (!tIn) {
        std::cerr << "ERROR: no syst_g4rw_proton_truth tree in CAF file.\n";
        fIn->Close();
        return;
    }

    std::cout << "Branches actually present in syst_g4rw_proton_truth:\n";
    for (auto* b : *tIn->GetListOfBranches()) {
        std::cout << "  " << b->GetName() << "\n";
    }

    // --- Required branches (already confirmed working elsewhere) ---
    if (!HasBranch(tIn, "cafFileIdx") || !HasBranch(tIn, "nu_id") || !HasBranch(tIn, "G4ID")) {
        std::cerr << "ERROR: one of the required branches (cafFileIdx, nu_id, G4ID) is missing.\n"
                   << "Cannot proceed -- check the branch list printed above.\n";
        fIn->Close();
        return;
    }

    int cafFileIdx, G4ID;
    long long nu_id;
    tIn->SetBranchAddress("cafFileIdx", &cafFileIdx);
    tIn->SetBranchAddress("nu_id",      &nu_id);
    tIn->SetBranchAddress("G4ID",       &G4ID);

    // --- Opportunistic branches ---
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

    if (!hasPdg)      std::cerr << "WARNING: no 'pdg' branch -- pdg column will be -999.\n";
    if (!hasCosTheta) std::cerr << "WARNING: no 'cosTheta' branch -- cosTheta column will be -999.\n";
    if (!hasLength)   std::cerr << "WARNING: no 'length' branch -- length column will be -999.\n";
    if (!hasWPatched) std::cerr << "WARNING: no 'w_patched' branch -- weight columns will be -999.\n"
                                 << "         (Fixed-array convention assumed: [0]=0.34,"
                                 << " [1]=0.67, [2]=1.33, [3]=1.66, [4]=2.00 -- same convention"
                                 << " as the pion tree. NOT independently verified here.\n";

    bool hasEndXYZ = HasBranch(tIn, "endX") && HasBranch(tIn, "endY") && HasBranch(tIn, "endZ");
    double endX = -9999, endY = -9999, endZ = -9999;
    if (hasEndXYZ) {
        tIn->SetBranchAddress("endX", &endX);
        tIn->SetBranchAddress("endY", &endY);
        tIn->SetBranchAddress("endZ", &endZ);
    } else {
        std::cerr << "WARNING: no endX/endY/endZ branches -- truth end-position columns will be -9999.\n";
    }


    // --- Reco proton tree join (for weight-vs-recoLength comparison) ---
    // syst_g4rw_proton_reco has NO weight of its own, only reco length --
    // but shares the same (nu_id, cafFileIdx, G4ID) keys as the truth
    // tree, so it can be joined to pull reco length alongside truth's
    // own weight for the same proton.
    struct RecoInfo { double length, endX, endY, endZ; };
    std::map<std::tuple<long long,int,int>, RecoInfo> recoInfoMap;
    TTree* tReco = (TTree*)fIn->Get("syst_g4rw_proton_reco");
    bool recoHasEndXYZ = tReco && tReco->GetBranch("endX") && tReco->GetBranch("endY") && tReco->GetBranch("endZ");
    if (tReco) {
        long long r_nu_id; int r_cafFileIdx, r_G4ID; double r_length;
        double r_endX = -9999, r_endY = -9999, r_endZ = -9999;
        tReco->SetBranchAddress("nu_id",      &r_nu_id);
        tReco->SetBranchAddress("cafFileIdx", &r_cafFileIdx);
        tReco->SetBranchAddress("G4ID",       &r_G4ID);
        tReco->SetBranchAddress("length",     &r_length);
        if (recoHasEndXYZ) {
            tReco->SetBranchAddress("endX", &r_endX);
            tReco->SetBranchAddress("endY", &r_endY);
            tReco->SetBranchAddress("endZ", &r_endZ);
        } else {
            std::cerr << "WARNING: no endX/endY/endZ branches in syst_g4rw_proton_reco -- recoEnd* columns will be -9999.\n";
        }
        for (Long64_t i = 0; i < tReco->GetEntries(); ++i) {
            tReco->GetEntry(i);
            recoInfoMap[{r_nu_id, r_cafFileIdx, r_G4ID}] = {r_length, r_endX, r_endY, r_endZ};
        }
        std::cout << "Loaded " << recoInfoMap.size() << " reco proton entries for length/endpos join.\n";
    } else {
        std::cerr << "WARNING: no syst_g4rw_proton_reco tree found -- recoLength/recoEnd* columns will be -999 for all rows.\n";
    }

    // --- Output tree ---
    TFile* fOut = new TFile(kOutRoot.c_str(), "RECREATE");
    TTree* tOut = new TTree("caf_proton_list", "Flat list of CAF-matched truth protons, for traceable candidate selection");

    int out_cafFileIdx, out_RunId, out_EventId, out_TrackId, out_pdg;
    double out_cosTheta, out_length;
    double out_w034, out_w067, out_w133, out_w166, out_w200;

    tOut->Branch("cafFileIdx", &out_cafFileIdx, "cafFileIdx/I");
    tOut->Branch("RunId",      &out_RunId,      "RunId/I");
    tOut->Branch("EventId",    &out_EventId,    "EventId/I");
    tOut->Branch("TrackId",    &out_TrackId,    "TrackId/I");
    tOut->Branch("pdg",        &out_pdg,        "pdg/I");
    tOut->Branch("cosTheta",   &out_cosTheta,   "cosTheta/D");
    tOut->Branch("length",     &out_length,     "length/D");
    tOut->Branch("w_0p34",     &out_w034,       "w_0p34/D");
    tOut->Branch("w_0p67",     &out_w067,       "w_0p67/D");
    tOut->Branch("w_1p33",     &out_w133,       "w_1p33/D");
    tOut->Branch("w_1p66",     &out_w166,       "w_1p66/D");
    tOut->Branch("w_2p00",     &out_w200,       "w_2p00/D");
    double out_recoLength; int out_hasReco;
    tOut->Branch("recoLength", &out_recoLength, "recoLength/D");
    tOut->Branch("hasReco",    &out_hasReco,    "hasReco/I");
    double out_endX, out_endY, out_endZ;
    tOut->Branch("endX", &out_endX, "endX/D");
    tOut->Branch("endY", &out_endY, "endY/D");
    tOut->Branch("endZ", &out_endZ, "endZ/D");
    double out_recoEndX, out_recoEndY, out_recoEndZ;
    tOut->Branch("recoEndX", &out_recoEndX, "recoEndX/D");
    tOut->Branch("recoEndY", &out_recoEndY, "recoEndY/D");
    tOut->Branch("recoEndZ", &out_recoEndZ, "recoEndZ/D");

    std::ofstream txt(kOutTxt);
    txt << std::setw(10) << "cafFileIdx" << std::setw(10) << "RunId" << std::setw(12) << "EventId"
        << std::setw(10) << "TrackId" << std::setw(8) << "pdg"
        << std::setw(12) << "cosTheta" << std::setw(12) << "length"
        << std::setw(10) << "w_0.34" << std::setw(10) << "w_0.67" << std::setw(10) << "w_1.33"
        << std::setw(10) << "w_1.66" << std::setw(10) << "w_2.00"
        << std::setw(12) << "recoLen" << std::setw(8) << "hasReco"
        << std::setw(10) << "endX" << std::setw(10) << "endY" << std::setw(10) << "endZ" << "\n";

    Long64_t N = tIn->GetEntries();
    std::cout << "\nDumping " << N << " entries from syst_g4rw_proton_truth ...\n";

    for (Long64_t i = 0; i < N; ++i) {
        tIn->GetEntry(i);

        out_cafFileIdx = cafFileIdx;
        out_RunId      = (int)(nu_id / 1000000);
        out_EventId    = (int)(nu_id % 1000000);
        out_TrackId    = G4ID;
        out_pdg        = pdg;
        out_cosTheta   = cosTheta;
        out_length     = length;
        out_w034 = w_patched[0]; out_w067 = w_patched[1]; out_w133 = w_patched[2];
        out_w166 = w_patched[3]; out_w200 = w_patched[4];

        out_endX = endX; out_endY = endY; out_endZ = endZ;

        out_recoLength = -999; out_hasReco = 0;
        out_recoEndX = -9999; out_recoEndY = -9999; out_recoEndZ = -9999;
        auto it = recoInfoMap.find({nu_id, cafFileIdx, G4ID});
        if (it != recoInfoMap.end()) {
            out_recoLength = it->second.length; out_hasReco = 1;
            out_recoEndX = it->second.endX; out_recoEndY = it->second.endY; out_recoEndZ = it->second.endZ;
        }

        tOut->Fill();

        txt << std::setw(10) << out_cafFileIdx << std::setw(10) << out_RunId << std::setw(12) << out_EventId
            << std::setw(10) << out_TrackId << std::setw(8) << out_pdg
            << std::setw(12) << std::fixed << std::setprecision(4) << out_cosTheta
            << std::setw(12) << std::setprecision(3) << out_length
            << std::setw(10) << std::setprecision(4) << out_w034 << std::setw(10) << out_w067
            << std::setw(10) << out_w133 << std::setw(10) << out_w166 << std::setw(10) << out_w200
            << std::setw(12) << out_recoLength << std::setw(8) << out_hasReco
            << std::setw(10) << std::setprecision(2) << out_endX
            << std::setw(10) << out_endY << std::setw(10) << out_endZ << "\n";
    }

    tOut->Write();
    fOut->Close();
    txt.close();
    fIn->Close();

    std::cout << "\nDone. Wrote " << N << " rows to:\n"
               << "  " << kOutRoot << "  (TTree 'caf_proton_list' -- fast to re-scan/cut with TTree::Draw)\n"
               << "  " << kOutTxt  << "  (plain text -- grep/sort by eye)\n"
               << "\nEvery row here is, by construction, already present in the CAF --\n"
               << "safe to pick candidates from THIS list going forward,\n"
               << "instead of the raw g4rwEvt weight file.\n";
}

int main() {
    dump_caf_proton_list_eq3p10();
    return 0;
}