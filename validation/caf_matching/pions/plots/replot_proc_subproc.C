// replot_proc_subproc.C
//
// Redraws the proc/subproc histogram directly from scan_all_proc_subproc.root
// -- no EDepSim files touched, no rescanning. Use this any time the plot
// styling needs a tweak (colors, labels, canvas size) without repeating
// the ~90-minute full-sample scan.
//
// Also builds weight (f=1.33) vs length 2D histograms, one per distinct
// (proc,subproc) last-point tag, by joining scan_all_proc_subproc.root's
// proc_subproc_by_track tree against caf_pion_list.root's own
// length/weight columns via (cafFileIdx, RunId, EventId, TrackId).
//
// Run with:
//   root -l -b -q replot_proc_subproc.C

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"

#include <map>
#include <string>
#include <tuple>

void plot_weight_vs_length_by_tag();   // forward declaration -- defined below, called from replot_proc_subproc()

void replot_proc_subproc() {

    gROOT->LoadMacro("protoDUNEStyle.C");
    gROOT->SetStyle("protoDUNEStyle");
    gROOT->ForceStyle();
    gSystem->mkdir("2d_plots", true);

    TFile* f = TFile::Open("scan_all_proc_subproc.root");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: cannot open scan_all_proc_subproc.root\n";
        return;
    }
    TH1D* h = (TH1D*)f->Get("h_proc_subproc");
    if (!h) {
        std::cerr << "ERROR: no h_proc_subproc histogram found in file.\n";
        return;
    }

    // Adjust styling here as needed -- e.g.:
    // h->SetFillColor(kOrange+1);

    TCanvas* c = new TCanvas("c_proc_subproc_replot", "", 800, 600);
    c->SetLeftMargin(0.13);
    c->SetBottomMargin(0.17);
    h->GetYaxis()->SetTitleOffset(1.3);
    h->GetXaxis()->SetTitle("Geant4 Process/Subprocess");
    h->Draw("HIST");

    TLatex tL; tL.SetNDC();
    tL.DrawLatex(0.13, 0.94, "#bf{DUNE:ND-LAr 2x2}");
    TLatex tR; tR.SetNDC(); tR.SetTextSize(0.035);
    tR.DrawLatex(0.55, 0.94, "Truth #pi^{+} (CAF), last-point tag");

    TLatex tLeg; tLeg.SetNDC(); tLeg.SetTextSize(0.03);
    tLeg.DrawLatex(0.55, 0.75, "4/121 = hadron inelastic");
    tLeg.DrawLatex(0.55, 0.70, "2/22 = EM/scintillation");
    tLeg.DrawLatex(0.55, 0.65, "6/201 = decay");

    c->SaveAs("2d_plots/proc_subproc_hist_replot.png");
    std::cout << "Wrote 2d_plots/proc_subproc_hist_replot.png\n";

    f->Close();

    plot_weight_vs_length_by_tag();
}

// ---------------------------------------------------------------------
// Weight (f=1.33) vs length, one 2D histogram per distinct (proc,subproc)
// last-point tag.
// ---------------------------------------------------------------------
void plot_weight_vs_length_by_tag() {

    TFile* fTags = TFile::Open("scan_all_proc_subproc.root");
    if (!fTags || fTags->IsZombie()) {
        std::cerr << "ERROR: cannot open scan_all_proc_subproc.root\n";
        return;
    }
    TTree* tTags = (TTree*)fTags->Get("proc_subproc_by_track");
    if (!tTags) {
        std::cerr << "ERROR: no proc_subproc_by_track tree found -- this file was written\n"
                   << "before the per-track tree was added. Re-run scan_all_proc_subproc.c++\n"
                   << "with that edit applied first.\n";
        fTags->Close();
        return;
    }

    int t_cafFileIdx, t_RunId, t_EventId, t_TrackId, t_lastProc, t_lastSubproc;
    tTags->SetBranchAddress("cafFileIdx", &t_cafFileIdx);
    tTags->SetBranchAddress("RunId",      &t_RunId);
    tTags->SetBranchAddress("EventId",    &t_EventId);
    tTags->SetBranchAddress("TrackId",    &t_TrackId);
    tTags->SetBranchAddress("lastProc",   &t_lastProc);
    tTags->SetBranchAddress("lastSubproc",&t_lastSubproc);

    // Build lookup: (cafFileIdx,RunId,EventId,TrackId) -> "proc/subproc" label
    std::map<std::tuple<int,int,int,int>, std::string> tagLookup;
    Long64_t nTags = tTags->GetEntries();
    for (Long64_t i = 0; i < nTags; ++i) {
        tTags->GetEntry(i);
        std::string label = std::to_string(t_lastProc) + "/" + std::to_string(t_lastSubproc);
        tagLookup[{t_cafFileIdx, t_RunId, t_EventId, t_TrackId}] = label;
    }
    std::cout << "Loaded " << tagLookup.size() << " track tags for the length/weight join.\n";
    fTags->Close();

    TFile* fCaf = TFile::Open("caf_pion_list.root");
    if (!fCaf || fCaf->IsZombie()) {
        std::cerr << "ERROR: cannot open caf_pion_list.root\n";
        return;
    }
    TTree* tCaf = (TTree*)fCaf->Get("caf_pion_list");
    if (!tCaf) {
        std::cerr << "ERROR: no caf_pion_list tree found.\n";
        fCaf->Close();
        return;
    }

    int c_cafFileIdx, c_RunId, c_EventId, c_TrackId, c_pdg;
    double c_length, c_w133;
    tCaf->SetBranchAddress("cafFileIdx", &c_cafFileIdx);
    tCaf->SetBranchAddress("RunId",      &c_RunId);
    tCaf->SetBranchAddress("EventId",    &c_EventId);
    tCaf->SetBranchAddress("TrackId",    &c_TrackId);
    tCaf->SetBranchAddress("pdg",        &c_pdg);
    tCaf->SetBranchAddress("length",     &c_length);
    tCaf->SetBranchAddress("w_1p33",     &c_w133);

    const int nBinsX = 80, nBinsY = 100;
    const double xMin = 0.0, xMax = 400.0, yMin = 0.0, yMax = 2.0;

    std::map<std::string, TH2D*> histsByTag;
    long nMatched = 0, nUnmatched = 0;

    Long64_t nCaf = tCaf->GetEntries();
    for (Long64_t i = 0; i < nCaf; ++i) {
        tCaf->GetEntry(i);
        if (c_pdg != 211) continue;

        auto it = tagLookup.find({c_cafFileIdx, c_RunId, c_EventId, c_TrackId});
        if (it == tagLookup.end()) { nUnmatched++; continue; }
        const std::string& label = it->second;

        if (histsByTag.find(label) == histsByTag.end()) {
            std::string hname = "h_wlen_tag_" + label;
            for (auto& ch : hname) if (ch == '/') ch = '_';   // sanitize for ROOT object name
            std::string title = ";Track Length [cm];G4RW Weight";
            TH2D* newHist = new TH2D(hname.c_str(), title.c_str(), nBinsX, xMin, xMax, nBinsY, yMin, yMax);
            newHist->SetDirectory(0);   // detach from fCaf so it survives fCaf->Close() below
            histsByTag[label] = newHist;
        }
        histsByTag[label]->Fill(c_length, c_w133);
        nMatched++;
    }
    fCaf->Close();

    std::cout << "Matched " << nMatched << " pi+ entries to a tag (skipped " << nUnmatched
               << " with no tag match).\n";

    gSystem->mkdir("2d_plots", true);
    for (auto& kv : histsByTag) {
        const std::string& label = kv.first;
        TH2D* hist = kv.second;

        // Sanitize ONCE, reuse for both the canvas name and the output
        // filename -- a raw "/" in a ROOT object name (e.g. TCanvas)
        // causes ROOT's internal object bookkeeping to misbehave and
        // segfault, which is exactly what happened when the canvas name
        // was left unsanitized in an earlier version of this script.
        std::string safeLabel = label;
        for (auto& ch : safeLabel) if (ch == '/') ch = '_';

        TCanvas* c2 = new TCanvas(("c_wlen_" + safeLabel).c_str(), "", 800, 600);
        c2->SetLeftMargin(0.14);
        c2->SetRightMargin(0.14);

        hist->GetYaxis()->SetTitleOffset(1.3);
        hist->Draw("COLZ");

        TLatex tL2; tL2.SetNDC();
        tL2.DrawLatex(0.13, 0.94, "#bf{DUNE:ND-LAr 2x2}");
        TLatex tR2; tR2.SetNDC(); tR2.SetTextSize(0.035);
        tR2.DrawLatex(0.50, 0.94, ("Truth #pi^{+} (CAF) (f=1.33) for " + label).c_str());
        
        std::string outPng = "2d_plots/weight_vs_length_tag_" + safeLabel + "_f1p33.png";
        c2->SaveAs(outPng.c_str());
        std::cout << "  [" << label << "] " << hist->GetEntries() << " entries -> " << outPng << "\n";
        delete c2;
    }

    std::cout << "\nDone.\n";
}