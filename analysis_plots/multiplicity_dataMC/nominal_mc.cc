#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TEfficiency.h"

#include <iostream>
#include <fstream>
#include <string>

#include "duneanaobj/StandardRecord/StandardRecord.h" // Ideally, this should be SRProxy.h, but there is an include error for that now.
// Alternatively, you can use SetBranchStatus function in TreeLoader, but it does not work for the common branch (to do)

double minTrkLength = 3;

long nRecoMatchedProtonAll          = 0; 
long nRecoProtonMatchedContained    = 0; 
long nRecoProtonMatchedNonContained = 0;

long nTrueMatchedProtonAll          = 0; 
long nTrueProtonMatchedContained    = 0; 
long nTrueProtonMatchedNonContained = 0;

long nRecoMatchedPionAll          = 0; 
long nRecoPionMatchedContained    = 0; 
long nRecoPionMatchedNonContained = 0;

long nTrueMatchedPionAll          = 0; 
long nTruePionMatchedContained    = 0; 
long nTruePionMatchedNonContained = 0;

bool requireMinerva = true;   // false = MINERvA OFF, true = ON


double signed3dDistance(
	double firstPoint1, double firstPoint2, double firstPoint3,
	double secondPoint1, double secondPoint2, double secondPoint3,
	TVector3 point
) {
	double denominator = sqrt(
		(secondPoint2 - firstPoint2) * (secondPoint2 - firstPoint2) +
		(secondPoint1 - firstPoint1) * (secondPoint1 - firstPoint1) +
		(secondPoint3 - firstPoint3) * (secondPoint3 - firstPoint3)
	);

	double X1 = point.X() - firstPoint1;
	double Y1 = point.Y() - firstPoint2;
	double Z1 = point.Z() - firstPoint3;

	double X2 = point.X() - secondPoint1;
	double Y2 = point.Y() - secondPoint2;
	double Z2 = point.Z() - secondPoint3;

	double numerator = (X1 * Y2 - Y1 * X2) - (X1 * Z2 - X2 * Z1) + (Y1 * Z2 - Z1 * Y2);

	return numerator / denominator;
}


int caf_plotter(std::string input_file_list, std::string output_rootfile, bool mcOnly) {

	double minLength = 0;
	int goodIntNum = 0;
	int badIntNum = 0;
	int goodIntInFidVol = 0;
	int badIntInFidVol = 0;
	int goodMINERvAMatch = 0;
	int totalMINERvAMatch = 0;
	double offset = 0;
			
	// Counters for MINERvA matched
	int count_match_total = 0, count_match_qe = 0, count_match_mec = 0, count_match_dis = 0;
	int count_match_res = 0, count_match_coh = 0, count_match_nc = 0, count_match_rock = 0, count_match_sec = 0;

	// Counters for MINERvA not matched
	int count_alt_total = 0, count_alt_qe = 0, count_alt_mec = 0, count_alt_dis = 0;
	int count_alt_res = 0, count_alt_coh = 0, count_alt_nc = 0, count_alt_rock = 0, count_alt_sec = 0;

	// Reconstructed distributions
	TH1D* part_energy_hist    = new TH1D("recpart_energy", "Reco particle energy in GeV", 100, 0, 1); // Distribution of reconstructed particle energy (GeV)
	TH1D* track_length        = new TH1D("track_length", "track_length", 100, 0, 100); // Distribution of reconstructed track lengths
	TH1D* part_multTrkOnly    = new TH1D("part_multTrkOnly", "part_multTrkOnly", 20, 0, 20); // Primary particle multiplicity with track length > threshold
	TH1D* part_mult           = new TH1D("part_mult", "part_mult", 20, 0, 20); // Total reconstructed primary particle multiplicity
	TH1D* bad_origin          = new TH1D("rock", "rock", 2, 0, 2); // Type of bad interaction origin: rock (1) or secondary (0)
	
	// Track multiplicity categorized by GENIE interaction mode
	TH1D* track_mult          = new TH1D("track_mult", "track_mult", 20, 0, 20); 
	TH1D* track_multMode      = new TH1D("track_multMode", "track_multMode", 20, 0, 20);
	TH1D* track_multQE        = new TH1D("track_multQE", "track_multQE", 20, 0, 20);
	TH1D* track_multMEC       = new TH1D("track_multMEC", "track_multMEC", 20, 0, 20);
	TH1D* track_multRES       = new TH1D("track_multRES", "track_multRES", 20, 0, 20);
	TH1D* track_multCOH       = new TH1D("track_multCOH", "track_multCOH", 20, 0, 20);
	TH1D* track_multDIS       = new TH1D("track_multDIS", "track_multDIS", 20, 0, 20);
	TH1D* track_multNC        = new TH1D("track_multNC", "track_multNC", 20, 0, 20);
	TH1D* track_multRock      = new TH1D("track_multRock", "track_multRock", 20, 0, 20); 			// Track multiplicity for interactions backtracked to rock
	TH1D* track_multSec       = new TH1D("track_multSec", "track_multSec", 20, 0, 20); 				// Track multiplicity for secondary (non-neutrino) interactions
	TH1D* track_multBad       = new TH1D("track_multBad", "track_multBad", 20, 0, 20); 				// Track multiplicity for poorly matched reco interactions
	TH1D* track_multGood      = new TH1D("track_multGood", "track_multGood", 20, 0, 20);			// Track multiplicity for well-matched reco interactions

TH1D* track_multGood_muon = new TH1D("track_multGood_muon", "track_multGood_muon", 20, 0, 20);
TH1D* track_multGood_pion = new TH1D("track_multGood_pion", "track_multGood_pion", 20, 0, 20);
TH1D* track_multGood_proton = new TH1D("track_multGood_proton", "track_multGood_proton", 20, 0, 20);
TH1D* track_multGood_kaon = new TH1D("track_multGood_kaon", "track_multGood_kaon", 20, 0, 20);


TH1D* hadron_mult         = new TH1D("hadron_mult", "hadron_mult", 20, 0, 20); 
TH1D* hadron_multMode     = new TH1D("hadron_multMode", "hadron_multMode", 20, 0, 20);
TH1D* hadron_multQE       = new TH1D("hadron_multQE", "hadron_multQE", 20, 0, 20);
TH1D* hadron_multMEC      = new TH1D("hadron_multMEC", "hadron_multMEC", 20, 0, 20);
TH1D* hadron_multRES      = new TH1D("hadron_multRES", "hadron_multRES", 20, 0, 20);
TH1D* hadron_multCOH      = new TH1D("hadron_multCOH", "hadron_multCOH", 20, 0, 20);
TH1D* hadron_multDIS      = new TH1D("hadron_multDIS", "hadron_multDIS", 20, 0, 20);
TH1D* hadron_multNC       = new TH1D("hadron_multNC", "hadron_multNC", 20, 0, 20);
TH1D* hadron_multRock     = new TH1D("hadron_multRock", "hadron_multRock", 20, 0, 20);
TH1D* hadron_multSec      = new TH1D("hadron_multSec", "hadron_multSec", 20, 0, 20);
TH1D* hadron_multBad      = new TH1D("hadron_multBad", "hadron_multBad", 20, 0, 20);
TH1D* hadron_multGood     = new TH1D("hadron_multGood", "hadron_multGood", 20, 0, 20);


	// Truth-level particle histograms (energy, length, cosTheta)
	// True Muon distributions
	TH1D* true_muon_energy = new TH1D("true_muon_energy", "True Muon Energy;E [GeV];Events", 100, 0, 2);
	TH1D* true_muon_length = new TH1D("true_muon_length", "True Muon Track Length;Length [cm];Events", 100, 0, 100);
	TH1D* true_muon_cosTheta = new TH1D("true_muon_cosTheta", "True Muon cos(#theta);cos(#theta);Events", 100, -1, 1);
	TH1D* true_muon_startX = new TH1D("true_muon_startX", "True Muon Start X;X [cm];Events", 70, -70, 70);
	TH1D* true_muon_startY = new TH1D("true_muon_startY", "True Muon Start Y;Y [cm];Events", 70, -70, 70);
	TH1D* true_muon_startZ = new TH1D("true_muon_startZ", "True Muon Start Z;Z [cm];Events", 70, -70, 70);
	TH2D* true_muon_length_vs_energy = new TH2D(
		"true_muon_length_vs_energy",
		"True Muon: Track Length vs Energy;Track Length [cm];Energy [GeV]",
		100, 0, 100,   // length binning same as true_muon_length
		100, 0, 2      // energy binning same as true_muon_energy
	);


	// True Charged Pion distributions (π+ and π-)
	TH1D* true_pion_energy = new TH1D("true_pion_energy", "True Charged Pion Energy;E [GeV];Events", 100, 0, 2);
	TH1D* true_pion_length = new TH1D("true_pion_length", "True Charged Pion Track Length;Length [cm];Events", 100, 0, 100);
	TH1D* true_pion_cosTheta = new TH1D("true_pion_cosTheta", "True Charged Pion cos(#theta);cos(#theta);Events", 100, -1, 1);
	TH1D* true_pion_startX = new TH1D("true_pion_startX", "True Charged Pion Start X;X [cm];Events", 70, -70, 70);
	TH1D* true_pion_startY = new TH1D("true_pion_startY", "True Charged Pion Start Y;Y [cm];Events", 70, -70, 70);
	TH1D* true_pion_startZ = new TH1D("true_pion_startZ", "True Charged Pion Start Z;Z [cm];Events", 70, -70, 70);
	TH1D* true_pion_endX = new TH1D("true_pion_endX", "True Charged Pion End X;X [cm];Events", 70, -70, 70);
	TH1D* true_pion_endY = new TH1D("true_pion_endY", "True Charged Pion End Y;Y [cm];Events", 70, -70, 70);
	TH1D* true_pion_endZ = new TH1D("true_pion_endZ", "True Charged Pion End Z;Z [cm];Events", 70, -70, 70);
	TH1D* true_pion_momentumX = new TH1D("true_pion_momentumX", "True Charged Pion Momentum X;P_{x} [GeV/c];Events", 100, 0, 2);
	TH1D* true_pion_momentumY = new TH1D("true_pion_momentumY", "True Charged Pion Momentum Y;P_{y} [GeV/c];Events", 100, 0, 2);
	TH1D* true_pion_momentumZ = new TH1D("true_pion_momentumZ", "True Charged Pion Momentum Z;P_{z} [GeV/c];Events", 100, 0, 2);
	TH2D* true_pion_length_vs_energy = new TH2D(
		"true_pion_length_vs_energy",
		"True Charged Pion: Track Length vs Energy;Track Length [cm];Energy [GeV]",
		100, 0, 100,   // length bins (same as true_pion_length)
		100, 0, 2      // energy bins (same as true_pion_energy)
	);
	TH2D* true_pion_startX_vs_startY = new TH2D(
		"true_pion_startX_vs_startY",
		"True Charged Pion: StartX vs StartY;Start X [cm];Start Y [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* true_pion_startZ_vs_startX = new TH2D(
		"true_pion_startZ_vs_startX",
		"True Charged Pion: StartZ vs StartX;Start Z [cm];Start X [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* true_pion_startZ_vs_startY = new TH2D(
		"true_pion_startZ_vs_startY",
		"True Charged Pion: StartZ vs StartY;Start Z [cm];Start Y [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* true_pion_startX_vs_endX = new TH2D(
		"true_pion_startX_vs_endX",
		"True Charged Pion: StartX vs EndX;Start X [cm];End X [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* true_pion_startY_vs_endY = new TH2D(
		"true_pion_startY_vs_endY",
		"True Charged Pion: StartY vs EndY;Start Y [cm];End Y [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* true_pion_startZ_vs_endZ = new TH2D(
		"true_pion_startZ_vs_endZ",
		"True Charged Pion: StartZ vs EndZ;Start Z [cm];End Z [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* true_pion_length_vs_cosTheta = new TH2D(
		"true_pion_length_vs_cosTheta",
		"True Charged Pion: Length vs cos(#theta);Track Length [cm];cos(#theta)",
		100, 0, 100,
		100, -1, 1
	);
	TH2D* true_pion_energy_vs_cosTheta = new TH2D(
		"true_pion_energy_vs_cosTheta",
		"True Charged Pion: Energy vs cos(#theta);Energy [GeV];cos(#theta)",
		100, 0, 2,
		100, -1, 1
	);
	TH2D* true_pion_dX_vs_length = new TH2D(
		"true_pion_dX_vs_length",
		"True Charged Pion: #DeltaX vs Length;Track Length [cm];#DeltaX [cm]",
		100, 0, 100,
		100, -140, 140
	);
	TH2D* true_pion_dY_vs_length = new TH2D(
		"true_pion_dY_vs_length",
		"True Charged Pion: #DeltaY vs Length;Track Length [cm];#DeltaY [cm]",
		100, 0, 100,
		100, -140, 140
	);
	TH2D* true_pion_dZ_vs_length = new TH2D(
		"true_pion_dZ_vs_length",
		"True Charged Pion: #DeltaZ vs Length;Track Length [cm];#DeltaZ [cm]",
		100, 0, 100,
		100, -140, 140
	);	
	TH2D* true_pion_endX_vs_endY = new TH2D(
		"true_pion_endX_vs_endY",
		"True Charged Pion: EndX vs EndY;End X [cm];End Y [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* true_pion_endZ_vs_endX = new TH2D(
		"true_pion_endZ_vs_endX",
		"True Charged Pion: EndZ vs EndX;End Z [cm];End X [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* true_pion_endZ_vs_endY = new TH2D(
		"true_pion_endZ_vs_endY",
		"True Charged Pion: EndZ vs EndY;End Z [cm];End Y [cm]",
		70, -70, 70,
		70, -70, 70
	);	

	
	// TH1D* true_pion_energy = new TH1D("true_pion_energy", "True Charged Pion Energy;E [GeV];Events", 100, 0, 2);
	// TH1D* true_pion_length = new TH1D("true_pion_length", "True Charged Pion Track Length;Length [cm];Events", 100, 0, 100);
	// TH1D* true_pion_cosTheta = new TH1D("true_pion_cosTheta", "True Charged Pion cos(#theta);cos(#theta);Events", 100, -1, 1);
	// TH1D* true_pion_startX = new TH1D("true_pion_startX", "True Charged Pion Start X;X [cm];Events", 70, -70, 70);
	// TH1D* true_pion_startY = new TH1D("true_pion_startY", "True Charged Pion Start Y;Y [cm];Events", 70, -70, 70);
	// TH1D* true_pion_startZ = new TH1D("true_pion_startZ", "True Charged Pion Start Z;Z [cm];Events", 70, -70, 70);
	// TH2D* true_pion_length_vs_energy = new TH2D(
		// "true_pion_length_vs_energy",
		// "True Charged Pion: Track Length vs Energy;Track Length [cm];Energy [GeV]",
		// 100, 0, 100,   // same as true_pion_length
		// 100, 0, 2      // same as true_pion_energy
	// );


	// True Proton distributions
	TH1D* true_proton_energy = new TH1D("true_proton_energy", "True Proton Energy;E [GeV];Events", 100, 0, 2);
	TH1D* true_proton_length = new TH1D("true_proton_length", "True Proton Track Length;Length [cm];Events", 100, 0, 100);
	TH1D* true_proton_cosTheta = new TH1D("true_proton_cosTheta", "True Proton cos(#theta);cos(#theta);Events", 100, -1, 1);
	TH1D* true_proton_startX = new TH1D("true_proton_startX", "True Proton Start X;X [cm];Events", 70, -70, 70);
	TH1D* true_proton_startY = new TH1D("true_proton_startY", "True Proton Start Y;Y [cm];Events", 70, -70, 70);
	TH1D* true_proton_startZ = new TH1D("true_proton_startZ", "True Proton Start Z;Z [cm];Events", 70, -70, 70);
	TH1D* true_proton_endX = new TH1D("true_proton_endX", "True Proton End X;X [cm];Events", 70, -70, 70);
	TH1D* true_proton_endY = new TH1D("true_proton_endY", "True Proton End Y;Y [cm];Events", 70, -70, 70);
	TH1D* true_proton_endZ = new TH1D("true_proton_endZ", "True Proton End Z;Z [cm];Events", 70, -70, 70);
	TH1D* true_proton_momentumX = new TH1D("true_proton_momentumX", "True Proton Momentum X;P_{x} [GeV/c];Events", 100, 0, 2);
	TH1D* true_proton_momentumY = new TH1D("true_proton_momentumY", "True Proton Momentum Y;P_{y} [GeV/c];Events", 100, 0, 2);
	TH1D* true_proton_momentumZ = new TH1D("true_proton_momentumZ", "True Proton Momentum Z;P_{z} [GeV/c];Events", 100, 0, 2);
	TH2D* true_proton_length_vs_energy = new TH2D(
		"true_proton_length_vs_energy",
		"True Proton: Track Length vs Energy;Track Length [cm];Energy [GeV]",
		100, 0, 100,   // length bins (same as true_proton_length)
		100, 0, 2      // energy bins (same as true_proton_energy)
	);
	TH2D* true_proton_startX_vs_startY = new TH2D(
		"true_proton_startX_vs_startY",
		"True Proton: StartX vs StartY;Start X [cm];Start Y [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* true_proton_startZ_vs_startX = new TH2D(
		"true_proton_startZ_vs_startX",
		"True Proton: StartZ vs StartX;Start Z [cm];Start X [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* true_proton_startZ_vs_startY = new TH2D(
		"true_proton_startZ_vs_startY",
		"True Proton: StartZ vs StartY;Start Z [cm];Start Y [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* true_proton_startX_vs_endX = new TH2D(
		"true_proton_startX_vs_endX",
		"True Proton: StartX vs EndX;Start X [cm];End X [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* true_proton_startY_vs_endY = new TH2D(
		"true_proton_startY_vs_endY",
		"True Proton: StartY vs EndY;Start Y [cm];End Y [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* true_proton_startZ_vs_endZ = new TH2D(
		"true_proton_startZ_vs_endZ",
		"True Proton: StartZ vs EndZ;Start Z [cm];End Z [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* true_proton_length_vs_cosTheta = new TH2D(
		"true_proton_length_vs_cosTheta",
		"True Proton: Length vs cos(#theta);Track Length [cm];cos(#theta)",
		100, 0, 100,
		100, -1, 1
	);
	TH2D* true_proton_energy_vs_cosTheta = new TH2D(
		"true_proton_energy_vs_cosTheta",
		"True Proton: Energy vs cos(#theta);Energy [GeV];cos(#theta)",
		100, 0, 2,
		100, -1, 1
	);
	TH2D* true_proton_dX_vs_length = new TH2D(
		"true_proton_dX_vs_length",
		"True Proton: #DeltaX vs Length;Track Length [cm];#DeltaX [cm]",
		100, 0, 100,
		100, -140, 140
	);
	TH2D* true_proton_dY_vs_length = new TH2D(
		"true_proton_dY_vs_length",
		"True Proton: #DeltaY vs Length;Track Length [cm];#DeltaY [cm]",
		100, 0, 100,
		100, -140, 140
	);
	TH2D* true_proton_dZ_vs_length = new TH2D(
		"true_proton_dZ_vs_length",
		"True Proton: #DeltaZ vs Length;Track Length [cm];#DeltaZ [cm]",
		100, 0, 100,
		100, -140, 140
	);	
	TH2D* true_proton_endX_vs_endY = new TH2D(
		"true_proton_endX_vs_endY",
		"True Proton: EndX vs EndY;End X [cm];End Y [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* true_proton_endZ_vs_endX = new TH2D(
		"true_proton_endZ_vs_endX",
		"True Proton: EndZ vs EndX;End Z [cm];End X [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* true_proton_endZ_vs_endY = new TH2D(
		"true_proton_endZ_vs_endY",
		"True Proton: EndZ vs EndY;End Z [cm];End Y [cm]",
		70, -70, 70,
		70, -70, 70
	);	
	
	// True Charged Kaon distributions (K+ and K-)
	TH1D* true_kaon_energy = new TH1D("true_kaon_energy", "True Charged Kaon Energy;E [GeV];Events", 100, 0, 2);
	TH1D* true_kaon_length = new TH1D("true_kaon_length", "True Charged Kaon Track Length;Length [cm];Events", 100, 0, 100);
	TH1D* true_kaon_cosTheta = new TH1D("true_kaon_cosTheta", "True Charged Kaon cos(#theta);cos(#theta);Events", 100, -1, 1);
	TH1D* true_kaon_startX = new TH1D("true_kaon_startX", "True Charged Kaon Start X;X [cm];Events", 70, -70, 70);
	TH1D* true_kaon_startY = new TH1D("true_kaon_startY", "True Charged Kaon Start Y;Y [cm];Events", 70, -70, 70);
	TH1D* true_kaon_startZ = new TH1D("true_kaon_startZ", "True Charged Kaon Start Z;Z [cm];Events", 70, -70, 70);

	// Truth-level neutrino energy histograms (NO selection cuts - all events)
	TH1F* true_nuE_all_nosel = new TH1F("true_nuE_all_nosel", "True Neutrino Energy (No Selection);E_{#nu} [GeV];Events", 100, 0, 20);

	TH1F* true_nuE_QE_nosel  = new TH1F("true_nuE_QE_nosel", "True Neutrino Energy (QE, No Selection);E_{#nu} [GeV];Events", 100, 0, 20);
	TH1F* true_nuE_RES_nosel = new TH1F("true_nuE_RES_nosel", "True Neutrino Energy (RES, No Selection);E_{#nu} [GeV];Events", 100, 0, 20);
	TH1F* true_nuE_DIS_nosel = new TH1F("true_nuE_DIS_nosel", "True Neutrino Energy (DIS, No Selection);E_{#nu} [GeV];Events", 100, 0, 20);
	TH1F* true_nuE_MEC_nosel = new TH1F("true_nuE_MEC_nosel", "True Neutrino Energy (MEC, No Selection);E_{#nu} [GeV];Events", 100, 0, 20);
	TH1F* true_nuE_COH_nosel = new TH1F("true_nuE_COH_nosel", "True Neutrino Energy (COH, No Selection);E_{#nu} [GeV];Events", 100, 0, 20);

// True track-like multiplicity by interaction mode (NO selection cuts)
	TH1D* true_multTrkOnly_QE_nosel     = new TH1D("true_multTrkOnly_QE_nosel", "true_multTrkOnly_QE_nosel", 20, 0, 20);
	TH1D* true_multTrkOnly_MEC_nosel    = new TH1D("true_multTrkOnly_MEC_nosel", "true_multTrkOnly_MEC_nosel", 20, 0, 20);
	TH1D* true_multTrkOnly_RES_nosel    = new TH1D("true_multTrkOnly_RES_nosel", "true_multTrkOnly_RES_nosel", 20, 0, 20);
	TH1D* true_multTrkOnly_DIS_nosel    = new TH1D("true_multTrkOnly_DIS_nosel", "true_multTrkOnly_DIS_nosel", 20, 0, 20);
	TH1D* true_multTrkOnly_COH_nosel    = new TH1D("true_multTrkOnly_COH_nosel", "true_multTrkOnly_COH_nosel", 20, 0, 20);
	TH1D* true_multTrkOnly_nosel        = new TH1D("true_multTrkOnly_nosel", "true_multTrkOnly_nosel", 20, 0, 20);

	// Truth-level neutrino energy histograms (after FV + CC cuts)
	TH1F* true_nuE_all = new TH1F("true_nuE_all", "True Neutrino Energy;E_{#nu} [GeV];Events", 100, 0, 20);

	TH1F* true_nuE_QE  = new TH1F("true_nuE_QE", "True Neutrino Energy (QE);E_{#nu} [GeV];Events", 100, 0, 20);
	TH1F* true_nuE_RES = new TH1F("true_nuE_RES", "True Neutrino Energy (RES);E_{#nu} [GeV];Events", 100, 0, 20);
	TH1F* true_nuE_DIS = new TH1F("true_nuE_DIS", "True Neutrino Energy (DIS);E_{#nu} [GeV];Events", 100, 0, 20);
	TH1F* true_nuE_MEC = new TH1F("true_nuE_MEC", "True Neutrino Energy (MEC);E_{#nu} [GeV];Events", 100, 0, 20);
	TH1F* true_nuE_COH = new TH1F("true_nuE_COH", "True Neutrino Energy (COH);E_{#nu} [GeV];Events", 100, 0, 20);

	// Truth-level neutrino vertices
	TH1F* true_vtxX = new TH1F("true_vtxX", "True Neutrino Vertex X;X [cm];Events", 70, -70, 70);
	TH1F* true_vtxY = new TH1F("true_vtxY", "True Neutrino Vertex Y;Y [cm];Events", 70, -70, 70);
	TH1F* true_vtxZ = new TH1F("true_vtxZ", "True Neutrino Vertex Z;Z [cm];Events", 70, -70, 70);

	// Reconstructed neutrino vertices
	TH1D* recoVtxX = new TH1D("recoVtxX", "Reco Neutrino Vertex X;X [cm];Events", 70, -70, 70);
	TH1D* recoVtxY = new TH1D("recoVtxY", "Reco Neutrino Vertex Y;Y [cm];Events", 70, -70, 70);
	TH1D* recoVtxZ = new TH1D("recoVtxZ", "Reco Neutrino Vertex Z;Z [cm];Events", 70, -70, 70);

	// Reconstructed neutrino vertices (Truth Matched)
	TH1D* truthMatchedRecoVtxX = new TH1D("truthMatchedRecoVtxX", "Truth (Matched) Reco Neutrino Vertex X;X [cm];Events", 70, -70, 70);
	TH1D* truthMatchedRecoVtxY = new TH1D("truthMatchedRecoVtxY", "Truth (Matched) Reco Neutrino Vertex Y;Y [cm];Events", 70, -70, 70);
	TH1D* truthMatchedRecoVtxZ = new TH1D("truthMatchedRecoVtxZ", "Truth (Matched) Reco Neutrino Vertex Z;Z [cm];Events", 70, -70, 70);


	// Truth neutrino energy for truth-matched reconstructed interactions
	TH1F* truthMatchedReco_nuE = new TH1F("truthMatchedReco_nuE", "True Nu Energy (Truth-Matched Reco);E_{#nu} [GeV];Events", 100, 0, 20);
	TH1F* truthMatchedReco_nuE_QE = new TH1F("truthMatchedReco_nuE_QE", "True Nu Energy (Truth-Matched Reco QE);E_{#nu} [GeV];Events", 100, 0, 20);
	TH1F* truthMatchedReco_nuE_RES = new TH1F("truthMatchedReco_nuE_RES", "True Nu Energy (Truth-Matched Reco RES);E_{#nu} [GeV];Events", 100, 0, 20);
	TH1F* truthMatchedReco_nuE_DIS = new TH1F("truthMatchedReco_nuE_DIS", "True Nu Energy (Truth-Matched Reco DIS);E_{#nu} [GeV];Events", 100, 0, 20);
	TH1F* truthMatchedReco_nuE_MEC = new TH1F("truthMatchedReco_nuE_MEC", "True Nu Energy (Truth-Matched Reco MEC);E_{#nu} [GeV];Events", 100, 0, 20);
	TH1F* truthMatchedReco_nuE_COH = new TH1F("truthMatchedReco_nuE_COH", "True Nu Energy (Truth-Matched Reco COH);E_{#nu} [GeV];Events", 100, 0, 20);

	// Classification of reconstructed tracks based on truth matching:
	// Bin 0: Correct match → Reco track backtracks to correct primary (muon, proton, pion, kaon) with sufficient overlap
	// Bin 1: Wrong type → Matched to truth, but PDG is unexpected (e.g., neutron, electron, or other)
	// Bin 2: No match → Either no backtracking or overlap below threshold; classified as unmatched
	// Bin 3: Secondary match → Reco track matches to a secondary particle (PartType::kSecondary)
	TH1D* trackCorrectness = new TH1D("track_correctness", "track_correctness", 4, 0, 4); 

	// Classification of reconstructed showers based on truth matching:
	// Bin 0: Correct match → Shower backtracks to primary electromagnetic particle (e⁻, γ, π⁰) with sufficient overlap
	// Bin 1: Wrong type → Matched to a particle like neutron or other non-EM source
	// Bin 2: No match → No usable truth match found (e.g., low overlap or unclassified)
	// Bin 3: Secondary match → Shower matches to a secondary particle (PartType::kSecondary)
	TH1D* showerCorrectness = new TH1D("shower_corrrectness", "shower_correctness", 4, 0, 4); 

	// Reconstructed interaction vertex
	TH1D* recoHistVertexY     = new TH1D("recoHistVertexY", "recoHistVertexY", 70, -338, -198);
	TH1D* recoHistVertexX     = new TH1D("recoHistVertexX", "recoHistVertexX", 70, -70, 70);
	TH1D* recoHistVertexZ     = new TH1D("recoHistVertexZ", "recoHistVertexZ", 70, 1230, 1370);

	// Difference between true and primary multiplicity
	TH1D* diffPip             = new TH1D("diffPip", "diffPip", 8, -4, 4);
	TH1D* diffPim             = new TH1D("diffPim", "diffPim", 8, -4, 4);
	TH1D* diffP               = new TH1D("diffP", "diffP", 8, -4, 4);
	TH1D* diffN               = new TH1D("diffN", "diff", 8, -4, 4);

	// Backtracked CC classification of interactions in argon, interactions in rock, and secondary interactions 
	TH1D* recoBacktrackCCAr   = new TH1D("recoBacktrackCCAr", "recoBacktrackCCAr", 2, 0, 2);
	TH1D* recoBacktrackCCRock = new TH1D("recoBacktrackCCRock", "recoBacktrackCCRock", 2, 0, 2);
	TH1D* recoBacktrackCCSec  = new TH1D("recoBacktrackCCSec", "recoBacktrackCCSec", 2, 0, 2);

	// True vs reco response matrices
	TH2D* responseMult        = new TH2D("responseMult", "responseMult", 15, 0, 15, 15, 0, 15);			// track multiplicity
	TH2D* responseCosL        = new TH2D("responseCosL", "responseCosL", 15, 0.85, 1, 15, 0.85, 1);		// muon cos(theta)
	TH2D* responseMultAlt        = new TH2D("responseMultAlt", "responseMultAlt", 15, 0, 15, 15, 0, 15);			// track multiplicity


	// Backtracked PDG
	TH1D* recoBacktrackPDGAr  = new TH1D("recoBacktrackPDGAr", "recoBacktrackPDGAr", 40, -20, 20);		// for Ar-matched events
	TH1D* recoBacktrackPDGRock= new TH1D("recoBacktrackPDGRock", "recoBacktrackPDGRock", 40, -20, 20);	// for rock events
	TH1D* recoBacktrackPDGSec = new TH1D("recoBacktrackPDGSec", "recoBacktrackPDGSec", 40, -20, 20);	// for secondary events

	// Ar-matched reco distributions
	TH1D* recoBacktrackElAr     = new TH1D("recoBacktrackElAr", "recoBacktrackElAr", 50, 0, 20);		// Reco muon energy (Ar-matched)
	TH1D* recoBacktrackCoslAr  = new TH1D("recoBacktrackCoslAr", "recoBackTrackCoslAr", 100, 0.85, 1);	// Reco muon cos(theta) (Ar-matched)
	
	// Truth particle matched to MINERvA
	TH1D* minervaMatchPDG      = new TH1D("minervaMatchPDG", "minervaMatchPDG", 6000, -3000, 3000);		// PDG 
	TH1D* minervaMatchE        = new TH1D("minervaMatchE", "minervaMatchE", 50, 0, 20);					// Energy
	TH1D* minervaMatchCos      = new TH1D("minervaMatchCos", "minervaMatchCos", 100, -1, 1);			// cos(theta)

	// Reco energy of matched pion, proton, kaon tracks
	TH1D* selPionE             = new TH1D("selPionE", "selPionE", 15, 0, 0.2);
	TH1D* selProtonE           = new TH1D("selProtonE", "selProtonE", 15, 0, 0.2);
	TH1D* selKaonE             = new TH1D("selKaonE", "selKaonE", 15, 0, 0.2);

	// Reco direction of matched pions
	TH1D* selPionDirX          = new TH1D("selPionDirX", "selPionDirX", 20, -1, 1);
	TH1D* selPionDirY          = new TH1D("selPionDirY", "selPionDirY", 20, -1, 1);
	TH1D* selPionDirZ          = new TH1D("selPionDirZ", "selPionDirZ", 20, -1, 1);

	// Reco direction of matched protons
	TH1D* selProtonDirX        = new TH1D("selProtonDirX", "selProtonDirX", 20, -1, 1);
	TH1D* selProtonDirY        = new TH1D("selProtonDirY", "selProtonDirY", 20, -1, 1);
	TH1D* selProtonDirZ        = new TH1D("selProtonDirZ", "selProtonDirZ", 20, -1, 1);

	// Reco direction of matched kaons
	TH1D* selKaonDirX          = new TH1D("selKaonDirX", "selKaonDirX", 20, -1, 1);
	TH1D* selKaonDirY          = new TH1D("selKaonDirY", "selKaonDirY", 20, -1, 1);
	TH1D* selKaonDirZ          = new TH1D("selKaonDirZ", "selKaonDirZ", 20, -1, 1);

	// Length of matched tracks 
	TH1D* selTrkLen            = new TH1D("selTrkLen", "selTrkLen", 20, 0, 10);				// all modes
	TH1D* selTrkLenNonQE       = new TH1D("selTrkLenNonQE", "selTrkLenNonQE", 20, 0, 10);	// non-QE
	TH1D* selTrkLenQE          = new TH1D("selTrkLenQE", "selTrkLenQE", 20, 0, 10);			// QE

	// Direction cosine of matched protons
	TH1D* selProtonXZ          = new TH1D("selProtonXZ", "selProtonXZ", 20, -1, 1);
	TH1D* selProtonYZ          = new TH1D("selProtonYZ", "selProtonYZ", 20, -1, 1);
	TH1D* selProtonXY          = new TH1D("selProtonXY", "selProtonXY", 20, -1, 1);

	// Direction cosine of matched pions
	TH1D* selPionXZ            = new TH1D("selPionXZ", "selPionXZ", 20, -1, 1);
	TH1D* selPionYZ            = new TH1D("selPionYZ", "selPionYZ", 20, -1, 1);
	TH1D* selPionXY            = new TH1D("selPionXY", "selPionXY", 20, -1, 1);

	// Direction cosine of matched kaons
	TH1D* selKaonXZ            = new TH1D("selKaonXZ", "selKaonXZ", 20, -1, 1);
	TH1D* selKaonYZ            = new TH1D("selKaonYZ", "selKaonYZ", 20, -1, 1);
	TH1D* selKaonXY            = new TH1D("selKaonXY", "selKaonXY", 20, -1, 1);

	// Reco vs truth PDG classification (μ, p, π, other)
	TH2D* confusionMatrix      = new TH2D("confusionMatrix", "confusionMatrix", 4, 0, 4, 4, 0, 4);
	TH2D* confusionMatrixAlt   = new TH2D("confusionMatrixAlt", "confusionMatrixAlt", 4, 0, 4, 4, 0, 4);

	TH2D* confusionMatrix_truthContained      = new TH2D("confusionMatrix_truthContained", "confusionMatrix_truthContained", 4, 0, 4, 4, 0, 4);
	TH2D* confusionMatrix_truthNotContained      = new TH2D("confusionMatrix_truthNotContained", "confusionMatrix_truthNotContained", 4, 0, 4, 4, 0, 4);


	// Number of true tracks in event
	TH1D* nPi                  = new TH1D("nPi", "nPi", 20, 0, 20);						// pions
	TH1D* nP                   = new TH1D("nP", "nP", 20, 0, 20);                       // protons
	TH1D* escapePi             = new TH1D("escapePi", "escapePi", 20, 0, 20);           // pions exiting detector

	// Track length of contained particles
	TH1D* containPiLen         = new TH1D("containPiLen", "containPiLen", 100, 0, 20);  // pions
	TH1D* containPLen          = new TH1D("containP", "containP", 100, 0, 20);          // protons

	// True KE in matched events
	TH1D* truePionEWithRecoInt        = new TH1D("truePionEWithRecoInt", "truePionEWithRecoInt", 15, 0, 0.2);		// pions
	TH1D* trueProtonEWithRecoInt      = new TH1D("trueProtonEWithRecoInt", "trueProtonEWithRecoInt", 15, 0, 0.2);   // protons
	TH1D* trueKaonEWithRecoInt        = new TH1D("trueKaonEWithRecoInt", "trueKaonEWithRecoInt", 15, 0, 0.2);       // kaons

	// Lengths of true tracks
	TH1D* trueTrkLenQE                = new TH1D("trueTrkLenQE", "trueTrkLenQE", 20, 0, 10);		// QE
	TH1D* trueTrkLenNonQE            = new TH1D("trueTrkLenNonQE", "trueTrkLenNonQE", 20, 0, 10);   // non-QE
	TH1D* trueTrkLen                 = new TH1D("trueTrkLen", "trueTrkLen", 20, 0, 10);             // All

	// True direction components of protons matched to reconstructed interactions
	TH1D* trueProtonWithRecoIntDirX  = new TH1D("trueProtonWithRecoIntDirX", "trueProtonWithRecoIntDirX", 20, -1, 1);
	TH1D* trueProtonWithRecoIntDirY  = new TH1D("trueProtonWithRecoIntDirY", "trueProtonWithRecoIntDirY", 20, -1, 1);
	TH1D* trueProtonWithRecoIntDirZ  = new TH1D("trueProtonWithRecoIntDirZ", "trueProtonWithRecoIntDirZ", 20, -1, 1);

	// True direction components of pions matched to reconstructed interactions
	TH1D* truePionWithRecoIntDirX    = new TH1D("truePionWithRecoIntDirX", "truePionWithRecoIntDirX", 20, -1, 1);
	TH1D* truePionWithRecoIntDirY    = new TH1D("truePionWithRecoIntDirY", "truePionWithRecoIntDirY", 20, -1, 1);
	TH1D* truePionWithRecoIntDirZ    = new TH1D("truePionWithRecoIntDirZ", "truePionWithRecoIntDirZ", 20, -1, 1);

	// True direction components of kaons matched to reconstructed interactions
	TH1D* trueKaonWithRecoIntDirX    = new TH1D("trueKaonWithRecoIntDirX", "trueKaonWithRecoIntDirX", 20, -1, 1);
	TH1D* trueKaonWithRecoIntDirY    = new TH1D("trueKaonWithRecoIntDirY", "trueKaonWithRecoIntDirY", 20, -1, 1);
	TH1D* trueKaonWithRecoIntDirZ    = new TH1D("trueKaonWithRecoIntDirZ", "trueKaonWithRecoIntDirZ", 20, -1, 1);

	// Angular correlations (cosine between axis pairs) for true protons in matched interactions
	TH1D* trueProtonWithRecoIntXY    = new TH1D("trueProtonWithRecoIntXY", "trueProtonWithRecoIntXY", 20, -1, 1);
	TH1D* trueProtonWithRecoIntXZ    = new TH1D("trueProtonWithRecoIntXZ", "trueProtonWithRecoIntXZ", 20, -1, 1);
	TH1D* trueProtonWithRecoIntYZ    = new TH1D("trueProtonWithRecoIntYZ", "trueProtonWithRecoIntYZ", 20, -1, 1);

	// Angular correlations (cosine between axis pairs) for true pions in matched interactions
	TH1D* truePionWithRecoIntXY		= new TH1D("truePionWithRecoIntXY", "truePionWithRecoIntXY", 20, -1, 1);
	TH1D* truePionWithRecoIntXZ		= new TH1D("truePionWithRecoIntXZ", "truePionWithRecoIntXZ", 20, -1, 1);
	TH1D* truePionWithRecoIntYZ		= new TH1D("truePionWithRecoIntYZ", "truePionWithRecoIntYZ", 20, -1, 1);

	// Angular correlations (cosine between axis pairs) for true kaons in matched interactions
	TH1D* trueKaonWithRecoIntXY		= new TH1D("trueKaonWithRecoIntXY", "trueKaonWithRecoIntXY", 20, -1, 1);
	TH1D* trueKaonWithRecoIntXZ		= new TH1D("trueKaonWithRecoIntXZ", "trueKaonWithRecoIntXZ", 20, -1, 1);
	TH1D* trueKaonWithRecoIntYZ		= new TH1D("trueKaonWithRecoIntYZ", "trueKaonWithRecoIntYZ", 20, -1, 1);

	// True particle multiplicity histograms (based on MC truth info)
	TH1D* true_mult					= new TH1D("true_mult", "true_mult", 20, 0, 20);					// All true primary particles
	TH1D* true_multTrkOnly			= new TH1D("true_multTrkOnly", "true_multTrkOnly", 20, 0, 20);		// True particles with track-like lengths
	TH1D* true_multGENIE			= new TH1D("true_multGENIE", "true_multGENIE", 20, 0, 20);			// From GENIE primaries only

// True track-like multiplicity by interaction mode
	TH1D* true_multTrkOnly_QE     = new TH1D("true_multTrkOnly_QE", "true_multTrkOnly_QE", 20, 0, 20);
	TH1D* true_multTrkOnly_MEC    = new TH1D("true_multTrkOnly_MEC", "true_multTrkOnly_MEC", 20, 0, 20);
	TH1D* true_multTrkOnly_RES    = new TH1D("true_multTrkOnly_RES", "true_multTrkOnly_RES", 20, 0, 20);
	TH1D* true_multTrkOnly_DIS    = new TH1D("true_multTrkOnly_DIS", "true_multTrkOnly_DIS", 20, 0, 20);
	TH1D* true_multTrkOnly_COH    = new TH1D("true_multTrkOnly_COH", "true_multTrkOnly_COH", 20, 0, 20);

TH1D* true_multTrkOnly_muon = new TH1D("true_multTrkOnly_muon", "true_multTrkOnly_muon", 20, 0, 20);
TH1D* true_multTrkOnly_pion = new TH1D("true_multTrkOnly_pion", "true_multTrkOnly_pion", 20, 0, 20);
TH1D* true_multTrkOnly_proton = new TH1D("true_multTrkOnly_proton", "true_multTrkOnly_proton", 20, 0, 20);
TH1D* true_multTrkOnly_kaon = new TH1D("true_multTrkOnly_kaon", "true_multTrkOnly_kaon", 20, 0, 20);

	// Reco-to-truth multiplicity response matrices
	TH2D* responseGenieToG4			= new TH2D("responseGenieToG4", "responseGenieToG4", 15, 0, 15, 15, 0, 15);			// GENIE multiplicity vs. track-only
	TH2D* recoVertex2DNoCuts		= new TH2D("recoVertex2DNoCuts", "recoVertex2DNoCuts", 70, -70, 70, 70, -70, 70);	// Reco vertex (X,Z) before cuts
	TH2D* recoVertex2D				= new TH2D("recoVertex2D", "recoVertex2D", 60, -60, 60, 60, -60, 60);				// Reco vertex (X,Z) after basic cuts

	// Matching-based histograms
	TH1D* matchTrue_mult			= new TH1D("matchTrue_mult", "matchTrue_mult", 20, 0, 20);					// Matched truth multiplicity
	TH1D* matchTrue_multTrkOnly		= new TH1D("matchTrue_multTrkOnly", "matchTrue_multTrkOnly", 20, 0, 20);	// Track-only version
	TH1D* matchTrue_multGENIE		= new TH1D("matchTrue_multGENIE", "matchTrue_multGENIE", 20, 0, 20);		// Matched GENIE multiplicity
	TH1D* matchHistEl				= new TH1D("matchHistEl", "matchHistEl", 50, 0, 20);						// Electron energy (matched)
	TH1D* matchHistCosl				= new TH1D("matchHistCosl", "matchHistCosl", 20, -1, 1);					// Muon cos(θ) from matched interaction

	// Muon kinematics (truth level)
	TH1D* histEl					= new TH1D("histEl", "histEl", 50, 0, 20);                          // Muon energy from truth
	TH1D* histCosl					= new TH1D("histCosl", "histCosl", 100, -1, 1);						// Muon cos(θ) from truth
	TH1D* trueDiffPosvsPDirZ		= new TH1D("trueDiffPosvsPDirZ", "trueDiffPosvsPDirZ", 20, -1, 1);	// Cos(θ) from displacement vs. momentum

	// Vertex diagnostics
	TH2D* recoVertex2DBadYZ			= new TH2D("recoVertex2DBadYZ", "recoVertex2DBadYZ", 70, -70, 70, 70, -70, 70);	// Y-Z projection for bad reco vertices
	TH2D* recoVertex2DBad			= new TH2D("recoVertex2DBad", "recoVertex2DBad", 70, -70, 70, 70, -70, 70);		// X-Z projection for bad reco vertices
	TH2D* trueVertex2DBad			= new TH2D("trueVertex2DBad", "trueVertex2DBad", 60, -60, 60, 60, -60, 60);		// True vertex X-Z for poorly reconstructed events


	// Other diagnostics
	TH1D* recoCosL					= new TH1D("recoCosL", "recoCosL", 100, -1, 1);  // Reco muon cos(θ)
	TH1D* longest_track				= new TH1D("longest_track", "longest_track", 100, 0, 100); // Length of longest reconstructed track in each interaction


	// Reco muon kinematics (truth matched)
	TH1D* recoMuonLength  = new TH1D("recoMuonLength", "Reco Muon Track Length", 100, 0, 100);
	TH1D* recoMuonCosL    = new TH1D("recoMuonCosL", "Reco Muon cos(θ)", 100, -1, 1);
	TH1D* recoMuonE       = new TH1D("recoMuonE", "Reco Muon Energy", 100, 0, 2);
	TH1D* recoMuonStartX  = new TH1D("recoMuonStartX", "Reco Muon Start X", 70, -70, 70);
	TH1D* recoMuonStartY  = new TH1D("recoMuonStartY", "Reco Muon Start Y", 70, -70, 70);
	TH1D* recoMuonStartZ  = new TH1D("recoMuonStartZ", "Reco Muon Start Z", 70, -70, 70);
	TH2D* recoMuonLength_vs_E = new TH2D(
		"recoMuonLength_vs_E",
		"Truth-Matched Reco Muon: Track Length vs Energy;Track Length [cm];Energy [GeV]",
		100, 0, 100,   // same binning as recoMuonLength
		100, 0, 2      // same binning as recoMuonE
	);

	TH1D* recoMuonLength_QE  = new TH1D("recoMuonLength_QE", "Reco Muon Length (QE)", 100, 0, 100);
	TH1D* recoMuonLength_RES = new TH1D("recoMuonLength_RES", "Reco Muon Length (RES)", 100, 0, 100);
	TH1D* recoMuonLength_DIS = new TH1D("recoMuonLength_DIS", "Reco Muon Length (DIS)", 100, 0, 100);
	TH1D* recoMuonLength_MEC = new TH1D("recoMuonLength_MEC", "Reco Muon Length (MEC)", 100, 0, 100);
	TH1D* recoMuonLength_COH = new TH1D("recoMuonLength_COH", "Reco Muon Length (COH)", 100, 0, 100);
	TH1D* recoMuonLength_Rock = new TH1D("recoMuonLength_Rock", "Reco Muon Length (Rock)", 100, 0, 100);
	TH1D* recoMuonLength_Sec = new TH1D("recoMuonLength_Sec", "Reco Muon Length (Secondary)", 100, 0, 100);

	// Reco pion kinematics (truth matched)
	TH1D* recoPionLength  = new TH1D("recoPionLength", "Reco Pion Track Length", 300, 0, 300); //100, 0, 100
	TH1D* recoPionCosL    = new TH1D("recoPionCosL", "Reco Pion cos(θ)", 100, -1, 1);
	TH1D* recoPionE       = new TH1D("recoPionE", "Reco Pion Energy", 100, 0, 2);
	TH1D* recoPionStartX  = new TH1D("recoPionStartX", "Reco Pion Start X", 70, -70, 70);
	TH1D* recoPionStartY  = new TH1D("recoPionStartY", "Reco Pion Start Y", 70, -70, 70);
	TH1D* recoPionStartZ  = new TH1D("recoPionStartZ", "Reco Pion Start Z", 70, -70, 70);

	TH1D* recoPionEndX = new TH1D("recoPionEndX", "Reco Pion End X;X [cm];Events", 70, -70, 70);
	TH1D* recoPionEndY = new TH1D("recoPionEndY", "Reco Pion End Y;Y [cm];Events", 70, -70, 70);
	TH1D* recoPionEndZ = new TH1D("recoPionEndZ", "Reco Pion End Z;Z [cm];Events", 70, -70, 70);
	TH1D* recoPionMomentumX = new TH1D("recoPionMomentumX", "Reco Pion Momentum X;P_{x} [GeV/c];Events", 100, 0, 2);
	TH1D* recoPionMomentumY = new TH1D("recoPionMomentumY", "Reco Pion Momentum Y;P_{y} [GeV/c];Events", 100, 0, 2);
	TH1D* recoPionMomentumZ = new TH1D("recoPionMomentumZ", "Reco Pion Momentum Z;P_{z} [GeV/c];Events", 100, 0, 2);
	TH2D* recoPionLength_vs_E = new TH2D(
		"recoPionLength_vs_E",
		"Truth-Matched Reco Pion: Track Length vs Energy;Track Length [cm];Energy [GeV]",
		100, 0, 100,   // same binning as recoPionLength
		100, 0, 2      // same binning as recoPionE
	);
	TH2D* recoPion_startX_vs_startY = new TH2D(
		"recoPion_startX_vs_startY",
		"Truth-Matched Reco Pion: StartX vs StartY;Start X [cm];Start Y [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* recoPion_startZ_vs_startX = new TH2D(
		"recoPion_startZ_vs_startX",
		"Truth-Matched Reco Pion: StartZ vs StartX;Start Z [cm];Start X [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* recoPion_startZ_vs_startY = new TH2D(
		"recoPion_startZ_vs_startY",
		"Truth-Matched Reco Pion: StartZ vs StartY;Start Z [cm];Start Y [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* recoPion_startX_vs_endX = new TH2D(
		"recoPion_startX_vs_endX",
		"Truth-Matched Reco Pion: StartX vs EndX;Start X [cm];End X [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* recoPion_startY_vs_endY = new TH2D(
		"recoPion_startY_vs_endY",
		"Truth-Matched Reco Pion: StartY vs EndY;Start Y [cm];End Y [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* recoPion_startZ_vs_endZ = new TH2D(
		"recoPion_startZ_vs_endZ",
		"Truth-Matched Reco Pion: StartZ vs EndZ;Start Z [cm];End Z [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* recoPion_length_vs_cosL = new TH2D(
		"recoPion_length_vs_cosL",
		"Truth-Matched Reco Pion: Length vs cos(#theta);Track Length [cm];cos(#theta)",
		100, 0, 100,
		100, -1, 1
	);
	TH2D* recoPion_energy_vs_cosL = new TH2D(
		"recoPion_energy_vs_cosL",
		"Truth-Matched Reco Pion: Energy vs cos(#theta);Energy [GeV];cos(#theta)",
		100, 0, 2,
		100, -1, 1
	);
	
	TH2D* recoPion_cosL_vs_startz = new TH2D(
		"recoPion_cosL_vs_startz",
		"Truth-Matched Pion: cos(#theta) vs StartZ;cos(#theta);Start Z [cm]",
		100, -1, 1,
		70, -70, 70
	);

	TH2D* recoPion_cosL_vs_endz = new TH2D(
		"recoPion_cosL_vs_endz",
		"Truth-Matched Pion: cos(#theta) vs EndZ;cos(#theta);End Z [cm]",
		100, -1, 1,
		70, -70, 70
	);	
	
	TH2D* recoPion_dX_vs_length = new TH2D(
		"recoPion_dX_vs_length",
		"Truth-Matched Reco Pion: #DeltaX vs Length;Track Length [cm];#DeltaX [cm]",
		100, 0, 100,
		100, -140, 140
	);
	TH2D* recoPion_dY_vs_length = new TH2D(
		"recoPion_dY_vs_length",
		"Truth-Matched Reco Pion: #DeltaY vs Length;Track Length [cm];#DeltaY [cm]",
		100, 0, 100,
		100, -140, 140
	);
	TH2D* recoPion_dZ_vs_length = new TH2D(
		"recoPion_dZ_vs_length",
		"Truth-Matched Reco Pion: #DeltaZ vs Length;Track Length [cm];#DeltaZ [cm]",
		100, 0, 100,
		100, -140, 140
	);
	TH2D* recoPion_endX_vs_endY = new TH2D(
		"recoPion_endX_vs_endY",
		"Truth-Matched Reco Pion: EndX vs EndY;End X [cm];End Y [cm]",
		70, -70, 70,
		70, -70, 70
	);

	TH2D* recoPion_endZ_vs_endX = new TH2D(
		"recoPion_endZ_vs_endX",
		"Truth-Matched Reco Pion: EndZ vs EndX;End Z [cm];End X [cm]",
		70, -70, 70,
		70, -70, 70
	);

	TH2D* recoPion_endZ_vs_endY = new TH2D(
		"recoPion_endZ_vs_endY",
		"Truth-Matched Reco Pion: EndZ vs EndY;End Z [cm];End Y [cm]",
		70, -70, 70,
		70, -70, 70
	);

	TH1D* recoPionCosAlpha = new TH1D(
		"recoPionCosAlpha",
		"Contained Matched Truth Pions;cos(#alpha);Number of Pions",
		50, -1.0, 1.0
	);

	TH2D* recoPionCosAlpha_vs_E = new TH2D(
		"recoPionCosAlpha_vs_E",
		"Contained Matched Truth Pions;Kinetic Energy [GeV];cos(#alpha)",
		50, 0.0, 2.0,
		50, -1.0, 1.0
	);

	TH2D* recoPionCosAlpha_vs_length = new TH2D(
		"recoPionCosAlpha_vs_length",
		"Contained Matched Truth Pions;Track Length [cm];cos(#alpha)",
		50, 0.0, 120.0,
		50, -1.0, 1.0
	);

	// Number of contained truth primary pions per interaction
	TH1D* primPionMultiplicity = new TH1D(
				"primPionMultiplicity", "Contained truth primary pion multiplicity;N_{prim p};Interactions",
				 20, 0, 20);
	
	// TH1D* recoPionLength  = new TH1D("recoPionLength", "Reco Pion Track Length", 100, 0, 100);
	// TH1D* recoPionCosL    = new TH1D("recoPionCosL", "Reco Pion cos(θ)", 100, -1, 1);
	// TH1D* recoPionE       = new TH1D("recoPionE", "Reco Pion Energy", 100, 0, 2);
	// TH1D* recoPionStartX  = new TH1D("recoPionStartX", "Reco Pion Start X", 70, -70, 70);
	// TH1D* recoPionStartY  = new TH1D("recoPionStartY", "Reco Pion Start Y", 70, -70, 70);
	// TH1D* recoPionStartZ  = new TH1D("recoPionStartZ", "Reco Pion Start Z", 70, -70, 70);
	// TH2D* recoPionLength_vs_E = new TH2D(
		// "recoPionLength_vs_E",
		// "Truth-Matched Reco Pion: Track Length vs Energy;Track Length [cm];Energy [GeV]",
		// 100, 0, 100,   // same binning as recoPionLength
		// 100, 0, 2      // same binning as recoPionE
	// );

	// Reco proton kinematics (truth matched)
	TH1D* recoProtonLength  = new TH1D("recoProtonLength", "Reco Proton Track Length", 300, 0, 300); // 100, 0, 100 // 300, 0, 300
	TH1D* recoProtonCosL    = new TH1D("recoProtonCosL", "Reco Proton cos(θ)", 100, -1, 1);
	TH1D* recoProtonE       = new TH1D("recoProtonE", "Reco Proton Energy", 100, 0, 2);
	TH1D* recoProtonStartX  = new TH1D("recoProtonStartX", "Reco Proton Start X", 70, -70, 70);
	TH1D* recoProtonStartY  = new TH1D("recoProtonStartY", "Reco Proton Start Y", 70, -70, 70);
	TH1D* recoProtonStartZ  = new TH1D("recoProtonStartZ", "Reco Proton Start Z", 70, -70, 70);

	TH1D* recoProtonEndX = new TH1D("recoProtonEndX", "Reco Proton End X;X [cm];Events", 200, -100, 300); //70, -70, 70
	TH1D* recoProtonEndY = new TH1D("recoProtonEndY", "Reco Proton End Y;Y [cm];Events", 200, -100, 300);
	TH1D* recoProtonEndZ = new TH1D("recoProtonEndZ", "Reco Proton End Z;Z [cm];Events", 200, -100, 300);
	TH1D* recoProtonMomentumX = new TH1D("recoProtonMomentumX", "Reco Proton Momentum X;P_{x} [GeV/c];Events", 100, 0, 2);
	TH1D* recoProtonMomentumY = new TH1D("recoProtonMomentumY", "Reco Proton Momentum Y;P_{y} [GeV/c];Events", 100, 0, 2);
	TH1D* recoProtonMomentumZ = new TH1D("recoProtonMomentumZ", "Reco Proton Momentum Z;P_{z} [GeV/c];Events", 100, 0, 2);
	TH2D* recoProtonLength_vs_E = new TH2D(
		"recoProtonLength_vs_E",
		"Truth-Matched Reco Proton: Track Length vs Energy;Track Length [cm];Energy [GeV]",
		300, 0, 300,   // same binning as recoProtonLength
		100, 0, 2      // same binning as recoProtonE
	);
	TH2D* recoProton_startX_vs_startY = new TH2D(
		"recoProton_startX_vs_startY",
		"Truth-Matched Reco Proton: StartX vs StartY;Start X [cm];Start Y [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* recoProton_startZ_vs_startX = new TH2D(
		"recoProton_startZ_vs_startX",
		"Truth-Matched Reco Proton: StartZ vs StartX;Start Z [cm];Start X [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* recoProton_startZ_vs_startY = new TH2D(
		"recoProton_startZ_vs_startY",
		"Truth-Matched Reco Proton: StartZ vs StartY;Start Z [cm];Start Y [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* recoProton_startX_vs_endX = new TH2D(
		"recoProton_startX_vs_endX",
		"Truth-Matched Reco Proton: StartX vs EndX;Start X [cm];End X [cm]",
		70, -70, 70,
		200, -100, 300
	);
	TH2D* recoProton_startY_vs_endY = new TH2D(
		"recoProton_startY_vs_endY",
		"Truth-Matched Reco Proton: StartY vs EndY;Start Y [cm];End Y [cm]",
		70, -70, 70,
		200, -100, 300
	);
	TH2D* recoProton_startZ_vs_endZ = new TH2D(
		"recoProton_startZ_vs_endZ",
		"Truth-Matched Reco Proton: StartZ vs EndZ;Start Z [cm];End Z [cm]",
		70, -70, 70,
		200, -100, 300
	);
	TH2D* recoProton_length_vs_cosL = new TH2D(
		"recoProton_length_vs_cosL",
		"Truth-Matched Reco Proton: Length vs cos(#theta);Track Length [cm];cos(#theta)",
		300, 0, 300,
		100, -1, 1
	);
	TH2D* recoProton_energy_vs_cosL = new TH2D(
		"recoProton_energy_vs_cosL",
		"Truth-Matched Reco Proton: Energy vs cos(#theta);Energy [GeV];cos(#theta)",
		100, 0, 2,
		100, -1, 1
	);
	
	TH2D* recoProton_cosL_vs_startz = new TH2D(
		"recoProton_cosL_vs_startz",
		"Truth-Matched Proton: cos(#theta) vs StartZ;cos(#theta);Start Z [cm]",
		100, -1, 1,
		70, -70, 70
	);

	TH2D* recoProton_cosL_vs_endz = new TH2D(
		"recoProton_cosL_vs_endz",
		"Truth-Matched Proton: cos(#theta) vs EndZ;cos(#theta);End Z [cm]",
		100, -1, 1,
//		70, -70, 70
		200, -100, 300    
	);	
	
	TH2D* recoProton_dX_vs_length = new TH2D(
		"recoProton_dX_vs_length",
		"Truth-Matched Reco Proton: #DeltaX vs Length;Track Length [cm];#DeltaX [cm]",
		100, 0, 100,
		100, -140, 140
	);
	TH2D* recoProton_dY_vs_length = new TH2D(
		"recoProton_dY_vs_length",
		"Truth-Matched Reco Proton: #DeltaY vs Length;Track Length [cm];#DeltaY [cm]",
		100, 0, 100,
		100, -140, 140
	);
	TH2D* recoProton_dZ_vs_length = new TH2D(
		"recoProton_dZ_vs_length",
		"Truth-Matched Reco Proton: #DeltaZ vs Length;Track Length [cm];#DeltaZ [cm]",
		100, 0, 100,
		100, -140, 140
	);
	TH2D* recoProton_endX_vs_endY = new TH2D(
		"recoProton_endX_vs_endY",
		"Truth-Matched Reco Proton: EndX vs EndY;End X [cm];End Y [cm]",
		200, -100, 300,
		200, -100, 300
	);

	TH2D* recoProton_endZ_vs_endX = new TH2D(
		"recoProton_endZ_vs_endX",
		"Truth-Matched Reco Proton: EndZ vs EndX;End Z [cm];End X [cm]",
		200, -100, 300,
		200, -100, 300
	);

	TH2D* recoProton_endZ_vs_endY = new TH2D(
		"recoProton_endZ_vs_endY",
		"Truth-Matched Reco Proton: EndZ vs EndY;End Z [cm];End Y [cm]",
		200, -100, 300,
		200, -100, 300
	);

	TH1D* recoProtonCosAlpha = new TH1D(
		"recoProtonCosAlpha",
		"Contained Matched Truth Protons;cos(#alpha);Number of Protons",
		50, -1.0, 1.0
	);

	TH2D* recoProtonCosAlpha_vs_E = new TH2D(
		"recoProtonCosAlpha_vs_E",
		"Contained Matched Truth Protons;Kinetic Energy [GeV];cos(#alpha)",
		50, 0.0, 2.0,
		50, -1.0, 1.0
	);

	TH2D* recoProtonCosAlpha_vs_length = new TH2D(
		"recoProtonCosAlpha_vs_length",
		"Contained Matched Truth Protons;Track Length [cm];cos(#alpha)",
		50, 0.0, 120.0,
		50, -1.0, 1.0
	);

	// Number of contained truth primary protons per interaction
	TH1D* primProtonMultiplicity = new TH1D(
				"primProtonMultiplicity", "Contained truth primary proton multiplicity;N_{prim p};Interactions",
				 20, 0, 20);


	TH1D* recoProtonLength_QE  = new TH1D("recoProtonLength_QE", "Reco Proton Length (QE)", 100, 0, 100);
	TH1D* recoProtonLength_RES = new TH1D("recoProtonLength_RES", "Reco Proton Length (RES)", 100, 0, 100);
	TH1D* recoProtonLength_DIS = new TH1D("recoProtonLength_DIS", "Reco Proton Length (DIS)", 100, 0, 100);
	TH1D* recoProtonLength_MEC = new TH1D("recoProtonLength_MEC", "Reco Proton Length (MEC)", 100, 0, 100);
	TH1D* recoProtonLength_COH = new TH1D("recoProtonLength_COH", "Reco Proton Length (COH)", 100, 0, 100);
	TH1D* recoProtonLength_Rock = new TH1D("recoProtonLength_Rock", "Reco Proton Length (Rock)", 100, 0, 100);
	TH1D* recoProtonLength_Sec = new TH1D("recoProtonLength_Sec", "Reco Proton Length (Secondary)", 100, 0, 100);


	// Reco Muon kinematics (all reco Muons, no truth match)
	TH1D* recoMuonLengthAll  = new TH1D("recoMuonLengthAll", "Reco Muon Track Length (All)", 100, 0, 100);
	TH1D* recoMuonCosLAll    = new TH1D("recoMuonCosLAll", "Reco Muon cos(#theta) (All)", 100, -1, 1);
	TH1D* recoMuonEAll       = new TH1D("recoMuonEAll", "Reco Muon Energy (All)", 100, 0, 2);
	TH1D* recoMuonStartXAll  = new TH1D("recoMuonStartXAll", "Reco Muon Start X (All)", 70, -70, 70);
	TH1D* recoMuonStartYAll  = new TH1D("recoMuonStartYAll", "Reco Muon Start Y (All)", 70, -70, 70);
	TH1D* recoMuonStartZAll  = new TH1D("recoMuonStartZAll", "Reco Muon Start Z (All)", 70, -70, 70);
	TH1D* recoMuonEndXAll  = new TH1D("recoMuonEndXAll", "Reco Muon End X (All)", 70, -70, 70);
	TH1D* recoMuonEndYAll  = new TH1D("recoMuonEndYAll", "Reco Muon End Y (All)", 70, -70, 70);
	TH1D* recoMuonEndZAll  = new TH1D("recoMuonEndZAll", "Reco Muon End Z (All)", 70, -70, 70);
	TH1D* recoMuonMomentumXAll  = new TH1D("recoMuonMomentumXAll", "Reco Muon Momentum X (All)", 100, 0, 2);
	TH1D* recoMuonMomentumYAll  = new TH1D("recoMuonMomentumYAll", "Reco Muon Momentum Y (All)", 100, 0, 2);
	TH1D* recoMuonMomentumZAll  = new TH1D("recoMuonMomentumZAll", "Reco Muon Momentum Z (All)", 100, 0, 2);
	TH2D* recoMuonLengthAll_vs_E = new TH2D(
		"recoMuonLengthAll_vs_E",
		"Reco Muon (All): Track Length vs Energy;Track Length [cm];Energy [GeV]",
		100, 0, 100,   // same binning as recoMuonLengthAll
		100, 0, 2      // same binning as recoMuonEAll
	);


	// Reco proton kinematics (all reco protons, no truth match)
	TH1D* recoProtonLengthAll  = new TH1D("recoProtonLengthAll", "Reco Proton Track Length (All)", 100, 0, 100);
	TH1D* recoProtonCosLAll    = new TH1D("recoProtonCosLAll", "Reco Proton cos(#theta) (All)", 100, -1, 1);
	TH1D* recoProtonEAll       = new TH1D("recoProtonEAll", "Reco Proton Energy (All)", 100, 0, 2);
	TH1D* recoProtonStartXAll  = new TH1D("recoProtonStartXAll", "Reco Proton Start X (All)", 70, -70, 70);
	TH1D* recoProtonStartYAll  = new TH1D("recoProtonStartYAll", "Reco Proton Start Y (All)", 70, -70, 70);
	TH1D* recoProtonStartZAll  = new TH1D("recoProtonStartZAll", "Reco Proton Start Z (All)", 70, -70, 70);
	TH1D* recoProtonEndXAll  = new TH1D("recoProtonEndXAll", "Reco Proton End X (All)", 70, -70, 70);
	TH1D* recoProtonEndYAll  = new TH1D("recoProtonEndYAll", "Reco Proton End Y (All)", 70, -70, 70);
	TH1D* recoProtonEndZAll  = new TH1D("recoProtonEndZAll", "Reco Proton End Z (All)", 70, -70, 70);
	TH1D* recoProtonMomentumXAll  = new TH1D("recoProtonMomentumXAll", "Reco Proton Momentum X (All)", 100, 0, 2);
	TH1D* recoProtonMomentumYAll  = new TH1D("recoProtonMomentumYAll", "Reco Proton Momentum Y (All)", 100, 0, 2);
	TH1D* recoProtonMomentumZAll  = new TH1D("recoProtonMomentumZAll", "Reco Proton Momentum Z (All)", 100, 0, 2);
	TH2D* recoProtonLength_vs_EAll = new TH2D(
		"recoProtonLength_vs_EAll",
		"Reco Proton (All): Track Length vs Energy;Track Length [cm];Energy [GeV]",
		100, 0, 100,   // same binning as recoProtonLengthAll
		100, 0, 2      // same binning as recoProtonEAll
	);
	TH2D* recoProton_startX_vs_startYAll = new TH2D(
		"recoProton_startX_vs_startYAll",
		"Reco Proton (All): StartX vs StartY;Start X [cm];Start Y [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* recoProton_startZ_vs_startXAll = new TH2D(
		"recoProton_startZ_vs_startXAll",
		"Reco Proton (All): StartZ vs StartX;Start Z [cm];Start X [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* recoProton_startZ_vs_startYAll = new TH2D(
		"recoProton_startZ_vs_startYAll",
		"Reco Proton (All): StartZ vs StartY;Start Z [cm];Start Y [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* recoProton_startX_vs_endXAll = new TH2D(
		"recoProton_startX_vs_endXAll",
		"Reco Proton (All): StartX vs EndX;Start X [cm];End X [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* recoProton_startY_vs_endYAll = new TH2D(
		"recoProton_startY_vs_endYAll",
		"Reco Proton (All): StartY vs EndY;Start Y [cm];End Y [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* recoProton_startZ_vs_endZAll = new TH2D(
		"recoProton_startZ_vs_endZAll",
		"Reco Proton (All): StartZ vs EndZ;Start Z [cm];End Z [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* recoProton_length_vs_cosLAll = new TH2D(
		"recoProton_length_vs_cosLAll",
		"Reco Proton (All): Length vs cos(#theta);Track Length [cm];cos(#theta)",
		100, 0, 100,
		100, -1, 1
	);
	TH2D* recoProton_energy_vs_cosLAll = new TH2D(
		"recoProton_energy_vs_cosLAll",
		"Reco Proton (All): Energy vs cos(#theta);Energy [GeV];cos(#theta)",
		100, 0, 2,
		100, -1, 1
	);
	TH2D* recoProton_dX_vs_lengthAll = new TH2D(
		"recoProton_dX_vs_lengthAll",
		"Reco Proton (All): #DeltaX vs Length;Track Length [cm];#DeltaX [cm]",
		100, 0, 100,
		100, -140, 140
	);
	TH2D* recoProton_dY_vs_lengthAll = new TH2D(
		"recoProton_dY_vs_lengthAll",
		"Reco Proton (All): #DeltaY vs Length;Track Length [cm];#DeltaY [cm]",
		100, 0, 100,
		100, -140, 140
	);
	TH2D* recoProton_dZ_vs_lengthAll = new TH2D(
		"recoProton_dZ_vs_lengthAll",
		"Reco Proton (All): #DeltaZ vs Length;Track Length [cm];#DeltaZ [cm]",
		100, 0, 100,
		100, -140, 140
	);
	TH2D* recoProton_endX_vs_endYAll = new TH2D(
		"recoProton_endX_vs_endYAll",
		"Reco Proton (All): EndX vs EndY;End X [cm];End Y [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* recoProton_endZ_vs_endXAll = new TH2D(
		"recoProton_endZ_vs_endXAll",
		"Reco Proton (All): EndZ vs EndX;End Z [cm];End X [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* recoProton_endZ_vs_endYAll = new TH2D(
		"recoProton_endZ_vs_endYAll",
		"Reco Proton (All): EndZ vs EndY;End Z [cm];End Y [cm]",
		70, -70, 70,
		70, -70, 70
	);



	// Reco pion kinematics (all reco pions, no truth match)
	TH1D* recoPionLengthAll  = new TH1D("recoPionLengthAll", "Reco Pion Track Length (All)", 100, 0, 100);
	TH1D* recoPionCosLAll    = new TH1D("recoPionCosLAll", "Reco Pion cos(#theta) (All)", 100, -1, 1);
	TH1D* recoPionEAll       = new TH1D("recoPionEAll", "Reco Pion Energy (All)", 100, 0, 2);
	TH1D* recoPionStartXAll  = new TH1D("recoPionStartXAll", "Reco Pion Start X (All)", 70, -70, 70);
	TH1D* recoPionStartYAll  = new TH1D("recoPionStartYAll", "Reco Pion Start Y (All)", 70, -70, 70);
	TH1D* recoPionStartZAll  = new TH1D("recoPionStartZAll", "Reco Pion Start Z (All)", 70, -70, 70);
	TH1D* recoPionEndXAll  = new TH1D("recoPionEndXAll", "Reco Pion End X (All)", 70, -70, 70);
	TH1D* recoPionEndYAll  = new TH1D("recoPionEndYAll", "Reco Pion End Y (All)", 70, -70, 70);
	TH1D* recoPionEndZAll  = new TH1D("recoPionEndZAll", "Reco Pion End Z (All)", 70, -70, 70);
	TH1D* recoPionMomentumXAll  = new TH1D("recoPionMomentumXAll", "Reco Pion Momentum X (All)", 100, 0, 2);
	TH1D* recoPionMomentumYAll  = new TH1D("recoPionMomentumYAll", "Reco Pion Momentum Y (All)", 100, 0, 2);
	TH1D* recoPionMomentumZAll  = new TH1D("recoPionMomentumZAll", "Reco Pion Momentum Z (All)", 100, 0, 2);
	TH2D* recoPionLength_vs_EAll = new TH2D(
		"recoPionLength_vs_EAll",
		"Reco Pion (All): Track Length vs Energy;Track Length [cm];Energy [GeV]",
		100, 0, 100,   // same binning as recoPionLengthAll
		100, 0, 2      // same binning as recoPionEAll
	);
	TH2D* recoPion_startX_vs_startYAll = new TH2D(
		"recoPion_startX_vs_startYAll",
		"Reco Pion (All): StartX vs StartY;Start X [cm];Start Y [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* recoPion_startZ_vs_startXAll = new TH2D(
		"recoPion_startZ_vs_startXAll",
		"Reco Pion (All): StartZ vs StartX;Start Z [cm];Start X [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* recoPion_startZ_vs_startYAll = new TH2D(
		"recoPion_startZ_vs_startYAll",
		"Reco Pion (All): StartZ vs StartY;Start Z [cm];Start Y [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* recoPion_startX_vs_endXAll = new TH2D(
		"recoPion_startX_vs_endXAll",
		"Reco Pion (All): StartX vs EndX;Start X [cm];End X [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* recoPion_startY_vs_endYAll = new TH2D(
		"recoPion_startY_vs_endYAll",
		"Reco Pion (All): StartY vs EndY;Start Y [cm];End Y [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* recoPion_startZ_vs_endZAll = new TH2D(
		"recoPion_startZ_vs_endZAll",
		"Reco Pion (All): StartZ vs EndZ;Start Z [cm];End Z [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* recoPion_length_vs_cosLAll = new TH2D(
		"recoPion_length_vs_cosLAll",
		"Reco Pion (All): Length vs cos(#theta);Track Length [cm];cos(#theta)",
		100, 0, 100,
		100, -1, 1
	);
	TH2D* recoPion_energy_vs_cosLAll = new TH2D(
		"recoPion_energy_vs_cosLAll",
		"Reco Pion (All): Energy vs cos(#theta);Energy [GeV];cos(#theta)",
		100, 0, 2,
		100, -1, 1
	);
	TH2D* recoPion_dX_vs_lengthAll = new TH2D(
		"recoPion_dX_vs_lengthAll",
		"Reco Pion (All): #DeltaX vs Length;Track Length [cm];#DeltaX [cm]",
		100, 0, 100,
		100, -140, 140
	);
	TH2D* recoPion_dY_vs_lengthAll = new TH2D(
		"recoPion_dY_vs_lengthAll",
		"Reco Pion (All): #DeltaY vs Length;Track Length [cm];#DeltaY [cm]",
		100, 0, 100,
		100, -140, 140
	);
	TH2D* recoPion_dZ_vs_lengthAll = new TH2D(
		"recoPion_dZ_vs_lengthAll",
		"Reco Pion (All): #DeltaZ vs Length;Track Length [cm];#DeltaZ [cm]",
		100, 0, 100,
		100, -140, 140
	);
	TH2D* recoPion_endX_vs_endYAll = new TH2D(
		"recoPion_endX_vs_endYAll",
		"Reco Pion (All): EndX vs EndY;End X [cm];End Y [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* recoPion_endZ_vs_endXAll = new TH2D(
		"recoPion_endZ_vs_endXAll",
		"Reco Pion (All): EndZ vs EndX;End Z [cm];End X [cm]",
		70, -70, 70,
		70, -70, 70
	);
	TH2D* recoPion_endZ_vs_endYAll = new TH2D(
		"recoPion_endZ_vs_endYAll",
		"Reco Pion (All): EndZ vs EndY;End Z [cm];End Y [cm]",
		70, -70, 70,
		70, -70, 70
	);


	// TH1D* recoPionLengthAll  = new TH1D("recoPionLengthAll", "Reco Pion Track Length (All)", 100, 0, 100);
	// TH1D* recoPionCosLAll    = new TH1D("recoPionCosLAll", "Reco Pion cos(#theta) (All)", 100, -1, 1);
	// TH1D* recoPionEAll       = new TH1D("recoPionEAll", "Reco Pion Energy (All)", 100, 0, 2);
	// TH1D* recoPionStartXAll  = new TH1D("recoPionStartXAll", "Reco Pion Start X (All)", 70, -70, 70);
	// TH1D* recoPionStartYAll  = new TH1D("recoPionStartYAll", "Reco Pion Start Y (All)", 70, -70, 70);
	// TH1D* recoPionStartZAll  = new TH1D("recoPionStartZAll", "Reco Pion Start Z (All)", 70, -70, 70);
	// TH1D* recoPionEndXAll  = new TH1D("recoPionEndXAll", "Reco Pion End X (All)", 70, -70, 70);
	// TH1D* recoPionEndYAll  = new TH1D("recoPionEndYAll", "Reco Pion End Y (All)", 70, -70, 70);
	// TH1D* recoPionEndZAll  = new TH1D("recoPionEndZAll", "Reco Pion End Z (All)", 70, -70, 70);
	// TH1D* recoPionMomentumXAll  = new TH1D("recoPionMomentumXAll", "Reco Pion Momentum X (All)", 100, 0, 2);
	// TH1D* recoPionMomentumYAll  = new TH1D("recoPionMomentumYAll", "Reco Pion Momentum Y (All)", 100, 0, 2);
	// TH1D* recoPionMomentumZAll  = new TH1D("recoPionMomentumZAll", "Reco Pion Momentum Z (All)", 100, 0, 2);
	// TH2D* recoPionLengthAll_vs_E = new TH2D(
		// "recoPionLengthAll_vs_E",
		// "Reco Pion (All): Track Length vs Energy;Track Length [cm];Energy [GeV]",
		// 100, 0, 100,   // same binning as recoPionLengthAll
		// 100, 0, 2      // same binning as recoPionEAll
	// );


	// Reco Kaon kinematics (all reco Kaons, no truth match)
	TH1D* recoKaonLengthAll  = new TH1D("recoKaonLengthAll", "Reco Kaon Track Length (All)", 100, 0, 100);
	TH1D* recoKaonCosLAll    = new TH1D("recoKaonCosLAll", "Reco Kaon cos(#theta) (All)", 100, -1, 1);
	TH1D* recoKaonEAll       = new TH1D("recoKaonEAll", "Reco Kaon Energy (All)", 100, 0, 2);
	TH1D* recoKaonStartXAll  = new TH1D("recoKaonStartXAll", "Reco Kaon Start X (All)", 70, -70, 70);
	TH1D* recoKaonStartYAll  = new TH1D("recoKaonStartYAll", "Reco Kaon Start Y (All)", 70, -70, 70);
	TH1D* recoKaonStartZAll  = new TH1D("recoKaonStartZAll", "Reco Kaon Start Z (All)", 70, -70, 70);
	TH1D* recoKaonEndXAll  = new TH1D("recoKaonEndXAll", "Reco Kaon End X (All)", 70, -70, 70);
	TH1D* recoKaonEndYAll  = new TH1D("recoKaonEndYAll", "Reco Kaon End Y (All)", 70, -70, 70);
	TH1D* recoKaonEndZAll  = new TH1D("recoKaonEndZAll", "Reco Kaon End Z (All)", 70, -70, 70);
	TH1D* recoKaonMomentumXAll  = new TH1D("recoKaonMomentumXAll", "Reco Kaon Momentum X (All)", 100, 0, 2);
	TH1D* recoKaonMomentumYAll  = new TH1D("recoKaonMomentumYAll", "Reco Kaon Momentum Y (All)", 100, 0, 2);
	TH1D* recoKaonMomentumZAll  = new TH1D("recoKaonMomentumZAll", "Reco Kaon Momentum Z (All)", 100, 0, 2);


// Reco (Matched)

TH1D* recoVtxX_matched = new TH1D("recoVtxX_matched", "Reco Nu Vertex X (Matched)", 70, -70, 70);
TH1D* recoVtxY_matched = new TH1D("recoVtxY_matched", "Reco Nu Vertex Y (Matched)", 70, -70, 70);
TH1D* recoVtxZ_matched = new TH1D("recoVtxZ_matched", "Reco Nu Vertex Z (Matched)", 70, -70, 70);

// ===========================================================
// =============== Reco (Matched) Proton Histograms ==========
// ===========================================================

// 1D distributions
TH1D* reco_matched_proton_length  = new TH1D("reco_matched_proton_length",  "Reco-Matched Proton Track Length", 300, 0, 300);
TH1D* reco_matched_proton_cosL    = new TH1D("reco_matched_proton_cosL",    "Reco-Matched Proton cos(#theta)", 100, -1, 1);
TH1D* reco_matched_proton_energy  = new TH1D("reco_matched_proton_energy",  "Reco-Matched Proton Energy", 100, 0, 2);

TH1D* reco_matched_proton_startX  = new TH1D("reco_matched_proton_startX", "Reco-Matched Proton Start X;Start X [cm];Events", 70, -70, 70);
TH1D* reco_matched_proton_startY  = new TH1D("reco_matched_proton_startY", "Reco-Matched Proton Start Y;Start Y [cm];Events", 70, -70, 70);
TH1D* reco_matched_proton_startZ  = new TH1D("reco_matched_proton_startZ", "Reco-Matched Proton Start Z;Start Z [cm];Events", 70, -70, 70);

TH1D* reco_matched_proton_endX    = new TH1D("reco_matched_proton_endX",   "Reco-Matched Proton End X;End X [cm];Events", 200, -100, 300); // 200, -100, 300
TH1D* reco_matched_proton_endY    = new TH1D("reco_matched_proton_endY",   "Reco-Matched Proton End Y;End Y [cm];Events", 200, -100, 300);
TH1D* reco_matched_proton_endZ    = new TH1D("reco_matched_proton_endZ",   "Reco-Matched Proton End Z;End Z [cm];Events", 200, -100, 300);

// Interaction-type breakdown of reco-matched proton length and cosL
TH1D* reco_matched_proton_length_QE  = new TH1D("reco_matched_proton_length_QE",  "Reco-Matched Proton Length (QE);Track Length [cm];Events",  300, 0, 300);
TH1D* reco_matched_proton_length_MEC = new TH1D("reco_matched_proton_length_MEC", "Reco-Matched Proton Length (MEC);Track Length [cm];Events", 300, 0, 300);
TH1D* reco_matched_proton_length_RES = new TH1D("reco_matched_proton_length_RES", "Reco-Matched Proton Length (RES);Track Length [cm];Events", 300, 0, 300);
TH1D* reco_matched_proton_length_DIS = new TH1D("reco_matched_proton_length_DIS", "Reco-Matched Proton Length (DIS);Track Length [cm];Events", 300, 0, 300);
TH1D* reco_matched_proton_length_COH = new TH1D("reco_matched_proton_length_COH", "Reco-Matched Proton Length (COH);Track Length [cm];Events", 300, 0, 300);
TH1D* reco_matched_proton_cosL_QE    = new TH1D("reco_matched_proton_cosL_QE",    "Reco-Matched Proton cos(#theta) (QE);cos(#theta);Events",   100, -1, 1);
TH1D* reco_matched_proton_cosL_MEC   = new TH1D("reco_matched_proton_cosL_MEC",   "Reco-Matched Proton cos(#theta) (MEC);cos(#theta);Events",  100, -1, 1);
TH1D* reco_matched_proton_cosL_RES   = new TH1D("reco_matched_proton_cosL_RES",   "Reco-Matched Proton cos(#theta) (RES);cos(#theta);Events",  100, -1, 1);
TH1D* reco_matched_proton_cosL_DIS   = new TH1D("reco_matched_proton_cosL_DIS",   "Reco-Matched Proton cos(#theta) (DIS);cos(#theta);Events",  100, -1, 1);
TH1D* reco_matched_proton_cosL_COH   = new TH1D("reco_matched_proton_cosL_COH",   "Reco-Matched Proton cos(#theta) (COH);cos(#theta);Events",  100, -1, 1);

// Particle-type breakdown
TH1D* reco_matched_proton_length_trueMuon   = new TH1D("reco_matched_proton_length_trueMuon",   "Reco-Matched Proton Length (True #mu);Track Length [cm];Events",  300, 0, 300);
TH1D* reco_matched_proton_length_trueProton = new TH1D("reco_matched_proton_length_trueProton", "Reco-Matched Proton Length (True p);Track Length [cm];Events",    300, 0, 300);
TH1D* reco_matched_proton_length_truePion   = new TH1D("reco_matched_proton_length_truePion",   "Reco-Matched Proton Length (True #pi);Track Length [cm];Events",  300, 0, 300);
TH1D* reco_matched_proton_length_trueOther  = new TH1D("reco_matched_proton_length_trueOther",  "Reco-Matched Proton Length (True Other);Track Length [cm];Events",300, 0, 300);
TH1D* reco_matched_proton_cosL_trueMuon     = new TH1D("reco_matched_proton_cosL_trueMuon",     "Reco-Matched Proton cos(#theta) (True #mu);cos(#theta);Events",   100, -1, 1);
TH1D* reco_matched_proton_cosL_trueProton   = new TH1D("reco_matched_proton_cosL_trueProton",   "Reco-Matched Proton cos(#theta) (True p);cos(#theta);Events",     100, -1, 1);
TH1D* reco_matched_proton_cosL_truePion     = new TH1D("reco_matched_proton_cosL_truePion",     "Reco-Matched Proton cos(#theta) (True #pi);cos(#theta);Events",   100, -1, 1);
TH1D* reco_matched_proton_cosL_trueOther    = new TH1D("reco_matched_proton_cosL_trueOther",    "Reco-Matched Proton cos(#theta) (True Other);cos(#theta);Events", 100, -1, 1);

TH1D* reco_matched_proton_momentumX = new TH1D("reco_matched_proton_momentumX",
    "Reco-Matched Proton Momentum X;P_{x} [GeV/c];Events", 100, 0, 2);
TH1D* reco_matched_proton_momentumY = new TH1D("reco_matched_proton_momentumY",
    "Reco-Matched Proton Momentum Y;P_{y} [GeV/c];Events", 100, 0, 2);
TH1D* reco_matched_proton_momentumZ = new TH1D("reco_matched_proton_momentumZ",
    "Reco-Matched Proton Momentum Z;P_{z} [GeV/c];Events", 100, 0, 2);


// 2D distributions

TH2D* reco_matched_proton_length_vs_energy = new TH2D(
    "reco_matched_proton_length_vs_energy",
    "Reco-Matched Proton: Track Length vs Energy;Track Length [cm];Energy [GeV]",
    300, 0, 300,
    100, 0, 2
);

TH2D* reco_matched_proton_startX_vs_startY = new TH2D(
    "reco_matched_proton_startX_vs_startY",
    "Reco-Matched Proton: StartX vs StartY;Start X [cm];Start Y [cm]",
    70, -70, 70,
    70, -70, 70
);

TH2D* reco_matched_proton_startZ_vs_startX = new TH2D(
    "reco_matched_proton_startZ_vs_startX",
    "Reco-Matched Proton: StartZ vs StartX;Start Z [cm];Start X [cm]",
    70, -70, 70,
    70, -70, 70
);

TH2D* reco_matched_proton_startZ_vs_startY = new TH2D(
    "reco_matched_proton_startZ_vs_startY",
    "Reco-Matched Proton: StartZ vs StartY;Start Z [cm];Start Y [cm]",
    70, -70, 70,
    70, -70, 70
);

TH2D* reco_matched_proton_startX_vs_endX = new TH2D(
    "reco_matched_proton_startX_vs_endX",
    "Reco-Matched Proton: StartX vs EndX;Start X [cm];End X [cm]",
    70, -70, 70,
    200, -100, 300
);

TH2D* reco_matched_proton_startY_vs_endY = new TH2D(
    "reco_matched_proton_startY_vs_endY",
    "Reco-Matched Proton: StartY vs EndY;Start Y [cm];End Y [cm]",
    70, -70, 70,
    200, -100, 300
);

TH2D* reco_matched_proton_startZ_vs_endZ = new TH2D(
    "reco_matched_proton_startZ_vs_endZ",
    "Reco-Matched Proton: StartZ vs EndZ;Start Z [cm];End Z [cm]",
    70, -70, 70,
    200, -100, 300
);

TH2D* reco_matched_proton_length_vs_cosL = new TH2D(
    "reco_matched_proton_length_vs_cosL",
    "Reco-Matched Proton: Length vs cos(#theta);Track Length [cm];cos(#theta)",
    300, 0, 300,
    100, -1, 1
);

TH2D* reco_matched_proton_cosL_vs_startZ = new TH2D(
    "reco_matched_proton_cosL_vs_startZ",
    "Reco-Matched Proton: cos(#theta) vs StartZ;cos(#theta);Start Z [cm]",
    100, -1, 1,
    70, -70, 70
);

TH2D* reco_matched_proton_cosL_vs_endZ = new TH2D(
    "reco_matched_proton_cosL_vs_endZ",
    "Reco-Matched Proton: cos(#theta) vs EndZ;cos(#theta);End Z [cm]",
    100, -1, 1,
    200, -100, 300
);

TH2D* reco_matched_proton_energy_vs_cosL = new TH2D(
    "reco_matched_proton_energy_vs_cosL",
    "Reco-Matched Proton: Energy vs cos(#theta);Energy [GeV];cos(#theta)",
    100, 0, 2,
    100, -1, 1
);

TH2D* reco_matched_proton_dX_vs_length = new TH2D(
    "reco_matched_proton_dX_vs_length",
    "Reco-Matched Proton: #DeltaX vs Length;Track Length [cm];#DeltaX [cm]",
    100, 0, 100,
    100, -140, 140
);

TH2D* reco_matched_proton_dY_vs_length = new TH2D(
    "reco_matched_proton_dY_vs_length",
    "Reco-Matched Proton: #DeltaY vs Length;Track Length [cm];#DeltaY [cm]",
    100, 0, 100,
    100, -140, 140
);

TH2D* reco_matched_proton_dZ_vs_length = new TH2D(
    "reco_matched_proton_dZ_vs_length",
    "Reco-Matched Proton: #DeltaZ vs Length;Track Length [cm];#DeltaZ [cm]",
    100, 0, 100,
    100, -140, 140
);

TH2D* reco_matched_proton_endX_vs_endY = new TH2D(
    "reco_matched_proton_endX_vs_endY",
    "Reco-Matched Proton: EndX vs EndY;End X [cm];End Y [cm]",
    200, -100, 300,
    200, -100, 300
);

TH2D* reco_matched_proton_endZ_vs_endX = new TH2D(
    "reco_matched_proton_endZ_vs_endX",
    "Reco-Matched Proton: EndZ vs EndX;End Z [cm];End X [cm]",
    200, -100, 300,
    200, -100, 300
);

TH2D* reco_matched_proton_endZ_vs_endY = new TH2D(
    "reco_matched_proton_endZ_vs_endY",
    "Reco-Matched Proton: EndZ vs EndY;End Z [cm];End Y [cm]",
    200, -100, 300,
    200, -100, 300
);

// Reco proton multiplicity per matched reco interaction
TH1D* recoProtonMultiplicity_matched =
	new TH1D(
		"recoProtonMultiplicity_matched",
		"Reco Proton Multiplicity (Matched Reco Interactions)",
		20, 0, 20
	);


// ===========================================================
// =============== Reco (Matched) Pion Histograms ==========
// ===========================================================

// 1D distributions
TH1D* reco_matched_pion_length  = new TH1D("reco_matched_pion_length",  "Reco-Matched Pion Track Length", 300, 0, 300); //100, 0, 100
TH1D* reco_matched_pion_cosL    = new TH1D("reco_matched_pion_cosL",    "Reco-Matched Pion cos(#theta)", 100, -1, 1);
TH1D* reco_matched_pion_energy  = new TH1D("reco_matched_pion_energy",  "Reco-Matched Pion Energy", 100, 0, 2);

// Reco pion kinetic energy from reconstructed momentum
TH1D* reco_matched_pion_energy_fromP =
	new TH1D("reco_matched_pion_energy_fromP",
	         "Reco matched pion KE from momentum;T_{#pi}^{kin} [GeV];Entries",
	         100, 0, 2);
TH2D* reco_matched_pion_length_vs_energy_fromP =
	new TH2D("reco_matched_pion_length_vs_energy_fromP",
	         "Reco matched pion length vs KE from momentum;Length [cm];T_{#pi}^{kin} [GeV]",
	         100, 0, 100,
	         100, 0, 2);

// (Optional but highly recommended diagnostic)
TH1D* reco_matched_pion_energy_diff =
	new TH1D("reco_matched_pion_energy_diff",
	         "Reco pion energy: stored - momentum-based;#DeltaT [GeV];Entries",
	         100, -1.0, 1.0);


TH1D* reco_matched_pion_startX  = new TH1D("reco_matched_pion_startX", "Reco-Matched Pion Start X;Start X [cm];Events", 70, -70, 70);
TH1D* reco_matched_pion_startY  = new TH1D("reco_matched_pion_startY", "Reco-Matched Pion Start Y;Start Y [cm];Events", 70, -70, 70);
TH1D* reco_matched_pion_startZ  = new TH1D("reco_matched_pion_startZ", "Reco-Matched Pion Start Z;Start Z [cm];Events", 70, -70, 70);

TH1D* reco_matched_pion_endX    = new TH1D("reco_matched_pion_endX",   "Reco-Matched Pion End X;End X [cm];Events", 70, -70, 70);
TH1D* reco_matched_pion_endY    = new TH1D("reco_matched_pion_endY",   "Reco-Matched Pion End Y;End Y [cm];Events", 70, -70, 70);
TH1D* reco_matched_pion_endZ    = new TH1D("reco_matched_pion_endZ",   "Reco-Matched Pion End Z;End Z [cm];Events", 70, -70, 70);

TH1D* reco_matched_pion_momentumX = new TH1D("reco_matched_pion_momentumX",
    "Reco-Matched Pion Momentum X;P_{x} [GeV/c];Events", 100, 0, 2);
TH1D* reco_matched_pion_momentumY = new TH1D("reco_matched_pion_momentumY",
    "Reco-Matched Pion Momentum Y;P_{y} [GeV/c];Events", 100, 0, 2);
TH1D* reco_matched_pion_momentumZ = new TH1D("reco_matched_pion_momentumZ",
    "Reco-Matched Pion Momentum Z;P_{z} [GeV/c];Events", 100, 0, 2);


// 2D distributions

TH2D* reco_matched_pion_length_vs_energy = new TH2D(
    "reco_matched_pion_length_vs_energy",
    "Reco-Matched Pion: Track Length vs Energy;Track Length [cm];Energy [GeV]",
    100, 0, 100,
    100, 0, 2
);

TH2D* reco_matched_pion_startX_vs_startY = new TH2D(
    "reco_matched_pion_startX_vs_startY",
    "Reco-Matched Pion: StartX vs StartY;Start X [cm];Start Y [cm]",
    70, -70, 70,
    70, -70, 70
);

TH2D* reco_matched_pion_startZ_vs_startX = new TH2D(
    "reco_matched_pion_startZ_vs_startX",
    "Reco-Matched Pion: StartZ vs StartX;Start Z [cm];Start X [cm]",
    70, -70, 70,
    70, -70, 70
);

TH2D* reco_matched_pion_startZ_vs_startY = new TH2D(
    "reco_matched_pion_startZ_vs_startY",
    "Reco-Matched Pion: StartZ vs StartY;Start Z [cm];Start Y [cm]",
    70, -70, 70,
    70, -70, 70
);

TH2D* reco_matched_pion_startX_vs_endX = new TH2D(
    "reco_matched_pion_startX_vs_endX",
    "Reco-Matched Pion: StartX vs EndX;Start X [cm];End X [cm]",
    70, -70, 70,
    70, -70, 70
);

TH2D* reco_matched_pion_startY_vs_endY = new TH2D(
    "reco_matched_pion_startY_vs_endY",
    "Reco-Matched Pion: StartY vs EndY;Start Y [cm];End Y [cm]",
    70, -70, 70,
    70, -70, 70
);

TH2D* reco_matched_pion_startZ_vs_endZ = new TH2D(
    "reco_matched_pion_startZ_vs_endZ",
    "Reco-Matched Pion: StartZ vs EndZ;Start Z [cm];End Z [cm]",
    70, -70, 70,
    70, -70, 70
);

TH2D* reco_matched_pion_length_vs_cosL = new TH2D(
    "reco_matched_pion_length_vs_cosL",
    "Reco-Matched Pion: Length vs cos(#theta);Track Length [cm];cos(#theta)",
    100, 0, 100,
    100, -1, 1
);

TH2D* reco_matched_pion_cosL_vs_startZ = new TH2D(
    "reco_matched_pion_cosL_vs_startZ",
    "Reco-Matched Pion: cos(#theta) vs StartZ;cos(#theta);Start Z [cm]",
    100, -1, 1,
    70, -70, 70
);

TH2D* reco_matched_pion_cosL_vs_endZ = new TH2D(
    "reco_matched_pion_cosL_vs_endZ",
    "Reco-Matched Pion: cos(#theta) vs EndZ;cos(#theta);End Z [cm]",
    100, -1, 1,
    70, -70, 70
);

TH2D* reco_matched_pion_energy_vs_cosL = new TH2D(
    "reco_matched_pion_energy_vs_cosL",
    "Reco-Matched Pion: Energy vs cos(#theta);Energy [GeV];cos(#theta)",
    100, 0, 2,
    100, -1, 1
);

TH2D* reco_matched_pion_dX_vs_length = new TH2D(
    "reco_matched_pion_dX_vs_length",
    "Reco-Matched Pion: #DeltaX vs Length;Track Length [cm];#DeltaX [cm]",
    100, 0, 100,
    100, -140, 140
);

TH2D* reco_matched_pion_dY_vs_length = new TH2D(
    "reco_matched_pion_dY_vs_length",
    "Reco-Matched Pion: #DeltaY vs Length;Track Length [cm];#DeltaY [cm]",
    100, 0, 100,
    100, -140, 140
);

TH2D* reco_matched_pion_dZ_vs_length = new TH2D(
    "reco_matched_pion_dZ_vs_length",
    "Reco-Matched Pion: #DeltaZ vs Length;Track Length [cm];#DeltaZ [cm]",
    100, 0, 100,
    100, -140, 140
);

TH2D* reco_matched_pion_endX_vs_endY = new TH2D(
    "reco_matched_pion_endX_vs_endY",
    "Reco-Matched Pion: EndX vs EndY;End X [cm];End Y [cm]",
    70, -70, 70,
    70, -70, 70
);

TH2D* reco_matched_pion_endZ_vs_endX = new TH2D(
    "reco_matched_pion_endZ_vs_endX",
    "Reco-Matched Pion: EndZ vs EndX;End Z [cm];End X [cm]",
    70, -70, 70,
    70, -70, 70
);

TH2D* reco_matched_pion_endZ_vs_endY = new TH2D(
    "reco_matched_pion_endZ_vs_endY",
    "Reco-Matched Pion: EndZ vs EndY;End Z [cm];End Y [cm]",
    70, -70, 70,
    70, -70, 70
);

// Reco Pion multiplicity per matched reco interaction
TH1D* recoPionMultiplicity_matched =
	new TH1D(
		"recoPionMultiplicity_matched",
		"Reco Pion Multiplicity (Matched Reco Interactions)",
		20, 0, 20
	);



	// Track multiplicity for events failing MINERvA match but with one long muon
	TH1D* track_multAltCut       = new TH1D("track_multAltCut", "Track multiplicity (Alt cut: long muon + X)", 20, 0, 20);
	TH1D* track_multAltCutQE     = new TH1D("track_multAltCutQE", "Track multiplicity (Alt Cut) - QE", 20, 0, 20);
	TH1D* track_multAltCutMEC    = new TH1D("track_multAltCutMEC", "Track multiplicity (Alt Cut) - MEC", 20, 0, 20);
	TH1D* track_multAltCutRES    = new TH1D("track_multAltCutRES", "Track multiplicity (Alt Cut) - RES", 20, 0, 20);
	TH1D* track_multAltCutDIS    = new TH1D("track_multAltCutDIS", "Track multiplicity (Alt Cut) - DIS", 20, 0, 20);
	TH1D* track_multAltCutCOH    = new TH1D("track_multAltCutCOH", "Track multiplicity (Alt Cut) - COH", 20, 0, 20);
	TH1D* track_multAltCutNC     = new TH1D("track_multAltCutNC", "Track multiplicity (Alt Cut) - NC", 20, 0, 20);
	TH1D* track_multAltCutRock      = new TH1D("track_multAltCutRock", "track_multAltCutRock", 20, 0, 20); 			
	TH1D* track_multAltCutSec       = new TH1D("track_multAltCutSec", "track_multAltCutSec", 20, 0, 20); 			

	// Muon track length for events failing MINERvA match but with one long muon
	TH1D* muon_lenAltCutQE   = new TH1D("muon_lenAltCutQE",   ";Muon Length (cm);Interactions", 100, 0, 180);
	TH1D* muon_lenAltCutMEC  = new TH1D("muon_lenAltCutMEC",  ";Muon Length (cm);Interactions", 100, 0, 180);
	TH1D* muon_lenAltCutRES  = new TH1D("muon_lenAltCutRES",  ";Muon Length (cm);Interactions", 100, 0, 180);
	TH1D* muon_lenAltCutDIS  = new TH1D("muon_lenAltCutDIS",  ";Muon Length (cm);Interactions", 100, 0, 180);
	TH1D* muon_lenAltCutCOH  = new TH1D("muon_lenAltCutCOH",  ";Muon Length (cm);Interactions", 100, 0, 180);
	TH1D* muon_lenAltCutNC   = new TH1D("muon_lenAltCutNC",   ";Muon Length (cm);Interactions", 100, 0, 180);
	TH1D* muon_lenAltCutRock = new TH1D("muon_lenAltCutRock", ";Muon Length (cm);Interactions", 100, 0, 180);
	TH1D* muon_lenAltCutSec  = new TH1D("muon_lenAltCutSec",  ";Muon Length (cm);Interactions", 100, 0, 180);

	// Muon cosine for events failing MINERvA match but with one long muon
	TH1D* muon_cosAltCutQE   = new TH1D("muon_cosAltCutQE",   ";Muon cos(#theta);Interactions", 40, -1.0, 1.0);
	TH1D* muon_cosAltCutMEC  = new TH1D("muon_cosAltCutMEC",  ";Muon cos(#theta);Interactions", 40, -1.0, 1.0);
	TH1D* muon_cosAltCutRES  = new TH1D("muon_cosAltCutRES",  ";Muon cos(#theta);Interactions", 40, -1.0, 1.0);
	TH1D* muon_cosAltCutDIS  = new TH1D("muon_cosAltCutDIS",  ";Muon cos(#theta);Interactions", 40, -1.0, 1.0);
	TH1D* muon_cosAltCutCOH  = new TH1D("muon_cosAltCutCOH",  ";Muon cos(#theta);Interactions", 40, -1.0, 1.0);
	TH1D* muon_cosAltCutNC   = new TH1D("muon_cosAltCutNC",   ";Muon cos(#theta);Interactions", 40, -1.0, 1.0);
	TH1D* muon_cosAltCutRock = new TH1D("muon_cosAltCutRock", ";Muon cos(#theta);Interactions", 40, -1.0, 1.0);
	TH1D* muon_cosAltCutSec  = new TH1D("muon_cosAltCutSec",  ";Muon cos(#theta);Interactions", 40, -1.0, 1.0);


	// Muon start X position for events failing MINERvA match but with one long muon
	TH1D* muon_startXAltCutQE   = new TH1D("muon_startXAltCutQE",   ";Muon Start X [cm];Interactions", 70, -70, 70);
	TH1D* muon_startXAltCutMEC  = new TH1D("muon_startXAltCutMEC",  ";Muon Start X [cm];Interactions", 70, -70, 70);
	TH1D* muon_startXAltCutRES  = new TH1D("muon_startXAltCutRES",  ";Muon Start X [cm];Interactions", 70, -70, 70);
	TH1D* muon_startXAltCutDIS  = new TH1D("muon_startXAltCutDIS",  ";Muon Start X [cm];Interactions", 70, -70, 70);
	TH1D* muon_startXAltCutCOH  = new TH1D("muon_startXAltCutCOH",  ";Muon Start X [cm];Interactions", 70, -70, 70);
	TH1D* muon_startXAltCutNC   = new TH1D("muon_startXAltCutNC",   ";Muon Start X [cm];Interactions", 70, -70, 70);
	TH1D* muon_startXAltCutRock = new TH1D("muon_startXAltCutRock", ";Muon Start X [cm];Interactions", 70, -70, 70);
	TH1D* muon_startXAltCutSec  = new TH1D("muon_startXAltCutSec",  ";Muon Start X [cm];Interactions", 70, -70, 70);

	// Muon start Y position for events failing MINERvA match but with one long muon
	TH1D* muon_startYAltCutQE   = new TH1D("muon_startYAltCutQE",   ";Muon Start Y [cm];Interactions", 70, -70, 70);
	TH1D* muon_startYAltCutMEC  = new TH1D("muon_startYAltCutMEC",  ";Muon Start Y [cm];Interactions", 70, -70, 70);
	TH1D* muon_startYAltCutRES  = new TH1D("muon_startYAltCutRES",  ";Muon Start Y [cm];Interactions", 70, -70, 70);
	TH1D* muon_startYAltCutDIS  = new TH1D("muon_startYAltCutDIS",  ";Muon Start Y [cm];Interactions", 70, -70, 70);
	TH1D* muon_startYAltCutCOH  = new TH1D("muon_startYAltCutCOH",  ";Muon Start Y [cm];Interactions", 70, -70, 70);
	TH1D* muon_startYAltCutNC   = new TH1D("muon_startYAltCutNC",   ";Muon Start Y [cm];Interactions", 70, -70, 70);
	TH1D* muon_startYAltCutRock = new TH1D("muon_startYAltCutRock", ";Muon Start Y [cm];Interactions", 70, -70, 70);
	TH1D* muon_startYAltCutSec  = new TH1D("muon_startYAltCutSec",  ";Muon Start Y [cm];Interactions", 70, -70, 70);

	// Muon start Z position for events failing MINERvA match but with one long muon
	TH1D* muon_startZAltCutQE   = new TH1D("muon_startZAltCutQE",   ";Muon Start Z [cm];Interactions", 70, -70, 70);
	TH1D* muon_startZAltCutMEC  = new TH1D("muon_startZAltCutMEC",  ";Muon Start Z [cm];Interactions", 70, -70, 70);
	TH1D* muon_startZAltCutRES  = new TH1D("muon_startZAltCutRES",  ";Muon Start Z [cm];Interactions", 70, -70, 70);
	TH1D* muon_startZAltCutDIS  = new TH1D("muon_startZAltCutDIS",  ";Muon Start Z [cm];Interactions", 70, -70, 70);
	TH1D* muon_startZAltCutCOH  = new TH1D("muon_startZAltCutCOH",  ";Muon Start Z [cm];Interactions", 70, -70, 70);
	TH1D* muon_startZAltCutNC   = new TH1D("muon_startZAltCutNC",   ";Muon Start Z [cm];Interactions", 70, -70, 70);
	TH1D* muon_startZAltCutRock = new TH1D("muon_startZAltCutRock", ";Muon Start Z [cm];Interactions", 70, -70, 70);
	TH1D* muon_startZAltCutSec  = new TH1D("muon_startZAltCutSec",  ";Muon Start Z [cm];Interactions", 70, -70, 70);


// Histograms for 1 muon + 1 proton
TH1F* true_nuE_muP = new TH1F("true_nuE_muP", "True Neutrino Energy (1#mu+1p);E_{#nu} [GeV];Events", 20, 0, 20);
TH1F* true_mult_muP = new TH1F("true_mult_muP", "True multiplicity (1#mu+1p only);Number of True Tracks", 10, 0, 10);
TH2F* h2_true_nuE_vs_muEndZ = new TH2F("h2_true_nuE_vs_muEndZ", "True #nu Energy vs Muon End Z (1#mu+1p);True E_{#nu} [GeV];Muon End Z [cm]", 20, 0, 20, 100, 0, 16000);


TH1F* true_mu_len_muP = new TH1F("true_mu_len_muP", "Muon Length (1#mu+1p only);Length [cm]", 2000, 0, 2000);
TH1F* true_mu_cos_muP = new TH1F("true_mu_cos_muP", "Muon cos(#theta_{z}) (1#mu+1p only);cos(#theta_{z})", 40, -1, 1);
TH1F* true_mu_startX_muP = new TH1F("true_mu_startX_muP", "Muon start X (1#mu+1p only);X [cm]", 70, -70, 70);
TH1F* true_mu_startY_muP = new TH1F("true_mu_startY_muP", "Muon start Y (1#mu+1p only);Y [cm]", 70, -70, 70);
TH1F* true_mu_startZ_muP = new TH1F("true_mu_startZ_muP", "Muon start Z (1#mu+1p only);Z [cm]", 70, -70, 70);
TH1F* true_mu_endX_muP = new TH1F("true_mu_endX_muP", "Muon end X (1#mu+1p only);X [cm]", 100, -1200, 1200);
TH1F* true_mu_endY_muP = new TH1F("true_mu_endY_muP", "Muon end Y (1#mu+1p only);Y [cm]", 100, -1800, 1800);
TH1F* true_mu_endZ_muP = new TH1F("true_mu_endZ_muP", "Muon end Z (1#mu+1p only);Z [cm]", 150, 0, 16000); 

TH1F* true_p_len_muP = new TH1F("true_p_len_muP", "Proton Length (1#mu+1p only);Length [cm]", 100, 0, 100);
TH1F* true_p_cos_muP = new TH1F("true_p_cos_muP", "Proton cos(#theta_{z}) (1#mu+1p only);cos(#theta_{z})", 40, -1, 1);
TH1F* true_p_startX_muP = new TH1F("true_p_startX_muP", "Proton start X (1#mu+1p only);X [cm]", 70, -70, 70);
TH1F* true_p_startY_muP = new TH1F("true_p_startY_muP", "Proton start Y (1#mu+1p only);Y [cm]", 70, -70, 70);
TH1F* true_p_startZ_muP = new TH1F("true_p_startZ_muP", "Proton start Z (1#mu+1p only);Z [cm]", 70, -70, 70);
TH1F* true_p_endX_muP = new TH1F("true_p_endX_muP", "Proton end X (1#mu+1p only);X [cm]", 100, -100, 600);
TH1F* true_p_endY_muP = new TH1F("true_p_endY_muP", "Proton end Y (1#mu+1p only);Y [cm]", 100, -300, 100);
TH1F* true_p_endZ_muP = new TH1F("true_p_endZ_muP", "Proton end Z (1#mu+1p only);Z [cm]", 100, -60, 550);


// Reconstructed (MINERvA-matched) histograms for 1 muon + 1 proton
TH1F* reco_nuE_muP = new TH1F("reco_nuE_muP", "Reco Neutrino Energy (1#mu+1p, matched);E_{#nu} [GeV];Events", 20, 0, 20);
TH1F* reco_mult_muP = new TH1F("reco_mult_muP", "Reco multiplicity (1#mu+1p only, matched);Number of Reco Tracks", 10, 0, 10);
TH2F* h2_reco_nuE_vs_muEndZ = new TH2F("h2_reco_nuE_vs_muEndZ", "Reco #nu Energy vs Muon End Z (1#mu+1p, matched);True E_{#nu} [GeV];Muon End Z [cm]", 20, 0, 20, 70, -70, 70);

TH1F* reco_mu_len_muP     = new TH1F("reco_mu_len_muP", "Reco Muon Length (1#mu+1p only, matched);Length [cm]", 100, 0, 200);
TH1F* reco_mu_cos_muP     = new TH1F("reco_mu_cos_muP", "Reco Muon cos(#theta_{z}) (1#mu+1p only, matched);cos(#theta_{z})", 100, -1, 1);
TH1F* reco_mu_startX_muP  = new TH1F("reco_mu_startX_muP", "Reco Muon start X (1#mu+1p only, matched);X [cm]", 70, -70, 70);
TH1F* reco_mu_startY_muP  = new TH1F("reco_mu_startY_muP", "Reco Muon start Y (1#mu+1p only, matched);Y [cm]", 70, -70, 70);
TH1F* reco_mu_startZ_muP  = new TH1F("reco_mu_startZ_muP", "Reco Muon start Z (1#mu+1p only, matched);Z [cm]", 70, -70, 70);
TH1F* reco_mu_endX_muP    = new TH1F("reco_mu_endX_muP", "Reco Muon end X (1#mu+1p only, matched);X [cm]", 70, -70, 70);
TH1F* reco_mu_endY_muP    = new TH1F("reco_mu_endY_muP", "Reco Muon end Y (1#mu+1p only, matched);Y [cm]", 70, -70, 70);
TH1F* reco_mu_endZ_muP    = new TH1F("reco_mu_endZ_muP", "Reco Muon end Z (1#mu+1p only, matched);Z [cm]", 70, -70, 70);

TH1F* reco_p_len_muP      = new TH1F("reco_p_len_muP", "Reco Proton Length (1#mu+1p only, matched);Length [cm]", 50, 0, 100);
TH1F* reco_p_cos_muP      = new TH1F("reco_p_cos_muP", "Reco Proton cos(#theta_{z}) (1#mu+1p only, matched);cos(#theta_{z})", 100, -1, 1);
TH1F* reco_p_startX_muP   = new TH1F("reco_p_startX_muP", "Reco Proton start X (1#mu+1p only, matched);X [cm]", 70, -70, 70);
TH1F* reco_p_startY_muP   = new TH1F("reco_p_startY_muP", "Reco Proton start Y (1#mu+1p only, matched);Y [cm]", 70, -70, 70);
TH1F* reco_p_startZ_muP   = new TH1F("reco_p_startZ_muP", "Reco Proton start Z (1#mu+1p only, matched);Z [cm]", 70, -70, 70);
TH1F* reco_p_endX_muP     = new TH1F("reco_p_endX_muP", "Reco Proton end X (1#mu+1p only, matched);X [cm]", 70, -70, 70);
TH1F* reco_p_endY_muP     = new TH1F("reco_p_endY_muP", "Reco Proton end Y (1#mu+1p only, matched);Y [cm]", 70, -70, 70);
TH1F* reco_p_endZ_muP     = new TH1F("reco_p_endZ_muP", "Reco Proton end Z (1#mu+1p only, matched);Z [cm]", 70, -70, 70);

	// Protons (Contained and Non-Contained)
	
	TH1D* trueProtonMatchedLengthContained =
	new TH1D("trueProtonMatchedLengthContained",
	         "Truth Proton Track Length (Contained);Track Length [cm];Events",
	         100, 0, 100);

	TH1D* trueProtonMatchedLengthNonContained =
		new TH1D("trueProtonMatchedLengthNonContained",
				 "Truth Proton Track Length (Non-contained);Track Length [cm];Events",
				 100, 0, 100);

	TH1D* trueProtonMatchedEContained =
		new TH1D("trueProtonMatchedEContained",
				 "Truth Proton Energy (Contained);Energy [GeV];Events",
				 100, 0, 2);

	TH1D* trueProtonMatchedENonContained =
		new TH1D("trueProtonMatchedENonContained",
				 "Truth Proton Energy (Non-contained);Energy [GeV];Events",
				 100, 0, 2);

	TH2D* trueProtonMatchedLength_vs_EContained =
		new TH2D("trueProtonMatchedLength_vs_EContained",
				 "Truth Proton (Contained): Track Length vs Energy;"
				 "Track Length [cm];Energy [GeV]",
				 100, 0, 100,
				 100, 0, 2);

	TH2D* trueProtonMatchedLength_vs_ENonContained =
		new TH2D("trueProtonMatchedLength_vs_ENonContained",
				 "Truth Proton (Non-contained): Track Length vs Energy;"
				 "Track Length [cm];Energy [GeV]",
				 100, 0, 100,
				 100, 0, 2);

	TH1D* recoProtonMatchedLengthContained =
	new TH1D("recoProtonMatchedLengthContained",
	         "Reco Proton Track Length (Contained);Track Length [cm];Events",
	         100, 0, 100);

	TH1D* recoProtonMatchedLengthNonContained =
		new TH1D("recoProtonMatchedLengthNonContained",
				 "Reco Proton Track Length (Non-contained);Track Length [cm];Events",
				 100, 0, 100);

	TH1D* recoProtonMatchedEContained =
		new TH1D("recoProtonMatchedEContained",
				 "Reco Proton Energy (Contained);Energy [GeV];Events",
				 100, 0, 2);

	TH1D* recoProtonMatchedENonContained =
		new TH1D("recoProtonMatchedENonContained",
				 "Reco Proton Energy (Non-contained);Energy [GeV];Events",
				 100, 0, 2);

	TH2D* recoProtonMatchedLength_vs_EContained =
		new TH2D("recoProtonMatchedLength_vs_EContained",
				 "Reco Proton (Contained): Track Length vs Energy;"
				 "Track Length [cm];Energy [GeV]",
				 100, 0, 100,
				 100, 0, 2);

	TH2D* recoProtonMatchedLength_vs_ENonContained =
		new TH2D("recoProtonMatchedLength_vs_ENonContained",
				 "Reco Proton (Non-contained): Track Length vs Energy;"
				 "Track Length [cm];Energy [GeV]",
				 100, 0, 100,
				 100, 0, 2);


	// Pions (Contained and Non-Contained)
	
	TH1D* truePionMatchedLengthContained =
	new TH1D("truePionMatchedLengthContained",
	         "Truth Pion Track Length (Contained);Track Length [cm];Events",
	         100, 0, 100);

	TH1D* truePionMatchedLengthNonContained =
		new TH1D("truePionMatchedLengthNonContained",
				 "Truth Pion Track Length (Non-contained);Track Length [cm];Events",
				 100, 0, 100);

	TH1D* truePionMatchedEContained =
		new TH1D("truePionMatchedEContained",
				 "Truth Pion Energy (Contained);Energy [GeV];Events",
				 100, 0, 2);

	TH1D* truePionMatchedENonContained =
		new TH1D("truePionMatchedENonContained",
				 "Truth Pion Energy (Non-contained);Energy [GeV];Events",
				 100, 0, 2);

	TH2D* truePionMatchedLength_vs_EContained =
		new TH2D("truePionMatchedLength_vs_EContained",
				 "Truth Pion (Contained): Track Length vs Energy;"
				 "Track Length [cm];Energy [GeV]",
				 100, 0, 100,
				 100, 0, 2);

	TH2D* truePionMatchedLength_vs_ENonContained =
		new TH2D("truePionMatchedLength_vs_ENonContained",
				 "Truth Pion (Non-contained): Track Length vs Energy;"
				 "Track Length [cm];Energy [GeV]",
				 100, 0, 100,
				 100, 0, 2);

	TH1D* recoPionMatchedLengthContained =
	new TH1D("recoPionMatchedLengthContained",
	         "Reco Pion Track Length (Contained);Track Length [cm];Events",
	         100, 0, 100);

	TH1D* recoPionMatchedLengthNonContained =
		new TH1D("recoPionMatchedLengthNonContained",
				 "Reco Pion Track Length (Non-contained);Track Length [cm];Events",
				 100, 0, 100);

	TH1D* recoPionMatchedEContained =
		new TH1D("recoPionMatchedEContained",
				 "Reco Pion Energy (Contained);Energy [GeV];Events",
				 100, 0, 2);

	TH1D* recoPionMatchedENonContained =
		new TH1D("recoPionMatchedENonContained",
				 "Reco Pion Energy (Non-contained);Energy [GeV];Events",
				 100, 0, 2);

	TH2D* recoPionMatchedLength_vs_EContained =
		new TH2D("recoPionMatchedLength_vs_EContained",
				 "Reco Pion (Contained): Track Length vs Energy;"
				 "Track Length [cm];Energy [GeV]",
				 100, 0, 100,
				 100, 0, 2);

	TH2D* recoPionMatchedLength_vs_ENonContained =
		new TH2D("recoPionMatchedLength_vs_ENonContained",
				 "Reco Pion (Non-contained): Track Length vs Energy;"
				 "Track Length [cm];Energy [GeV]",
				 100, 0, 100,
				 100, 0, 2);

TH1D* recoProtonTrueID = new TH1D("recoProtonTrueID", 
    "True PDG of Reco Proton Candidates (after cuts);True Particle ID;Count", 
    4, 0, 4);

	// Give an input list
	std::ifstream caf_list(input_file_list.c_str());

	// Check if input list is present
	if(!caf_list.is_open()){
		std::cerr << Form("File %s not found", input_file_list.c_str()) << std::endl;
		return 1;
	}

	// Add files to CAF chain from input list
	std::string tmp;
	TChain* caf_chain = new TChain("cafTree");

	while(caf_list >> tmp){
		caf_chain->Add(tmp.c_str());
		std::cout << Form("Adding File %s", tmp.c_str()) << std::endl;
	}

	// Check if CAF tree is present
	if(!caf_chain){
		std::cerr << Form("There is no tree in %s", tmp.c_str()) << std::endl;
		return 1;
	}

	long Nentries = caf_chain->GetEntries();
	std::cout << Form("Total number of spills = %ld", Nentries) << std::endl;

	// Define Standard Record and link it to the CAF tree branch "rec"
	auto sr = new caf::StandardRecord;
	caf_chain->SetBranchAddress("rec", &sr);

	bool skipEvent = false;

	// for nusyst start
	// --- Output TTree for nusyst reweighting ---
	int out_genieIdx = -1;
	int out_trackMult = -1;
	int out_hadronMult = -1;
	double out_trueEnu = -1;
	int out_cafFileIdx = -1;

	TTree* syst_tree = new TTree("syst_reco_events", "Reco-matched events for nusyst reweighting");
	syst_tree->Branch("genieIdx", &out_genieIdx, "genieIdx/I");
	syst_tree->Branch("trackMult", &out_trackMult, "trackMult/I");
	syst_tree->Branch("hadronMult", &out_hadronMult, "hadronMult/I");
	syst_tree->Branch("trueEnu", &out_trueEnu, "trueEnu/D");
	syst_tree->Branch("cafFileIdx", &out_cafFileIdx, "cafFileIdx/I");

	// --- Additional output TTrees for nusyst reweighting ---
	int out_trueMultTrkOnly = -1;
	TTree* syst_tree_trueMult = new TTree("syst_true_mult_events", "Truth-selected events for nusyst reweighting");
	syst_tree_trueMult->Branch("genieIdx", &out_genieIdx, "genieIdx/I");
	syst_tree_trueMult->Branch("trueMultTrkOnly", &out_trueMultTrkOnly, "trueMultTrkOnly/I");
	syst_tree_trueMult->Branch("cafFileIdx", &out_cafFileIdx, "cafFileIdx/I");

	double out_recoProtonCosL = -999;
	TTree* syst_tree_recoProtonCosL = new TTree("syst_reco_proton_cosL", "Reco matched proton cosL for nusyst reweighting");
	syst_tree_recoProtonCosL->Branch("genieIdx", &out_genieIdx, "genieIdx/I");
	syst_tree_recoProtonCosL->Branch("recoProtonCosL", &out_recoProtonCosL, "recoProtonCosL/D");
	syst_tree_recoProtonCosL->Branch("cafFileIdx", &out_cafFileIdx, "cafFileIdx/I");

    double out_recoProtonLength = -999;
	TTree* syst_tree_recoProtonLength = new TTree("syst_reco_proton_Length", "Reco matched proton Length for nusyst reweighting");
	syst_tree_recoProtonLength->Branch("genieIdx", &out_genieIdx, "genieIdx/I");
	syst_tree_recoProtonLength->Branch("recoProtonLength", &out_recoProtonLength, "recoProtonLength/D");
	syst_tree_recoProtonLength->Branch("cafFileIdx", &out_cafFileIdx, "cafFileIdx/I");
	// for nusyst end
	

// G4RW syst tree — true pions from matched interactions
TTree* syst_tree_g4rw_pion = new TTree("syst_g4rw_pion", "G4RW true pion syst");
long   g4rw_nu_id    = -1;
int    g4rw_fileIdx  = -1;
int    g4rw_G4ID     = -1;
double g4rw_cosTheta = -9;
syst_tree_g4rw_pion->Branch("nu_id",      &g4rw_nu_id,    "nu_id/L");
syst_tree_g4rw_pion->Branch("cafFileIdx", &g4rw_fileIdx,  "cafFileIdx/I");
syst_tree_g4rw_pion->Branch("G4ID",       &g4rw_G4ID,     "G4ID/I");
syst_tree_g4rw_pion->Branch("cosTheta",   &g4rw_cosTheta, "cosTheta/D");
// G4RW end

// G4RW syst tree — reco pions from matched interactions
TTree* syst_tree_g4rw_pion_reco = new TTree("syst_g4rw_pion_reco", "G4RW reco pion syst");
double g4rw_reco_cosTheta = -9;
syst_tree_g4rw_pion_reco->Branch("nu_id",      &g4rw_nu_id,         "nu_id/L");
syst_tree_g4rw_pion_reco->Branch("cafFileIdx", &g4rw_fileIdx,       "cafFileIdx/I");
syst_tree_g4rw_pion_reco->Branch("G4ID",       &g4rw_G4ID,          "G4ID/I");
syst_tree_g4rw_pion_reco->Branch("cosTheta",   &g4rw_reco_cosTheta, "cosTheta/D");
// G4RW end

// G4RW proton truth tree
TTree* syst_tree_g4rw_proton = new TTree("syst_g4rw_proton", "G4RW truth proton syst");
double g4rw_proton_cosTheta = -9;
double g4rw_proton_length   = -9;
syst_tree_g4rw_proton->Branch("nu_id",      &g4rw_nu_id,           "nu_id/L");
syst_tree_g4rw_proton->Branch("cafFileIdx", &g4rw_fileIdx,         "cafFileIdx/I");
syst_tree_g4rw_proton->Branch("G4ID",       &g4rw_G4ID,            "G4ID/I");
syst_tree_g4rw_proton->Branch("cosTheta",   &g4rw_proton_cosTheta, "cosTheta/D");
syst_tree_g4rw_proton->Branch("length",     &g4rw_proton_length,   "length/D");

// G4RW proton reco tree
TTree* syst_tree_g4rw_proton_reco = new TTree("syst_g4rw_proton_reco", "G4RW reco proton syst");
double g4rw_proton_reco_cosTheta = -9;
double g4rw_proton_reco_length   = -9;
int    g4rw_truePDG = 0;   // NEW: signed truth PDG of the matched truth particle
syst_tree_g4rw_proton_reco->Branch("nu_id",      &g4rw_nu_id,                "nu_id/L");
syst_tree_g4rw_proton_reco->Branch("cafFileIdx", &g4rw_fileIdx,              "cafFileIdx/I");
syst_tree_g4rw_proton_reco->Branch("G4ID",       &g4rw_G4ID,                 "G4ID/I");
syst_tree_g4rw_proton_reco->Branch("truePDG",    &g4rw_truePDG,              "truePDG/I");   // NEW
syst_tree_g4rw_proton_reco->Branch("cosTheta",   &g4rw_proton_reco_cosTheta, "cosTheta/D");
syst_tree_g4rw_proton_reco->Branch("length",     &g4rw_proton_reco_length,   "length/D");
// G4RW end

// G4RW syst tree — ALL multiplicity-loop species (mu, pi+/pi-, proton, kaon),
// built from the SAME selection that fills track_multGood (5cm vertex-
// consistency, length > minTrkLength, primary==true) — unlike the proton/
// pion trees above, which are scoped to a separate, stricter downstream
// block (active containment cut, length > 3cm) and don't match this study.
TTree* syst_tree_g4rw_mult_reco = new TTree("syst_g4rw_mult_reco", "G4RW reco syst for track multiplicity (all species)");
int    g4rw_mult_G4ID     = -1;
int    g4rw_mult_truePDG  = 0;      // signed truth PDG -> species + sign downstream
double g4rw_mult_cosTheta = -9;
double g4rw_mult_length   = -9;
syst_tree_g4rw_mult_reco->Branch("nu_id",      &g4rw_nu_id,         "nu_id/L");
syst_tree_g4rw_mult_reco->Branch("cafFileIdx", &g4rw_fileIdx,       "cafFileIdx/I");
syst_tree_g4rw_mult_reco->Branch("G4ID",       &g4rw_mult_G4ID,     "G4ID/I");
syst_tree_g4rw_mult_reco->Branch("truePDG",    &g4rw_mult_truePDG,  "truePDG/I");
syst_tree_g4rw_mult_reco->Branch("cosTheta",   &g4rw_mult_cosTheta, "cosTheta/D");
syst_tree_g4rw_mult_reco->Branch("length",     &g4rw_mult_length,   "length/D");
// G4RW mult end


	// Loop over all entries (spills) in the CAF chain
	for(long n = 0; n < Nentries; ++n) { 
		
		// Initialize variables for the true neutrino vertex position
		double trueVtxX = -9999; 
		double trueVtxY = -999; 
		double trueVtxZ =- 9999;
        skipEvent = false;

		if (n % 10000 == 0) std::cout << Form("Processing trigger %ld of %ld", n, Nentries) << std::endl;

		caf_chain->GetEntry(n);	// Get spill from tree
		
		// for nusyst start
		// Determine which CAF file this entry belongs to
		out_cafFileIdx = caf_chain->GetTreeNumber();
		// for nusyst end


		bool hasANeutrino = false;

		double mnvOffsetX = -10;
		double mnvOffsetY = 5;

		// If running on Monte Carlo only, remove MINERvA vertex offset
		// (since truth vertices are already in DUNE coordinates)
		if (mcOnly) {
			mnvOffsetX = 0; 
			mnvOffsetY = 0;
											  
			// Loop over all true neutrino interactions in the event (MC truth info)
			for(long unsigned ntrue = 0; ntrue < sr->mc.nu.size(); ntrue++) {
				auto vertex = sr->mc.nu[ntrue].vtx;
				trueVtxX = vertex.x; 
				trueVtxY = vertex.y; 
				trueVtxZ = vertex.z;
      
				auto truePrimary = sr->mc.nu[ntrue].prim;
				
				double Enu = sr->mc.nu[ntrue].E;
				int mode   = sr->mc.nu[ntrue].mode;  // GENIE mode code

				// Fil neutrino energy without any selection
				true_nuE_all_nosel->Fill(Enu);

				if (mode == 1 || mode == 1001) true_nuE_QE_nosel->Fill(Enu);
				else if (mode == 10)           true_nuE_MEC_nosel->Fill(Enu);
				else if (mode == 3)            true_nuE_DIS_nosel->Fill(Enu);
				else if (mode == 4)            true_nuE_RES_nosel->Fill(Enu);
				else if (mode == 5)            true_nuE_COH_nosel->Fill(Enu);


				// Fill no-selection histograms FIRST (before any cuts)
				int truePartTrkOnly_temp = 0;    // Temporary counter for no-selection
				
				for(int k = 0; k < sr->mc.nu[ntrue].prim.size(); k++) {
					int pdg_temp = abs(sr->mc.nu[ntrue].prim[k].pdg);
					if (pdg_temp == 13 || pdg_temp == 211 || pdg_temp == 2212 || pdg_temp == 321) {
						truePartTrkOnly_temp++;
					}
				}
				
				if (truePartTrkOnly_temp > 0) {
					true_multTrkOnly_nosel->Fill(truePartTrkOnly_temp);
					
					int mode_temp = sr->mc.nu[ntrue].mode;
					if (mode_temp == 1 || mode_temp == 1001) {
						true_multTrkOnly_QE_nosel->Fill(truePartTrkOnly_temp);
					}
					else if (mode_temp == 10) {
						true_multTrkOnly_MEC_nosel->Fill(truePartTrkOnly_temp);
					}
					else if (mode_temp == 4) {
						true_multTrkOnly_RES_nosel->Fill(truePartTrkOnly_temp);
					}
					else if (mode_temp == 3) {
						true_multTrkOnly_DIS_nosel->Fill(truePartTrkOnly_temp);
					}
					else if (mode_temp == 5) {
						true_multTrkOnly_COH_nosel->Fill(truePartTrkOnly_temp);
					}
				}


				int truePart = 0; 			// Total number of primary particles (all types)
				int truePartTrkOnly = 0;    // Number of primary particles that are trackable (e.g., charged)
				int truePartNoG4 = 0;       // Number of particles not propagated by GEANT4    
				
				int trueMuonCount = 0;
				int truePionCount = 0;
				int trueProtonCount = 0;
				int trueKaonCount = 0;
				
				if (sr->mc.nu[ntrue].targetPDG != 1000180400) continue;	// skip if target is not Argon-40
				if (sr->mc.nu[ntrue].iscc == false) continue;			// skip if not charged current
				if (abs(sr->mc.nu[ntrue].pdg) != 14) continue;			// skip if not muon neutrino

				true_nuE_all->Fill(Enu);
				if (mode == 1 || mode == 1001) true_nuE_QE->Fill(Enu);
				else if (mode == 10)           true_nuE_MEC->Fill(Enu);
				else if (mode == 3)            true_nuE_DIS->Fill(Enu);
				else if (mode == 4)            true_nuE_RES->Fill(Enu);
				else if (mode == 5)            true_nuE_COH->Fill(Enu);

				// Apply true neutrino vertex cuts to stay within fiducial volume
				if (abs(sr->mc.nu[ntrue].vtx.x) > 59 || abs(sr->mc.nu[ntrue].vtx.x) < 5) continue;
				if (abs(sr->mc.nu[ntrue].vtx.z) > 59.5 || abs(sr->mc.nu[ntrue].vtx.z) < 5) continue;
				if (abs(sr->mc.nu[ntrue].vtx.y) > 57) continue;
				
				true_vtxX->Fill(sr->mc.nu[ntrue].vtx.x);
				true_vtxY->Fill(sr->mc.nu[ntrue].vtx.y);
				true_vtxZ->Fill(sr->mc.nu[ntrue].vtx.z);
				
				// Fill neutrino energy histograms after all truth-level cuts
				// double Enu = sr->mc.nu[ntrue].E;
				// int mode   = sr->mc.nu[ntrue].mode;  // GENIE mode code

				// true_nuE_all->Fill(Enu);

				// if (mode == 1 || mode == 1001) true_nuE_QE->Fill(Enu);
				// else if (mode == 10)           true_nuE_MEC->Fill(Enu);
				// else if (mode == 3)            true_nuE_DIS->Fill(Enu);
				// else if (mode == 4)            true_nuE_RES->Fill(Enu);
				// else if (mode == 5)            true_nuE_COH->Fill(Enu);

				// Get GENIE truth-level counts of final state particles
				int npip      = sr->mc.nu[ntrue].npip;
				int npim      = sr->mc.nu[ntrue].npim;
				int np        = sr->mc.nu[ntrue].nproton;
				int nneutrons = sr->mc.nu[ntrue].nneutron;

				// Will be used to count *primary* versions of these particles (from GENIE prim[])
				int npipPrimaries      = 0;
				int npimPrimaries      = 0;
				int npPrimaries        = 0;
				int nneutronsPrimaries = 0;
				
				// Loop over GENIE primary particles and count particles
				for(int k = 0; k < sr->mc.nu[ntrue].prim.size(); k++) {
					if(sr->mc.nu[ntrue].prim[k].pdg ==  211)  npipPrimaries++;
					if(sr->mc.nu[ntrue].prim[k].pdg == 2212)  npPrimaries++;
					if(sr->mc.nu[ntrue].prim[k].pdg == -211)  npimPrimaries++;
					if(sr->mc.nu[ntrue].prim[k].pdg == 2112)  nneutronsPrimaries++;
				}

				// Fill histograms with the difference between GENIE primaries and total final-state counts
				diffPip->Fill(npipPrimaries - npip);
				diffPim->Fill(npimPrimaries - npim);
				diffP->Fill(npPrimaries - np);
				diffN->Fill(nneutronsPrimaries - nneutrons);

				int nProton = 0;		// Number of mc truth protons
				int nPion = 0;			// Number of mc truth charged pions
				int escapingPions = 0;	// Number of mc truth pions escaping detector volume

				double Elep = -999;  // Muon energy
				double cosL = -999;  // Muon cosθ
				
				// Loop through GENIE true primary particles
				for (long unsigned primaries = 0; primaries < truePrimary.size(); primaries++) {
					int pdg = truePrimary[primaries].pdg;
					double E = truePrimary[primaries].p.E;

					double px = truePrimary[primaries].p.px;
					double py = truePrimary[primaries].p.py;
					double pz = truePrimary[primaries].p.pz;

					double totP = TMath::Sqrt(px * px + py * py + pz * pz);
					
					// Filter: consider only muons, protons, charged pions, and charged kaons
					if ((abs(pdg) == 13 || abs(pdg) == 2212 || abs(pdg) == 211 || abs(pdg) == 321)) {
						truePartNoG4++;		// Increment count of true primary particles not from Geant4

						if (abs(pdg) == 211) nPion++;
						if (abs(pdg) == 2212) nProton++;

						auto start_pos = sr->mc.nu[ntrue].prim[primaries].start_pos;
						auto end_pos = sr->mc.nu[ntrue].prim[primaries].end_pos;
						
						// Optional check for NaN (commented out)
						// if (std::isnan(start_pos.z)) { continue; }

						auto p = sr->mc.nu[ntrue].prim[primaries].p;

						double dX = end_pos.x - start_pos.x;
						double dY = end_pos.y - start_pos.y;
						double dZ = end_pos.z - start_pos.z;

						// double dX = abs(end_pos.x - start_pos.x);
						// double dY = abs(end_pos.y - start_pos.y);
						// double dZ = abs(end_pos.z - start_pos.z);

						double length = TMath::Sqrt(dX * dX + dY * dY + dZ * dZ);
						double lengthP = TMath::Sqrt(p.px * p.px + p.py * p.py + p.pz * p.pz);
						
						double cosTheta = (totP > 0) ? pz / totP : 0;
						
						 if (abs(pdg) == 13) { // Muon
							true_muon_energy->Fill(E);
							true_muon_length->Fill(length);
							true_muon_cosTheta->Fill(cosTheta);
							true_muon_startX->Fill(start_pos.x);
							true_muon_startY->Fill(start_pos.y);
							true_muon_startZ->Fill(start_pos.z);
							true_muon_length_vs_energy->Fill(length,E);

						}
						else if (abs(pdg) == 211) { // Charged pion
							true_pion_energy->Fill(E - 0.139570);
							true_pion_length->Fill(length);
							true_pion_cosTheta->Fill(cosTheta);
							true_pion_startX->Fill(start_pos.x);
							true_pion_startY->Fill(start_pos.y);
							true_pion_startZ->Fill(start_pos.z);
							true_pion_endX->Fill(end_pos.x);
							true_pion_endY->Fill(end_pos.y);
							true_pion_endZ->Fill(end_pos.z);
							true_pion_momentumX->Fill(p.px);
							true_pion_momentumY->Fill(p.py);
							true_pion_momentumZ->Fill(p.pz);
							true_pion_length_vs_energy->Fill(length, E - 0.139570);
							true_pion_startX_vs_startY->Fill(start_pos.x, start_pos.y);
							true_pion_startZ_vs_startX->Fill(start_pos.z, start_pos.x);
							true_pion_startZ_vs_startY->Fill(start_pos.z, start_pos.y);
							true_pion_startX_vs_endX->Fill(start_pos.x, end_pos.x);
							true_pion_startY_vs_endY->Fill(start_pos.y, end_pos.y);
							true_pion_startZ_vs_endZ->Fill(start_pos.z, end_pos.z);							
							true_pion_endX_vs_endY->Fill(end_pos.x, end_pos.y);
							true_pion_endZ_vs_endX->Fill(end_pos.z, end_pos.x);
							true_pion_endZ_vs_endY->Fill(end_pos.z, end_pos.y);
							true_pion_length_vs_cosTheta->Fill(length, cosTheta);
							true_pion_energy_vs_cosTheta->Fill(E - 0.139570, cosTheta);
							true_pion_dX_vs_length->Fill(length, dX);
							true_pion_dY_vs_length->Fill(length, dY);
							true_pion_dZ_vs_length->Fill(length, dZ);

							// true_pion_energy->Fill(E);
							// true_pion_length->Fill(length);
							// true_pion_cosTheta->Fill(cosTheta);
							// true_pion_startX->Fill(start_pos.x);
							// true_pion_startY->Fill(start_pos.y);
							// true_pion_startZ->Fill(start_pos.z);
							// true_pion_length_vs_energy->Fill(length,E);
						}
						else if (pdg == 2212) { // Proton
							true_proton_energy->Fill(E - 0.938272);
							true_proton_length->Fill(length);
							true_proton_cosTheta->Fill(cosTheta);
							true_proton_startX->Fill(start_pos.x);
							true_proton_startY->Fill(start_pos.y);
							true_proton_startZ->Fill(start_pos.z);
							true_proton_endX->Fill(end_pos.x);
							true_proton_endY->Fill(end_pos.y);
							true_proton_endZ->Fill(end_pos.z);
							true_proton_momentumX->Fill(p.px);
							true_proton_momentumY->Fill(p.py);
							true_proton_momentumZ->Fill(p.pz);
							true_proton_length_vs_energy->Fill(length, E - 0.938272);
							true_proton_startX_vs_startY->Fill(start_pos.x, start_pos.y);
							true_proton_startZ_vs_startX->Fill(start_pos.z, start_pos.x);
							true_proton_startZ_vs_startY->Fill(start_pos.z, start_pos.y);
							true_proton_startX_vs_endX->Fill(start_pos.x, end_pos.x);
							true_proton_startY_vs_endY->Fill(start_pos.y, end_pos.y);
							true_proton_startZ_vs_endZ->Fill(start_pos.z, end_pos.z);							
							true_proton_endX_vs_endY->Fill(end_pos.x, end_pos.y);
							true_proton_endZ_vs_endX->Fill(end_pos.z, end_pos.x);
							true_proton_endZ_vs_endY->Fill(end_pos.z, end_pos.y);
							true_proton_length_vs_cosTheta->Fill(length, cosTheta);
							true_proton_energy_vs_cosTheta->Fill(E - 0.938272, cosTheta);
							true_proton_dX_vs_length->Fill(length, dX);
							true_proton_dY_vs_length->Fill(length, dY);
							true_proton_dZ_vs_length->Fill(length, dZ);

						}
						else if (abs(pdg) == 321) { // Charged kaon
							true_kaon_energy->Fill(E);
							true_kaon_length->Fill(length);
							true_kaon_cosTheta->Fill(cosTheta);
							true_kaon_startX->Fill(start_pos.x);
							true_kaon_startY->Fill(start_pos.y);
							true_kaon_startZ->Fill(start_pos.z);
						}

						// Check if the interaction is CC muon-neutrino and the particle is a muon
						if (sr->mc.nu[ntrue].iscc && abs(sr->mc.nu[ntrue].pdg) == 14 && abs(pdg) == 13) {

							Elep = sr->mc.nu[ntrue].prim[primaries].p.E;
							cosL = dZ / length;

							histEl->Fill(Elep);
							histCosl->Fill(cosL);

							// std::cout << n << "," << ntrue << "," << primaries << "," << dX << "," << dY << "," << dZ << "," << p.pz / lengthP << std::endl;

							// Fill histogram of the difference between position-based and momentum-based cos_theta
							trueDiffPosvsPDirZ->Fill(cosL - p.pz / lengthP);
						}
						
						// If the pion's end z or x position is outside detector bounds, count it as escaping
						if ((abs(end_pos.z) > 60 || abs(end_pos.x) > 60) && abs(pdg) == 211) {
							escapingPions++;
						}
						// If pion is contained in x and z bounds, fill its track length
						else if (abs(pdg) == 211) {
							containPiLen->Fill(length);
						}

						// If proton is fully contained in x and z, fill its length
						if (abs(end_pos.z) < 60 && abs(end_pos.x) < 60 && abs(pdg) == 2212) {
							containPLen->Fill(length);
						}

						// if (cosL < 0.9) continue;

						// Only count particles (of accepted types) that pass minimum track length threshold
						if (length > minTrkLength) {
							truePartTrkOnly++;
							
							if (abs(pdg) == 13) trueMuonCount++;
							else if (abs(pdg) == 211) truePionCount++;
							else if (abs(pdg) == 2212) trueProtonCount++;
							else if (abs(pdg) == 321) trueKaonCount++;
						}
					}
					
					truePart++;  // Count as a valid GENIE-level primary particle
				}
				
				// if (cosL < 0.9 && Elep < 1) continue;


// Extra logic: store start position and length for 1 mu + 1 p truth events
int nMu = 0; int nProton_muP = 0;
int nOther = 0;
for (auto const& part : truePrimary) {
	if (abs(part.pdg) == 13) ++nMu;
	else if (abs(part.pdg) == 2212) ++nProton_muP;
	else ++nOther;
}
if (nMu == 1 && nProton_muP == 1 && nOther == 0) {
// if (nProton == 1 && truePart == 2) {
	double mu_len = -1, mu_cos = -2, mu_startX = -999, mu_startY = -999, mu_startZ = -999, mu_endX = -999, mu_endY = -999, mu_endZ = -999;
	double p_len = -1, p_cos = -2, p_startX = -999, p_startY = -999, p_startZ = -999, p_endX = -999, p_endY = -999, p_endZ = -999;

	for (size_t i = 0; i < truePrimary.size(); ++i) {
		auto const& part = truePrimary[i];
		if (abs(part.pdg) == 13) {
			auto start = part.start_pos;
			auto end   = part.end_pos;
			double dx = end.x - start.x;
			double dy = end.y - start.y;
			double dz = end.z - start.z;
			double len = std::sqrt(dx*dx + dy*dy + dz*dz);
			double cosz = len > 0 ? dz / len : -2;
			mu_len = len;
			mu_cos = cosz;
			mu_startX = start.x;
			mu_startY = start.y;
			mu_startZ = start.z;
			mu_endX = end.x;
			mu_endY = end.y;
			mu_endZ = end.z;

			std::cout << "mu_len: " << mu_len << ", mu_endX: " << mu_endX << ", mu_endY: " << mu_endY << ", mu_endZ: " << mu_endZ << std::endl;



		}
		if (abs(part.pdg) == 2212) {
			auto start = part.start_pos;
			auto end   = part.end_pos;
			double dx = end.x - start.x;
			double dy = end.y - start.y;
			double dz = end.z - start.z;
			double len = std::sqrt(dx*dx + dy*dy + dz*dz);
			double cosz = len > 0 ? dz / len : -2;
			p_len = len;
			p_cos = cosz;
			p_startX = start.x;
			p_startY = start.y;
			p_startZ = start.z;
			p_endX = end.x;
			p_endY = end.y;
			p_endZ = end.z;
			
			std::cout << "p_len: " << p_len << ", p_endX: " << p_endX << ", p_endY: " << p_endY << ", p_endZ: " << p_endZ << std::endl;

		}
	}
					
	double trueE = sr->mc.nu[ntrue].E;
	true_nuE_muP->Fill(trueE);
	std::cout << "true_nuE: " << trueE << std::endl;
					
	true_mult_muP->Fill(nMu + nProton_muP);  

	true_mu_len_muP->Fill(mu_len);
	true_mu_cos_muP->Fill(mu_cos);
	true_mu_startX_muP->Fill(mu_startX);
	true_mu_startY_muP->Fill(mu_startY);
	true_mu_startZ_muP->Fill(mu_startZ);
	true_mu_endX_muP->Fill(mu_endX);
	true_mu_endY_muP->Fill(mu_endY);
	true_mu_endZ_muP->Fill(mu_endZ);
	
	true_p_len_muP->Fill(p_len);
	true_p_cos_muP->Fill(p_cos);
	true_p_startX_muP->Fill(p_startX);
	true_p_startY_muP->Fill(p_startY);
	true_p_startZ_muP->Fill(p_startZ);
	true_p_endX_muP->Fill(p_endX);
	true_p_endY_muP->Fill(p_endY);
	true_p_endZ_muP->Fill(p_endZ);
	
	h2_true_nuE_vs_muEndZ->Fill(trueE, mu_endZ);
}


				escapePi->Fill(escapingPions);
				nPi->Fill(nPion);
				nP->Fill(nProton);
				true_mult->Fill(truePart);

				// if (truePartTrkOnly > 0)
					// true_multTrkOnly->Fill(truePartTrkOnly);
				if (truePartTrkOnly > 0) true_multTrkOnly->Fill(truePartTrkOnly);
				if (trueMuonCount > 0) true_multTrkOnly_muon->Fill(trueMuonCount);
				if (truePionCount > 0) true_multTrkOnly_pion->Fill(truePionCount);
				if (trueProtonCount > 0) true_multTrkOnly_proton->Fill(trueProtonCount);
				if (trueKaonCount > 0) true_multTrkOnly_kaon->Fill(trueKaonCount);

				// for nusyst start
				if (truePartTrkOnly > 0) {
					out_genieIdx = sr->mc.nu[ntrue].genieIdx;
					out_trueMultTrkOnly = truePartTrkOnly;
					out_cafFileIdx = caf_chain->GetTreeNumber();
					syst_tree_trueMult->Fill();
				}
				// for nusyst end

				// Fill by interaction mode
				if (mode == 1 || mode == 1001) {
					if (truePartTrkOnly > 0) true_multTrkOnly_QE->Fill(truePartTrkOnly);
				}
				else if (mode == 10) {
					if (truePartTrkOnly > 0) true_multTrkOnly_MEC->Fill(truePartTrkOnly);
				}
				else if (mode == 4) {
					if (truePartTrkOnly > 0) true_multTrkOnly_RES->Fill(truePartTrkOnly);
				}
				else if (mode == 3) {
					if (truePartTrkOnly > 0) true_multTrkOnly_DIS->Fill(truePartTrkOnly);
				}
				else if (mode == 5) {
					if (truePartTrkOnly > 0) true_multTrkOnly_COH->Fill(truePartTrkOnly);
				}

				true_multGENIE->Fill(truePartNoG4);						// Fill GENIE-level multiplicity (ignoring G4 containment)
				responseGenieToG4->Fill(truePart, truePartTrkOnly);
				hasANeutrino = true;
			}
		}
		
		int rock = 0;
		bool goodInteraction = false;

		std::vector<int> trueInteractionIndex;
		std::vector<int> primaryTrkIndex;
	
		// Loop over all reconstructed interactions
		for (long unsigned nixn = 0;  nixn < sr->common.ixn.dlp.size();  nixn++) {
			
			bool oneContained = false;
			bool oneNotContained = false;
			goodInteraction = false;

			double biggestMatch = -999;
			int biggestMatchIndex = -999;
			double maxDotProductDS = -999;
			double maxDotProductUS = -999;

			int maxEventPar = -999;
			int maxEventTyp = -9999;
			int maxEventIxn = -999;

			double dirZExiting = -999;
			double startZMuonCand = -999;
			
			double recoVertexX = sr->common.ixn.dlp[nixn].vtx.x;
			double recoVertexY = sr->common.ixn.dlp[nixn].vtx.y;
			double recoVertexZ = sr->common.ixn.dlp[nixn].vtx.z;

			if (/*abs(abs(sr->common.ixn.dlp[nixn].vtx.x)-33)<1 || */ 
				abs(sr->common.ixn.dlp[nixn].vtx.x) > 59 ||
				abs(sr->common.ixn.dlp[nixn].vtx.x) < 5  ||
				abs(sr->common.ixn.dlp[nixn].vtx.y) > 57 ||
				abs(sr->common.ixn.dlp[nixn].vtx.z) < 5  ||
				abs(sr->common.ixn.dlp[nixn].vtx.z) > 59.5)
				continue;

			recoVertex2DNoCuts->Fill(sr->common.ixn.dlp[nixn].vtx.x, sr->common.ixn.dlp[nixn].vtx.z);
        

			if (mcOnly){
       
				// Loop over truth-matching information for this reconstructed interaction
				for (int ntruth = 0; ntruth < sr->common.ixn.dlp[nixn].truth.size(); ntruth++) {

					// Check if this GENIE interaction has a higher overlap score than the previous best
					if (biggestMatch < sr->common.ixn.dlp[nixn].truthOverlap.at(ntruth)) {
						biggestMatch = sr->common.ixn.dlp[nixn].truthOverlap.at(ntruth);
						biggestMatchIndex = sr->common.ixn.dlp[nixn].truth.at(ntruth);
					}
				}
			
				if (sr->mc.nu[biggestMatchIndex].id > 1E9) {
					rock = 1;
				}
				
				trueVtxX = sr->mc.nu[biggestMatchIndex].vtx.x;
				trueVtxY = sr->mc.nu[biggestMatchIndex].vtx.y;
				trueVtxZ = sr->mc.nu[biggestMatchIndex].vtx.z;

				if (abs(sr->mc.nu[biggestMatchIndex].pdg) == 14 &&
					abs(sr->mc.nu[biggestMatchIndex].iscc) == true &&
					sr->mc.nu[biggestMatchIndex].targetPDG == 1000180400 &&
					abs(recoVertexX - trueVtxX) < 5 &&
					abs(recoVertexY - trueVtxY) < 5 &&
					abs(recoVertexZ - trueVtxZ) < 5 &&
					abs(trueVtxX) > 5 &&
					abs(trueVtxX) < 59 &&
					abs(trueVtxZ) > 5 &&
					abs(trueVtxZ) < 59.5 &&
					abs(trueVtxY) < 57)
				{
					goodInteraction = true;
				}
				
				recoHistVertexY->Fill(sr->common.ixn.dlp[nixn].vtx.y);
				recoHistVertexX->Fill(sr->common.ixn.dlp[nixn].vtx.x);
				recoHistVertexZ->Fill(sr->common.ixn.dlp[nixn].vtx.z);
				
				if (goodInteraction == true)
					goodIntNum++;
				else
					badIntNum++;

				trueInteractionIndex.push_back(biggestMatchIndex);

				if (goodInteraction == true)
					goodIntInFidVol++;
				else
					badIntInFidVol++;
			} 
			
			int partMult = 0;
			int partMultTrkOnly = 0;
			bool hasAtLeastOneProton = false;
			bool hasAtLeastOneMuon = false;
			bool hasAtLeastOne = false;

			// std::cout<<"Interaction: "<<nixn<<std::endl;

			int correctTrack = 2;
			int correctShower = 2;

			double longestTrk = -9999;

			int trackMult = 0;
			int trackMultExit = 0;
			int minervaTracks = 0;
			int minervaThrough = 0;
			
			int hadronMult = 0;
			
			int muonMult_good = 0;
			int pionMult_good = 0;
			int protonMult_good = 0;
			int kaonMult_good = 0;
			
			// Loop over reconstructed particles
			for (long unsigned npart = 0; npart < sr->common.ixn.dlp[nixn].part.dlp.size(); npart++) { 
			
				if (!sr->common.ixn.dlp[nixn].part.dlp[npart].primary) continue;

				int pdg = sr->common.ixn.dlp[nixn].part.dlp[npart].pdg;

				if (abs(pdg) == 2212 || abs(pdg) == 13 || abs(pdg) == 211 || abs(pdg) == 321) {
					if (sr->common.ixn.dlp[nixn].part.dlp[npart].primary == true)
						trackCorrectness->Fill(correctTrack);
				

					auto start_pos = sr->common.ixn.dlp[nixn].part.dlp[npart].start;
					auto end_pos   = sr->common.ixn.dlp[nixn].part.dlp[npart].end;

					double diffVertexdZ = abs(start_pos.z - sr->common.ixn.dlp[nixn].vtx.z);
					double diffVertexdX = abs(start_pos.x - sr->common.ixn.dlp[nixn].vtx.x);
					double diffVertexdY = abs(start_pos.y - sr->common.ixn.dlp[nixn].vtx.y);

					double diffVertex = TMath::Sqrt(
						diffVertexdZ * diffVertexdZ +
						diffVertexdY * diffVertexdY +
						diffVertexdX * diffVertexdX
					);

					if (diffVertex > 5) continue;
					
					double dX = (end_pos.x - start_pos.x);
					double dY = (end_pos.y - start_pos.y);
					double dZ = (end_pos.z - start_pos.z);

					double length = TMath::Sqrt(dX * dX + dY * dY + dZ * dZ);

					double dirX = dX / length;
					double dirY = dY / length;
					double dirZ = dZ / length;
					
					if (dirZ < 0) {
						dirZ = -dirZ;
						dirX = -dirX;
						dirY = -dirY;
						auto temp = start_pos;
						end_pos = start_pos;
						end_pos = temp;
					}

					if (std::isnan(start_pos.z)) length = -999;
					
					if (length > longestTrk) longestTrk = length;
					
					if (sr->common.ixn.dlp[nixn].part.dlp[npart].primary == true && length > minTrkLength) {
						partMult++;
						trackMult++;
						
						if (abs(pdg) == 13) {
							muonMult_good++;
						} else if (abs(pdg) == 211) {
							pionMult_good++;
						} else if (abs(pdg) == 2212) {
							protonMult_good++;
						} else if (abs(pdg) == 321) {
							kaonMult_good++;
						}
						
						// if (abs(pdg) != 13) {
						if (abs(pdg) == 2212 || abs(pdg) == 211 || abs(pdg) == 321) {
							hadronMult++;
						}

						if ((start_pos.z > 60) || (end_pos.z > 60)) trackMultExit++;

						int maxPartMinerva    = -999;
						int maxTypeMinerva    = -999;
						int maxIxnMinerva     = -999;
						int maxPartMinervaUS  = -999;
						int maxTypeMinervaUS  = -999;
						
						if ((abs(start_pos.z) > 58 || abs(end_pos.z) > 58)) {
							int minervaPass = 0;
							
							double dotProductDS    = -999;
							double deltaExtrapYUS  = -999;
							double deltaExtrapY    = -999;
							double dotProductUS    = -999;
							double deltaExtrapX    = -999;
							double deltaExtrapXUS  = -999;
							
							for (int i = 0; i < sr->nd.minerva.ixn.size(); i++) {

								for (int j = 0; j < sr->nd.minerva.ixn[i].ntracks; j++) {
									double dir_z   = sr->nd.minerva.ixn[i].tracks[j].dir.z;

									double end_z   = sr->nd.minerva.ixn[i].tracks[j].end.z;
									double start_z = sr->nd.minerva.ixn[i].tracks[j].start.z;

									double end_x   = sr->nd.minerva.ixn[i].tracks[j].end.x;
									double start_x = sr->nd.minerva.ixn[i].tracks[j].start.x;

									double end_y   = sr->nd.minerva.ixn[i].tracks[j].end.y;
									double start_y = sr->nd.minerva.ixn[i].tracks[j].start.y;
	
									if (start_z > 0 && ((start_pos.z > 60) || (end_pos.z > 60))) {
										
										int truthPart = sr->nd.minerva.ixn[i].tracks[j].truth[0].part;

										double dXMnv = sr->nd.minerva.ixn[i].tracks[j].end.x - sr->nd.minerva.ixn[i].tracks[j].start.x;
										double dYMnv = sr->nd.minerva.ixn[i].tracks[j].end.y - sr->nd.minerva.ixn[i].tracks[j].start.y;
										double dZMnv = sr->nd.minerva.ixn[i].tracks[j].end.z - sr->nd.minerva.ixn[i].tracks[j].start.z;

										double lengthMinerva = TMath::Sqrt(dXMnv * dXMnv + dYMnv * dYMnv + dZMnv * dZMnv);

										if (lengthMinerva < 10) continue;
										
										double dirXMinerva = dXMnv / lengthMinerva;
										double dirYMinerva = dYMnv / lengthMinerva;
										double dirZMinerva = dZMnv / lengthMinerva;

										double dotProduct = dirXMinerva * dirX + dirYMinerva * dirY + dirZ * dirZMinerva;

										double extrapdZ = start_z - end_pos.z;

										double extrapY = dirY / dirZ * extrapdZ + end_pos.y - start_y;
										double extrapX = dirX / dirZ * extrapdZ + end_pos.x - start_x;

										double diffExtrap = TMath::Sqrt(TMath::Power(extrapY - start_y, 2));

										if (dotProductDS < dotProduct &&
											abs(extrapY - mnvOffsetY) < 15 &&
											abs(TMath::ATan(dirXMinerva / dirZMinerva) - TMath::ATan(dirX / dirZ)) < 0.06 &&
											abs(TMath::ATan(dirYMinerva / dirZMinerva) - TMath::ATan(dirY / dirZ)) < 0.06 &&
											abs(extrapX - mnvOffsetX) < 15)
										{
											dotProductDS = dotProduct;
											deltaExtrapY = extrapY;
											deltaExtrapX = extrapX;
											dirZExiting  = dirZ;
											
											if (mcOnly) {
												maxPartMinerva = sr->nd.minerva.ixn[i].tracks[j].truth[0].part;
												maxTypeMinerva = sr->nd.minerva.ixn[i].tracks[j].truth[0].type;
												maxIxnMinerva  = sr->nd.minerva.ixn[i].tracks[j].truth[0].ixn;
											}

											// if (end_z > 300) { minervaPass = 1; }
											// if (dirZExiting < dirZ) { dirZExiting = dirZ; }
										}
									}
									
									if (start_z < 0 && end_z > 0 &&
										((start_pos.z < -58 && end_pos.z > 58) ||
										 (start_pos.z > 58 && end_pos.z < -58))) 
									{
										int truthPart = sr->nd.minerva.ixn[i].tracks[j].truth[0].part;

										double dXMnv = sr->nd.minerva.ixn[i].tracks[j].end.x - sr->nd.minerva.ixn[i].tracks[j].start.x;
										double dYMnv = sr->nd.minerva.ixn[i].tracks[j].end.y - sr->nd.minerva.ixn[i].tracks[j].start.y;
										double dZMnv = sr->nd.minerva.ixn[i].tracks[j].end.z - sr->nd.minerva.ixn[i].tracks[j].start.z;

										double lengthMinerva = TMath::Sqrt(dXMnv * dXMnv + dYMnv * dYMnv + dZMnv * dZMnv);

										double dirXMinerva = dXMnv / lengthMinerva;
										double dirYMinerva = dYMnv / lengthMinerva;
										double dirZMinerva = dZMnv / lengthMinerva;

										double dotProduct = dirXMinerva * dirX + dirYMinerva * dirY + dirZ * dirZMinerva;

										double extrapdZUS = end_z - end_pos.z;

										double extrapYUS = dirY / dirZ * extrapdZUS + end_pos.y - end_y;
										double extrapXUS = dirX / dirZ * extrapdZUS + end_pos.x - end_x;

										// double diffExtrap = TMath::Sqrt(TMath::Power(extrapY - end_y, 2));
										
										if (dotProductUS < dotProduct &&
										abs(extrapYUS - mnvOffsetY) < 15 &&
										abs(extrapXUS - mnvOffsetX) < 15 &&
										abs(TMath::ATan(dirXMinerva / dirZMinerva) - TMath::ATan(dirX / dirZ)) < 0.06 &&
										abs(TMath::ATan(dirYMinerva / dirZMinerva) - TMath::ATan(dirY / dirZ)) < 0.06)
											dotProductUS = dotProduct;
											
											deltaExtrapYUS = extrapYUS;
											deltaExtrapXUS = extrapXUS;
											
										if (mcOnly) {
											maxPartMinervaUS = sr->nd.minerva.ixn[i].tracks[j].truth[0].part;
											maxTypeMinervaUS = sr->nd.minerva.ixn[i].tracks[j].truth[0].type;
										}
									}
								}
							} // Minerva
							
							if (dotProductDS > maxDotProductDS) {
								maxDotProductDS = dotProductDS;
								maxEventPar     = maxPartMinerva;
								maxEventTyp     = maxTypeMinerva;
								maxEventIxn     = maxIxnMinerva;
							}
							
							if (dotProductDS > 0.99) {
								minervaTracks++;

								if (minervaPass == 1) {
									minervaThrough++;
									startZMuonCand = start_pos.z;

									if (start_pos.z > end_pos.z)
										startZMuonCand = end_pos.z;
								}
							}
							
							if (dotProductUS > maxDotProductUS)
								maxDotProductUS = dotProductUS;

							if (dotProductUS > 0.99) {
								/*
								if (maxPartNumber == maxPartMinervaUS && maxTypeMinervaUS == maxTypeNumber) {
									deltaYGoodUS->Fill(deltaExtrapYUS);
									deltaXGoodUS->Fill(deltaExtrapXUS);
									goodMINERvAMatchUS++;
								} else {
									deltaYBadUS->Fill(deltaExtrapYUS);
									deltaXBadUS->Fill(deltaExtrapXUS);
								}
								totalMINERvAMatchUS++;
								*/
							}

							// oneNotContained = true;
						} // exiting particles
						
						//else oneContained = true;
						
					} // primary of particle greater than 5
				} 
			} // particles
			
			// if (/*minervaThrough < 1 ||*/ maxDotProductDS < 0.99) continue;
			bool passesMinerva = (maxDotProductDS > 0.99);

			if (!requireMinerva) passesMinerva = true;

			recoCosL->Fill(dirZExiting);
			recoVertex2D->Fill(sr->common.ixn.dlp[nixn].vtx.x, sr->common.ixn.dlp[nixn].vtx.z);
			
			auto recoVtx = sr->common.ixn.dlp[nixn].vtx;
			recoVtxX->Fill(recoVtx.x);
			recoVtxY->Fill(recoVtx.y);
			recoVtxZ->Fill(recoVtx.z);


			if (mcOnly) {
				
				if (maxEventTyp == 1) {
					minervaMatchPDG->Fill(sr->mc.nu[maxEventIxn].prim[maxEventPar].pdg);
					int minervaPDG = sr->mc.nu[maxEventIxn].prim[maxEventPar].pdg;
					
					if (abs(minervaPDG) == 13) {
						minervaMatchE->Fill(sr->mc.nu[maxEventIxn].prim[maxEventPar].p.E);
						
						auto start_posMnv = sr->mc.nu[maxEventIxn].prim[maxEventPar].start_pos;
						auto end_posMnv   = sr->mc.nu[maxEventIxn].prim[maxEventPar].end_pos;

						double dX = end_posMnv.x - start_posMnv.x;
						double dY = end_posMnv.y - start_posMnv.y;
						double dZ = end_posMnv.z - start_posMnv.z;

						auto p = sr->mc.nu[maxEventIxn].prim[maxEventPar].p;

						double length = TMath::Sqrt(p.px * p.px + p.py * p.py + p.pz * p.pz);
						double cosL   = p.pz / length;

						if (abs(startZMuonCand - start_posMnv.z) < 5)
							minervaMatchCos->Fill(cosL);

						// std::cout << startZMuonCand << "," << start_posMnv.z << "," << cosL << std::endl;
						// std::cout << biggestMatchIndex << "," << maxEventIxn << std::endl;

					}
				}              
				
				// std::cout << "counting protons and kaons" << std::endl;

if (sr->mc.nu[biggestMatchIndex].prim.size() < 500 &&
    sr->mc.nu[biggestMatchIndex].targetPDG == 1000180400 &&
    passesMinerva)
{
    // Always fill reco vertices for selected interactions
    // auto recoVtx = sr->common.ixn.dlp[nixn].vtx;
    // recoVtxX->Fill(recoVtx.x);
    // recoVtxY->Fill(recoVtx.y);
    // recoVtxZ->Fill(recoVtx.z);

    // Only fill truth vertices when the interaction is “goodInteraction”
    if (goodInteraction) {
		
		// Fill reco-of-matched
		auto recoVtx = sr->common.ixn.dlp[nixn].vtx;
		recoVtxX_matched->Fill(recoVtx.x);
		recoVtxY_matched->Fill(recoVtx.y);
		recoVtxZ_matched->Fill(recoVtx.z);
		
        truthMatchedRecoVtxX->Fill(sr->mc.nu[biggestMatchIndex].vtx.x);
        truthMatchedRecoVtxY->Fill(sr->mc.nu[biggestMatchIndex].vtx.y);
        truthMatchedRecoVtxZ->Fill(sr->mc.nu[biggestMatchIndex].vtx.z);
    }
}


				// if (passesMinerva && goodInteraction) {
				if (sr->mc.nu[biggestMatchIndex].prim.size() < 500 &&
					sr->mc.nu[biggestMatchIndex].targetPDG == 1000180400)
				{
					if (passesMinerva && goodInteraction) {
						


						
						// auto recoVtx = sr->common.ixn.dlp[nixn].vtx;
						// recoVtxX->Fill(recoVtx.x);
						// recoVtxY->Fill(recoVtx.y);
						// recoVtxZ->Fill(recoVtx.z);
						
						int cc   = sr->mc.nu[biggestMatchIndex].iscc;
						int mode = sr->mc.nu[biggestMatchIndex].mode;
						double Enu = sr->mc.nu[biggestMatchIndex].E;
						
						// truthMatchedRecoVtxX->Fill(sr->mc.nu[biggestMatchIndex].vtx.x);
						// truthMatchedRecoVtxY->Fill(sr->mc.nu[biggestMatchIndex].vtx.y);
						// truthMatchedRecoVtxZ->Fill(sr->mc.nu[biggestMatchIndex].vtx.z);
						
						
						// Fill truth neutrino energy for truth-matched reco interactions
						truthMatchedReco_nuE->Fill(Enu);

						if (cc == 1) {
							if (mode == 1 || mode == 1001) {
								truthMatchedReco_nuE_QE->Fill(Enu);
							}
							if (mode == 10) {
								truthMatchedReco_nuE_MEC->Fill(Enu);
							}
							if (mode == 3) {
								truthMatchedReco_nuE_DIS->Fill(Enu);
							}
							if (mode == 4) {
								truthMatchedReco_nuE_RES->Fill(Enu);
							}
							if (mode == 5) {
								truthMatchedReco_nuE_COH->Fill(Enu);
							}
						}
						
						for (int primaries = 0; primaries < sr->mc.nu[biggestMatchIndex].prim.size(); primaries++) {

							auto start_pos = sr->mc.nu[biggestMatchIndex].prim[primaries].start_pos;
							auto end_pos   = sr->mc.nu[biggestMatchIndex].prim[primaries].end_pos;

							// if (std::isnan(start_pos.z)) continue;
							// std::cout << "Got a position" << std::endl;

							auto p = sr->mc.nu[biggestMatchIndex].prim[primaries].p;

							double dX = p.px; // (end_pos.x - start_pos.x);
							double dY = p.py; // (end_pos.y - start_pos.y);
							double dZ = p.pz; // (end_pos.z - start_pos.z);

							double length = TMath::Sqrt(dX * dX + dY * dY + dZ * dZ);

							double dirX = dX / length;
							double dirY = dY / length;
							double dirZ = dZ / length;

							double dXLen = end_pos.x - start_pos.x;
							double dYLen = end_pos.y - start_pos.y;
							double dZLen = end_pos.z - start_pos.z;

							double lengthPos = TMath::Sqrt(dXLen * dXLen + dYLen * dYLen + dZLen * dZLen);

							if (!std::isnan(start_pos.z) &&
								(abs(sr->mc.nu[biggestMatchIndex].prim[primaries].pdg) == 211 /*|| abs(...) == 13*/ ||
								 abs(sr->mc.nu[biggestMatchIndex].prim[primaries].pdg) == 2212))
							{
								trueTrkLen->Fill(lengthPos);

								int cc   = sr->mc.nu[biggestMatchIndex].iscc;
								int mode = sr->mc.nu[biggestMatchIndex].mode;

								if (cc == 1 && (mode == 1 || mode == 1001 || mode == 10))
									trueTrkLenQE->Fill(lengthPos);

								if (cc == 1 && (mode == 3 || mode == 4))
									trueTrkLenNonQE->Fill(lengthPos);
							}
							
							if (abs(sr->mc.nu[biggestMatchIndex].prim[primaries].pdg) == 2212) {
								trueProtonEWithRecoInt->Fill(sr->mc.nu[biggestMatchIndex].prim[primaries].p.E - 0.938272);
								trueProtonWithRecoIntDirZ->Fill(dirZ);
								trueProtonWithRecoIntDirX->Fill(dirX);
								trueProtonWithRecoIntDirY->Fill(dirY);                                
							}
							
							if (abs(sr->mc.nu[biggestMatchIndex].prim[primaries].pdg) == 211) {
								truePionEWithRecoInt->Fill(sr->mc.nu[biggestMatchIndex].prim[primaries].p.E - 0.13957);
								truePionWithRecoIntDirZ->Fill(dirZ);
								truePionWithRecoIntDirX->Fill(dirX);
								truePionWithRecoIntDirY->Fill(dirY);

                                // for G4RW start
                                g4rw_nu_id    = (long)sr->mc.nu[biggestMatchIndex].id;
                                g4rw_fileIdx  = caf_chain->GetTreeNumber();
                                g4rw_G4ID     = (int)sr->mc.nu[biggestMatchIndex].prim[primaries].G4ID;
                                g4rw_cosTheta = dirZ;  // pz/|p|, already computed above
                                syst_tree_g4rw_pion->Fill();
                                // for G4RW end
							}
							
							if (abs(sr->mc.nu[biggestMatchIndex].prim[primaries].pdg) == 321) {
								trueKaonEWithRecoInt->Fill(sr->mc.nu[biggestMatchIndex].prim[primaries].p.E - 0.493677);
								trueKaonWithRecoIntDirZ->Fill(dirZ);
								trueKaonWithRecoIntDirX->Fill(dirX);
								trueKaonWithRecoIntDirY->Fill(dirY);
							}
						}
					}
				}
				
				
				// ==========================================================
				// Matched-interaction multiplicity counters
				// (DECLARE ONCE PER INTERACTION — THIS IS THE KEY POINT)
				// ==========================================================
				int nPrimProtons_true    = 0;
				int nRecoProtons_matched = 0;				

				int nPrimPions_true    = 0;
				int nRecoPions_matched = 0;	
						
				for (long unsigned npart = 0; npart < sr->common.ixn.dlp[nixn].part.dlp.size(); npart++) { // loop over particles
					// std::cout << "Containment variable: " << sr->common.ixn.dlp[nixn].part.dlp[npart].contained << std::endl;
					

					if (!sr->common.ixn.dlp[nixn].part.dlp[npart].primary) continue;

					int pdg = sr->common.ixn.dlp[nixn].part.dlp[npart].pdg;

					auto start_pos = sr->common.ixn.dlp[nixn].part.dlp[npart].start;
					auto end_pos   = sr->common.ixn.dlp[nixn].part.dlp[npart].end;

					auto momentum   = sr->common.ixn.dlp[nixn].part.dlp[npart].p;

					double dX = end_pos.x - start_pos.x;
					double dY = end_pos.y - start_pos.y;
					double dZ = end_pos.z - start_pos.z;

					double length = TMath::Sqrt(dX * dX + dY * dY + dZ * dZ);

					double dirX = dX / length;
					double dirY = dY / length;
					double dirZ = dZ / length;
					
					if (abs(pdg) == 13) {
						auto start_pos = sr->common.ixn.dlp[nixn].part.dlp[npart].start;
						auto end_pos   = sr->common.ixn.dlp[nixn].part.dlp[npart].end;
						
						auto momentum   = sr->common.ixn.dlp[nixn].part.dlp[npart].p;

						double dX = end_pos.x - start_pos.x;
						double dY = end_pos.y - start_pos.y;
						double dZ = end_pos.z - start_pos.z;

						double length = TMath::Sqrt(dX * dX + dY * dY + dZ * dZ);
						double dirZ = dZ / length;

						// auto p       = sr->mc.nu[interactionNumber].prim[partNumber].p;

						auto energy = sr->common.ixn.dlp[nixn].part.dlp[npart].E;

						recoMuonLengthAll->Fill(length);
						recoMuonCosLAll->Fill(dirZ);
						recoMuonEAll->Fill(energy);
						recoMuonStartXAll->Fill(start_pos.x);
						recoMuonStartYAll->Fill(start_pos.y);
						recoMuonStartZAll->Fill(start_pos.z);
						recoMuonEndXAll->Fill(end_pos.x);
						recoMuonEndYAll->Fill(end_pos.y);
						recoMuonEndZAll->Fill(end_pos.z);
						recoMuonMomentumXAll->Fill(momentum.x);
						recoMuonMomentumYAll->Fill(momentum.y);
						recoMuonMomentumZAll->Fill(momentum.z);
						recoMuonLengthAll_vs_E->Fill(length, energy);
					}					
					
					if (abs(pdg) == 2212) {
						auto start_pos = sr->common.ixn.dlp[nixn].part.dlp[npart].start;
						auto end_pos   = sr->common.ixn.dlp[nixn].part.dlp[npart].end;
						
						auto momentum   = sr->common.ixn.dlp[nixn].part.dlp[npart].p;

						double dX = end_pos.x - start_pos.x;
						double dY = end_pos.y - start_pos.y;
						double dZ = end_pos.z - start_pos.z;

						double length = TMath::Sqrt(dX * dX + dY * dY + dZ * dZ);
						double dirZ = dZ / length;

						// auto p       = sr->mc.nu[interactionNumber].prim[partNumber].p;

						auto energy = sr->common.ixn.dlp[nixn].part.dlp[npart].E;

						recoProtonLengthAll->Fill(length);
						recoProtonCosLAll->Fill(dirZ);
						recoProtonEAll->Fill(energy);
						recoProtonStartXAll->Fill(start_pos.x);
						recoProtonStartYAll->Fill(start_pos.y);
						recoProtonStartZAll->Fill(start_pos.z);
						recoProtonEndXAll->Fill(end_pos.x);
						recoProtonEndYAll->Fill(end_pos.y);
						recoProtonEndZAll->Fill(end_pos.z);
						recoProtonMomentumXAll->Fill(momentum.x);
						recoProtonMomentumYAll->Fill(momentum.y);
						recoProtonMomentumZAll->Fill(momentum.z);
						recoProtonLength_vs_EAll->Fill(length, energy);
						recoProton_startX_vs_startYAll->Fill(start_pos.x, start_pos.y);
						recoProton_startZ_vs_startXAll->Fill(start_pos.z, start_pos.x);
						recoProton_startZ_vs_startYAll->Fill(start_pos.z, start_pos.y);
						recoProton_startX_vs_endXAll->Fill(start_pos.x, end_pos.x);
						recoProton_startY_vs_endYAll->Fill(start_pos.y, end_pos.y);
						recoProton_startZ_vs_endZAll->Fill(start_pos.z, end_pos.z);
						recoProton_endX_vs_endYAll->Fill(end_pos.x, end_pos.y);
						recoProton_endZ_vs_endXAll->Fill(end_pos.z, end_pos.x);
						recoProton_endZ_vs_endYAll->Fill(end_pos.z, end_pos.y);						
						recoProton_length_vs_cosLAll->Fill(length, dirZ);
						recoProton_energy_vs_cosLAll->Fill(energy, dirZ);
						recoProton_dX_vs_lengthAll->Fill(length, dX);
						recoProton_dY_vs_lengthAll->Fill(length, dY);
						recoProton_dZ_vs_lengthAll->Fill(length, dZ);


					}

					if (abs(pdg) == 211) {
						auto start_pos = sr->common.ixn.dlp[nixn].part.dlp[npart].start;
						auto end_pos   = sr->common.ixn.dlp[nixn].part.dlp[npart].end;
						
						auto momentum   = sr->common.ixn.dlp[nixn].part.dlp[npart].p;

						double dX = end_pos.x - start_pos.x;
						double dY = end_pos.y - start_pos.y;
						double dZ = end_pos.z - start_pos.z;

						double length = TMath::Sqrt(dX * dX + dY * dY + dZ * dZ);
						double dirZ = dZ / length;

						// auto p       = sr->mc.nu[interactionNumber].prim[partNumber].p;

						auto energy = sr->common.ixn.dlp[nixn].part.dlp[npart].E;

						recoPionLengthAll->Fill(length);
						recoPionCosLAll->Fill(dirZ);
						recoPionEAll->Fill(energy);
						recoPionStartXAll->Fill(start_pos.x);
						recoPionStartYAll->Fill(start_pos.y);
						recoPionStartZAll->Fill(start_pos.z);
						recoPionEndXAll->Fill(end_pos.x);
						recoPionEndYAll->Fill(end_pos.y);
						recoPionEndZAll->Fill(end_pos.z);
						recoPionMomentumXAll->Fill(momentum.x);
						recoPionMomentumYAll->Fill(momentum.y);
						recoPionMomentumZAll->Fill(momentum.z);
						recoPionLength_vs_EAll->Fill(length, energy);
						recoPion_startX_vs_startYAll->Fill(start_pos.x, start_pos.y);
						recoPion_startZ_vs_startXAll->Fill(start_pos.z, start_pos.x);
						recoPion_startZ_vs_startYAll->Fill(start_pos.z, start_pos.y);
						recoPion_startX_vs_endXAll->Fill(start_pos.x, end_pos.x);
						recoPion_startY_vs_endYAll->Fill(start_pos.y, end_pos.y);
						recoPion_startZ_vs_endZAll->Fill(start_pos.z, end_pos.z);
						recoPion_endX_vs_endYAll->Fill(end_pos.x, end_pos.y);
						recoPion_endZ_vs_endXAll->Fill(end_pos.z, end_pos.x);
						recoPion_endZ_vs_endYAll->Fill(end_pos.z, end_pos.y);						
						recoPion_length_vs_cosLAll->Fill(length, dirZ);
						recoPion_energy_vs_cosLAll->Fill(energy, dirZ);
						recoPion_dX_vs_lengthAll->Fill(length, dX);
						recoPion_dY_vs_lengthAll->Fill(length, dY);
						recoPion_dZ_vs_lengthAll->Fill(length, dZ);

						// recoPionLengthAll->Fill(length);
						// recoPionCosLAll->Fill(dirZ);
						// recoPionEAll->Fill(energy);
						// recoPionStartXAll->Fill(start_pos.x);
						// recoPionStartYAll->Fill(start_pos.y);
						// recoPionStartZAll->Fill(start_pos.z);
						// recoPionEndXAll->Fill(end_pos.x);
						// recoPionEndYAll->Fill(end_pos.y);
						// recoPionEndZAll->Fill(end_pos.z);
						// recoPionMomentumXAll->Fill(momentum.x);
						// recoPionMomentumYAll->Fill(momentum.y);
						// recoPionMomentumZAll->Fill(momentum.z);
						// recoPionLengthAll_vs_E->Fill(length, energy);
					}


					if (abs(pdg) == 321) {
						auto start_pos = sr->common.ixn.dlp[nixn].part.dlp[npart].start;
						auto end_pos   = sr->common.ixn.dlp[nixn].part.dlp[npart].end;

						auto momentum   = sr->common.ixn.dlp[nixn].part.dlp[npart].p;

						double dX = end_pos.x - start_pos.x;
						double dY = end_pos.y - start_pos.y;
						double dZ = end_pos.z - start_pos.z;

						double length = TMath::Sqrt(dX * dX + dY * dY + dZ * dZ);
						double dirZ = dZ / length;

						// auto p       = sr->mc.nu[interactionNumber].prim[partNumber].p;

						auto energy = sr->common.ixn.dlp[nixn].part.dlp[npart].E;

						recoKaonLengthAll->Fill(length);
						recoKaonCosLAll->Fill(dirZ);
						recoKaonEAll->Fill(energy);
						recoKaonStartXAll->Fill(start_pos.x);
						recoKaonStartYAll->Fill(start_pos.y);
						recoKaonStartZAll->Fill(start_pos.z);
						recoKaonEndXAll->Fill(end_pos.x);
						recoKaonEndYAll->Fill(end_pos.y);
						recoKaonEndZAll->Fill(end_pos.z);
						recoKaonMomentumXAll->Fill(momentum.x);
						recoKaonMomentumYAll->Fill(momentum.y);
						recoKaonMomentumZAll->Fill(momentum.z);
					}
					
					if (abs(pdg) == 22 || abs(pdg) == 11 || abs(pdg) == 111) {
						auto truthSize = sr->common.ixn.dlp[nixn].part.dlp[npart].truth.size();
						double maxPartTruthOverlap = 0;
						
						for (int backTrack = 0; backTrack < truthSize; backTrack++) {
							int parType = sr->common.ixn.dlp[nixn].part.dlp[npart].truth[backTrack].type;
							int partNumber = sr->common.ixn.dlp[nixn].part.dlp[npart].truth[backTrack].part;
							int interactionNumber = sr->common.ixn.dlp[nixn].part.dlp[npart].truth[backTrack].ixn;

							double partTruthOverlap = sr->common.ixn.dlp[nixn].part.dlp[npart].truthOverlap[backTrack];

							if (maxPartTruthOverlap < partTruthOverlap) {
								maxPartTruthOverlap = partTruthOverlap;
							}
						}
						
						for (int backTrack = 0; backTrack < truthSize; backTrack++) {
							int parType = sr->common.ixn.dlp[nixn].part.dlp[npart].truth[backTrack].type;
							int partNumber = sr->common.ixn.dlp[nixn].part.dlp[npart].truth[backTrack].part;
							int interactionNumber = sr->common.ixn.dlp[nixn].part.dlp[npart].truth[backTrack].ixn;
							
							double partTruthOverlap = sr->common.ixn.dlp[nixn].part.dlp[npart].truthOverlap[backTrack];
							
							int backtracked = -9999;
							
							if (partTruthOverlap > maxPartTruthOverlap) {
								
								if (parType < 3) {
									backtracked = sr->mc.nu[interactionNumber].prim[partNumber].pdg;
									primaryTrkIndex.push_back(partNumber);
								}
								else
									backtracked = sr->mc.nu[interactionNumber].sec[partNumber].pdg;

								correctShower = 0;
								
								if (parType > 2)
									correctShower = 3;
								else if (backtracked == 11 || backtracked == 22 || backtracked == 111 || backtracked == 2112) {
									correctShower = 1;
									
								}
							}
							else backtracked = -1;
						}
						
						if (sr->common.ixn.dlp[nixn].part.dlp[npart].primary==true) showerCorrectness->Fill(correctShower);
					}
					
					
					int maxPartNumber = -999;
					int maxTypeNumber = -999;	
					
					if (passesMinerva && goodInteraction) {
						
						
						// === Reco proton kinematics for matched interactions (no truth match required) ===
						if (abs(pdg) == 2212) {
							// Reco start/end from DLp
							auto r_start = sr->common.ixn.dlp[nixn].part.dlp[npart].start;
							auto r_end   = sr->common.ixn.dlp[nixn].part.dlp[npart].end;

							double r_energy = sr->common.ixn.dlp[nixn].part.dlp[npart].E;

							// Reco geometric differences
							double r_dX = r_end.x - r_start.x;
							double r_dY = r_end.y - r_start.y;
							double r_dZ = r_end.z - r_start.z;

							// Reco track length
							double r_length = TMath::Sqrt(r_dX*r_dX + r_dY*r_dY + r_dZ*r_dZ);

							// Reco cos(theta_Z)
							double r_cosZ = (r_length > 0) ? r_dZ / r_length : -999.0;

							// Momentum from DLp
							auto r_p = sr->common.ixn.dlp[nixn].part.dlp[npart].p;
														
							// if (sr->common.ixn.dlp[nixn].part.dlp[npart].contained == true) continue;			// skip if not contained

							nRecoMatchedProtonAll++;

							bool contained_geom =
								(r_start.x > -55 && r_start.x < 55) &&
								(r_start.y > -55 && r_start.y < 55) &&
								(r_start.z > -55 && r_start.z < 55) &&
								(r_end.x   > -55 && r_end.x   < 55) &&
								(r_end.y   > -55 && r_end.y   < 55) &&
								(r_end.z   > -55 && r_end.z   < 55);

							// bool contained_geom =	sr->common.ixn.dlp[nixn].part.dlp[npart].contained;

							if (!contained_geom) continue;   // skip non-contained

                            // if (r_cosZ <= 0.8) continue; // keep protons with cos>0.8
                            // if (r_cosZ > 0.8) continue; // keep protons with cos<=0.8    

                            // if (r_start.z <= 0.0) continue;  // keep only protons starting at z > 0

                            // if (r_length < 30.0 || r_length > 70.0) continue; // select protons with 30–70 cm track length


                            // if (length <= 2.0) continue;   // require proton track length > 2 cm

                            if (length <= 3.0) continue;   // require proton track length > 3 cm
                            // if (length <= 3.20) continue;   // require proton track length > 2 cm

                            // if (r_end.z < -7.0 || r_end.z > 7.0) continue;  // keep only ±7 cm around z = 0
                            // if (r_end.z >= -5.0 && r_end.z <= 5.0) continue; // reject |z_end| ≤ 5 cm → keep outside 
                            if (r_end.z >= -5.0 && r_end.z <= 0.0) continue; //  reject only −z side (−5 to 0 cm):
                            if (r_cosZ <= -0.95) continue;  // rejects cosθ ≤ −0.95, keeps cosθ > −0.95

                            // if (r_cosZ < -0.2 || r_cosZ > 0.2) continue; // Keep only −0.2 ≤ cosZ ≤ +0.2  (transverse protons)
                            // if (r_cosZ < -0.05 || r_cosZ > 0.05) continue; // Keep only −0.05 ≤ cosZ ≤ +0.05  (transverse protons)

                            
							// if (contained_geom) continue;    // skip contained

							nRecoProtons_matched++;

                            // After nRecoProtons_matched++, add:
                            int truePDGBin = 3; // default: other
                            // find best truth match
                            double maxOvlp = 0;
                            int backtrackedPDG = -1;
                            for (unsigned bt = 0; bt < sr->common.ixn.dlp[nixn].part.dlp[npart].truth.size(); ++bt) {
                                double ovlp = sr->common.ixn.dlp[nixn].part.dlp[npart].truthOverlap[bt];
                                if (ovlp > maxOvlp) {
                                    maxOvlp = ovlp;
                                    int parType = sr->common.ixn.dlp[nixn].part.dlp[npart].truth[bt].type;
                                    int partNum = sr->common.ixn.dlp[nixn].part.dlp[npart].truth[bt].part;
                                    int ixnNum  = sr->common.ixn.dlp[nixn].part.dlp[npart].truth[bt].ixn;
                                    if (parType < 3)
                                        backtrackedPDG = sr->mc.nu[ixnNum].prim[partNum].pdg;
                                    else
                                        backtrackedPDG = sr->mc.nu[ixnNum].sec[partNum].pdg;
                                }
                            }
                            if      (abs(backtrackedPDG) == 13)   truePDGBin = 0; // muon
                            else if (abs(backtrackedPDG) == 2212) truePDGBin = 1; // proton
                            else if (abs(backtrackedPDG) == 211)  truePDGBin = 2; // pion
                            
                            recoProtonTrueID->Fill(truePDGBin); // new 1D histogram: what are reco protons really?

                            // G4RW reco proton fill
                            {
                                int reco_proton_G4ID = -1;
                                int reco_proton_truePDG = 0;   // NEW: signed PDG of the same winning truth candidate
                                double reco_proton_maxOverlap = 0;
                                auto reco_proton_truthSize = sr->common.ixn.dlp[nixn].part.dlp[npart].truth.size();
                                for (int bt = 0; bt < (int)reco_proton_truthSize; bt++) {
                                    int parType    = sr->common.ixn.dlp[nixn].part.dlp[npart].truth[bt].type;
                                    int partNumber = sr->common.ixn.dlp[nixn].part.dlp[npart].truth[bt].part;
                                    double overlap = sr->common.ixn.dlp[nixn].part.dlp[npart].truthOverlap[bt];
                                    if (overlap > reco_proton_maxOverlap && parType < 3) {
                                        reco_proton_maxOverlap = overlap;
                                        reco_proton_G4ID       = sr->mc.nu[biggestMatchIndex].prim[partNumber].G4ID;
                                        reco_proton_truePDG    = sr->mc.nu[biggestMatchIndex].prim[partNumber].pdg;   // NEW
                                    }
                                }
                                g4rw_nu_id               = (long long)sr->mc.nu[biggestMatchIndex].id;
                                g4rw_fileIdx             = caf_chain->GetTreeNumber();
                                g4rw_G4ID                = reco_proton_G4ID;
                                g4rw_truePDG             = reco_proton_truePDG;   // NEW: save it as its own branch
                                g4rw_proton_reco_cosTheta = r_cosZ;
                                g4rw_proton_reco_length   = r_length;
                                syst_tree_g4rw_proton_reco->Fill();
                            }
                            // g4rw end

							// Fill 1D histograms
							reco_matched_proton_length->Fill(r_length);
							reco_matched_proton_cosL->Fill(r_cosZ);
							reco_matched_proton_energy->Fill(r_energy);

                            // Interaction-type fills
                            {
                                int cc_mode   = sr->mc.nu[biggestMatchIndex].iscc;
                                int mode_mode = sr->mc.nu[biggestMatchIndex].mode;
                                if (cc_mode == 1) {
                                    if      (mode_mode == 1 || mode_mode == 1001) { reco_matched_proton_length_QE->Fill(r_length);  reco_matched_proton_cosL_QE->Fill(r_cosZ);  }
                                    else if (mode_mode == 10)                     { reco_matched_proton_length_MEC->Fill(r_length); reco_matched_proton_cosL_MEC->Fill(r_cosZ); }
                                    else if (mode_mode == 4)                      { reco_matched_proton_length_RES->Fill(r_length); reco_matched_proton_cosL_RES->Fill(r_cosZ); }
                                    else if (mode_mode == 3)                      { reco_matched_proton_length_DIS->Fill(r_length); reco_matched_proton_cosL_DIS->Fill(r_cosZ); }
                                    else if (mode_mode == 5)                      { reco_matched_proton_length_COH->Fill(r_length); reco_matched_proton_cosL_COH->Fill(r_cosZ); }
                                }
                            }

                            // Particle-type fills (truePDGBin already computed above: 0=mu, 1=p, 2=pi, 3=other)
                            if      (truePDGBin == 0) { reco_matched_proton_length_trueMuon->Fill(r_length);   reco_matched_proton_cosL_trueMuon->Fill(r_cosZ);   }
                            else if (truePDGBin == 1) { reco_matched_proton_length_trueProton->Fill(r_length); reco_matched_proton_cosL_trueProton->Fill(r_cosZ); }
                            else if (truePDGBin == 2) { reco_matched_proton_length_truePion->Fill(r_length);   reco_matched_proton_cosL_truePion->Fill(r_cosZ);   }
                            else                      { reco_matched_proton_length_trueOther->Fill(r_length);  reco_matched_proton_cosL_trueOther->Fill(r_cosZ);  }
                            
							// for nusyst start
							out_genieIdx = sr->mc.nu[biggestMatchIndex].genieIdx;
							out_recoProtonCosL = r_cosZ;
                            out_recoProtonLength = r_length;
							out_cafFileIdx = caf_chain->GetTreeNumber();
							syst_tree_recoProtonCosL->Fill();
                            syst_tree_recoProtonLength->Fill();
							// for nusyst end

							reco_matched_proton_startX->Fill(r_start.x);
							reco_matched_proton_startY->Fill(r_start.y);
							reco_matched_proton_startZ->Fill(r_start.z);

							reco_matched_proton_endX->Fill(r_end.x);
							reco_matched_proton_endY->Fill(r_end.y);
							reco_matched_proton_endZ->Fill(r_end.z);

							reco_matched_proton_momentumX->Fill(r_p.x);
							reco_matched_proton_momentumY->Fill(r_p.y);
							reco_matched_proton_momentumZ->Fill(r_p.z);

							// Fill 2D correlations
							reco_matched_proton_length_vs_energy->Fill(r_length, r_energy);

							reco_matched_proton_startX_vs_startY->Fill(r_start.x, r_start.y);
							reco_matched_proton_startZ_vs_startX->Fill(r_start.z, r_start.x);
							reco_matched_proton_startZ_vs_startY->Fill(r_start.z, r_start.y);

							reco_matched_proton_startX_vs_endX->Fill(r_start.x, r_end.x);
							reco_matched_proton_startY_vs_endY->Fill(r_start.y, r_end.y);
							reco_matched_proton_startZ_vs_endZ->Fill(r_start.z, r_end.z);

							reco_matched_proton_length_vs_cosL->Fill(r_length, r_cosZ);
							reco_matched_proton_energy_vs_cosL->Fill(r_energy, r_cosZ);

							reco_matched_proton_dX_vs_length->Fill(r_length, r_dX);
							reco_matched_proton_dY_vs_length->Fill(r_length, r_dY);
							reco_matched_proton_dZ_vs_length->Fill(r_length, r_dZ);

							reco_matched_proton_endX_vs_endY->Fill(r_end.x, r_end.y);
							reco_matched_proton_endZ_vs_endX->Fill(r_end.z, r_end.x);
							reco_matched_proton_endZ_vs_endY->Fill(r_end.z, r_end.y);
							
							reco_matched_proton_cosL_vs_startZ->Fill(r_cosZ, r_start.z);
							reco_matched_proton_cosL_vs_endZ->Fill(r_cosZ, r_end.z);
							
							if (contained_geom) {
								nRecoProtonMatchedContained++;
								
								recoProtonMatchedLengthContained->Fill(r_length);
								recoProtonMatchedEContained->Fill(r_energy);
								recoProtonMatchedLength_vs_EContained->Fill(r_length, r_energy);
							}

							if (!contained_geom) {
								nRecoProtonMatchedNonContained++;
								
								recoProtonMatchedLengthNonContained->Fill(r_length);
								recoProtonMatchedENonContained->Fill(r_energy);
								recoProtonMatchedLength_vs_ENonContained->Fill(r_length, r_energy);
							}
							
						}
						

						
						
												// === Reco proton kinematics for matched interactions (no truth match required - NOT SURE ABOUT THIS. I THINK THESE ARE FROM MATHCED INTERACTIONS) ===
						if (abs(pdg) == 211) {
							// Reco start/end from DLp
							auto r_start = sr->common.ixn.dlp[nixn].part.dlp[npart].start;
							auto r_end   = sr->common.ixn.dlp[nixn].part.dlp[npart].end;

							double r_energy = sr->common.ixn.dlp[nixn].part.dlp[npart].E;

							// Reco geometric differences
							double r_dX = r_end.x - r_start.x;
							double r_dY = r_end.y - r_start.y;
							double r_dZ = r_end.z - r_start.z;

							// Reco track length
							double r_length = TMath::Sqrt(r_dX*r_dX + r_dY*r_dY + r_dZ*r_dZ);

							// Reco cos(theta_Z)
							double r_cosZ = (r_length > 0) ? r_dZ / r_length : -999.0;

							// Momentum from DLp
							auto r_p = sr->common.ixn.dlp[nixn].part.dlp[npart].p;
														
							// if (sr->common.ixn.dlp[nixn].part.dlp[npart].contained == true) continue;			// skip if not contained

							// Charged pion mass [GeV]
							const double mpi = 0.139570;

							// Reco momentum magnitude [GeV]
							const double r_p2   = r_p.x*r_p.x + r_p.y*r_p.y + r_p.z*r_p.z;
							const double r_pMag = (r_p2 > 0.0) ? std::sqrt(r_p2) : 0.0;

							// Reco kinetic energy from momentum [GeV]
							const double r_energy_fromP = (r_p2 > 0.0) ? (std::sqrt(r_p2 + mpi*mpi) - mpi) : -999.0;


							nRecoMatchedPionAll++;

							bool contained_geom =
								(r_start.x > -55 && r_start.x < 55) &&
								(r_start.y > -55 && r_start.y < 55) &&
								(r_start.z > -55 && r_start.z < 55) &&
								(r_end.x   > -55 && r_end.x   < 55) &&
								(r_end.y   > -55 && r_end.y   < 55) &&
								(r_end.z   > -55 && r_end.z   < 55);

							// bool contained_geom =	sr->common.ixn.dlp[nixn].part.dlp[npart].contained;

							// if (!contained_geom) continue;   // skip non-contained

                            // if (r_length > 20.0) continue; // reject pions longer than 20 cm
                            // if (r_length <= 20.0) continue; // reject pions shorter than 20 cm
							
							// if (contained_geom) continue;    // skip contained

							nRecoPions_matched++;

                            // for G4RW reco start
int reco_matched_G4ID = -1;
double reco_maxOverlap = 0;
auto reco_truthSize = sr->common.ixn.dlp[nixn].part.dlp[npart].truth.size();
for (int bt = 0; bt < (int)reco_truthSize; bt++) {
    int parType    = sr->common.ixn.dlp[nixn].part.dlp[npart].truth[bt].type;
    int partNumber = sr->common.ixn.dlp[nixn].part.dlp[npart].truth[bt].part;
    double overlap = sr->common.ixn.dlp[nixn].part.dlp[npart].truthOverlap[bt];
    if (overlap > reco_maxOverlap && parType < 3) {
        reco_maxOverlap   = overlap;
        reco_matched_G4ID = sr->mc.nu[biggestMatchIndex].prim[partNumber].G4ID;
    }
}
g4rw_nu_id         = (long long)sr->mc.nu[biggestMatchIndex].id;
g4rw_fileIdx       = caf_chain->GetTreeNumber();
g4rw_G4ID          = reco_matched_G4ID;
g4rw_reco_cosTheta = r_cosZ;
syst_tree_g4rw_pion_reco->Fill();
// for G4RW reco end

							// Fill 1D histograms
							reco_matched_pion_length->Fill(r_length);
							reco_matched_pion_cosL->Fill(r_cosZ);
							reco_matched_pion_energy->Fill(r_energy);
							reco_matched_pion_energy_fromP->Fill(r_energy_fromP);



							reco_matched_pion_length_vs_energy_fromP->Fill(r_length, r_energy_fromP);
							reco_matched_pion_energy_diff->Fill(r_energy - r_energy_fromP);

							reco_matched_pion_startX->Fill(r_start.x);
							reco_matched_pion_startY->Fill(r_start.y);
							reco_matched_pion_startZ->Fill(r_start.z);

							reco_matched_pion_endX->Fill(r_end.x);
							reco_matched_pion_endY->Fill(r_end.y);
							reco_matched_pion_endZ->Fill(r_end.z);

							reco_matched_pion_momentumX->Fill(r_p.x);
							reco_matched_pion_momentumY->Fill(r_p.y);
							reco_matched_pion_momentumZ->Fill(r_p.z);

							// Fill 2D correlations
							reco_matched_pion_length_vs_energy->Fill(r_length, r_energy);

							reco_matched_pion_startX_vs_startY->Fill(r_start.x, r_start.y);
							reco_matched_pion_startZ_vs_startX->Fill(r_start.z, r_start.x);
							reco_matched_pion_startZ_vs_startY->Fill(r_start.z, r_start.y);

							reco_matched_pion_startX_vs_endX->Fill(r_start.x, r_end.x);
							reco_matched_pion_startY_vs_endY->Fill(r_start.y, r_end.y);
							reco_matched_pion_startZ_vs_endZ->Fill(r_start.z, r_end.z);

							reco_matched_pion_length_vs_cosL->Fill(r_length, r_cosZ);
							reco_matched_pion_energy_vs_cosL->Fill(r_energy, r_cosZ);

							reco_matched_pion_dX_vs_length->Fill(r_length, r_dX);
							reco_matched_pion_dY_vs_length->Fill(r_length, r_dY);
							reco_matched_pion_dZ_vs_length->Fill(r_length, r_dZ);

							reco_matched_pion_endX_vs_endY->Fill(r_end.x, r_end.y);
							reco_matched_pion_endZ_vs_endX->Fill(r_end.z, r_end.x);
							reco_matched_pion_endZ_vs_endY->Fill(r_end.z, r_end.y);
							
							reco_matched_pion_cosL_vs_startZ->Fill(r_cosZ, r_start.z);
							reco_matched_pion_cosL_vs_endZ->Fill(r_cosZ, r_end.z);
							
							if (contained_geom) {
								nRecoPionMatchedContained++;
								
								recoPionMatchedLengthContained->Fill(r_length);
								recoPionMatchedEContained->Fill(r_energy);
								recoPionMatchedLength_vs_EContained->Fill(r_length, r_energy);
							}

							if (!contained_geom) {
								nRecoPionMatchedNonContained++;
								
								recoPionMatchedLengthNonContained->Fill(r_length);
								recoPionMatchedENonContained->Fill(r_energy);
								recoPionMatchedLength_vs_ENonContained->Fill(r_length, r_energy);
							}
							
						}

						
						
						
						if (abs(pdg) == 2212 || abs(pdg) == 13 || abs(pdg) == 211 || abs(pdg) == 321) {
							auto truthSize = sr->common.ixn.dlp[nixn].part.dlp[npart].truth.size();
							double maxPartTruthOverlap = 0;
							
							for (int backTrack = 0; backTrack < truthSize; backTrack++) {
								int parType = sr->common.ixn.dlp[nixn].part.dlp[npart].truth[backTrack].type;
								int partNumber = sr->common.ixn.dlp[nixn].part.dlp[npart].truth[backTrack].part;
								int interactionNumber = sr->common.ixn.dlp[nixn].part.dlp[npart].truth[backTrack].ixn;

								double partTruthOverlap = sr->common.ixn.dlp[nixn].part.dlp[npart].truthOverlap[backTrack];

								// std::cout << "Truth has overlap of: " << partTruthOverlap << std::endl;

								if (maxPartTruthOverlap < partTruthOverlap) {
									maxPartTruthOverlap = partTruthOverlap;
								}
							}
							
							for (int backTrack=0; backTrack<truthSize; backTrack++){
								int parType=sr->common.ixn.dlp[nixn].part.dlp[npart].truth[backTrack].type;
								int partNumber=sr->common.ixn.dlp[nixn].part.dlp[npart].truth[backTrack].part;
								int interactionNumber=sr->common.ixn.dlp[nixn].part.dlp[npart].truth[backTrack].ixn;
								double  partTruthOverlap=sr->common.ixn.dlp[nixn].part.dlp[npart].truthOverlap[backTrack];
								
								if (interactionNumber != biggestMatchIndex) continue;
								
								if (partTruthOverlap < maxPartTruthOverlap) continue;


								int backtracked=-9999;
								
								if (partTruthOverlap > 0.5) {
									maxPartNumber = partNumber;
									maxTypeNumber = parType;
									
									if (parType < 3)
										backtracked = sr->mc.nu[interactionNumber].prim[partNumber].pdg;
									else
										backtracked = sr->mc.nu[interactionNumber].sec[partNumber].pdg;

									if (parType < 3) {
										int pdgNumber = 0;

										if (abs(pdg) == 13) {
											pdgNumber = 0;
										} else if (abs(pdg) == 2212) {
											pdgNumber = 1;
										} else if (abs(pdg) == 211) {
											pdgNumber = 2;
										} else {
											pdgNumber = 3;
										}

										int backtrackedPDG = 0;
										
										if (abs(backtracked) == 13)
											backtrackedPDG = 0;
										else if (abs(backtracked) == 2212) {
											backtrackedPDG = 1;
										} else if (abs(backtracked) == 211) {
											backtrackedPDG = 2;
										} else {
											backtrackedPDG = 3;
										}

										confusionMatrix->Fill(pdgNumber, backtrackedPDG);


										auto start_pos = sr->mc.nu[interactionNumber].prim[partNumber].start_pos;
										// if (std::isnan(start_pos.z)) continue;

										auto p       = sr->mc.nu[interactionNumber].prim[partNumber].p;
										auto end_pos = sr->mc.nu[interactionNumber].prim[partNumber].end_pos;

										double dX = p.px; // (end_pos.x - start_pos.x);
										double dY = p.py; // (end_pos.y - start_pos.y);
										double dZ = p.pz; // (end_pos.z - start_pos.z);

										double length = TMath::Sqrt(dX * dX + dY * dY + dZ * dZ);

										double dXLen = end_pos.x - start_pos.x;
										double dYLen = end_pos.y - start_pos.y;
										double dZLen = end_pos.z - start_pos.z;

										double lengthPos = TMath::Sqrt(dXLen * dXLen + dYLen * dYLen + dZLen * dZLen);

										double dirX = dX / length;
										double dirY = dY / length;
										double dirZ = dZ / length;
										
										// --- Direction–chord alignment (cos alpha) ---
										double px = p.px;
										double py = p.py;
										double pz = p.pz;

										double pMag = TMath::Sqrt(px*px + py*py + pz*pz);
										double lMag = lengthPos;

										double cosAlpha = -999.0;
										if (pMag > 0 && lMag > 0) {
											cosAlpha = (dXLen*px + dYLen*py + dZLen*pz) / (lMag * pMag);
										}
										
										int cc   = sr->mc.nu[biggestMatchIndex].iscc;
										int mode = sr->mc.nu[biggestMatchIndex].mode;
										
										bool truth_contained_geom =
											(start_pos.x > -55 && start_pos.x < 55) &&
											(start_pos.y > -55 && start_pos.y < 55) &&
											(start_pos.z > -55 && start_pos.z < 55) &&
											(end_pos.x   > -55 && end_pos.x   < 55) &&
											(end_pos.y   > -55 && end_pos.y   < 55) &&
											(end_pos.z   > -55 && end_pos.z   < 55);

										

										// bool truth_contained_geom =
											// (start_pos.x > -59.0 && start_pos.x < 59.0) &&
											// (end_pos.x   > -59.0 && end_pos.x   < 59.0) &&
											// (start_pos.y > -57.0 && start_pos.y < 57.0) &&
											// (end_pos.y   > -57.0 && end_pos.y   < 57.0) &&
											// (start_pos.z > -59.5 && start_pos.z < 59.5) &&
											// (end_pos.z   > -59.5 && end_pos.z   < 59.5);

											// if (truth_contained_geom)
												// confusionMatrix_truthContained->Fill(pdgNumber, backtrackedPDG);
											// else
												// confusionMatrix_truthNotContained->Fill(pdgNumber, backtrackedPDG);


										// if (!truth_contained_geom) continue;	// skip non-contained

										// if (truth_contained_geom) continue;	// skip contained

										// if ( (p.E - 0.938272) > 0.5) continue;


										if (!std::isnan(start_pos.z) &&
											(backtracked == 2212 /*|| abs(backtracked) == 13*/ || abs(backtracked) == 211) &&
											parType < 3)
										{
											selTrkLen->Fill(lengthPos);

											if (cc == 1 && (mode == 1 || mode == 1001 || mode == 10))
												selTrkLenQE->Fill(lengthPos);

											if (cc == 1 && (mode == 3 || mode == 4))
												selTrkLenNonQE->Fill(lengthPos);
										}
										

										// if (!(passesMinerva && goodInteraction)) continue;  // Only count multiplicity for MINERvA-passed, truth-matched (good) interactions
	
										if (sr->mc.nu[biggestMatchIndex].targetPDG == 1000180400 &&
											backtracked == 2212 && parType < 3)
										{

                                            // if (lengthPos > 80.0) continue;   // require proton track length < 80 cm

											
											nTrueMatchedProtonAll++;
											nPrimProtons_true++;

                                            // G4RW truth proton fill
                                            g4rw_nu_id           = (long long)sr->mc.nu[biggestMatchIndex].id;
                                            g4rw_fileIdx         = caf_chain->GetTreeNumber();
                                            g4rw_G4ID            = sr->mc.nu[interactionNumber].prim[partNumber].G4ID;
                                            g4rw_proton_cosTheta = dirZ;
                                            g4rw_proton_length   = lengthPos;
                                            syst_tree_g4rw_proton->Fill();
											
											
											selProtonE->Fill(sr->mc.nu[interactionNumber].prim[partNumber].p.E - 0.938272);
											selProtonDirX->Fill(dirX);
											selProtonDirZ->Fill(dirZ);
											selProtonDirY->Fill(dirY);
											
											recoProtonLength->Fill(lengthPos);      // geometric length from start to end
											recoProtonCosL->Fill(dirZ);             // cos(θ) = Z-component of direction
											recoProtonE->Fill(p.E - 0.938272); // or .p.E
											recoProtonStartX->Fill(start_pos.x);
											recoProtonStartY->Fill(start_pos.y);
											recoProtonStartZ->Fill(start_pos.z);
											recoProtonEndX->Fill(end_pos.x);
											recoProtonEndY->Fill(end_pos.y);
											recoProtonEndZ->Fill(end_pos.z);
											recoProtonMomentumX->Fill(p.px);
											recoProtonMomentumY->Fill(p.py);
											recoProtonMomentumZ->Fill(p.pz);
											recoProtonLength_vs_E->Fill(lengthPos, p.E - 0.938272);
											recoProton_startX_vs_startY->Fill(start_pos.x, start_pos.y);
											recoProton_startZ_vs_startX->Fill(start_pos.z, start_pos.x);
											recoProton_startZ_vs_startY->Fill(start_pos.z, start_pos.y);
											recoProton_startX_vs_endX->Fill(start_pos.x, end_pos.x);
											recoProton_startY_vs_endY->Fill(start_pos.y, end_pos.y);
											recoProton_startZ_vs_endZ->Fill(start_pos.z, end_pos.z);
											recoProton_endX_vs_endY->Fill(end_pos.x, end_pos.y);
											recoProton_endZ_vs_endX->Fill(end_pos.z, end_pos.x);
											recoProton_endZ_vs_endY->Fill(end_pos.z, end_pos.y);
											recoProton_length_vs_cosL->Fill(lengthPos, dirZ);
											recoProton_energy_vs_cosL->Fill(p.E - 0.938272, dirZ);
											recoProton_dX_vs_length->Fill(lengthPos, dXLen);
											recoProton_dY_vs_length->Fill(lengthPos, dYLen);
											recoProton_dZ_vs_length->Fill(lengthPos, dZLen);
											
											recoProton_cosL_vs_startz->Fill(dirZ, start_pos.z);
											recoProton_cosL_vs_endz->Fill(dirZ, end_pos.z);
											
											recoProtonCosAlpha->Fill(cosAlpha);
											recoProtonCosAlpha_vs_E->Fill(p.E - 0.938272, cosAlpha);
											recoProtonCosAlpha_vs_length->Fill(lengthPos, cosAlpha);
											
											if (truth_contained_geom) {
												nTrueProtonMatchedContained++;
												
												trueProtonMatchedLengthContained->Fill(lengthPos);
												trueProtonMatchedEContained->Fill(p.E - 0.938272);
												trueProtonMatchedLength_vs_EContained->Fill(lengthPos, p.E - 0.938272);
											}

											if (!truth_contained_geom) {
												nTrueProtonMatchedNonContained++;
												
												trueProtonMatchedLengthNonContained->Fill(lengthPos);
												trueProtonMatchedENonContained->Fill(p.E - 0.938272);
												trueProtonMatchedLength_vs_ENonContained->Fill(lengthPos, p.E - 0.938272);
											}
											
											
											if (!goodInteraction) {
												if (rock == 1)
													recoProtonLength_Rock->Fill(lengthPos);
												else
													recoProtonLength_Sec->Fill(lengthPos);  // <--- Secondary interaction
											}
											else {
												if (cc == 1 && mode == 1)
													recoProtonLength_QE->Fill(lengthPos);
												else if (cc == 1 && mode == 10)
													recoProtonLength_MEC->Fill(lengthPos);
												else if (cc == 1 && mode == 4)
													recoProtonLength_RES->Fill(lengthPos);
												else if (cc == 1 && mode == 3)
													recoProtonLength_DIS->Fill(lengthPos);
												else if (cc == 1 && mode == 5)
													recoProtonLength_COH->Fill(lengthPos);
											}										
										}
										
										// This commented part is yet to fix
										if (sr->mc.nu[biggestMatchIndex].targetPDG == 1000180400 &&
											backtracked == 13 && parType < 3)
										{
											
											recoMuonLength->Fill(lengthPos);      // geometric length from start to end
											recoMuonCosL->Fill(dirZ);             // cos(θ) = Z-component of direction
											recoMuonE->Fill(p.E);
											recoMuonStartX->Fill(start_pos.x);
											recoMuonStartY->Fill(start_pos.y);
											recoMuonStartZ->Fill(start_pos.z);
											recoMuonLength_vs_E->Fill(lengthPos, p.E);

											if (!goodInteraction) {
												if (rock == 1)
													recoMuonLength_Rock->Fill(lengthPos);
												else
													recoMuonLength_Sec->Fill(lengthPos);  // <--- Secondary interaction
											}
											else {
												if (cc == 1 && mode == 1)
													recoMuonLength_QE->Fill(lengthPos);
												else if (cc == 1 && mode == 10)
													recoMuonLength_MEC->Fill(lengthPos);
												else if (cc == 1 && mode == 4)
													recoMuonLength_RES->Fill(lengthPos);
												else if (cc == 1 && mode == 3)
													recoMuonLength_DIS->Fill(lengthPos);
												else if (cc == 1 && mode == 5)
													recoMuonLength_COH->Fill(lengthPos);
											}
										}
										
										

										// if (passesMinerva && goodInteraction) {
											if (sr->mc.nu[biggestMatchIndex].targetPDG == 1000180400 && abs(backtracked) == 211 && parType < 3)
											{
												
												nTrueMatchedPionAll++; 
												nPrimPions_true++;
												
												selPionE->Fill(sr->mc.nu[interactionNumber].prim[partNumber].p.E - 0.13957039);
												selPionDirX->Fill(dirX);
												selPionDirZ->Fill(dirZ);
												selPionDirY->Fill(dirY);
												
												recoPionLength->Fill(lengthPos);      // geometric length from start to end
												recoPionCosL->Fill(dirZ);             // cos(θ) = Z-component of direction
												recoPionE->Fill(p.E - 0.139570); // or .p.E
												recoPionStartX->Fill(start_pos.x);
												recoPionStartY->Fill(start_pos.y);
												recoPionStartZ->Fill(start_pos.z);
												recoPionEndX->Fill(end_pos.x);
												recoPionEndY->Fill(end_pos.y);
												recoPionEndZ->Fill(end_pos.z);
												recoPionMomentumX->Fill(p.px);
												recoPionMomentumY->Fill(p.py);
												recoPionMomentumZ->Fill(p.pz);
												recoPionLength_vs_E->Fill(lengthPos, p.E - 0.139570);
												recoPion_startX_vs_startY->Fill(start_pos.x, start_pos.y);
												recoPion_startZ_vs_startX->Fill(start_pos.z, start_pos.x);
												recoPion_startZ_vs_startY->Fill(start_pos.z, start_pos.y);
												recoPion_startX_vs_endX->Fill(start_pos.x, end_pos.x);
												recoPion_startY_vs_endY->Fill(start_pos.y, end_pos.y);
												recoPion_startZ_vs_endZ->Fill(start_pos.z, end_pos.z);
												recoPion_endX_vs_endY->Fill(end_pos.x, end_pos.y);
												recoPion_endZ_vs_endX->Fill(end_pos.z, end_pos.x);
												recoPion_endZ_vs_endY->Fill(end_pos.z, end_pos.y);
												recoPion_length_vs_cosL->Fill(lengthPos, dirZ);
												recoPion_energy_vs_cosL->Fill(p.E - 0.139570, dirZ);
												recoPion_dX_vs_length->Fill(lengthPos, dXLen);
												recoPion_dY_vs_length->Fill(lengthPos, dYLen);
												recoPion_dZ_vs_length->Fill(lengthPos, dZLen);
												
												recoPion_cosL_vs_startz->Fill(dirZ, start_pos.z);
												recoPion_cosL_vs_endz->Fill(dirZ, end_pos.z);
												
												recoPionCosAlpha->Fill(cosAlpha);
												recoPionCosAlpha_vs_E->Fill(p.E - 0.139570, cosAlpha);
												recoPionCosAlpha_vs_length->Fill(lengthPos, cosAlpha);
												
												if (truth_contained_geom) {
													nTruePionMatchedContained++;
													
													truePionMatchedLengthContained->Fill(lengthPos);
													truePionMatchedEContained->Fill(p.E - 0.139570);
													truePionMatchedLength_vs_EContained->Fill(lengthPos, p.E - 0.139570);
												}

												if (!truth_contained_geom) {
													nTruePionMatchedNonContained++;
													
													truePionMatchedLengthNonContained->Fill(lengthPos);
													truePionMatchedENonContained->Fill(p.E - 0.139570);
													truePionMatchedLength_vs_ENonContained->Fill(lengthPos, p.E - 0.139570);
												}
												
												
												// recoPionLength->Fill(lengthPos);      // geometric length from start to end
												// recoPionCosL->Fill(dirZ);             // cos(θ) = Z-component of direction
												// recoPionE->Fill(p.E); // or .p.E
												// recoPionStartX->Fill(start_pos.x);
												// recoPionStartY->Fill(start_pos.y);
												// recoPionStartZ->Fill(start_pos.z);
												// recoPionLength_vs_E->Fill(lengthPos, p.E);

											}

											if (sr->mc.nu[biggestMatchIndex].targetPDG == 1000180400 && abs(backtracked) == 321 && parType < 3)
											{
												selKaonE->Fill(sr->mc.nu[interactionNumber].prim[partNumber].p.E - 0.493677);
												selKaonDirX->Fill(dirX);
												selKaonDirZ->Fill(dirZ);
												selKaonDirY->Fill(dirY);
											}
										// }
									}
									
									correctTrack = 0;

									// std::cout << "Backtrackilng for " << parType << " : " << pdg << "," << backtracked << std::endl;

									if (parType > 2)
										correctTrack = 3;
									else if (backtracked == 11 || backtracked == 22 || backtracked == 111 || backtracked == 2112) {
										correctTrack = 1;
									}
								}
								
								if (partTruthOverlap < 0.5)
									correctTrack = 2;
								else
									backtracked = -1;							
							}
						}
						
						// if (nPrimProtons_true > 0) {
							// primProtonMultiplicity->Fill(nPrimProtons_true);
						// }
						
						// if (nTrueMatchedProtonAll > 0) {
							// primProtonMultiplicity->Fill(nTrueMatchedProtonAll);
						// }
						
						
						
						// if (nRecoProtons_matched > 0) {
							// recoProtonMultiplicity_matched->Fill(nRecoProtons_matched);
						// }
						
						// if (nPrimPions_true > 0) {
							// primPionMultiplicity->Fill(nPrimPions_true);
						// }
						
						// if (nRecoPions_matched > 0) {
							// recoPionMultiplicity_matched->Fill(nRecoPions_matched);
						// }
						
					}
				}
				
				// AFTER particle loop — once per interaction
				if (nRecoProtons_matched > 0) {
					recoProtonMultiplicity_matched->Fill(nRecoProtons_matched);
				}
				
				if (nPrimProtons_true > 0) {
					primProtonMultiplicity->Fill(nPrimProtons_true);
				}
				
				if (nPrimPions_true > 0) {
					primPionMultiplicity->Fill(nPrimPions_true);
				}
				
				if (nRecoPions_matched > 0) {
					recoPionMultiplicity_matched->Fill(nRecoPions_matched);
				}
				
				int trueMatchMult = 0;

				if (trackMult > 0 && trackMultExit > 0 /*&& oneContained && oneNotContained*/) {
					// recoCosL->Fill(dirZExiting);

					if (goodInteraction) {
						recoBacktrackCCAr->Fill(sr->mc.nu[biggestMatchIndex].iscc);
						recoBacktrackPDGAr->Fill(sr->mc.nu[biggestMatchIndex].pdg);
						
						if (sr->mc.nu[biggestMatchIndex].iscc && abs(sr->mc.nu[biggestMatchIndex].pdg) == 14) {
							
							for (int primaries = 0; primaries < sr->mc.nu[biggestMatchIndex].prim.size(); primaries++) {
								
								double Elep = sr->mc.nu[biggestMatchIndex].prim[primaries].p.E;

								auto start_pos = sr->mc.nu[biggestMatchIndex].prim[primaries].start_pos;
								auto end_pos   = sr->mc.nu[biggestMatchIndex].prim[primaries].end_pos;
								auto p         = sr->mc.nu[biggestMatchIndex].prim[primaries].p;

								// if (std::isnan(start_pos.z)) continue;

								double dX = end_pos.x - start_pos.x;
								double dY = end_pos.y - start_pos.y;
								double dZ = end_pos.z - start_pos.z;

								double length     = TMath::Sqrt(p.px * p.px + p.py * p.py + p.pz * p.pz);
								double lengthTrk  = TMath::Sqrt(dX * dX + dY * dY + dZ * dZ);
								double cosL       = p.pz / length;

								int truePDG = sr->mc.nu[biggestMatchIndex].prim[primaries].pdg;
								
								if ((abs(truePDG) == 13 || abs(truePDG) == 211 || abs(truePDG) == 2212) && lengthTrk > minTrkLength)
									trueMatchMult++;

								// if (goodInteraction)
								// std::cout << biggestMatch << "," << sr->common.ixn.dlp[nixn].vtx.z << ","
								          // << sr->mc.nu[biggestMatchIndex].vtx.z << "," << start_pos.z << ","
								          // << end_pos.z << std::endl;

								if (abs(sr->mc.nu[biggestMatchIndex].prim[primaries].pdg) != 13)
									continue;

								recoBacktrackElAr->Fill(Elep);
								recoBacktrackCoslAr->Fill(cosL);
								responseCosL->Fill(dirZExiting, cosL);

						    }
						   
						    responseMult->Fill(trackMult,trueMatchMult); 
					    }
					}
				}
				
				if (!goodInteraction && sr->mc.nu[biggestMatchIndex].id > 1E9) {
					recoBacktrackCCRock->Fill(sr->mc.nu[biggestMatchIndex].iscc);
					recoBacktrackPDGRock->Fill(sr->mc.nu[biggestMatchIndex].pdg);
				}

				if (!goodInteraction && sr->mc.nu[biggestMatchIndex].id < 1E9) {
					recoBacktrackCCSec->Fill(sr->mc.nu[biggestMatchIndex].iscc);
					recoBacktrackPDGSec->Fill(sr->mc.nu[biggestMatchIndex].pdg);
				}


				// This block looks for FV interactions in which muon is matched to MINERvA.
				if (passesMinerva && goodInteraction) {
					track_multGood->Fill(trackMult);
					hadron_multGood->Fill(hadronMult);
					
					// for nusyst start
					// Save for nusyst reweighting
					out_genieIdx = sr->mc.nu[biggestMatchIndex].genieIdx;
					out_trackMult = trackMult;
					out_hadronMult = hadronMult;
					out_trueEnu = sr->mc.nu[biggestMatchIndex].E;
					syst_tree->Fill();

					// for nusyst end


                    // G4RW reco fill for multiplicity study — re-walk the SAME particles
                    // counted into trackMult above, same selection, NO containment cut,
                    // NO length override (minTrkLength is already 2 from the global setting).
                    for (long unsigned npart2 = 0; npart2 < sr->common.ixn.dlp[nixn].part.dlp.size(); npart2++) {
                    
                        if (!sr->common.ixn.dlp[nixn].part.dlp[npart2].primary) continue;
                    
                        int pdg2 = sr->common.ixn.dlp[nixn].part.dlp[npart2].pdg;
                        if (!(abs(pdg2) == 2212 || abs(pdg2) == 13 || abs(pdg2) == 211 || abs(pdg2) == 321)) continue;
                    
                        auto start_pos2 = sr->common.ixn.dlp[nixn].part.dlp[npart2].start;
                        auto end_pos2   = sr->common.ixn.dlp[nixn].part.dlp[npart2].end;
                    
                        double dVz = abs(start_pos2.z - sr->common.ixn.dlp[nixn].vtx.z);
                        double dVx = abs(start_pos2.x - sr->common.ixn.dlp[nixn].vtx.x);
                        double dVy = abs(start_pos2.y - sr->common.ixn.dlp[nixn].vtx.y);
                        double diffVertex2 = TMath::Sqrt(dVz*dVz + dVy*dVy + dVx*dVx);
                        if (diffVertex2 > 5) continue;
                    
                        double dX2 = end_pos2.x - start_pos2.x;
                        double dY2 = end_pos2.y - start_pos2.y;
                        double dZ2 = end_pos2.z - start_pos2.z;
                        double length2 = TMath::Sqrt(dX2*dX2 + dY2*dY2 + dZ2*dZ2);
                        double cosZ2    = (length2 > 0) ? dZ2 / length2 : -999.0;
                    
                        if (length2 <= minTrkLength) continue;
                    
                        // Truth-match this reco particle -- same pattern as the existing
                        // proton/pion G4RW blocks further down in this file.
                        int    matched_G4ID    = -1;
                        int    matched_truePDG = 0;
                        double maxOverlap2     = 0;
                        auto   truthSize2      = sr->common.ixn.dlp[nixn].part.dlp[npart2].truth.size();
                        for (int bt = 0; bt < (int)truthSize2; bt++) {
                            int parType    = sr->common.ixn.dlp[nixn].part.dlp[npart2].truth[bt].type;
                            int partNumber = sr->common.ixn.dlp[nixn].part.dlp[npart2].truth[bt].part;
                            double overlap = sr->common.ixn.dlp[nixn].part.dlp[npart2].truthOverlap[bt];
                            if (overlap > maxOverlap2 && parType < 3) {
                                maxOverlap2     = overlap;
                                matched_G4ID    = sr->mc.nu[biggestMatchIndex].prim[partNumber].G4ID;
                                matched_truePDG = sr->mc.nu[biggestMatchIndex].prim[partNumber].pdg;
                            }
                        }
                    
                        g4rw_nu_id         = (long long)sr->mc.nu[biggestMatchIndex].id;
                        g4rw_fileIdx       = caf_chain->GetTreeNumber();
                        g4rw_mult_G4ID     = matched_G4ID;
                        g4rw_mult_truePDG  = matched_truePDG;
                        g4rw_mult_cosTheta = cosZ2;
                        g4rw_mult_length   = length2;
                        syst_tree_g4rw_mult_reco->Fill();
                    }
                    // G4RW mult reco end
					
					// track_multGood_muon->Fill(muonMult_good);
					// track_multGood_pion->Fill(pionMult_good);      // Remove: if (pionMult_good > 0)
					// track_multGood_proton->Fill(protonMult_good);  // Remove: if (protonMult_good > 0)
					// track_multGood_kaon->Fill(kaonMult_good);      // Remove: if (kaonMult_good > 0)

if (muonMult_good > 0) track_multGood_muon->Fill(muonMult_good);      
if (pionMult_good > 0) track_multGood_pion->Fill(pionMult_good);      
if (protonMult_good > 0) track_multGood_proton->Fill(protonMult_good);    
if (kaonMult_good > 0) track_multGood_kaon->Fill(kaonMult_good);

					// if (muonMult_good > 0) track_multGood_muon->Fill(track_multGood_muon);						
					// if (pionMult_good > 0) track_multGood_pion->Fill(track_multGood_pion);
					// if (protonMult_good > 0) track_multGood_proton->Fill(hadronMult);
					// if (kaonMult_good > 0) track_multGood_kaon->Fill(hadronMult);
					
					// if (muonMult_good > 0) track_multGood_muon->Fill(muonMult_good);
					// if (pionMult_good > 0) track_multGood_pion->Fill(pionMult_good);
					// if (protonMult_good > 0) track_multGood_proton->Fill(protonMult_good);
					// if (kaonMult_good > 0) track_multGood_kaon->Fill(kaonMult_good);

					int cc   = sr->mc.nu[biggestMatchIndex].iscc;
					int mode = sr->mc.nu[biggestMatchIndex].mode;

					track_multMode->Fill(mode);
					hadron_multMode->Fill(mode);

					if (cc == 1) {
						if (mode == 1 || mode == 1001) {
							track_multQE->Fill(trackMult);
							hadron_multQE->Fill(hadronMult);
						}
						if (mode == 10) {
							track_multMEC->Fill(trackMult);
							hadron_multMEC->Fill(hadronMult);
						}
						if (mode == 3) {
							track_multDIS->Fill(trackMult);
							hadron_multDIS->Fill(hadronMult);
						}
						if (mode == 4) {
							track_multRES->Fill(trackMult);
							hadron_multRES->Fill(hadronMult);
						}
						if (mode == 5) {
							track_multCOH->Fill(trackMult);
							hadron_multCOH->Fill(hadronMult);
						}
					} else {
						track_multNC->Fill(trackMult);
						hadron_multNC->Fill(hadronMult);
					}

					// if (cc == 1) {
					// if (mode == 1 || mode == 1001) 
							// track_multQE->Fill(trackMult);
							// hadron_multQE->Fill(hadronMult);
					// if (mode == 10)
							// track_multMEC->Fill(trackMult);
							// hadron_multMEC->Fill(hadronMult);
						// if (mode == 3)
							// track_multDIS->Fill(trackMult);
							// hadron_multDIS->Fill(hadronMult);
						// if (mode == 4)
							// track_multRES->Fill(trackMult);
							// hadron_multRES->Fill(hadronMult);
						// if (mode == 5)
							// track_multCOH->Fill(trackMult);
							// hadron_multCOH->Fill(hadronMult);
					// } else {
						// track_multNC->Fill(trackMult);
						// hadron_multNC->Fill(hadronMult);
					// }
					
					
					
					if (cc == 1) {
						if (mode == 1 || mode == 1001) {
							++count_match_qe;
							++count_match_total;
						}
						if (mode == 10) {
							++count_match_mec;
							++count_match_total;
						}
						if (mode == 3) {
							++count_match_dis;
							++count_match_total;
						}
						if (mode == 4) {
							++count_match_res;
							++count_match_total;
						}
						if (mode == 5) {
							++count_match_coh;
							++count_match_total;
						}
					} else {
						++count_match_nc;
						++count_match_total;
					}
				}
				
								
				// This block looks for rock and secondary interactions in which muon is matched to MINERvA.				
				if (passesMinerva && !goodInteraction) {
					if (rock == 1) {
						track_multRock->Fill(trackMult);
						hadron_multRock->Fill(hadronMult);
						++count_match_rock;
						++count_match_total;
					} else {
						track_multSec->Fill(trackMult);
						hadron_multSec->Fill(hadronMult);
						++count_match_sec;
						++count_match_total;
					}


					bad_origin->Fill(rock);
					track_multBad->Fill(trackMult);
					hadron_multBad->Fill(hadronMult);


					recoVertex2DBad->Fill(sr->common.ixn.dlp[nixn].vtx.x, sr->common.ixn.dlp[nixn].vtx.z);
					trueVertex2DBad->Fill(trueVtxX, trueVtxZ);
					recoVertex2DBadYZ->Fill(sr->common.ixn.dlp[nixn].vtx.y, sr->common.ixn.dlp[nixn].vtx.z);
				}
				
				
			
				// if (passesMinerva && goodInteraction && cc == 1 && (mode == 1 || mode == 1001)) {
				if (passesMinerva && goodInteraction) {

					int cc   = sr->mc.nu[biggestMatchIndex].iscc;
					int mode = sr->mc.nu[biggestMatchIndex].mode;	

					if (cc == 1 && (mode == 1 || mode == 1001)) {

						int nMu = 0, nProton = 0;

						// Muon kinematics
						double mu_len = -999, mu_cos = -999;
						double mu_startX = -999, mu_startY = -999, mu_startZ = -999;
						double mu_endX = -999, mu_endY = -999, mu_endZ = -999;

						// Proton kinematics
						double p_len = -999, p_cos = -999;
						double p_startX = -999, p_startY = -999, p_startZ = -999;
						double p_endX = -999, p_endY = -999, p_endZ = -999;

						for (size_t npart = 0; npart < sr->common.ixn.dlp[nixn].part.dlp.size(); ++npart) {
							const auto& part = sr->common.ixn.dlp[nixn].part.dlp[npart];
							if (!part.primary || part.truth.empty()) continue;

							if (part.truth[0].type > 2 || part.truthOverlap[0] < 0.5) continue;
							int pdg_backtracked = sr->mc.nu[biggestMatchIndex].prim[part.truth[0].part].pdg;

							auto start = part.start;
							auto end   = part.end;

							double dx = end.x - start.x, dy = end.y - start.y, dz = end.z - start.z;
							double len = std::sqrt(dx*dx + dy*dy + dz*dz);
							double cosz = len > 0 ? dz / len : -2;

							if (abs(pdg_backtracked) == 13) {
								++nMu;
								mu_len = len; 
								mu_cos = cosz; 
								mu_startX = start.x; 
								mu_startY = start.y; 
								mu_startZ = start.z;
								mu_endX = end.x; 
								mu_endY = end.y;
								mu_endZ = end.z;
							}
							
							if (abs(pdg_backtracked) == 2212) {
								++nProton;
								p_len = len; 
								p_cos = cosz;
								p_startX = start.x; 
								p_startY = start.y; 
								p_startZ = start.z;
								p_endX = end.x; 
								p_endY = end.y; 
								p_endZ = end.z;
							}
						}

						if (nMu == 1 && nProton == 1) {
							reco_nuE_muP->Fill(sr->mc.nu[biggestMatchIndex].E);
							reco_mult_muP->Fill(nMu + nProton);
							h2_reco_nuE_vs_muEndZ->Fill(sr->mc.nu[biggestMatchIndex].E, mu_endZ);

							reco_mu_len_muP->Fill(mu_len);
							reco_mu_cos_muP->Fill(mu_cos);
							reco_mu_startX_muP->Fill(mu_startX);
							reco_mu_startY_muP->Fill(mu_startY);
							reco_mu_startZ_muP->Fill(mu_startZ);
							reco_mu_endX_muP->Fill(mu_endX);
							reco_mu_endY_muP->Fill(mu_endY);
							reco_mu_endZ_muP->Fill(mu_endZ);

							reco_p_len_muP->Fill(p_len);
							reco_p_cos_muP->Fill(p_cos);
							reco_p_startX_muP->Fill(p_startX);
							reco_p_startY_muP->Fill(p_startY);
							reco_p_startZ_muP->Fill(p_startZ);
							reco_p_endX_muP->Fill(p_endX);
							reco_p_endY_muP->Fill(p_endY);
							reco_p_endZ_muP->Fill(p_endZ);
							
							
							// Write to CSV
							static std::ofstream out_csv("reco_1mu1p_data.csv", std::ios::app);
							if (out_csv.tellp() == 0) {
								out_csv << "nuE,mu_len,mu_cos,mu_startX,mu_startY,mu_startZ,mu_endX,mu_endY,mu_endZ,"
									<< "p_len,p_cos,p_startX,p_startY,p_startZ,p_endX,p_endY,p_endZ\n";
							}

							out_csv << sr->mc.nu[biggestMatchIndex].E << ","
									<< mu_len << "," << mu_cos << "," << mu_startX << "," << mu_startY << "," << mu_startZ << ","
									<< mu_endX << "," << mu_endY << "," << mu_endZ << ","
									<< p_len << "," << p_cos << "," << p_startX << "," << p_startY << "," << p_startZ << ","
									<< p_endX << "," << p_endY << "," << p_endZ << "\n";
						}
					}
				}


				
				// This block looks for FV interactions in which muon is not matched to MINERvA.
				if (!passesMinerva && goodInteraction) {
					int nMuons = 0;
					double muon_len = -1;
					double muon_cos = -2;  // invalid default
					double start_x = -999, start_y = -999, start_z = -999;

					for (unsigned npart = 0; npart < sr->common.ixn.dlp[nixn].part.dlp.size(); ++npart) {
						auto const& part = sr->common.ixn.dlp[nixn].part.dlp[npart];
						if (!part.primary) continue;

						// ===== Fill confusion matrix for muon/proton/pion primaries =====
						int pdgNumber = 3;
						int backtrackedPDG = 3;
						int backtracked = -9999;
						double maxOverlap = 0;

						for (unsigned bt = 0; bt < part.truth.size(); ++bt) {
							int parType = part.truth[bt].type;
							int partNumber = part.truth[bt].part;
							int ixn = part.truth[bt].ixn;
							double overlap = part.truthOverlap[bt];

							if (overlap > maxOverlap) {
								maxOverlap = overlap;

								if (parType < 3)
									backtracked = sr->mc.nu[ixn].prim[partNumber].pdg;
								else
									backtracked = sr->mc.nu[ixn].sec[partNumber].pdg;
							}
						}

						if (abs(part.pdg) == 13) pdgNumber = 0;
						else if (abs(part.pdg) == 2212) pdgNumber = 1;
						else if (abs(part.pdg) == 211)  pdgNumber = 2;

						if (abs(backtracked) == 13)      backtrackedPDG = 0;
						else if (abs(backtracked) == 2212) backtrackedPDG = 1;
						else if (abs(backtracked) == 211)  backtrackedPDG = 2;

						confusionMatrixAlt->Fill(pdgNumber, backtrackedPDG); // includes rock + non-rock interactions
						
						// If only want to consider non-rock interactions
						// if (rock != 1) {  // Only non-rock interactions
							// confusionMatrixAlt->Fill(pdgNumber, backtrackedPDG);
						// }

						// ===== Proceed to long muon logic only for muons =====
						if (abs(part.pdg) != 13) continue;

						auto start = part.start;
						auto end   = part.end;

						double dx = end.x - start.x;
						double dy = end.y - start.y;
						double dz = end.z - start.z;

						double len = sqrt(dx*dx + dy*dy + dz*dz);
						double cosz = len > 0 ? dz / len : -2;

						if (len > 10.0) {
							++nMuons;
							muon_len = len;
							muon_cos = cosz;
							start_x = start.x;
							start_y = start.y;
							start_z = start.z;
						}
					}

					// if (nMuons > 0) {
					if (nMuons == 1) {
						track_multAltCut->Fill(trackMult);

						int cc   = sr->mc.nu[biggestMatchIndex].iscc;
						int mode = sr->mc.nu[biggestMatchIndex].mode;

						if (cc == 1) {
							if (mode == 1 || mode == 1001) {
								muon_lenAltCutQE->Fill(muon_len);
								muon_cosAltCutQE->Fill(muon_cos);
								muon_startXAltCutQE->Fill(start_x);
								muon_startYAltCutQE->Fill(start_y);
								muon_startZAltCutQE->Fill(start_z);
							}
							else if (mode == 10) {
								muon_lenAltCutMEC->Fill(muon_len);
								muon_cosAltCutMEC->Fill(muon_cos);
								muon_startXAltCutMEC->Fill(start_x);
								muon_startYAltCutMEC->Fill(start_y);
								muon_startZAltCutMEC->Fill(start_z);
							}
							else if (mode == 3) {
								muon_lenAltCutDIS->Fill(muon_len);
								muon_cosAltCutDIS->Fill(muon_cos);
								muon_startXAltCutDIS->Fill(start_x);
								muon_startYAltCutDIS->Fill(start_y);
								muon_startZAltCutDIS->Fill(start_z);
							}
							else if (mode == 4) {
								muon_lenAltCutRES->Fill(muon_len);
								muon_cosAltCutRES->Fill(muon_cos);
								muon_startXAltCutRES->Fill(start_x);
								muon_startYAltCutRES->Fill(start_y);
								muon_startZAltCutRES->Fill(start_z);
							}
							else if (mode == 5) {
								muon_lenAltCutCOH->Fill(muon_len);
								muon_cosAltCutCOH->Fill(muon_cos);
								muon_startXAltCutCOH->Fill(start_x);
								muon_startYAltCutCOH->Fill(start_y);
								muon_startZAltCutCOH->Fill(start_z);
							}
						} else {
							muon_lenAltCutNC->Fill(muon_len);
							muon_cosAltCutNC->Fill(muon_cos);
							muon_startXAltCutNC->Fill(start_x);
							muon_startYAltCutNC->Fill(start_y);
							muon_startZAltCutNC->Fill(start_z);
						}

						if (cc == 1) {
							if (mode == 1 || mode == 1001) track_multAltCutQE->Fill(trackMult);
							if (mode == 10)               track_multAltCutMEC->Fill(trackMult);
							if (mode == 3)                track_multAltCutDIS->Fill(trackMult);
							if (mode == 4)                track_multAltCutRES->Fill(trackMult);
							if (mode == 5)                track_multAltCutCOH->Fill(trackMult);
						} else {
							track_multAltCutNC->Fill(trackMult);
						}

						responseMultAlt->Fill(trackMult, trueMatchMult); // Fill GENIE in response matrix

						if (cc == 1) {
							if (mode == 1 || mode == 1001) { ++count_alt_qe; ++count_alt_total; }
							if (mode == 10)               { ++count_alt_mec; ++count_alt_total; }
							if (mode == 3)                { ++count_alt_dis; ++count_alt_total; }
							if (mode == 4)                { ++count_alt_res; ++count_alt_total; }
							if (mode == 5)                { ++count_alt_coh; ++count_alt_total; }
						} else {
							++count_alt_nc;
							++count_alt_total;
						}
					}
				}

				// This block looks for rock and secondary interactions in which muon is not matched to MINERvA.
				if (!passesMinerva && goodInteraction) {
					if (rock == 1) {
						bool long_muon_found = false;
						double muon_len = -1;
						double muon_cos = -2;
						double start_x = -999, start_y = -999, start_z = -999;

						for (unsigned npart = 0; npart < sr->common.ixn.dlp[nixn].part.dlp.size(); ++npart) {
							auto const& part = sr->common.ixn.dlp[nixn].part.dlp[npart];
							if (!part.primary) continue;
							if (abs(part.pdg) != 13) continue;

							auto start = part.start;
							auto end   = part.end;

							double dx = end.x - start.x;
							double dy = end.y - start.y;
							double dz = end.z - start.z;
							double len = sqrt(dx*dx + dy*dy + dz*dz);
							double cosz = len > 0 ? dz / len : -2;

							if (len > 10.0) {
								long_muon_found = true;
								muon_len = len;
								muon_cos = cosz;
								start_x = start.x;
								start_y = start.y;
								start_z = start.z;
								break;
							}
						}
					
						responseMultAlt->Fill(trackMult, trueMatchMult);  // Fill ROCK in response matrix

						if (long_muon_found) {
							track_multAltCutRock->Fill(trackMult);
							++count_alt_rock;
							++count_alt_total;
						} else {
							muon_lenAltCutRock->Fill(muon_len);
							muon_cosAltCutRock->Fill(muon_cos);
							muon_startXAltCutRock->Fill(start_x);
							muon_startYAltCutRock->Fill(start_y);
							muon_startZAltCutRock->Fill(start_z);
						}

					} else {
						bool long_muon_found = false;
						double muon_len = -1;
						double muon_cos = -2;
						double start_x = -999, start_y = -999, start_z = -999;

						for (unsigned npart = 0; npart < sr->common.ixn.dlp[nixn].part.dlp.size(); ++npart) {
							auto const& part = sr->common.ixn.dlp[nixn].part.dlp[npart];
							if (!part.primary) continue;
							if (abs(part.pdg) != 13) continue;

							auto start = part.start;
							auto end   = part.end;

							double dx = end.x - start.x;
							double dy = end.y - start.y;
							double dz = end.z - start.z;
							double len = sqrt(dx*dx + dy*dy + dz*dz);
							double cosz = len > 0 ? dz / len : -2;

							if (len > 10.0) {
								long_muon_found = true;
								muon_len = len;
								muon_cos = cosz;
								start_x = start.x;
								start_y = start.y;
								start_z = start.z;
								break;
							}
						}

					responseMultAlt->Fill(trackMult, trueMatchMult);  // Fill Secondary in response matrix 
					
					if (long_muon_found) {
							track_multAltCutSec->Fill(trackMult);
							++count_alt_sec;
							++count_alt_total;
					} else {
							muon_lenAltCutSec->Fill(muon_len);
							muon_cosAltCutSec->Fill(muon_cos);
							muon_startXAltCutSec->Fill(start_x);
							muon_startYAltCutSec->Fill(start_y);
							muon_startZAltCutSec->Fill(start_z);
						}
					}
				}

 
			}
			
			if (trackMult > 0 && trackMultExit > 0) {
				track_mult->Fill(trackMult);
				part_mult->Fill(partMult);
			}
 
		} // end for interaction	
	} // end for spills


std::cout << "\n========================================\n";
std::cout << "Reco–truth matched reco proton summary\n";
std::cout << "----------------------------------------\n";
std::cout << "Total matched reco protons        : " << nRecoMatchedProtonAll << "\n";
std::cout << "Contained matched reco protons    : " << nRecoProtonMatchedContained << "\n";
std::cout << "Non-contained matched reco protons: " << nRecoProtonMatchedNonContained << "\n";

if (nRecoMatchedProtonAll > 0) {
	std::cout << "Contained fraction                : "
	          << double(nRecoProtonMatchedContained) / nRecoMatchedProtonAll << "\n";
	std::cout << "Non-contained fraction            : "
	          << double(nRecoProtonMatchedNonContained) / nRecoMatchedProtonAll << "\n";
}

if (nRecoMatchedProtonAll !=
	nRecoProtonMatchedContained + nRecoProtonMatchedNonContained) {
	std::cerr << "[WARNING] Reco matched proton containment counting mismatch!"
	          << std::endl;
}

std::cout << "========================================\n\n";


std::cout << "\n========================================\n";
std::cout << "Reco–truth matched true proton summary\n";
std::cout << "----------------------------------------\n";
std::cout << "Total matched true protons        : " << nTrueMatchedProtonAll << "\n";
std::cout << "Contained matched true protons    : " << nTrueProtonMatchedContained << "\n";
std::cout << "Non-contained matched true protons: " << nTrueProtonMatchedNonContained << "\n";

if (nTrueMatchedProtonAll > 0) {
	std::cout << "Contained fraction                : "
	          << double(nTrueProtonMatchedContained) / nTrueMatchedProtonAll << "\n";
	std::cout << "Non-contained fraction            : "
	          << double(nTrueProtonMatchedNonContained) / nTrueMatchedProtonAll << "\n";
}

if (nTrueMatchedProtonAll !=
	nTrueProtonMatchedContained + nTrueProtonMatchedNonContained) {
	std::cerr << "[WARNING] True matched proton containment counting mismatch!"
	          << std::endl;
}

std::cout << "========================================\n\n";


double nMuon   = recoProtonTrueID->GetBinContent(1);
double nProton = recoProtonTrueID->GetBinContent(2);
double nPion   = recoProtonTrueID->GetBinContent(3);
double nOther  = recoProtonTrueID->GetBinContent(4);
double nTotal  = nMuon + nProton + nPion + nOther;

std::cout << "\n========================================================" << std::endl;
std::cout << "=== Reco Proton Candidate Truth Composition ===" << std::endl;
std::cout << "Total reco proton candidates (after cuts): " << nTotal << std::endl;
std::cout << "  True Muons  (bin 0): " << nMuon   << " (" << 100.0*nMuon/nTotal   << "%)" << std::endl;
std::cout << "  True Protons(bin 1): " << nProton << " (" << 100.0*nProton/nTotal << "%)" << std::endl;
std::cout << "  True Pions  (bin 2): " << nPion   << " (" << 100.0*nPion/nTotal   << "%)" << std::endl;
std::cout << "  True Other  (bin 3): " << nOther  << " (" << 100.0*nOther/nTotal  << "%)" << std::endl;
std::cout << "========================================================\n" << std::endl;


std::cout << "\n========================================\n";
std::cout << "Reco–truth matched reco pion summary\n";
std::cout << "----------------------------------------\n";
std::cout << "Total matched reco pions        : " << nRecoMatchedPionAll << "\n";
std::cout << "Contained matched reco pions    : " << nRecoPionMatchedContained << "\n";
std::cout << "Non-contained matched reco pions: " << nRecoPionMatchedNonContained << "\n";

if (nRecoMatchedPionAll > 0) {
	std::cout << "Contained fraction                : "
	          << double(nRecoPionMatchedContained) / nRecoMatchedPionAll << "\n";
	std::cout << "Non-contained fraction            : "
	          << double(nRecoPionMatchedNonContained) / nRecoMatchedPionAll << "\n";
}

if (nRecoMatchedPionAll !=
	nRecoPionMatchedContained + nRecoPionMatchedNonContained) {
	std::cerr << "[WARNING] Reco matched pion containment counting mismatch!"
	          << std::endl;
}

std::cout << "========================================\n\n";


std::cout << "\n========================================\n";
std::cout << "Reco–truth matched true pion summary\n";
std::cout << "----------------------------------------\n";
std::cout << "Total matched true pions        : " << nTrueMatchedPionAll << "\n";
std::cout << "Contained matched true pions    : " << nTruePionMatchedContained << "\n";
std::cout << "Non-contained matched true pions: " << nTruePionMatchedNonContained << "\n";

if (nTrueMatchedPionAll > 0) {
	std::cout << "Contained fraction                : "
	          << double(nTruePionMatchedContained) / nTrueMatchedPionAll << "\n";
	std::cout << "Non-contained fraction            : "
	          << double(nTruePionMatchedNonContained) / nTrueMatchedPionAll << "\n";
}

if (nTrueMatchedPionAll !=
	nTruePionMatchedContained + nTruePionMatchedNonContained) {
	std::cerr << "[WARNING] True matched pion containment counting mismatch!"
	          << std::endl;
}

std::cout << "========================================\n\n";





	std::cout << "\n========================================================" << std::endl;


	std::cout << "Total number of truth-level neutrino interactions in the CAF file: "
			  << true_mult->GetEntries() << std::endl;

	std::cout << "Number of GENIE signal interactions (excluding NC and unknown modes) over total truth interactions: "
			  << goodIntNum << " / " << true_mult->GetEntries()
			  << " = " << double(goodIntNum) / true_mult->GetEntries() << std::endl;

	std::cout << "Signal purity among GENIE-classified interactions (i.e., signal / (signal + background)): "
			  << float(goodIntNum) / float(goodIntNum + badIntNum) << std::endl;

	std::cout << "Number of GENIE signal interactions within fiducial volume over total truth interactions: "
			  << goodIntInFidVol << " / " << true_mult->GetEntries()
			  << " = " << double(goodIntInFidVol) / true_mult->GetEntries() << std::endl;

	std::cout << "Signal purity in fiducial volume (i.e., fiducial signal / (fiducial signal + fiducial background)): "
			  << float(goodIntInFidVol) / float(goodIntInFidVol + badIntInFidVol) << std::endl;

	std::cout << "Number of reconstructed events with a MINERvA-matched muon passing nominal selection cuts: "
			  << track_multGood->GetEntries() << " / " << track_mult->GetEntries()
			  << " = " << double(track_multGood->GetEntries()) / true_mult->GetEntries() << std::endl;

	std::cout << "Fraction of MINERvA-matched events among all reconstructed events: "
			  << double(track_multGood->GetEntries()) / track_mult->GetEntries() << std::endl;

	std::cout << "MINERvA-matched events with two or more reconstructed tracks: "
			  << track_multGood->Integral(3, 50) << " / " << track_mult->Integral(3, 50)
			  << " = " << double(track_multGood->Integral(3, 50)) / true_mult->Integral(3, 50) << std::endl;

	std::cout << "Fraction of MINERvA-matched events among reconstructed events with ≥ 2 tracks: "
			  << double(track_multGood->Integral(3, 50)) / track_mult->Integral(3, 50) << std::endl;

	std::cout << "Fraction of total CAF entries that result in at least one reconstructed interaction (regardless of MINERvA match): "
          << track_mult->GetEntries() << " / " << Nentries << " = "
          << double(track_mult->GetEntries()) / Nentries << std::endl;


	std::cout << "\n===== Muon MATCH to MINERvA =====" << std::endl;
	std::cout << "QE  / Total = " << count_match_qe  << " / " << count_match_total << " = " << (count_match_total > 0 ? (float)count_match_qe  / count_match_total : 0) << std::endl;
	std::cout << "MEC / Total = " << count_match_mec << " / " << count_match_total << " = " << (count_match_total > 0 ? (float)count_match_mec / count_match_total : 0) << std::endl;
	std::cout << "DIS / Total = " << count_match_dis << " / " << count_match_total << " = " << (count_match_total > 0 ? (float)count_match_dis / count_match_total : 0) << std::endl;
	std::cout << "RES / Total = " << count_match_res << " / " << count_match_total << " = " << (count_match_total > 0 ? (float)count_match_res / count_match_total : 0) << std::endl;
	std::cout << "COH / Total = " << count_match_coh << " / " << count_match_total << " = " << (count_match_total > 0 ? (float)count_match_coh / count_match_total : 0) << std::endl;
	std::cout << "NC  / Total = " << count_match_nc  << " / " << count_match_total << " = " << (count_match_total > 0 ? (float)count_match_nc  / count_match_total : 0) << std::endl;
	std::cout << "Rock  / Total = " << count_match_rock  << " / " << count_match_total << " = " << (count_match_total > 0 ? (float)count_match_rock  / count_match_total : 0) << std::endl;
	std::cout << "Sec  / Total = " << count_match_sec  << " / " << count_match_total << " = " << (count_match_total > 0 ? (float)count_match_sec  / count_match_total : 0) << std::endl;


	std::cout << "\n===== Muon NOT MATCH to MINERvA =====" << std::endl;
	std::cout << "QE  / Total = " << count_alt_qe  << " / " << count_alt_total << " = " << (count_alt_total > 0 ? (float)count_alt_qe  / count_alt_total : 0) << std::endl;
	std::cout << "MEC / Total = " << count_alt_mec << " / " << count_alt_total << " = " << (count_alt_total > 0 ? (float)count_alt_mec / count_alt_total : 0) << std::endl;
	std::cout << "DIS / Total = " << count_alt_dis << " / " << count_alt_total << " = " << (count_alt_total > 0 ? (float)count_alt_dis / count_alt_total : 0) << std::endl;
	std::cout << "RES / Total = " << count_alt_res << " / " << count_alt_total << " = " << (count_alt_total > 0 ? (float)count_alt_res / count_alt_total : 0) << std::endl;
	std::cout << "COH / Total = " << count_alt_coh << " / " << count_alt_total << " = " << (count_alt_total > 0 ? (float)count_alt_coh / count_alt_total : 0) << std::endl;
	std::cout << "NC  / Total = " << count_alt_nc  << " / " << count_alt_total << " = " << (count_alt_total > 0 ? (float)count_alt_nc  / count_alt_total : 0) << std::endl;
	std::cout << "Rock  / Total = " << count_alt_rock  << " / " << count_alt_total << " = " << (count_alt_total > 0 ? (float)count_alt_rock  / count_alt_total : 0) << std::endl;
	std::cout << "Sec  / Total = " << count_alt_sec  << " / " << count_alt_total << " = " << (count_alt_total > 0 ? (float)count_alt_sec  / count_alt_total : 0) << std::endl;


	// std::cout << "Minerva Match Purity: " << goodMINERvAMatch << "/" << totalMINERvAMatch << std::endl;

	// Create output file and write your histograms
	TFile* caf_out_file = new TFile(output_rootfile.c_str(), "recreate");

	true_nuE_all->Write();
	true_nuE_QE->Write();
	true_nuE_RES->Write();
	true_nuE_DIS->Write();
	true_nuE_MEC->Write();
	true_nuE_COH->Write();
	
	true_nuE_all_nosel->Write();
	true_nuE_QE_nosel->Write();
	true_nuE_RES_nosel->Write();
	true_nuE_DIS_nosel->Write();
	true_nuE_MEC_nosel->Write();
	true_nuE_COH_nosel->Write();
	
	true_multTrkOnly_nosel->Write();
	true_multTrkOnly_QE_nosel->Write();
	true_multTrkOnly_MEC_nosel->Write();
	true_multTrkOnly_RES_nosel->Write();
	true_multTrkOnly_DIS_nosel->Write();
	true_multTrkOnly_COH_nosel->Write();
	
	true_vtxX->Write();
	true_vtxY->Write();
	true_vtxZ->Write();
	
	recoVtxX->Write();
	recoVtxY->Write();
	recoVtxZ->Write();
	
	truthMatchedRecoVtxX->Write();
	truthMatchedRecoVtxY->Write();
	truthMatchedRecoVtxZ->Write();
	
	recoVtxX_matched->Write();
	recoVtxY_matched->Write();
	recoVtxZ_matched->Write();
	
	truthMatchedReco_nuE->Write();
	truthMatchedReco_nuE_QE->Write();
	truthMatchedReco_nuE_RES->Write();
	truthMatchedReco_nuE_DIS->Write();
	truthMatchedReco_nuE_MEC->Write();
	truthMatchedReco_nuE_COH->Write();
	
	true_muon_energy->Write();
	true_muon_length->Write();
	true_muon_cosTheta->Write();
	true_muon_startX->Write();
	true_muon_startY->Write();
	true_muon_startZ->Write();
	true_muon_length_vs_energy->Write();
	
	// true_pion_energy->Write();
	// true_pion_length->Write();
	// true_pion_cosTheta->Write();
	// true_pion_startX->Write();
	// true_pion_startY->Write();
	// true_pion_startZ->Write();
	// true_pion_length_vs_energy->Write();

	true_pion_energy->Write();
	true_pion_length->Write();
	true_pion_cosTheta->Write();
	true_pion_startX->Write();
	true_pion_startY->Write();
	true_pion_startZ->Write();
	true_pion_endX->Write();
	true_pion_endY->Write();
	true_pion_endZ->Write();
	true_pion_momentumX->Write();
	true_pion_momentumY->Write();
	true_pion_momentumZ->Write();
	true_pion_length_vs_energy->Write();	
	true_pion_startX_vs_startY->Write();
	true_pion_startZ_vs_startX->Write();
	true_pion_startZ_vs_startY->Write();
	true_pion_startX_vs_endX->Write();
	true_pion_startY_vs_endY->Write();
	true_pion_startZ_vs_endZ->Write();	
	true_pion_endX_vs_endY->Write();
	true_pion_endZ_vs_endX->Write();
	true_pion_endZ_vs_endY->Write();
	true_pion_length_vs_cosTheta->Write();
	true_pion_energy_vs_cosTheta->Write();	
	true_pion_dX_vs_length->Write();
	true_pion_dY_vs_length->Write();
	true_pion_dZ_vs_length->Write();

	true_proton_energy->Write();
	true_proton_length->Write();
	true_proton_cosTheta->Write();
	true_proton_startX->Write();
	true_proton_startY->Write();
	true_proton_startZ->Write();
	true_proton_endX->Write();
	true_proton_endY->Write();
	true_proton_endZ->Write();
	true_proton_momentumX->Write();
	true_proton_momentumY->Write();
	true_proton_momentumZ->Write();
	true_proton_length_vs_energy->Write();	
	true_proton_startX_vs_startY->Write();
	true_proton_startZ_vs_startX->Write();
	true_proton_startZ_vs_startY->Write();
	true_proton_startX_vs_endX->Write();
	true_proton_startY_vs_endY->Write();
	true_proton_startZ_vs_endZ->Write();	
	true_proton_endX_vs_endY->Write();
	true_proton_endZ_vs_endX->Write();
	true_proton_endZ_vs_endY->Write();
	true_proton_length_vs_cosTheta->Write();
	true_proton_energy_vs_cosTheta->Write();	
	true_proton_dX_vs_length->Write();
	true_proton_dY_vs_length->Write();
	true_proton_dZ_vs_length->Write();
	
	true_kaon_energy->Write();
	true_kaon_length->Write();
	true_kaon_cosTheta->Write();
	true_kaon_startX->Write();
	true_kaon_startY->Write();
	true_kaon_startZ->Write();

	matchHistCosl->Write();
	matchHistEl->Write();
	histCosl->Write();
	histEl->Write();
	part_mult->Write();
	part_multTrkOnly->Write();
	part_energy_hist->Write();
	track_mult->Write();
	track_multGood->Write();
	track_multBad->Write();
	track_length->Write();
	track_multNC->Write();
	track_multMode->Write();
	track_multQE->Write();
	track_multDIS->Write();
	track_multMEC->Write();
	track_multRES->Write();
	track_multCOH->Write();
	track_multRock->Write();
	track_multSec->Write();
	minervaMatchPDG->Write();
	minervaMatchCos->Write();
	minervaMatchE->Write();
	responseCosL->Write();

	track_multGood_muon->Write();
	track_multGood_pion->Write();
	track_multGood_proton->Write();
	track_multGood_kaon->Write();

	hadron_mult->Write();
	hadron_multMode->Write();
	hadron_multQE->Write();
	hadron_multMEC->Write();
	hadron_multRES->Write();
	hadron_multCOH->Write();
	hadron_multDIS->Write();
	hadron_multNC->Write();
	hadron_multRock->Write();
	hadron_multSec->Write();
	hadron_multBad->Write();
	hadron_multGood->Write();


	recoBacktrackCCAr->Write();
	recoBacktrackCCRock->Write();
	recoBacktrackCCSec->Write();

	recoBacktrackPDGSec->Write();
	recoBacktrackPDGAr->Write();
	recoBacktrackPDGRock->Write();

	recoBacktrackElAr->Write();
	recoBacktrackCoslAr->Write();

	true_mult->Write();
	true_multTrkOnly->Write();
	true_multGENIE->Write();
	
	true_multTrkOnly_QE->Write();
	true_multTrkOnly_MEC->Write();
	true_multTrkOnly_RES->Write();
	true_multTrkOnly_DIS->Write();
	true_multTrkOnly_COH->Write();
	
	true_multTrkOnly_muon->Write();
	true_multTrkOnly_pion->Write();
	true_multTrkOnly_proton->Write();
	true_multTrkOnly_kaon->Write();

	recoHistVertexX->Write();
	recoHistVertexZ->Write();
	recoHistVertexY->Write();

	matchTrue_mult->Write();
	matchTrue_multTrkOnly->Write();
	matchTrue_multGENIE->Write();
	bad_origin->Write();

	recoVertex2D->Write();
	recoVertex2DNoCuts->Write();
	longest_track->Write();
	showerCorrectness->Write();
	trackCorrectness->Write();

	selPionE->Write();
	selKaonE->Write();
	selProtonE->Write();

	truePionEWithRecoInt->Write();
	trueProtonEWithRecoInt->Write();
	trueKaonEWithRecoInt->Write();

	diffPip->Write();
	diffPim->Write();
	diffP->Write();
	diffN->Write();

	recoVertex2DBadYZ->Write();
	recoVertex2DBad->Write();
	trueVertex2DBad->Write();

	auto protonKEEff = new TEfficiency(*selProtonE, *trueProtonEWithRecoInt);
	auto kaonKEEff = new TEfficiency(*selKaonE, *trueKaonEWithRecoInt);
	auto pionKEEff = new TEfficiency(*selPionE, *truePionEWithRecoInt);


	auto protonDirX = new TEfficiency(*selProtonDirX, *trueProtonWithRecoIntDirX);
	auto protonDirY = new TEfficiency(*selProtonDirY, *trueProtonWithRecoIntDirY);
	auto protonDirZ = new TEfficiency(*selProtonDirZ, *trueProtonWithRecoIntDirZ);

	auto pionDirX = new TEfficiency(*selPionDirX, *truePionWithRecoIntDirX);
	auto pionDirY = new TEfficiency(*selPionDirY, *truePionWithRecoIntDirY);
	auto pionDirZ = new TEfficiency(*selPionDirZ, *truePionWithRecoIntDirZ);

	auto kaonDirX = new TEfficiency(*selKaonDirX, *trueKaonWithRecoIntDirX);
	auto kaonDirY = new TEfficiency(*selKaonDirY, *trueKaonWithRecoIntDirY);
	auto kaonDirZ = new TEfficiency(*selKaonDirZ, *trueKaonWithRecoIntDirZ);

	auto effTrkLen = new TEfficiency(*selTrkLen, *trueTrkLen);
	auto effTrkLenQE = new TEfficiency(*selTrkLenQE, *trueTrkLenQE);
	auto effTrkLenNonQE = new TEfficiency(*selTrkLenNonQE, *trueTrkLenNonQE);
	
	auto multEff = new TEfficiency(*matchTrue_multTrkOnly, *true_multTrkOnly);

	effTrkLen->Write();
	effTrkLenQE->Write();
	effTrkLenNonQE->Write();

	protonDirX->Write();
	protonDirY->Write();
	protonDirZ->Write();

	pionDirX->Write();
	pionDirY->Write();
	pionDirZ->Write();

	kaonDirX->Write();
	kaonDirY->Write();
	kaonDirZ->Write();

	protonKEEff->Write();
	kaonKEEff->Write();
	pionKEEff->Write();
	
	multEff->Write();

	nP->Write();
	nPi->Write();
	escapePi->Write();
	containPiLen->Write();
	containPLen->Write();

	responseMult->Write();
	responseMultAlt->Write();
	responseGenieToG4->Write();
	confusionMatrix->Write();
	confusionMatrixAlt->Write();
	
	confusionMatrix_truthContained->Write();
	confusionMatrix_truthNotContained->Write();
	
	recoCosL->Write();
	trueDiffPosvsPDirZ->Write();

	TH1D* track_multShape = (TH1D*)track_mult->Clone("track_multShape");
	track_multShape->Scale(1.f / track_mult->GetEntries());
	track_multShape->Write();


	recoProtonLength->Write();
	recoProtonCosL->Write();
	recoProtonE->Write();
	recoProtonStartX->Write();
	recoProtonStartY->Write();
	recoProtonStartZ->Write();
	recoProtonEndX->Write();
	recoProtonEndY->Write();
	recoProtonEndZ->Write();
	recoProtonMomentumX->Write();
	recoProtonMomentumY->Write();
	recoProtonMomentumZ->Write();
	recoProtonLength_vs_E->Write();
	recoProton_startX_vs_startY->Write();
	recoProton_startZ_vs_startX->Write();
	recoProton_startZ_vs_startY->Write();
	recoProton_startX_vs_endX->Write();
	recoProton_startY_vs_endY->Write();
	recoProton_startZ_vs_endZ->Write();	
	recoProton_endX_vs_endY->Write();
	recoProton_endZ_vs_endX->Write();
	recoProton_endZ_vs_endY->Write();
	recoProton_length_vs_cosL->Write();
	recoProton_energy_vs_cosL->Write();
	recoProton_dX_vs_length->Write();
	recoProton_dY_vs_length->Write();
	recoProton_dZ_vs_length->Write();
	recoProtonCosAlpha->Write();
	recoProtonCosAlpha_vs_E->Write();
	recoProtonCosAlpha_vs_length->Write();
	primProtonMultiplicity->Write();
	recoProton_cosL_vs_startz->Write();
	recoProton_cosL_vs_endz->Write();

	recoPionLength->Write();
	recoPionCosL->Write();
	recoPionE->Write();
	recoPionStartX->Write();
	recoPionStartY->Write();
	recoPionStartZ->Write();
	recoPionEndX->Write();
	recoPionEndY->Write();
	recoPionEndZ->Write();
	recoPionMomentumX->Write();
	recoPionMomentumY->Write();
	recoPionMomentumZ->Write();
	recoPionLength_vs_E->Write();
	recoPion_startX_vs_startY->Write();
	recoPion_startZ_vs_startX->Write();
	recoPion_startZ_vs_startY->Write();
	recoPion_startX_vs_endX->Write();
	recoPion_startY_vs_endY->Write();
	recoPion_startZ_vs_endZ->Write();	
	recoPion_endX_vs_endY->Write();
	recoPion_endZ_vs_endX->Write();
	recoPion_endZ_vs_endY->Write();
	recoPion_length_vs_cosL->Write();
	recoPion_energy_vs_cosL->Write();
	recoPion_dX_vs_length->Write();
	recoPion_dY_vs_length->Write();
	recoPion_dZ_vs_length->Write();
	recoPionCosAlpha->Write();
	recoPionCosAlpha_vs_E->Write();
	recoPionCosAlpha_vs_length->Write();
	primPionMultiplicity->Write();
	recoPion_cosL_vs_startz->Write();
	recoPion_cosL_vs_endz->Write();
	
	// recoPionLength->Write();
	// recoPionCosL->Write();
	// recoPionE->Write();
	// recoPionStartX->Write();
	// recoPionStartY->Write();
	// recoPionStartZ->Write();
	// recoPionLength_vs_E->Write();


	recoProtonLength_QE->Write();
	recoProtonLength_RES->Write();
	recoProtonLength_DIS->Write();
	recoProtonLength_MEC->Write();
	recoProtonLength_COH->Write();
	recoProtonLength_Rock->Write();
	recoProtonLength_Sec->Write();


	recoMuonLength->Write();
	recoMuonCosL->Write();
	recoMuonE->Write();
	recoMuonStartX->Write();
	recoMuonStartY->Write();
	recoMuonStartZ->Write();
	recoMuonLength_vs_E->Write();
	
	recoMuonLength_QE->Write();
	recoMuonLength_RES->Write();
	recoMuonLength_DIS->Write();
	recoMuonLength_MEC->Write();
	recoMuonLength_COH->Write();
	recoMuonLength_Rock->Write();
	recoMuonLength_Sec->Write();

	recoMuonLengthAll->Write();
	recoMuonCosLAll->Write();
	recoMuonEAll->Write();
	recoMuonStartXAll->Write();
	recoMuonStartYAll->Write();
	recoMuonStartZAll->Write();
	recoMuonEndXAll->Write();
	recoMuonEndYAll->Write();
	recoMuonEndZAll->Write();
	recoMuonMomentumXAll->Write();
	recoMuonMomentumYAll->Write();
	recoMuonMomentumZAll->Write();
	recoMuonLengthAll_vs_E->Write();

	recoProtonLengthAll->Write();
	recoProtonCosLAll->Write();
	recoProtonEAll->Write();
	recoProtonStartXAll->Write();
	recoProtonStartYAll->Write();
	recoProtonStartZAll->Write();
	recoProtonEndXAll->Write();
	recoProtonEndYAll->Write();
	recoProtonEndZAll->Write();
	recoProtonMomentumXAll->Write();
	recoProtonMomentumYAll->Write();
	recoProtonMomentumZAll->Write();
	recoProtonLength_vs_EAll->Write();
	recoProton_startX_vs_startYAll->Write();
	recoProton_startZ_vs_startXAll->Write();
	recoProton_startZ_vs_startYAll->Write();
	recoProton_startX_vs_endXAll->Write();
	recoProton_startY_vs_endYAll->Write();
	recoProton_startZ_vs_endZAll->Write();	
	recoProton_endX_vs_endYAll->Write();
	recoProton_endZ_vs_endXAll->Write();
	recoProton_endZ_vs_endYAll->Write();
	recoProton_length_vs_cosLAll->Write();
	recoProton_energy_vs_cosLAll->Write();
	recoProton_dX_vs_lengthAll->Write();
	recoProton_dY_vs_lengthAll->Write();
	recoProton_dZ_vs_lengthAll->Write();

	recoPionLengthAll->Write();
	recoPionCosLAll->Write();
	recoPionEAll->Write();
	recoPionStartXAll->Write();
	recoPionStartYAll->Write();
	recoPionStartZAll->Write();
	recoPionEndXAll->Write();
	recoPionEndYAll->Write();
	recoPionEndZAll->Write();
	recoPionMomentumXAll->Write();
	recoPionMomentumYAll->Write();
	recoPionMomentumZAll->Write();
	recoPionLength_vs_EAll->Write();
	recoPion_startX_vs_startYAll->Write();
	recoPion_startZ_vs_startXAll->Write();
	recoPion_startZ_vs_startYAll->Write();
	recoPion_startX_vs_endXAll->Write();
	recoPion_startY_vs_endYAll->Write();
	recoPion_startZ_vs_endZAll->Write();	
	recoPion_endX_vs_endYAll->Write();
	recoPion_endZ_vs_endXAll->Write();
	recoPion_endZ_vs_endYAll->Write();
	recoPion_length_vs_cosLAll->Write();
	recoPion_energy_vs_cosLAll->Write();
	recoPion_dX_vs_lengthAll->Write();
	recoPion_dY_vs_lengthAll->Write();
	recoPion_dZ_vs_lengthAll->Write();

	// recoPionLengthAll->Write();
	// recoPionCosLAll->Write();
	// recoPionEAll->Write();
	// recoPionStartXAll->Write();
	// recoPionStartYAll->Write();
	// recoPionStartZAll->Write();
	// recoPionEndXAll->Write();
	// recoPionEndYAll->Write();
	// recoPionEndZAll->Write();
	// recoPionMomentumXAll->Write();
	// recoPionMomentumYAll->Write();
	// recoPionMomentumZAll->Write();
	// recoPionLengthAll_vs_E->Write();

	recoKaonLengthAll->Write();
	recoKaonCosLAll->Write();
	recoKaonEAll->Write();
	recoKaonStartXAll->Write();
	recoKaonStartYAll->Write();
	recoKaonStartZAll->Write();
	recoKaonEndXAll->Write();
	recoKaonEndYAll->Write();
	recoKaonEndZAll->Write();
	recoKaonMomentumXAll->Write();
	recoKaonMomentumYAll->Write();
	recoKaonMomentumZAll->Write();


// ===========================================================
// =============== Write Reco (Matched) Proton Histograms =====
// ===========================================================

// 1D
reco_matched_proton_length->Write();
reco_matched_proton_cosL->Write();
reco_matched_proton_energy->Write();

// Interaction-type
reco_matched_proton_length_QE->Write();
reco_matched_proton_length_MEC->Write();
reco_matched_proton_length_RES->Write();
reco_matched_proton_length_DIS->Write();
reco_matched_proton_length_COH->Write();
reco_matched_proton_cosL_QE->Write();
reco_matched_proton_cosL_MEC->Write();
reco_matched_proton_cosL_RES->Write();
reco_matched_proton_cosL_DIS->Write();
reco_matched_proton_cosL_COH->Write();

// Particle-type
reco_matched_proton_length_trueMuon->Write();
reco_matched_proton_length_trueProton->Write();
reco_matched_proton_length_truePion->Write();
reco_matched_proton_length_trueOther->Write();
reco_matched_proton_cosL_trueMuon->Write();
reco_matched_proton_cosL_trueProton->Write();
reco_matched_proton_cosL_truePion->Write();
reco_matched_proton_cosL_trueOther->Write();

reco_matched_proton_startX->Write();
reco_matched_proton_startY->Write();
reco_matched_proton_startZ->Write();

reco_matched_proton_endX->Write();
reco_matched_proton_endY->Write();
reco_matched_proton_endZ->Write();

reco_matched_proton_momentumX->Write();
reco_matched_proton_momentumY->Write();
reco_matched_proton_momentumZ->Write();

// 2D
reco_matched_proton_length_vs_energy->Write();
reco_matched_proton_startX_vs_startY->Write();
reco_matched_proton_startZ_vs_startX->Write();
reco_matched_proton_startZ_vs_startY->Write();

reco_matched_proton_startX_vs_endX->Write();
reco_matched_proton_startY_vs_endY->Write();
reco_matched_proton_startZ_vs_endZ->Write();

reco_matched_proton_length_vs_cosL->Write();
reco_matched_proton_energy_vs_cosL->Write();

reco_matched_proton_dX_vs_length->Write();
reco_matched_proton_dY_vs_length->Write();
reco_matched_proton_dZ_vs_length->Write();

reco_matched_proton_endX_vs_endY->Write();
reco_matched_proton_endZ_vs_endX->Write();
reco_matched_proton_endZ_vs_endY->Write();

reco_matched_proton_cosL_vs_startZ->Write();
reco_matched_proton_cosL_vs_endZ->Write();

recoProtonMultiplicity_matched->Write();

recoProtonTrueID->Write();


// ===========================================================
// =============== Write Reco (Matched) Pion Histograms =====
// ===========================================================

// 1D
reco_matched_pion_length->Write();
reco_matched_pion_cosL->Write();
reco_matched_pion_energy->Write();
reco_matched_pion_energy_fromP->Write();
reco_matched_pion_length_vs_energy_fromP->Write();
reco_matched_pion_energy_diff->Write(); 

reco_matched_pion_startX->Write();
reco_matched_pion_startY->Write();
reco_matched_pion_startZ->Write();

reco_matched_pion_endX->Write();
reco_matched_pion_endY->Write();
reco_matched_pion_endZ->Write();

reco_matched_pion_momentumX->Write();
reco_matched_pion_momentumY->Write();
reco_matched_pion_momentumZ->Write();

// 2D
reco_matched_pion_length_vs_energy->Write();
reco_matched_pion_startX_vs_startY->Write();
reco_matched_pion_startZ_vs_startX->Write();
reco_matched_pion_startZ_vs_startY->Write();

reco_matched_pion_startX_vs_endX->Write();
reco_matched_pion_startY_vs_endY->Write();
reco_matched_pion_startZ_vs_endZ->Write();

reco_matched_pion_length_vs_cosL->Write();
reco_matched_pion_energy_vs_cosL->Write();

reco_matched_pion_dX_vs_length->Write();
reco_matched_pion_dY_vs_length->Write();
reco_matched_pion_dZ_vs_length->Write();

reco_matched_pion_endX_vs_endY->Write();
reco_matched_pion_endZ_vs_endX->Write();
reco_matched_pion_endZ_vs_endY->Write();

recoPionMultiplicity_matched->Write();

reco_matched_pion_cosL_vs_startZ->Write();
reco_matched_pion_cosL_vs_endZ->Write();


	track_multAltCut->Write();
	track_multAltCutQE->Write();
	track_multAltCutMEC->Write();
	track_multAltCutRES->Write();
	track_multAltCutDIS->Write();
	track_multAltCutCOH->Write();
	track_multAltCutNC->Write();
	track_multAltCutRock->Write();
	track_multAltCutSec->Write();

	muon_lenAltCutQE->Write();
	muon_lenAltCutMEC->Write();
	muon_lenAltCutDIS->Write();
	muon_lenAltCutRES->Write();
	muon_lenAltCutCOH->Write();
	muon_lenAltCutNC->Write();
	muon_lenAltCutRock->Write();
	muon_lenAltCutSec->Write();


	muon_cosAltCutQE->Write();
	muon_cosAltCutMEC->Write();
	muon_cosAltCutDIS->Write();
	muon_cosAltCutRES->Write();
	muon_cosAltCutCOH->Write();
	muon_cosAltCutNC->Write();
	muon_cosAltCutRock->Write();
	muon_cosAltCutSec->Write();


	muon_startXAltCutQE->Write();
	muon_startXAltCutMEC->Write();
	muon_startXAltCutDIS->Write();
	muon_startXAltCutRES->Write();
	muon_startXAltCutCOH->Write();
	muon_startXAltCutNC->Write();
	muon_startXAltCutRock->Write();
	muon_startXAltCutSec->Write();

	muon_startYAltCutQE->Write();
	muon_startYAltCutMEC->Write();
	muon_startYAltCutDIS->Write();
	muon_startYAltCutRES->Write();
	muon_startYAltCutCOH->Write();
	muon_startYAltCutNC->Write();
	muon_startYAltCutRock->Write();
	muon_startYAltCutSec->Write();

	muon_startZAltCutQE->Write();
	muon_startZAltCutMEC->Write();
	muon_startZAltCutDIS->Write();
	muon_startZAltCutRES->Write();
	muon_startZAltCutCOH->Write();
	muon_startZAltCutNC->Write();
	muon_startZAltCutRock->Write();
	muon_startZAltCutSec->Write();


	true_mult_muP->Write(); 
	true_nuE_muP->Write();	

	true_mu_len_muP->Write();
	true_mu_cos_muP->Write();
	true_mu_startX_muP->Write();
	true_mu_startY_muP->Write();
	true_mu_startZ_muP->Write();
	true_mu_endX_muP->Write();
	true_mu_endY_muP->Write();
	true_mu_endZ_muP->Write();
	
	true_p_len_muP->Write();
	true_p_cos_muP->Write();
	true_p_startX_muP->Write();
	true_p_startY_muP->Write();
	true_p_startZ_muP->Write();
	true_p_endX_muP->Write();
	true_p_endY_muP->Write();
	true_p_endZ_muP->Write();
	
	h2_true_nuE_vs_muEndZ->Write();
	
	
reco_nuE_muP->Write();
reco_mult_muP->Write();
h2_reco_nuE_vs_muEndZ->Write();

reco_mu_len_muP->Write();
reco_mu_cos_muP->Write();
reco_mu_startX_muP->Write();
reco_mu_startY_muP->Write();
reco_mu_startZ_muP->Write();
reco_mu_endX_muP->Write();
reco_mu_endY_muP->Write();
reco_mu_endZ_muP->Write();

reco_p_len_muP->Write();
reco_p_cos_muP->Write();
reco_p_startX_muP->Write();
reco_p_startY_muP->Write();
reco_p_startZ_muP->Write();
reco_p_endX_muP->Write();
reco_p_endY_muP->Write();
reco_p_endZ_muP->Write();	



// ==========================
// Truth-matched Protons — Contained
// ==========================
trueProtonMatchedLengthContained->Write();
trueProtonMatchedEContained->Write();
trueProtonMatchedLength_vs_EContained->Write();

// ==========================
// Truth-matched Protons — Non-contained
// ==========================
trueProtonMatchedLengthNonContained->Write();
trueProtonMatchedENonContained->Write();
trueProtonMatchedLength_vs_ENonContained->Write();


// ==========================
// Reco-matched Protons — Contained
// ==========================
recoProtonMatchedLengthContained->Write();
recoProtonMatchedEContained->Write();
recoProtonMatchedLength_vs_EContained->Write();

// ==========================
// Reco-matched Protons — Non-contained
// ==========================
recoProtonMatchedLengthNonContained->Write();
recoProtonMatchedENonContained->Write();
recoProtonMatchedLength_vs_ENonContained->Write();


	// for nusyst start
	syst_tree_trueMult->Write();
	syst_tree_recoProtonCosL->Write();
    syst_tree_recoProtonLength->Write();
	syst_tree->Write();
	// for nusyst end

// for g4rw start
syst_tree_g4rw_pion->Write();
syst_tree_g4rw_pion_reco->Write();

syst_tree_g4rw_proton->Write();
syst_tree_g4rw_proton_reco->Write();

syst_tree_g4rw_mult_reco->Write();
	// for g4rw start



	caf_out_file->Close();

	return 1;

}


int main(int argc, char** argv) {

	if (argc != 4) {
		std::cout << "\n USAGE: " << argv[0] << " input_caf_file_list output_root_file mcOnly_flag (0 or 1)\n" << std::endl;
		return 1;
	}

	std::string input_file_list = argv[1];
	std::string output_rootfile = argv[2];
	std::string mcOnlyString = argv[3];

	bool mcOnly = true;
	if (mcOnlyString == "0") mcOnly = false;

	std::cout << mcOnly << "," << argv[3] << std::endl;

	caf_plotter(input_file_list, output_rootfile, mcOnly);

	return 0;
}
