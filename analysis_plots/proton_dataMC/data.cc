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

long nRecoProtonAll          = 0;
long nRecoProtonContained    = 0;
long nRecoProtonNonContained = 0;

long nRecoPionAll          = 0;
long nRecoPionContained    = 0;
long nRecoPionNonContained = 0;

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

    TH1D* track_mult = new TH1D("track_mult", "track_mult", 15, 0, 15);

			
	// Reconstructed distributions
	TH1D* part_energy_hist    = new TH1D("recpart_energy", "Reco particle energy in GeV", 100, 0, 1); // Distribution of reconstructed particle energy (GeV)
	TH1D* track_length        = new TH1D("track_length", "track_length", 100, 0, 100); // Distribution of reconstructed track lengths
	TH1D* part_multTrkOnly    = new TH1D("part_multTrkOnly", "part_multTrkOnly", 20, 0, 20); // Primary particle multiplicity with track length > threshold
	TH1D* part_mult           = new TH1D("part_mult", "part_mult", 20, 0, 20); // Total reconstructed primary particle multiplicity
	

	// Reconstructed neutrino vertices
	TH1D* recoVtxX = new TH1D("recoVtxX", "Reco Neutrino Vertex X;X [cm];Events", 70, -70, 70);
	TH1D* recoVtxY = new TH1D("recoVtxY", "Reco Neutrino Vertex Y;Y [cm];Events", 70, -70, 70);
	TH1D* recoVtxZ = new TH1D("recoVtxZ", "Reco Neutrino Vertex Z;Z [cm];Events", 70, -70, 70);
	
	TH2D* recoVertex2DNoCuts		= new TH2D("recoVertex2DNoCuts", "recoVertex2DNoCuts", 70, -70, 70, 70, -70, 70);	// Reco vertex (X,Z) before cuts
	TH2D* recoVertex2D				= new TH2D("recoVertex2D", "recoVertex2D", 60, -60, 60, 60, -60, 60);				// Reco vertex (X,Z) after basic cuts


	// Other diagnostics
	TH1D* recoCosL					= new TH1D("recoCosL", "recoCosL", 100, -1, 1);  // Reco muon cos(θ)
	TH1D* longest_track				= new TH1D("longest_track", "longest_track", 100, 0, 100); // Length of longest reconstructed track in each interaction


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
	TH1D* recoProtonLengthAll  = new TH1D("recoProtonLengthAll", "Reco Proton Track Length (All)", 300, 0, 300);
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
	TH1D* recoPionLengthAll  = new TH1D("recoPionLengthAll", "Reco Pion Track Length (All)", 300, 0, 300);
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


	// Protons
	
	TH1D* recoProtonLengthContained =
	new TH1D("recoProtonLengthContained",
	         "Reco Proton Track Length (Contained);Track Length [cm];Events",
	         100, 0, 100);

	TH1D* recoProtonLengthNonContained =
		new TH1D("recoProtonLengthNonContained",
				 "Reco Proton Track Length (Non-contained);Track Length [cm];Events",
				 100, 0, 100);

	TH1D* recoProtonEContained =
		new TH1D("recoProtonEContained",
				 "Reco Proton Energy (Contained);Energy [GeV];Events",
				 100, 0, 2);

	TH1D* recoProtonENonContained =
		new TH1D("recoProtonENonContained",
				 "Reco Proton Energy (Non-contained);Energy [GeV];Events",
				 100, 0, 2);

	TH2D* recoProtonLength_vs_EContained =
		new TH2D("recoProtonLength_vs_EContained",
				 "Reco Proton (Contained): Track Length vs Energy;"
				 "Track Length [cm];Energy [GeV]",
				 100, 0, 100,
				 100, 0, 2);

	TH2D* recoProtonLength_vs_ENonContained =
		new TH2D("recoProtonLength_vs_ENonContained",
				 "Reco Proton (Non-contained): Track Length vs Energy;"
				 "Track Length [cm];Energy [GeV]",
				 100, 0, 100,
				 100, 0, 2);

	// Pions
	
	TH1D* recoPionLengthContained =
	new TH1D("recoPionLengthContained",
	         "Reco Pion Track Length (Contained);Track Length [cm];Events",
	         100, 0, 100);

	TH1D* recoPionLengthNonContained =
		new TH1D("recoPionLengthNonContained",
				 "Reco Pion Track Length (Non-contained);Track Length [cm];Events",
				 100, 0, 100);

	TH1D* recoPionEContained =
		new TH1D("recoPionEContained",
				 "Reco Pion Energy (Contained);Energy [GeV];Events",
				 100, 0, 2);

	TH1D* recoPionENonContained =
		new TH1D("recoPionENonContained",
				 "Reco Pion Energy (Non-contained);Energy [GeV];Events",
				 100, 0, 2);

	TH2D* recoPionLength_vs_EContained =
		new TH2D("recoPionLength_vs_EContained",
				 "Reco Pion (Contained): Track Length vs Energy;"
				 "Track Length [cm];Energy [GeV]",
				 100, 0, 100,
				 100, 0, 2);

	TH2D* recoPionLength_vs_ENonContained =
		new TH2D("recoPionLength_vs_ENonContained",
				 "Reco Pion (Non-contained): Track Length vs Energy;"
				 "Track Length [cm];Energy [GeV]",
				 100, 0, 100,
				 100, 0, 2);

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



	// Loop over all entries (spills) in the CAF chain
	for(long n = 0; n < Nentries; ++n) { 
		
		// Initialize variables for the true neutrino vertex position
		double trueVtxX = -9999; 
		double trueVtxY = -999; 
		double trueVtxZ =- 9999;
        skipEvent = false;

		if (n % 10000 == 0) std::cout << Form("Processing trigger %ld of %ld", n, Nentries) << std::endl;

		caf_chain->GetEntry(n);	// Get spill from tree

		bool hasANeutrino = false;

		double mnvOffsetX = -10;
		double mnvOffsetY = 5;


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
					// if (sr->common.ixn.dlp[nixn].part.dlp[npart].primary == true)				

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
				} // track PDGs 
			} // particles
			
			// if (/*minervaThrough < 1 ||*/ maxDotProductDS < 0.99) continue;
			bool passesMinerva = (maxDotProductDS > 0.99);
			
			// Require MINERvA-matched muon for reco interaction selection
			// if (!passesMinerva) continue;
            track_mult->Fill(trackMult);


			recoCosL->Fill(dirZExiting);
			recoVertex2D->Fill(sr->common.ixn.dlp[nixn].vtx.x, sr->common.ixn.dlp[nixn].vtx.z);
			
			auto recoVtx = sr->common.ixn.dlp[nixn].vtx;
			recoVtxX->Fill(recoVtx.x);
			recoVtxY->Fill(recoVtx.y);
			recoVtxZ->Fill(recoVtx.z);


				int nPrimProtons_true    = 0;
				int nRecoProtons_matched = 0;				

						
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

                    // if (length <= 2.0) continue;   // require proton track length > 2 cm
                    // if (end_pos.z >= -5.0 && end_pos.z <= 5.0) continue; // reject |z_end| ≤ 5 cm → keep outside 

                    
					if (abs(pdg) == 13) {
						auto start_pos = sr->common.ixn.dlp[nixn].part.dlp[npart].start;
						auto end_pos   = sr->common.ixn.dlp[nixn].part.dlp[npart].end;
						
						auto momentum   = sr->common.ixn.dlp[nixn].part.dlp[npart].p;

						double dX = end_pos.x - start_pos.x;
						double dY = end_pos.y - start_pos.y;
						double dZ = end_pos.z - start_pos.z;

						double length = TMath::Sqrt(dX * dX + dY * dY + dZ * dZ);
						double dirZ = dZ / length;

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
						
						bool truth_contained_geom =
							(start_pos.x > -55 && start_pos.x < 55) &&
							(start_pos.y > -55 && start_pos.y < 55) &&
							(start_pos.z > -55 && start_pos.z < 55) &&
							(end_pos.x   > -55 && end_pos.x   < 55) &&
							(end_pos.y   > -55 && end_pos.y   < 55) &&
							(end_pos.z   > -55 && end_pos.z   < 55);
							
						// bool truth_contained_geom =	sr->common.ixn.dlp[nixn].part.dlp[npart].contained;
	

						if (!truth_contained_geom) continue;	// skip non-contained

						// if (truth_contained_geom) continue;	// skip contained
						
						auto momentum   = sr->common.ixn.dlp[nixn].part.dlp[npart].p;

						double dX = end_pos.x - start_pos.x;
						double dY = end_pos.y - start_pos.y;
						double dZ = end_pos.z - start_pos.z;

						double length = TMath::Sqrt(dX * dX + dY * dY + dZ * dZ);
						double dirZ = dZ / length;

                        // if (length <= 2.0) continue;   // require proton track length > 2 cm
                        if (length <= 3.0) continue;   // require proton track length > 3 cm

                        // if (length > 5.0) continue;   // keep only protons with L < 10 cm
                        // if (length < 2.0 || length > 5.0) continue;   // keep 2 cm ≤ L ≤ 5 cm

                        // if (end_pos.z >= -5.0 && end_pos.z <= 5.0) continue; // reject |z_end| ≤ 5 cm → keep outside 
						if (end_pos.z >= -5.0 && end_pos.z <= 0.0) continue; //  reject only −z side (−5 to 0 cm):

                        // if (dirZ > -0.95) continue;   // keep only backward-going protons: cosθ ∈ [−1, −0.9]
                        if (dirZ <= -0.95) continue;  // rejects cosθ ≤ −0.95, keeps cosθ > −0.95

                        // Count all reco protons
						nRecoProtonAll++;
                        
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
						
						// if (truth_contained_geom) {
						// 	nRecoProtonContained++;
							
						// 	recoProtonLengthContained->Fill(length);
						// 	recoProtonEContained->Fill(energy);
						// 	recoProtonLength_vs_EContained->Fill(length, energy);
						// }

						// if (!truth_contained_geom) {
						// 	nRecoProtonNonContained++;
							
						// 	recoProtonLengthNonContained->Fill(length);
						// 	recoProtonENonContained->Fill(energy);
						// 	recoProtonLength_vs_ENonContained->Fill(length, energy);
						// }

						
					}

					if (abs(pdg) == 211) {
						auto start_pos = sr->common.ixn.dlp[nixn].part.dlp[npart].start;
						auto end_pos   = sr->common.ixn.dlp[nixn].part.dlp[npart].end;
												
						bool truth_contained_geom =
							(start_pos.x > -55 && start_pos.x < 55) &&
							(start_pos.y > -55 && start_pos.y < 55) &&
							(start_pos.z > -55 && start_pos.z < 55) &&
							(end_pos.x   > -55 && end_pos.x   < 55) &&
							(end_pos.y   > -55 && end_pos.y   < 55) &&
							(end_pos.z   > -55 && end_pos.z   < 55);

                        if (!truth_contained_geom) continue;	// skip non-contained

						
						auto momentum   = sr->common.ixn.dlp[nixn].part.dlp[npart].p;

						double dX = end_pos.x - start_pos.x;
						double dY = end_pos.y - start_pos.y;
						double dZ = end_pos.z - start_pos.z;

						double length = TMath::Sqrt(dX * dX + dY * dY + dZ * dZ);
						double dirZ = dZ / length;

						auto energy = sr->common.ixn.dlp[nixn].part.dlp[npart].E;

                         if (length <= 2.0) continue;   // require proton track length > 2 cm
                        // if (end_pos.z >= -5.0 && end_pos.z <= 5.0) continue; // reject |z_end| ≤ 5 cm → keep outside 
						if (end_pos.z >= -5.0 && end_pos.z <= 0.0) continue; //  reject only −z side (−5 to 0 cm):

                  		// Count all reco protons
						nRecoPionAll++;

                        
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
						
						// if (truth_contained_geom) {
						// 	nRecoPionContained++;
							
						// 	recoPionLengthContained->Fill(length);
						// 	recoPionEContained->Fill(energy);
						// 	recoPionLength_vs_EContained->Fill(length, energy);
						// }

						// if (!truth_contained_geom) {
						// 	nRecoPionNonContained++;
							
						// 	recoPionLengthNonContained->Fill(length);
						// 	recoPionENonContained->Fill(energy);
						// 	recoPionLength_vs_ENonContained->Fill(length, energy);
						// }


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
					
				}
		} // end for interaction	
	} // end for spills




	// std::cout << "Minerva Match Purity: " << goodMINERvAMatch << "/" << totalMINERvAMatch << std::endl;

	// Create output file and write your histograms
	TFile* caf_out_file = new TFile(output_rootfile.c_str(), "recreate");


	std::cout << "\n========================================\n";
	std::cout << "Reco Proton summary\n";
	std::cout << "----------------------------------------\n";
	std::cout << "Total reco Protons        : " << nRecoProtonAll << "\n";
	std::cout << "Contained reco Protons    : " << nRecoProtonContained << "\n";
	std::cout << "Non-contained reco Protons: " << nRecoProtonNonContained << "\n";

	if (nRecoProtonAll > 0) {
		std::cout << "Contained fraction        : "
				  << double(nRecoProtonContained) / nRecoProtonAll << "\n";
		std::cout << "Non-contained fraction    : "
				  << double(nRecoProtonNonContained) / nRecoProtonAll << "\n";
	}
	
	if (nRecoProtonAll != nRecoProtonContained + nRecoProtonNonContained) {
		std::cerr << "[WARNING] Proton containment counting mismatch!" << std::endl;
	}

	std::cout << "========================================\n\n";
	
	std::cout << "\n========================================\n";
	std::cout << "Reco Pion summary\n";
	std::cout << "----------------------------------------\n";
	std::cout << "Total reco Pions        : " << nRecoPionAll << "\n";
	std::cout << "Contained reco Pions    : " << nRecoPionContained << "\n";
	std::cout << "Non-contained reco Pions: " << nRecoPionNonContained << "\n";

	if (nRecoPionAll > 0) {
		std::cout << "Contained fraction        : "
				  << double(nRecoPionContained) / nRecoPionAll << "\n";
		std::cout << "Non-contained fraction    : "
				  << double(nRecoPionNonContained) / nRecoPionAll << "\n";
	}
	
	if (nRecoPionAll != nRecoPionContained + nRecoPionNonContained) {
		std::cerr << "[WARNING] Pion containment counting mismatch!" << std::endl;
	}

	std::cout << "========================================\n\n";	
		
	recoVtxX->Write();
	recoVtxY->Write();
	recoVtxZ->Write();
	
	part_mult->Write();
	part_multTrkOnly->Write();
	part_energy_hist->Write();

	track_length->Write();

	recoVertex2D->Write();
	recoVertex2DNoCuts->Write();
	longest_track->Write();
	
	recoCosL->Write();

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

	// ==========================
	// Reco Proton — Contained
	// ==========================
	recoProtonLengthContained->Write();
	recoProtonEContained->Write();
	recoProtonLength_vs_EContained->Write();

	// ==========================
	// Reco Proton — Non-contained
	// ==========================
	recoProtonLengthNonContained->Write();
	recoProtonENonContained->Write();
	recoProtonLength_vs_ENonContained->Write();

	// ==========================
	// Reco Pion — Contained
	// ==========================
	recoPionLengthContained->Write();
	recoPionEContained->Write();
	recoPionLength_vs_EContained->Write();

	// ==========================
	// Reco Pion — Non-contained
	// ==========================
	recoPionLengthNonContained->Write();
	recoPionENonContained->Write();
	recoPionLength_vs_ENonContained->Write();

track_mult->Write();
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
