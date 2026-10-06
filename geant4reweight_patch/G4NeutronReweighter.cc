#include "geant4reweight/ReweightBase/G4NeutronReweighter.hh"
#include "geant4reweight/ReweightBase/G4ReweightTraj.hh"

#include "Geant4/G4Neutron.hh"

G4NeutronReweighter::G4NeutronReweighter(
    TFile * FSInput,
    const std::map<std::string, TH1D*> &FSScales,
    const fhicl::ParameterSet & material_pars,
    G4ReweightManager * rw_manager,
    TH1D * inputElasticBiasHist, bool fix)
  : G4Reweighter(FSInput, FSScales, material_pars, rw_manager,
                 {"0n0p", "0n1p", "1n0p", "1n1p", "Other"},
                 inputElasticBiasHist, fix) {

  part_def = neutron->Definition();
  fInelastic = "neutronInelastic";
  SetupProcesses();
}

std::string G4NeutronReweighter::GetInteractionSubtype(
    const G4ReweightTraj & theTraj) {
  int nProtons  = theTraj.HasChild(2212).size();
  int nNeutrons = theTraj.HasChild(2112).size();
  if      (nNeutrons == 0 && nProtons == 0) return "0n0p";
  else if (nNeutrons == 0 && nProtons == 1) return "0n1p";
  else if (nNeutrons == 1 && nProtons == 0) return "1n0p";
  else if (nNeutrons == 1 && nProtons == 1) return "1n1p";
  else                                       return "Other";
}

G4NeutronReweighter::~G4NeutronReweighter(){}