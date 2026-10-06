#ifndef G4Reweighter_h
#define G4Reweighter_h

#include "fhiclcpp/ParameterSet.h"

class G4ReweightManager;
class G4ReweightTraj;

class G4Box;
class G4CascadeDetectorConstruction;
//class G4CascadePhysicsList;
class G4DynamicParticle;
class G4HadronElasticProcess;
class G4HadronInelasticProcess;
class G4LogicalVolume;
class G4Material;
class G4ParticleDefinition;
class G4PionPlus;
class G4RunManager;
class G4Step;
class G4StepPoint;
class G4Track;
class G4VPhysicalVolume;

class TFile;
class TGraph;
class TH1D;

#include <map>
#include <set>
#include <string>
#include <vector>

class G4Reweighter{
  public:

    G4Reweighter(){};
    G4Reweighter(TFile *, const std::map<std::string, TH1D*> &,
                 const fhicl::ParameterSet & material_pars,
                 G4ReweightManager * rw_manager,
                 std::vector<std::string> the_ints,
                 TH1D * inputElasticBiasHist = nullptr, bool fix = false);
    virtual ~G4Reweighter();

    double GetWeight(const G4ReweightTraj * theTraj);
    virtual std::string GetInteractionSubtype(const G4ReweightTraj &);

    void SetMomentum(double p);
    double GetNominalMFP(double p);
    double GetBiasedMFP(double p);
    double GetNominalElasticMFP(double p);
    double GetBiasedElasticMFP(double p);
    double GetInelasticBias(double p);
    double GetElasticBias(double p);
    double GetExclusiveFactor(double p, std::string cut);

    double GetInelasticXSec(double p, bool mb_units = false);
    double GetExclusiveXSec(double p, std::string cut, bool mb_units = false);
    double GetElasticXSec(double p, bool mb_units = false);

    void SetNewHists(const std::map< std::string, TH1D* > &FSScales);
    void SetNewElasticHists(TH1D * inputElasticBiasHist);

    void SetInelasticPreBias(double bias) {fInelasticPreBias = bias;};
    void SetElasticPreBias(double bias) {fElasticPreBias = bias;};

    // Process names treated as the paper's case 4, "some other interaction"
    // (eq. 3.11). A step whose chosen process is in this set gets an
    // interaction weight with a scale factor of 1, not a survival weight.
    // Defaults to the two processes the paper folds into sigma_X: Coulomb
    // scattering and decay. Steps recording no interaction (e.g. "default",
    // "Transportation") must NOT be listed here.
    const std::set<std::string> & GetOtherInteractions() const {
        return fOtherInteractions;};
    void SetOtherInteractions(const std::set<std::string> & procs) {
        fOtherInteractions = procs;};

  protected:

    bool fix_total = false;

    std::map<std::string, TGraph *> exclusiveFracs;
    std::map<std::string, TH1D *> inelScales;

    fhicl::ParameterSet MaterialParameters;
    G4ReweightManager * RWManager;
    TH1D * elasticBias;

    // These should be set in the constructor of the actual reweighter you use (e.g. G4PiPlusReweighter/G4PiMinusReweighter/G4ProtonReweighter)
    std::string fInelastic/* = "pi+Inelastic"*/;

    G4RunManager * rm;
    G4Track * testTrack;
    G4Step * testStep;
    G4StepPoint * testPoint;
    G4Material * testMaterial;
    G4Box * solidWorld;
    G4LogicalVolume * logicWorld;
    G4VPhysicalVolume * physWorld;
    //G4CascadeDetectorConstruction * detector;
    //G4CascadePhysicsList * physList;

    G4PionPlus  * piplus;
    G4ParticleDefinition * part_def;
    G4DynamicParticle * dynamic_part;

    G4HadronElasticProcess * elastic_proc;
    G4HadronInelasticProcess * inelastic_proc;

    void SetupProcesses();
    //void SetupWorld();
    void SetupParticle();

    // --- eq. 3.9 / 3.10 / 3.11 helpers -------------------------------------
    // sigma  = sigma_I + sigma_E          (sigma_X is taken as 0, see .cc)
    double GetNominalTotalXSec(double p);
    // sigma' = sigma'_I + sigma'_E
    double GetBiasedTotalXSec(double p);
    // dL*(sigma - sigma'), the eq. 3.9 exponent for a non-interacting step
    double GetSurvivalExponentTerm(double p, double stepLength);
    // (sigma/sigma')*(1 - e^{-dL sigma'})/(1 - e^{-dL sigma}), the factor
    // shared by eq. 3.10 and eq. 3.11
    double GetInteractionRateFactor(double p, double stepLength);
    bool IsOtherInteraction(const std::string & proc) const;

    std::set<std::string> fOtherInteractions = {"Decay", "CoulombScat"};

    double fInelasticPreBias = 1.;
    double fElasticPreBias = 1.;

};

#endif
