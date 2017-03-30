#include "globals.hh"
#include "LYSimPhysicsList.hh"
#include "LYSimPhysicsListMessenger.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4ProcessManager.hh"

#include "G4Cerenkov.hh"
#include "LYSimScintillation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpWLS.hh"

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"

// EmPhys
#include "G4EmParameters.hh"
#include "G4LossTableManager.hh"
#include "G4EmSaturation.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"

using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LYSimPhysicsList::LYSimPhysicsList() : G4VUserPhysicsList()
{
    pMessenger = new LYSimPhysicsListMessenger(this);

    defaultCutValue = 0.1*nm;
    cutForGamma     = 1*nm;
    cutForElectron  = 10*nm;
    cutForPositron  = 10*nm;
    cutForProton    = 1*um;

    theCerenkovProcess           = NULL;
    theScintillationProcess      = NULL;
    theAbsorptionProcess         = NULL;
    theRayleighScatteringProcess = NULL;
    theMieHGScatteringProcess    = NULL;
    theBoundaryProcess           = NULL;
    theWLSProcess                = NULL;

    SetVerboseLevel(0);

    // EM physics
    emName = G4String("emlivermore");
    emPhysicsList = new G4EmLivermorePhysics;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LYSimPhysicsList::~LYSimPhysicsList() {
    delete emPhysicsList;
    delete pMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LYSimPhysicsList::ConstructParticle()
{
    // In this method, static member functions should be called
    // for all particles which you want to use.
    // This ensures that objects of these particle types will be
    // created in the program.

    ConstructBosons();
    ConstructLeptons();
    ConstructMesons();
    ConstructBaryons();
    ConstructIons();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LYSimPhysicsList::ConstructBosons()
{
    // bosons, pseudo-particles, ...
    G4BosonConstructor bConstructor;
    bConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LYSimPhysicsList::ConstructLeptons()
{
    // leptons
    G4LeptonConstructor lConstructor;
    lConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LYSimPhysicsList::ConstructMesons()
{
    //  mesons
    G4MesonConstructor mConstructor;
    mConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LYSimPhysicsList::ConstructBaryons()
{
    //  barions
    G4BaryonConstructor bConstructor;
    bConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LYSimPhysicsList::ConstructIons()
{
    // ions
    G4IonConstructor iConstructor;
    iConstructor.ConstructParticle();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LYSimPhysicsList::ConstructProcess()
{
    AddTransportation();
    ConstructGeneral();
    ConstructEM();
    ConstructHad();
    ConstructOp();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Decay.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LYSimPhysicsList::ConstructGeneral()
{
    // Add Decay Process
    G4Decay* theDecayProcess = new G4Decay();
    auto theParticleIterator=GetParticleIterator();
    theParticleIterator->reset();
    while( (*theParticleIterator)() ){
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        if (theDecayProcess->IsApplicable(*particle)) {
            pmanager ->AddProcess(theDecayProcess);
            // set ordering for PostStepDoIt and AtRestDoIt
            pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
            pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4LivermoreRayleighModel.hh"

// e-
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4LivermoreIonisationModel.hh"
#include "G4LivermoreBremsstrahlungModel.hh"


// e+
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

// muon
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuMultipleScattering.hh"

// others
#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4NuclearStopping.hh"
#include "G4hImpactIonisation.hh"
#include "G4IonParametrisedLossModel.hh"

#include "G4UAtomicDeexcitation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LYSimPhysicsList::ConstructEM()
{
    // Update EM options
    emPhysicsList->ConstructProcess();
    G4EmParameters* param = G4EmParameters::Instance();
    param->SetMinEnergy(1.*eV);
    param->SetMaxEnergy(1.*GeV);
    param->SetNumberOfBinsPerDecade(100);
    param->SetMscStepLimitType(fMinimal);
    param->SetFluo(true);
    param->SetPixe(true);
    param->SetAuger(true);

    // Deexcitation
    G4LossTableManager* ltman = G4LossTableManager::Instance();
    G4VAtomDeexcitation* ad = ltman->AtomDeexcitation();
    if(!ad) {
        G4UAtomicDeexcitation* de = new G4UAtomicDeexcitation();
        de->SetFluo(true);
        de->SetAuger(true);
        de->SetPIXE(true);
        ltman->SetAtomDeexcitation(de);
    }

//     G4UAtomicDeexcitation* de = new G4UAtomicDeexcitation();
//     de->SetFluo(true);
//     de->SetAuger(true);
//     de->SetPIXE(true);
//     G4LossTableManager::Instance()->SetAtomDeexcitation(de);


    auto theParticleIterator=GetParticleIterator();
    theParticleIterator->reset();
    while( (*theParticleIterator)() ){
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();

        if (particleName == "gamma") {
            // gamma
            // Construct processes for gamma
            pmanager->AddDiscreteProcess(new G4RayleighScattering());

//             pmanager->AddDiscreteProcess(new G4GammaConversion());
//             pmanager->AddDiscreteProcess(new G4ComptonScattering());
//             pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());
            // Livermore
            G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
            thePhotoElectricEffect->SetEmModel(new G4LivermorePhotoElectricModel());
            pmanager->AddDiscreteProcess(thePhotoElectricEffect);

            G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
            theComptonScattering->SetEmModel(new G4LivermoreComptonModel());
            pmanager->AddDiscreteProcess(theComptonScattering);

            G4GammaConversion* theGammaConversion = new G4GammaConversion();
            theGammaConversion->SetEmModel(new G4LivermoreGammaConversionModel());
            pmanager->AddDiscreteProcess(theGammaConversion);

        } else if (particleName == "e-") {
            //electron
            // Construct processes for electron
            G4eIonisation* eIonisation = new G4eIonisation();
            eIonisation->SetStepFunction(0.2, 100*um);
            pmanager->AddProcess(new G4eMultipleScattering(),-1, 1, 1);
            //pmanager->AddProcess(new G4eIonisation(),        -1, 2, 2);
            pmanager->AddProcess(eIonisation,                -1, 2, 2);
            pmanager->AddProcess(new G4eBremsstrahlung(),    -1, 3, 3);

        } else if (particleName == "e+") {
            //positron
            // Construct processes for positron
            pmanager->AddProcess(new G4eMultipleScattering(),-1, 1, 1);
            G4eIonisation* eIonisation = new G4eIonisation();
            eIonisation->SetStepFunction(0.2, 100*um);
            pmanager->AddProcess(eIonisation,                -1, 2, 2);
            //pmanager->AddProcess(new G4eIonisation(),        -1, 2, 2);
            pmanager->AddProcess(new G4eBremsstrahlung(),    -1, 3, 3);
            pmanager->AddProcess(new G4eplusAnnihilation(),   0,-1, 4);

        } else if(particleName == "mu+" ||
                  particleName == "mu-") {
            //muon
            // Construct processes for muon
            pmanager->AddProcess(new G4MuMultipleScattering(),-1, 1, 1);
            pmanager->AddProcess(new G4MuIonisation(),        -1, 2, 2);
            pmanager->AddProcess(new G4MuBremsstrahlung(),    -1, 3, 3);
            pmanager->AddProcess(new G4MuPairProduction(),    -1, 4, 4);

        } else if (//particleName == "alpha" ||
                   particleName == "proton"||
                   particleName == "anti_proton") {

            // Instantiate the G4hImpactIonisation process
            G4hImpactIonisation* hIonisation = new G4hImpactIonisation();
            G4cout << "[LYSim] Add Hadron Impact Ionisation and PIXE to " << particleName << G4endl;
            // Select the cross section models to be applied for K, L and M shell vacancy creation
            // (here the ECPSSR model is selected for K, L and M shell; one can mix and match
            // different models for each shell)
            hIonisation->SetPixeCrossSectionK("ecpssr");
            hIonisation->SetPixeCrossSectionL("ecpssr");
            hIonisation->SetPixeCrossSectionM("ecpssr");

            // Register the process with the processManager associated with alpha
            pmanager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
            pmanager->AddProcess(hIonisation,                 -1, 2, 2);

        } else if(particleName == "alpha" ||
                  particleName == "He3") {

            //G4ionIonisation* ionIoni = new G4ionIonisation();
            //ionIoni->SetStepFunction(0.1, 0.01*um);

            G4hImpactIonisation* hIonisation = new G4hImpactIonisation();
            G4cout << "[LYSim] Add Hadron Impact Ionisation and PIXE to " << particleName << G4endl;
            // Select the cross section models to be applied for K, L and M shell vacancy creation
            // (here the ECPSSR model is selected for K, L and M shell; one can mix and match
            // different models for each shell)
            hIonisation->SetPixeCrossSectionK("ecpssr");
            hIonisation->SetPixeCrossSectionL("ecpssr");
            hIonisation->SetPixeCrossSectionM("ecpssr");
            hIonisation->SetStepFunction(0.01, 0.01*um);

            //pmanager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
            //pmanager->AddProcess(ionIoni,                     -1, 2, 2);
            pmanager->AddProcess(hIonisation,                 -1, 2, 2);
            //pmanager->AddProcess(new G4NuclearStopping(),     -1, 3,-1);

        } else if(particleName == "GenericIon") {

            G4ionIonisation* ionIoni = new G4ionIonisation();
            ionIoni->SetEmModel(new G4IonParametrisedLossModel());//ICRU 73
            ionIoni->SetStepFunction(0.1, 1*um);
            pmanager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
            pmanager->AddProcess(ionIoni,                     -1, 2, 2);
            pmanager->AddProcess(new G4NuclearStopping(),     -1, 3,-1);


        } else {
            if ((particle->GetPDGCharge() != 0.0) &&
                (particle->GetParticleName() != "chargedgeantino")) {
                // all others charged particles except geantino
                pmanager->AddProcess(new G4hMultipleScattering(),-1, 1, 1);
                pmanager->AddProcess(new G4hIonisation(),        -1, 2, 2);
            }
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LYSimPhysicsList::ConstructOp()
{
    theCerenkovProcess           = new G4Cerenkov("Cerenkov");
    theScintillationProcess      = new LYSimScintillation("Scintillation");
    theAbsorptionProcess         = new G4OpAbsorption();
    theRayleighScatteringProcess = new G4OpRayleigh();
    theMieHGScatteringProcess    = new G4OpMieHG();
    theBoundaryProcess           = new G4OpBoundaryProcess();
    theWLSProcess                = new G4OpWLS();

    //theCerenkovProcess->DumpPhysicsTable();
    //theScintillationProcess->DumpPhysicsTable();
    //theRayleighScatteringProcess->DumpPhysicsTable();

    //Set verbose level to 0
    SetVerbose(0);

    // Cerenkov process options
    theCerenkovProcess->SetMaxNumPhotonsPerStep(200);
    theCerenkovProcess->SetMaxBetaChangePerStep(10.0);
    theCerenkovProcess->SetTrackSecondariesFirst(true);

    // Note: either SetScintillationByParticleType or use Birks's correction
    // Scintillation process options
    theScintillationProcess->SetFiniteRiseTime(true);
    theScintillationProcess->SetTrackSecondariesFirst(true);

    // Scintillation process options for alpha:
//     G4Scintillation* theScintProcessAlpha = new G4Scintillation("ScintillationA");
//     theScintProcessAlpha->SetFiniteRiseTime(true);
//     theScintProcessAlpha->SetTrackSecondariesFirst(true);
//     theScintProcessAlpha->SetScintillationYieldFactor(1.1);//default: 1.0
//     //theScintProcessAlpha->SetScintillationExcitationRatio(1.0);//default: 1.0
//     //theScintProcessAlpha->DumpPhysicsTable();

    // Use Birks Correction in the Scintillation process
    G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
    theScintillationProcess->AddSaturation(emSaturation);
    //theScintProcessAlpha->AddSaturation(emSaturation);
    //theScintillationProcess->SetScintillationByParticleType(true);
    //theScintProcessAlpha->SetScintillationByParticleType(true);

    theWLSProcess->UseTimeProfile("delta");

    auto theParticleIterator=GetParticleIterator();
    theParticleIterator->reset();
    while( (*theParticleIterator)() ){
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();
        if (theCerenkovProcess->IsApplicable(*particle)) {
            pmanager->AddProcess(theCerenkovProcess);
            pmanager->SetProcessOrdering(theCerenkovProcess,idxPostStep);
        }
        if (theScintillationProcess->IsApplicable(*particle)) {
//             if(particle->GetParticleName() == "alpha") {
//                 pmanager->AddProcess(theScintProcessAlpha);
//                 pmanager->SetProcessOrderingToLast(theScintProcessAlpha, idxAtRest);
//                 pmanager->SetProcessOrderingToLast(theScintProcessAlpha, idxPostStep);
//             }
//             else {
            pmanager->AddProcess(theScintillationProcess);
            pmanager->SetProcessOrderingToLast(theScintillationProcess, idxAtRest);
            pmanager->SetProcessOrderingToLast(theScintillationProcess, idxPostStep);
//             }
        }
        if (particleName == "opticalphoton") {
            G4cout << "[LYSim] AddDiscreteProcess to OpticalPhoton " << G4endl;
            pmanager->AddDiscreteProcess(theAbsorptionProcess);
            pmanager->AddDiscreteProcess(theRayleighScatteringProcess);
            pmanager->AddDiscreteProcess(theMieHGScatteringProcess);
            pmanager->AddDiscreteProcess(theBoundaryProcess);
            pmanager->AddDiscreteProcess(theWLSProcess);
        }
    }
}


void LYSimPhysicsList::ConstructIdealOp()
//Use instead of ConstructOp to only activate WLS and Boundary processes.
{
    theCerenkovProcess           = new G4Cerenkov("Cerenkov");
    theScintillationProcess      = new LYSimScintillation("Scintillation");
    theBoundaryProcess           = new G4OpBoundaryProcess();
    theWLSProcess                = new G4OpWLS();

    //theCerenkovProcess->DumpPhysicsTable();
    //theScintillationProcess->DumpPhysicsTable();

    //Verbose level set to 0
    SetVerbose(0);

    theCerenkovProcess->SetMaxNumPhotonsPerStep(20);
    theCerenkovProcess->SetMaxBetaChangePerStep(10.0);
    theCerenkovProcess->SetTrackSecondariesFirst(true);

    theScintillationProcess->SetScintillationYieldFactor(1.);
    theScintillationProcess->SetTrackSecondariesFirst(true);
    //theScintillationProcess->SetScintillationByParticleType(true);

    // Use Birks Correction in the Scintillation process
    G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
    theScintillationProcess->AddSaturation(emSaturation);

    theWLSProcess->UseTimeProfile("delta");

    auto theParticleIterator=GetParticleIterator();
    theParticleIterator->reset();
    while( (*theParticleIterator)() ){
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();
        if (theCerenkovProcess->IsApplicable(*particle)) {
            pmanager->AddProcess(theCerenkovProcess);
            pmanager->SetProcessOrdering(theCerenkovProcess,idxPostStep);
        }
        if (theScintillationProcess->IsApplicable(*particle)) {
            pmanager->AddProcess(theScintillationProcess);
            pmanager->SetProcessOrderingToLast(theScintillationProcess, idxAtRest);
            pmanager->SetProcessOrderingToLast(theScintillationProcess, idxPostStep);
        }
        if (particleName == "opticalphoton") {
            G4cout << "[LYSim] AddDiscreteProcess to OpticalPhoton " << G4endl;
            pmanager->AddDiscreteProcess(theBoundaryProcess);
            pmanager->AddDiscreteProcess(theWLSProcess);
        }
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Hadronic processes ////////////////////////////////////////////////////////

// Elastic processes:
#include "G4HadronElasticProcess.hh"
#include "G4ChipsElasticModel.hh"
#include "G4ElasticHadrNucleusHE.hh"

// Inelastic processes:
#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4AntiNeutronInelasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"

// High energy FTFP model and Bertini cascade
#include "G4FTFModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4TheoFSGenerator.hh"
#include "G4CascadeInterface.hh"

// Cross sections
#include "G4VCrossSectionDataSet.hh"
#include "G4CrossSectionDataSetRegistry.hh"

#include "G4CrossSectionElastic.hh"
#include "G4BGGPionElasticXS.hh"
#include "G4AntiNuclElastic.hh"

#include "G4CrossSectionInelastic.hh"
#include "G4PiNuclearCrossSection.hh"
#include "G4CrossSectionPairGG.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"
#include "G4ComponentGGNuclNuclXsc.hh"

#include "G4HadronElastic.hh"
#include "G4HadronCaptureProcess.hh"

// Neutron high-precision models: <20 MeV
#include "G4ParticleHPElastic.hh"
#include "G4ParticleHPElasticData.hh"
#include "G4ParticleHPCapture.hh"
#include "G4ParticleHPCaptureData.hh"
#include "G4ParticleHPInelastic.hh"
#include "G4ParticleHPInelasticData.hh"

// Stopping processes
#include "G4PiMinusAbsorptionBertini.hh"
#include "G4KaonMinusAbsorptionBertini.hh"
#include "G4AntiProtonAbsorptionFritiof.hh"

void LYSimPhysicsList::ConstructHad()
{
  //Elastic models
  const G4double elastic_elimitPi = 1.0*GeV;

  G4HadronElastic* elastic_lhep0 = new G4HadronElastic();
  G4HadronElastic* elastic_lhep1 = new G4HadronElastic();
  elastic_lhep1->SetMaxEnergy( elastic_elimitPi );
  G4ChipsElasticModel* elastic_chip = new G4ChipsElasticModel();
  G4ElasticHadrNucleusHE* elastic_he = new G4ElasticHadrNucleusHE();
  elastic_he->SetMinEnergy( elastic_elimitPi );

  
  // Inelastic scattering
  const G4double theFTFMin0 =    0.0*GeV;
  const G4double theFTFMin1 =    4.0*GeV;
  const G4double theFTFMax =   100.0*TeV;
  const G4double theBERTMin0 =   0.0*GeV;
  const G4double theBERTMin1 =  19.0*MeV;
  const G4double theBERTMax =    5.0*GeV;
  const G4double theHPMin =      0.0*GeV;
  const G4double theHPMax =     20.0*MeV;

  G4FTFModel* theStringModel = new G4FTFModel;
  G4ExcitedStringDecay* theStringDecay = new G4ExcitedStringDecay( new G4LundStringFragmentation );
  theStringModel->SetFragmentationModel( theStringDecay );
  G4PreCompoundModel* thePreEquilib = new G4PreCompoundModel( new G4ExcitationHandler );
  G4GeneratorPrecompoundInterface* theCascade = new G4GeneratorPrecompoundInterface( thePreEquilib );

  G4TheoFSGenerator* theFTFModel0 = new G4TheoFSGenerator( "FTFP" );
  theFTFModel0->SetHighEnergyGenerator( theStringModel );
  theFTFModel0->SetTransport( theCascade );
  theFTFModel0->SetMinEnergy( theFTFMin0 );
  theFTFModel0->SetMaxEnergy( theFTFMax );

  G4TheoFSGenerator* theFTFModel1 = new G4TheoFSGenerator( "FTFP" );
  theFTFModel1->SetHighEnergyGenerator( theStringModel );
  theFTFModel1->SetTransport( theCascade );
  theFTFModel1->SetMinEnergy( theFTFMin1 );
  theFTFModel1->SetMaxEnergy( theFTFMax );

  G4CascadeInterface* theBERTModel0 = new G4CascadeInterface;
  theBERTModel0->SetMinEnergy( theBERTMin0 );
  theBERTModel0->SetMaxEnergy( theBERTMax );

  G4CascadeInterface* theBERTModel1 = new G4CascadeInterface;
  theBERTModel1->SetMinEnergy( theBERTMin1 );
  theBERTModel1->SetMaxEnergy( theBERTMax );

  G4VCrossSectionDataSet* thePiData = new G4CrossSectionPairGG( new G4PiNuclearCrossSection, 91*GeV );
  G4VCrossSectionDataSet* theAntiNucleonData = new G4CrossSectionInelastic( new G4ComponentAntiNuclNuclearXS );
  G4ComponentGGNuclNuclXsc* ggNuclNuclXsec = new G4ComponentGGNuclNuclXsc();
  G4VCrossSectionDataSet* theGGNuclNuclData = new G4CrossSectionInelastic(ggNuclNuclXsec);

  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while ((*particleIterator)())
    {
      G4ParticleDefinition* particle = particleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();

      if (particleName == "pi+") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( new G4BGGPionElasticXS( particle ) );
          theElasticProcess->RegisterMe( elastic_lhep1 );
          theElasticProcess->RegisterMe( elastic_he );
	  pmanager->AddDiscreteProcess( theElasticProcess );
	  //Inelastic scattering
	  G4PionPlusInelasticProcess* theInelasticProcess = 
	    new G4PionPlusInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( thePiData );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	} 

      else if (particleName == "pi-") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( new G4BGGPionElasticXS( particle ) );
          theElasticProcess->RegisterMe( elastic_lhep1 );
          theElasticProcess->RegisterMe( elastic_he );
	  pmanager->AddDiscreteProcess( theElasticProcess );
	  //Inelastic scattering
	  G4PionMinusInelasticProcess* theInelasticProcess = 
	    new G4PionMinusInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( thePiData );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );	  
	  //Absorption
	  pmanager->AddRestProcess(new G4PiMinusAbsorptionBertini, ordDefault);
	}
      
      else if (particleName == "kaon+") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering	
	  G4KaonPlusInelasticProcess* theInelasticProcess = 
	    new G4KaonPlusInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
					   GetCrossSectionDataSet(G4ChipsKaonPlusInelasticXS::Default_Name()));
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}
      
      else if (particleName == "kaon0S") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering	 
	  G4KaonZeroSInelasticProcess* theInelasticProcess = 
	    new G4KaonZeroSInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
					   GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name()));
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );	  
	}

      else if (particleName == "kaon0L") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
	  // Inelastic scattering
	  G4KaonZeroLInelasticProcess* theInelasticProcess = 
	    new G4KaonZeroLInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
					   GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name()));
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );	  
	}

      else if (particleName == "kaon-") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4KaonMinusInelasticProcess* theInelasticProcess = 
	    new G4KaonMinusInelasticProcess("inelastic");	
          theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
					   GetCrossSectionDataSet(G4ChipsKaonMinusInelasticXS::Default_Name()));
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	  pmanager->AddRestProcess(new G4KaonMinusAbsorptionBertini, ordDefault);
	}

      else if (particleName == "proton") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->
					GetCrossSectionDataSet(G4ChipsProtonElasticXS::Default_Name()));
          theElasticProcess->RegisterMe( elastic_chip );
	  pmanager->AddDiscreteProcess( theElasticProcess );
	  // Inelastic scattering
	  G4ProtonInelasticProcess* theInelasticProcess = 
	    new G4ProtonInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( new G4BGGNucleonInelasticXS( G4Proton::Proton() ) );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}
      else if (particleName == "anti_proton") 
	{
	  // Elastic scattering
          const G4double elastic_elimitAntiNuc = 100.0*MeV;
          G4AntiNuclElastic* elastic_anuc = new G4AntiNuclElastic();
          elastic_anuc->SetMinEnergy( elastic_elimitAntiNuc );
          G4CrossSectionElastic* elastic_anucxs = new G4CrossSectionElastic( elastic_anuc->GetComponentCrossSection() );
          G4HadronElastic* elastic_lhep2 = new G4HadronElastic();
          elastic_lhep2->SetMaxEnergy( elastic_elimitAntiNuc );
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( elastic_anucxs );
          theElasticProcess->RegisterMe( elastic_lhep2 );
          theElasticProcess->RegisterMe( elastic_anuc );
	  pmanager->AddDiscreteProcess( theElasticProcess );
	  // Inelastic scattering
	  G4AntiProtonInelasticProcess* theInelasticProcess = 
	    new G4AntiProtonInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( theAntiNucleonData );
	  theInelasticProcess->RegisterMe( theFTFModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	  // Absorption
	  pmanager->AddRestProcess(new G4AntiProtonAbsorptionFritiof, ordDefault);
	}

      else if (particleName == "neutron") {
	// elastic scattering
	G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
        theElasticProcess->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsNeutronElasticXS::Default_Name()));
        G4HadronElastic* elastic_neutronChipsModel = new G4ChipsElasticModel();
	elastic_neutronChipsModel->SetMinEnergy( 19.0*MeV );
        theElasticProcess->RegisterMe( elastic_neutronChipsModel );
	G4ParticleHPElastic* theElasticNeutronHP = new G4ParticleHPElastic;
        theElasticNeutronHP->SetMinEnergy( theHPMin );
        theElasticNeutronHP->SetMaxEnergy( theHPMax );
	theElasticProcess->RegisterMe( theElasticNeutronHP );
	theElasticProcess->AddDataSet( new G4ParticleHPElasticData );
	pmanager->AddDiscreteProcess( theElasticProcess );
	// inelastic scattering		
	G4NeutronInelasticProcess* theInelasticProcess =
	  new G4NeutronInelasticProcess("inelastic");
	theInelasticProcess->AddDataSet( new G4BGGNucleonInelasticXS( G4Neutron::Neutron() ) );
	theInelasticProcess->RegisterMe( theFTFModel1 );
        theInelasticProcess->RegisterMe( theBERTModel1 );
	G4ParticleHPInelastic* theNeutronInelasticHPModel = new G4ParticleHPInelastic;
        theNeutronInelasticHPModel->SetMinEnergy( theHPMin );
        theNeutronInelasticHPModel->SetMaxEnergy( theHPMax );
	theInelasticProcess->RegisterMe( theNeutronInelasticHPModel );
	theInelasticProcess->AddDataSet( new G4ParticleHPInelasticData );
	pmanager->AddDiscreteProcess(theInelasticProcess);
	// capture
	G4HadronCaptureProcess* theCaptureProcess =
	  new G4HadronCaptureProcess;
	G4ParticleHPCapture* theLENeutronCaptureModel = new G4ParticleHPCapture;
	theLENeutronCaptureModel->SetMinEnergy(theHPMin);
	theLENeutronCaptureModel->SetMaxEnergy(theHPMax);
	theCaptureProcess->RegisterMe(theLENeutronCaptureModel);
	theCaptureProcess->AddDataSet( new G4ParticleHPCaptureData);
	pmanager->AddDiscreteProcess(theCaptureProcess);

      }
      else if (particleName == "anti_neutron") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering (include annihilation on-fly)
	  G4AntiNeutronInelasticProcess* theInelasticProcess = 
	    new G4AntiNeutronInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( theAntiNucleonData );
	  theInelasticProcess->RegisterMe( theFTFModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );	  
	}

      else if (particleName == "deuteron") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4DeuteronInelasticProcess* theInelasticProcess = 
	    new G4DeuteronInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( theGGNuclNuclData );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}
      
      else if (particleName == "triton") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4TritonInelasticProcess* theInelasticProcess = 
	    new G4TritonInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( theGGNuclNuclData );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}
      else if (particleName == "alpha") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4AlphaInelasticProcess* theInelasticProcess = 
	    new G4AlphaInelasticProcess("inelastic");	 
          theInelasticProcess->AddDataSet( theGGNuclNuclData );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}

    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LYSimPhysicsList::SetVerbose(G4int verbose)
{
    theCerenkovProcess->SetVerboseLevel(verbose);
    theScintillationProcess->SetVerboseLevel(verbose);
    theAbsorptionProcess->SetVerboseLevel(verbose);
    theRayleighScatteringProcess->SetVerboseLevel(verbose);
    theMieHGScatteringProcess->SetVerboseLevel(verbose);
    theBoundaryProcess->SetVerboseLevel(verbose);
    theWLSProcess->SetVerboseLevel(2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LYSimPhysicsList::SetNbOfPhotonsCerenkov(G4int MaxNumber)
{  
    theCerenkovProcess->SetMaxNumPhotonsPerStep(MaxNumber);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LYSimPhysicsList::SetCuts()
{
    //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
    //   the default cut value for all particle types
    // 
    SetCutsWithDefault();

    //special for low energy physics
    G4double lowlimit=100.*eV;
    G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowlimit,100.*MeV);

    // specified
    SetCutValue(cutForGamma,    "gamma");
    SetCutValue(cutForElectron, "e-");
    SetCutValue(cutForPositron, "e+");
    SetCutValue(cutForProton,   "proton");

    //temp disabled cut value display
    //     if (verboseLevel>0)
    DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LYSimPhysicsList::AddPhysicsList(const G4String& name)
{
  if (verboseLevel>-1) {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }

  if (name == emName) return;

  if (name == "emlivermore") {

    emName = name;
    delete emPhysicsList;
    emPhysicsList = new G4EmLivermorePhysics;

  } else if (name == "emstandard") {
    emName = name;
    delete emPhysicsList;
    emPhysicsList = new G4EmStandardPhysics();

  } else if (name == "emstandard_opt1") {

    emName = name;
    delete emPhysicsList;
    emPhysicsList = new G4EmStandardPhysics_option1();

  } else if (name == "emstandard_opt2") {

    emName = name;
    delete emPhysicsList;
    emPhysicsList = new G4EmStandardPhysics_option2();

  } else if (name == "emstandard_opt3") {

    emName = name;
    delete emPhysicsList;
    emPhysicsList = new G4EmStandardPhysics_option3();
    
  } else if (name == "empenelope"){
    emName = name;
    delete emPhysicsList;
    emPhysicsList = new G4EmPenelopePhysics();

  } else {

    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined"
           << G4endl;
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LYSimPhysicsList::SetCutForGamma(G4double cut)
{
  cutForGamma = cut;
  SetParticleCuts(cutForGamma, G4Gamma::Gamma());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LYSimPhysicsList::SetCutForElectron(G4double cut)
{
  cutForElectron = cut;
  SetParticleCuts(cutForElectron, G4Electron::Electron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LYSimPhysicsList::SetCutForPositron(G4double cut)
{
  cutForPositron = cut;
  SetParticleCuts(cutForPositron, G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LYSimPhysicsList::SetCutForProton(G4double cut)
{
  cutForProton = cut;
  SetParticleCuts(cutForProton, G4Proton::Proton());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
