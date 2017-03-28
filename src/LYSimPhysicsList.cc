#include "globals.hh"
#include "LYSimPhysicsList.hh"
#include "LYSimPhysicsListMessenger.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4ProcessManager.hh"
#include "G4EmParameters.hh"

#include "G4Cerenkov.hh"
#include "LYSimScintillation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpWLS.hh"

#include "G4LossTableManager.hh"
#include "G4EmSaturation.hh"

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"

// EmPhys
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

    defaultCutValue = 0.0001*nm;
    cutForGamma     = 1*nm;
    cutForElectron  = 1*nm;
    cutForPositron  = 1*nm;
    cutForProton    = 1*nm;

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

#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

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
    param->SetMinEnergy(100*eV);
    param->SetMaxEnergy(100*MeV);

    auto theParticleIterator=GetParticleIterator();
    theParticleIterator->reset();
    while( (*theParticleIterator)() ){
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();

        if (particleName == "gamma") {
            // gamma
            // Construct processes for gamma
            pmanager->AddDiscreteProcess(new G4GammaConversion());
            pmanager->AddDiscreteProcess(new G4ComptonScattering());
            pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());

        } else if (particleName == "e-") {
            //electron
            // Construct processes for electron
            G4eIonisation* eIonisation = new G4eIonisation();
            eIonisation->SetStepFunction(0.1, 1*um);
            pmanager->AddProcess(new G4eMultipleScattering(),-1, 1, 1);
            //pmanager->AddProcess(new G4eIonisation(),        -1, 2, 2);
            pmanager->AddProcess(eIonisation,                -1, 2, 2);
            pmanager->AddProcess(new G4eBremsstrahlung(),    -1, 3, 3);

        } else if (particleName == "e+") {
            //positron
            // Construct processes for positron
            pmanager->AddProcess(new G4eMultipleScattering(),-1, 1, 1);
            pmanager->AddProcess(new G4eIonisation(),        -1, 2, 2);
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

        } else if( particleName == "GenericIon" ) {

            G4ionIonisation* ionIoni = new G4ionIonisation();
            ionIoni->SetEmModel(new G4IonParametrisedLossModel());
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
    // Deexcitation
    //
    G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
    de->SetFluo(true);
    de->SetAuger(true);
    de->SetPIXE(true);
    G4LossTableManager::Instance()->SetAtomDeexcitation(de);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LYSimPhysicsList::ConstructOp()
{
    theCerenkovProcess           = new G4Cerenkov("Cerenkov");
    theComptonScatteringProcess  = new G4ComptonScattering();
    theScintillationProcess      = new LYSimScintillation("Scintillation");
    theAbsorptionProcess         = new G4OpAbsorption();
    theRayleighScatteringProcess = new G4OpRayleigh();
    theMieHGScatteringProcess    = new G4OpMieHG();
    theBoundaryProcess           = new G4OpBoundaryProcess();
    theWLSProcess                = new G4OpWLS();

    //theCerenkovProcess->DumpPhysicsTable();
    theScintillationProcess->DumpPhysicsTable();
    //theRayleighScatteringProcess->DumpPhysicsTable();

    //Set verbose level to 0
    SetVerbose(0);

    theCerenkovProcess->SetMaxNumPhotonsPerStep(200);
    theCerenkovProcess->SetMaxBetaChangePerStep(10.0);
    theCerenkovProcess->SetTrackSecondariesFirst(true);

    theScintillationProcess->SetScintillationByParticleType(true);
    theScintillationProcess->SetScintillationYieldFactor(1.);
    theScintillationProcess->SetScintillationExcitationRatio(1.);
    theScintillationProcess->SetTrackSecondariesFirst(true);

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
            pmanager->AddDiscreteProcess(theAbsorptionProcess);
            pmanager->AddDiscreteProcess(theComptonScatteringProcess);
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
    theScintillationProcess->SetScintillationByParticleType(true);

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
