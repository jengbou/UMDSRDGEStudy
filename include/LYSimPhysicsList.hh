#ifndef LYSimPhysicsList_h
#define LYSimPhysicsList_h

#include "globals.hh"
#include "G4VUserPhysicsList.hh"

class G4Cerenkov;
class G4ComptonScattering;
class LYSimScintillation;
class G4OpAbsorption;
class G4OpRayleigh;
class G4OpMieHG;
class G4OpBoundaryProcess;
class G4OpWLS;
class G4VPhysicsConstructor;
class LYSimPhysicsListMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class LYSimPhysicsList : public G4VUserPhysicsList
{
public:
    LYSimPhysicsList();
    ~LYSimPhysicsList();

    void ConstructParticle();
    void ConstructProcess();

    void SetCuts();
    void AddPhysicsList(const G4String& name);
    void SetCutForGamma(G4double);
    void SetCutForElectron(G4double);
    void SetCutForPositron(G4double);
    void SetCutForProton(G4double);

    //these methods Construct particles
    void ConstructBosons();
    void ConstructLeptons();
    void ConstructMesons();
    void ConstructBaryons();
    void ConstructIons();

    //these methods Construct physics processes and register them
    void ConstructGeneral();
    void ConstructEM();
    void ConstructHad();
    void ConstructOp();
    void ConstructIdealOp();

    //for the Messenger 
    void SetVerbose(G4int);
    void SetNbOfPhotonsCerenkov(G4int);

    void SetHadProc(G4bool boolHad) {hadProcess = boolHad;}

private:
    G4Cerenkov*          theCerenkovProcess;
    LYSimScintillation*  theScintillationProcess;
    G4OpAbsorption*      theAbsorptionProcess;
    G4OpRayleigh*        theRayleighScatteringProcess;
    G4OpMieHG*           theMieHGScatteringProcess;
    G4OpBoundaryProcess* theBoundaryProcess;
    G4OpWLS*             theWLSProcess;

    // for EM physics options
    G4String emName;
    G4VPhysicsConstructor*  emPhysicsList;
    G4double cutForGamma;
    G4double cutForElectron;
    G4double cutForPositron;
    G4double cutForProton;

    G4bool hadProcess;

    LYSimPhysicsListMessenger* pMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /* LYSimPhysicsList_h */
