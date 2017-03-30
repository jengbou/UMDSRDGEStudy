#include "Analysis.hh"
#include "AnalysisMessenger.hh"
#include "G4UnitsTable.hh"
#include "G4SDManager.hh"
#include "G4PhysicalConstants.hh"
#include "LYSimDetectorConstruction.hh"
#include "G4ios.hh"
#include "g4root.hh"

using namespace std;
using namespace CLHEP;

//ROOT Stuff
//#include "TProfile.h"
//#include "TFile.h"

Analysis* Analysis::singleton = 0;

//Constructor
Analysis::Analysis()
{
    fMessenger = new AnalysisMessenger(this);
    // //Delete previous contents of output file.
    // outputfile.open(fOutputFileName.c_str(), ofstream::out | ofstream::trunc);
    // outputfile.close();
    //Instantiate the analysis manager
    G4AnalysisManager::Instance();

}

Analysis::~Analysis()
{
    if(fMessenger) delete fMessenger;
}


void Analysis::PrepareNewEvent(const G4Event* /*anEvent*/)
{
}

void Analysis::EndOfEvent(const G4Event* anEvent)
{
    G4AnalysisManager* man = G4AnalysisManager::Instance();

    G4String hitCollName = "PMTHitsCollection";
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    static G4int hitCollID = -1;
    if ( hitCollID < 0 )
        hitCollID = SDman->GetCollectionID( hitCollName );

    G4HCofThisEvent* hitsCollections = 0;
    hitsCollections = anEvent->GetHCofThisEvent();

    LYSimPMTHitsCollection* hits = 0;
    if ( hitsCollections )
    {
        hits = static_cast<LYSimPMTHitsCollection*> ( hitsCollections->GetHC(hitCollID) );
    }
    else
    {
        G4cerr << "hitsCollection not found" << G4endl;
        return;
    }

    G4double EventEnergy = 0;
    G4int EventPhotonCount = 0;

    G4double nHits = hits->entries();
    G4double firstHitTime = 0;
    for (G4int i=0; i<nHits; i++)
    {
        G4double HitEnergy = (*hits)[i]->GetEnergy();
        G4double HitTime = (*hits)[i]->GetTime();
        if (i==0) firstHitTime = HitTime;

        if (HitTime - firstHitTime > 200.0)
            G4cout << "[LYSim] Late signal: delta Time [ns] = " << HitTime - firstHitTime << G4endl;

        EventEnergy += HitEnergy;
        EventPhotonCount = (*hits)[i]->GetPhotonCount();
        PhotonCount += EventPhotonCount;
        HitCount++;
        man->FillH1(1,HitEnergy/eV);

        if (anEvent->GetEventID()%100 == 0 && nHits > 0) {
            if (i==0) G4cout << "[LYSim] Number of hits in PMT: " << nHits << G4endl;
            G4cout << "[LYSim] hit[" << i << "] in PMT energy = "
                   << std::setprecision(4)
                   << std::fixed
                   << std::setw(6) << HitEnergy/eV << " [eV]; time = "
                   << std::setw(6) << HitTime/ns << " [ns]"
                   << G4endl << std::resetiosflags(std::ios::fixed);;
        }
    }
    man->FillH1(2,HitCount);
    man->FillH1(3,EventEnergy/eV);
}


void Analysis::PrepareNewRun(const G4Run* /*aRun*/ )
{
    //Reset variables relative to the run
    PhotonCount = 0;
    HitCount = 0;
}


void Analysis::EndOfRun(const G4Run* aRun)
{
    outputfile.open(fOutputFileName.c_str(), ofstream::out | ofstream::app);
    G4double detEff = (PhotonCount > 0 ? (G4double)HitCount/(G4double)PhotonCount: 0.0);
    G4cout << "Efficiency in this run is " << detEff  << G4endl;
    if (outputfile.is_open())
    {
        outputfile << "#Mu_tile [cm^-1]\tMu_fiber [cm^-1]\tEfficiency" << G4endl;
        outputfile << inducedMuTile << "\t" << inducedMuFiber << "\t" << detEff << G4endl;
    }
    else
    {
        G4cout << "Output file not open" << G4endl;
        G4cout << "#Mu_tile [cm^-1]\tMu_fiber [cm^-1]\tEfficiency" << G4endl;
        G4cout << inducedMuTile << "\t" << inducedMuFiber << "\t" << detEff << G4endl;
    }
    outputfile.close();

    // Save histograms 
    G4AnalysisManager* man = G4AnalysisManager::Instance();
    man->Write();
    man->CloseFile();

}
