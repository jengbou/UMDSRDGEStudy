//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file LYSim/src/LYSimPMTSD.cc
/// \brief Implementation of the LYSimPMTSD class
//
//
//
#include "Analysis.hh"

#include "LYSimPMTSD.hh"
#include "LYSimPMTHit.hh"

#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"

LYSimPMTSD::LYSimPMTSD(G4String name)
    : G4VSensitiveDetector(name),fPMTHitsCollection(0)
{
    collectionName.insert("PMTHitsCollection");
}

LYSimPMTSD::~LYSimPMTSD() { }

void LYSimPMTSD::Initialize(G4HCofThisEvent*)
{
    fPMTHitsCollection = new LYSimPMTHitsCollection(GetName(),collectionName[0]);
}

G4bool LYSimPMTSD::ProcessHits(G4Step* aStep, G4TouchableHistory* )
{
    return false; //ProcessHits should not be called
}

G4bool LYSimPMTSD::ProcessHits_constStep(const G4Step* aStep, G4TouchableHistory* ROhist)
{
    G4Track* theTrack = aStep->GetTrack();

    //need to know if this is an optical photon
    if(theTrack->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) return false;

    LYSimPMTHit* hit=new LYSimPMTHit();

    hit->AddEnergy(theTrack->GetTotalEnergy());
    hit->IncPhotonCount(); //increment hit for the selected pmt
    hit->SetTime(aStep->GetPostStepPoint()->GetGlobalTime());

    fPMTHitsCollection->insert(hit);

    return true;

}

void LYSimPMTSD::EndOfEvent(G4HCofThisEvent* HCE){
    //Store collection with event and keep ID
    static G4int HCID = -1;
    if (HCID<0) HCID = GetCollectionID(0);
    HCE->AddHitsCollection( HCID, fPMTHitsCollection );

    G4int nHits = fPMTHitsCollection->entries();
    if (verboseLevel>=1) {
        G4cout << "[LYSim] PMT collection: " << nHits << " hits" << G4endl;
        if (verboseLevel>=2)
            fPMTHitsCollection->PrintAllHits();
    }
    fPMTHitsCollection->PrintAllHits();
}

void LYSimPMTSD::clear(){ }

void LYSimPMTSD::DrawAll(){ }

void LYSimPMTSD::PrintAll(){ }
