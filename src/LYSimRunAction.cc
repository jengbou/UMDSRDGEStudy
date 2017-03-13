#include "Analysis.hh"

#include "LYSimRunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"

#include "LYSimDetectorConstruction.hh"
#include <iostream>
#include "g4root.hh"

LYSimRunAction::LYSimRunAction(LYSimDetectorConstruction* ipDetectorConstruction, std::string outfile)
{
//     runMessenger = new LYSimRunActionMessenger(this);
    pDetectorConstruction = ipDetectorConstruction;
    G4AnalysisManager::Instance();
    outFileName = outfile;
}

LYSimRunAction::~LYSimRunAction()
{
//     delete runMessenger;
}

void LYSimRunAction::BeginOfRunAction(const G4Run* aRun)
{
    Analysis::GetInstance()->PrepareNewRun(aRun);
    G4AnalysisManager* man = G4AnalysisManager::Instance();
    G4cout << "Output filename: " << outFileName << G4endl;
    man->OpenFile(outFileName.c_str());
    man->SetFirstHistoId(1);

    // Create histogram(s)
    man->CreateH1("h1","Optical photons energy (eV)", //histoID,histo name 
                  500,0.,5.); //bins' number, xmin, xmax
    man->CreateH1("h2","Number of Detected Photons",
                  40,0.,40.); //bins' number, xmin, xmax
    man->CreateH1("h3","Total optical photons energy (eV)",
                  100,0.,5.); //bins' number, xmin, xmax

}

void LYSimRunAction::EndOfRunAction(const G4Run* aRun)
{
    Analysis::GetInstance()->EndOfRun(aRun);
    // Save histograms 
    G4AnalysisManager* man = G4AnalysisManager::Instance();
    man->Write();
    man->CloseFile();

}
