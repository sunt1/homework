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
//
/// \file B1EventAction.cc
/// \brief Implementation of the B1EventAction class

#include <G4SIunits.hh>
#include "B1EventAction.hh"
#include "B1RunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "B1Analysis.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::B1EventAction(B1RunAction *runAction)
        : G4UserEventAction(),
          fRunAction(runAction),
          fEdep(0.) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::~B1EventAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::BeginOfEventAction(const G4Event *) {
    fEdep = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::EndOfEventAction(const G4Event *event) {
    // accumulate statistics in run action
    auto eventID = event->GetEventID();
    fRunAction->AddEdep(fEdep);
    // get analysis manager
    auto analysisManager = G4AnalysisManager::Instance();

    // Print per event
    //
    G4cout << "---> End of event: " << eventID << G4endl;
    for (G4int i = 0; i < fTime.size(); i++) {
        G4ThreeVector pos, mnt_truth;
        G4double t, v, beta, Pz, Pz_truth;
        t = fTime.at(i);
        pos = fPos.at(i);
        mnt_truth = fMnt_truth.at(i);
        v = pos.z() * 1.e-3 / (t * 1.e-9);
        beta = v / (3.e8);
        Pz = 938.272 * beta / sqrt(1 - pow(beta, 2));
        Pz_truth = mnt_truth.getZ();
        // fill histograms
        analysisManager->FillH1(0, Pz);
        analysisManager->FillH1(1, Pz_truth);
//        G4cout << "velocity: " << v << " beta: " << beta << G4endl;
//        G4cout << "Pz:" << Pz << G4endl;
//        G4cout << "truth Pz: " << Pz_truth << G4endl;
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
