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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// Delage et al. PDB4DNA: implementation of DNA geometry from the Protein Data
//                  Bank (PDB) description for Geant4-DNA Monte-Carlo
//                  simulations (submitted to Comput. Phys. Commun.)
// The Geant4-DNA web site is available at http://geant4-dna.org

#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "G4EventManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4SteppingManager.hh"
#include "G4SystemOfUnits.hh"


#include "Analysis.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <algorithm>
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIdirectory.hh"



EventAction::EventAction():G4UserEventAction(){

  fThresEdepForSSB=8.22*eV;
  fThresDistForDSB=10;
  fTotalEnergyDeposit=0;
  fpEventMessenger = new EventActionMessenger(this);
}

EventAction::~EventAction(){
  delete fpEventMessenger;
}

void EventAction::BeginOfEventAction( const G4Event*){

  fTotalEnergyDeposit=0.;
  fEdepStrand1.clear();
  fEdepStrand2.clear();
}


void EventAction::EndOfEventAction( const G4Event*){

  G4int sb[2] = {0,0};
  ComputeStrandBreaks(sb);
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  if ( fTotalEnergyDeposit>0. ){

    analysisManager->FillH1(1,fTotalEnergyDeposit);
    analysisManager->FillH2(1,fTotalEnergyDeposit,1);

  }
  if ( sb[0]>0 ){
    analysisManager->FillH1(2,sb[0]);
  }
  if ( sb[1]>0 ){
    analysisManager->FillH1(3,sb[1]);
  }
}

void EventAction::ComputeStrandBreaks(G4int* sb){

  G4int ssb1=0;
  G4int ssb2=0;
  G4int dsb=0;

  G4int nucl1;
  G4int nucl2;
  G4double edep1;
  G4double edep2;

  while ( !fEdepStrand1.empty() )
  {
    nucl1 = fEdepStrand1.begin()->first;
    edep1 = fEdepStrand1.begin()->second;
    fEdepStrand1.erase( fEdepStrand1.begin() );

    if ( edep1 >= fThresEdepForSSB/eV )
    {
      ssb1++;
    }


    if ( !fEdepStrand2.empty() )
    {
      do
      {
        nucl2 = fEdepStrand2.begin()->first;
        edep2 = fEdepStrand2.begin()->second;
        if ( edep2 >= fThresEdepForSSB/eV )
        {
          ssb2++;
        }
        fEdepStrand2.erase( fEdepStrand2.begin() );
      } while ( ((nucl1-nucl2)>fThresDistForDSB) && (!fEdepStrand2.empty()) );


      if ( nucl2-nucl1 > fThresDistForDSB )
      {
        fEdepStrand2[nucl2]=edep2;
        if ( edep2 >= fThresEdepForSSB/eV )
        {
          ssb2--;
        }
      }

      if ( std::abs(nucl2-nucl1) <= fThresDistForDSB )
      {
        if ( ( edep2 >= fThresEdepForSSB/eV ) &&
            ( edep1 >= fThresEdepForSSB/eV ) )
        {
          ssb1--;
          ssb2--;
          dsb++;
        }
      }
    }
  }

  while ( !fEdepStrand1.empty() )
  {
    nucl1 = fEdepStrand1.begin()->first;
    edep1 = fEdepStrand1.begin()->second;
    if ( edep1 >= fThresEdepForSSB/eV )
    {
      ssb1++;
    }
    fEdepStrand1.erase( fEdepStrand1.begin() );
  }

  while ( !fEdepStrand2.empty() )
  {
    nucl2 = fEdepStrand2.begin()->first;
    edep2 = fEdepStrand2.begin()->second;
    if ( edep2 >= fThresEdepForSSB/eV )
    {
      ssb2++;
    }
    fEdepStrand2.erase( fEdepStrand2.begin() );
  }

  sb[0]=ssb1+ssb2;
  sb[1]=dsb;
}



EventActionMessenger::EventActionMessenger(EventAction* EvAct)
:G4UImessenger(),fpEventAction(EvAct)
{
  fpPDBDir = new G4UIdirectory("/PDB4DNA/");
  fpPDBDir->SetGuidance("commands specific to this example");

  fpThresEdepCmd = new G4UIcmdWithADoubleAndUnit(
      "/PDB4DNA/event/setEnergyThres",
      this);
  fpThresEdepCmd->SetGuidance("Set energy threshold for SSB");
  fpThresEdepCmd->SetParameterName("EnergyThres",false);
  fpThresEdepCmd->SetRange("EnergyThres>0");
  fpThresEdepCmd->AvailableForStates(G4State_Idle);

  fpThresDistCmd = new G4UIcmdWithAnInteger("/PDB4DNA/event/setDistanceThres",
                                            this);
  fpThresDistCmd->SetGuidance("Set distance threshold for DSB");
  fpThresDistCmd->SetParameterName("DistanceThres",false);
  fpThresDistCmd->SetRange("DistanceThres>0");
  fpThresDistCmd->AvailableForStates(G4State_Idle);
}



EventActionMessenger::~EventActionMessenger()
{
  delete fpThresEdepCmd;
  delete fpThresDistCmd;
  delete fpPDBDir;
}


void EventActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if(command==fpThresEdepCmd)
  {
    fpEventAction->SetEnergyThresForSSB(fpThresEdepCmd->GetNewDoubleValue(
        newValue));
  }

  if(command==fpThresDistCmd)
  {
    fpEventAction->SetDistanceThresForDSB(fpThresDistCmd->GetNewIntValue(
        newValue));
  }

}


SteppingAction::SteppingAction()
:G4UserSteppingAction(),RunInitObserver(),fpEventAction(0),fpDetector(0){
}


SteppingAction::~SteppingAction(){
}


void SteppingAction::Initialize(){
  fpEventAction = (EventAction*) G4EventManager::GetEventManager()->
      GetUserEventAction();
      printf("stepping action \n");
  fpDetector = (DetectorConstruction*)G4RunManager::GetRunManager()->
      GetUserDetectorConstruction();
}


void SteppingAction::UserSteppingAction(const G4Step* theStep){
  if(theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()!=
      "Transportation")
  {
 
    G4double x = theStep->GetPreStepPoint()->GetPosition().x()/nanometer;
    G4double y = theStep->GetPreStepPoint()->GetPosition().y()/nanometer;
    G4double z = theStep->GetPreStepPoint()->GetPosition().z()/nanometer;
    G4double edepStep = theStep->GetTotalEnergyDeposit()/eV;

    G4LogicalVolume* targetVolume =
        G4LogicalVolumeStore::GetInstance()->GetVolume("BoundingLV");
    G4LogicalVolume* theVolume =
        theStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume();

    if ((edepStep > 0.) && (theVolume==targetVolume))    {

      fpEventAction->AddEdepEvent(edepStep);
      if (fpDetector->GetBarycenterList()==NULL){
        G4cout << "Barycenter list is null!!!" << G4endl;
      }
      else
      {
        CheckAndProcessDNAHit(x,y,z,edepStep);
      }
    }
  }
}


G4bool SteppingAction::CheckAndProcessDNAHit(G4double x,G4double y, G4double z,
    G4double edepStep)
{
  int numStrand=0;
  int numNucl=0;
  printf("Checking and processing Step action \n");
  int intResidue=-1; // 0 for Phospat, 1 for Sugar, 2 for Base
  unsigned short int hit = (fpDetector->GetPDBlib()).ComputeMatchEdepDNA(
      fpDetector->GetBarycenterList(),
      fpDetector->GetMoleculeList(),
      x*10., y*10., z*10.,// x10 => angstrom<->nm
      numStrand, numNucl, intResidue);

  if (hit==1)
  {
    if ((intResidue==0)||(intResidue==1)) //Edep in Phosphate or Sugar
    {
      fpEventAction->AddEdepToNucleotide(numStrand,numNucl,edepStep);
      return true;
    }
    else
    {
      return false;
    }
  }
  else
  {
    return false;
  }
}


