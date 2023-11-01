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


#ifndef EventAction_h
#define EventAction_h 1

#include "globals.hh"
#include "G4UserEventAction.hh"
#include "G4UImessenger.hh"
#include <map>

class EventActionMessenger;

class EventAction : public G4UserEventAction
{
public:
  EventAction();
  ~EventAction();

public:
  virtual void BeginOfEventAction(const G4Event*);
  virtual void EndOfEventAction(const G4Event*);

  void AddEdepEvent(G4double edep)
  {
    fTotalEnergyDeposit += edep;
  };
  G4double GetEdepEvent()
  {
    return fTotalEnergyDeposit;
  };

  void AddEdepToNucleotide(G4int numStrand,G4int numNucl,G4double edep)
  {
    if(numStrand==1)
    {
      fEdepStrand1[numNucl]+=edep;
    }
    else{
      fEdepStrand2[numNucl]+=edep;
    }
  }

  void SetEnergyThresForSSB(G4double val)
  {
    fThresEdepForSSB=val;
  };
  void SetDistanceThresForDSB(G4int val)
  {fThresDistForDSB=val;
  };

private:

  G4double fTotalEnergyDeposit;

  std::map<G4int,G4double>  fEdepStrand1;
  std::map<G4int,G4double>  fEdepStrand2;

  G4double fThresEdepForSSB;
  G4int fThresDistForDSB;

  EventActionMessenger* fpEventMessenger;

  void ComputeStrandBreaks(G4int*);
};

#endif

#ifndef EventActionMessenger_h
#define EventActionMessenger_h 1


class EventAction;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIdirectory;

class EventActionMessenger: public G4UImessenger
{
public:
  EventActionMessenger(EventAction*);
  ~EventActionMessenger();

  virtual void SetNewValue(G4UIcommand*, G4String);

private:
  EventAction*                  fpEventAction;
  G4UIdirectory*                fpPDBDir;
  G4UIcmdWithADoubleAndUnit*    fpThresEdepCmd;
  G4UIcmdWithAnInteger*         fpThresDistCmd;
};

#endif



#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "RunAction.hh"

class DetectorConstruction;

class SteppingAction : public G4UserSteppingAction, public RunInitObserver
{
public:
  SteppingAction();
  ~SteppingAction();

  virtual void UserSteppingAction(const G4Step*);
  virtual void Initialize();

private:
  G4bool CheckAndProcessDNAHit(
      G4double x,G4double y, G4double z,
      G4double edep);
  EventAction* fpEventAction;
  DetectorConstruction* fpDetector;
};

#endif

