
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// Delage et al. PDB4DNA: implementation of DNA geometry from the Protein Data
//                  Bank (PDB) description for Geant4-DNA Monte-Carlo
//                  simulations (submitted to Comput. Phys. Commun.)
// The Geant4-DNA web site is available at http://geant4-dna.org

#include "globals.hh"

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction();
  ~PrimaryGeneratorAction();

public:
  virtual void GeneratePrimaries(G4Event*);

private:
  G4ParticleGun* fpParticleGun;
};

#endif


#ifndef ActionInitialization_h
#define ActionInitialization_h 1

#include "DetectorConstruction.hh"
#include "G4VUserActionInitialization.hh"

class ActionInitialization : public G4VUserActionInitialization
{
public:
  ActionInitialization();
  virtual ~ActionInitialization();

  virtual void BuildForMaster() const;
  virtual void Build() const;
};

#endif

#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"

class G4VPhysicsConstructor;

class PhysicsList: public G4VModularPhysicsList
{
public:
  PhysicsList();
  virtual ~PhysicsList();
};

#endif


