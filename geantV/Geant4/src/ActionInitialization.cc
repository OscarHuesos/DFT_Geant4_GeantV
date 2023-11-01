

// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// Delage et al. PDB4DNA: implementation of DNA geometry from the Protein Data
//                  Bank (PDB) description for Geant4-DNA Monte-Carlo
//                  simulations (submitted to Comput. Phys. Commun.)
// The Geant4-DNA web site is available at http://geant4-dna.org

#include "ActionInitialization.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"

#include "G4EmDNAPhysics.hh"
#include "G4SystemOfUnits.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(){
  G4int n_particle = 1;
  fpParticleGun  = new G4ParticleGun(n_particle);

  G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle("e-");
  fpParticleGun->SetParticleDefinition(particle);
  fpParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fpParticleGun->SetParticleEnergy(0.1*MeV);
  fpParticleGun->SetParticlePosition(G4ThreeVector(0.*nm,0.*nm,0.*nm));
}



PrimaryGeneratorAction::~PrimaryGeneratorAction(){
  delete fpParticleGun;
}


void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent){

  G4double boundXHalfLength = 0;
  G4double boundYHalfLength = 0;
  G4double boundZHalfLength = 0;

  G4VPhysicalVolume* boundPV
  = G4PhysicalVolumeStore::GetInstance()->GetVolume("boundingPV");

  G4ThreeVector boundPos;
  if ( boundPV ){
    boundPos = boundPV->GetTranslation();
  }

  G4LogicalVolume* boundLV
  = G4LogicalVolumeStore::GetInstance()->GetVolume("BoundingLV");

  G4Box* boundBox = 0;
  if ( boundLV ){
    boundBox = dynamic_cast< G4Box*>(boundLV->GetSolid());
  }

  if ( boundBox ){
    boundXHalfLength = boundBox->GetXHalfLength();
    boundYHalfLength = boundBox->GetYHalfLength();
    boundZHalfLength = boundBox->GetZHalfLength();


    G4double radius = std::sqrt(boundXHalfLength*boundXHalfLength+
    boundYHalfLength*boundYHalfLength+
    boundZHalfLength*boundZHalfLength);
    G4double cosTheta = 2*G4UniformRand()-1;
    G4double sinTheta = std::sqrt(1.-cosTheta*cosTheta);
    G4double phi = twopi*G4UniformRand();

    G4ThreeVector positionStart(boundPos.x()+radius*sinTheta*std::cos(phi),
    boundPos.y()+radius*sinTheta*std::sin(phi),
    boundPos.z()+radius*cosTheta);

    fpParticleGun->SetParticlePosition(positionStart);

    G4ThreeVector positionDir(
        boundPos.x()+boundXHalfLength*(2*G4UniformRand()-1),
        boundPos.y()+boundYHalfLength*(2*G4UniformRand()-1),
        boundPos.z()+boundZHalfLength*(2*G4UniformRand()-1));

    fpParticleGun->SetParticleMomentumDirection(
        (positionDir-positionStart).unit());
  }  else{
    G4cerr << "Bounding volume not found." << G4endl;
    G4cerr << "Default particle kinematic used" << G4endl;
  }

  fpParticleGun->GeneratePrimaryVertex(anEvent);
}



ActionInitialization::ActionInitialization(): G4VUserActionInitialization(){}

ActionInitialization::~ActionInitialization(){}

void ActionInitialization::BuildForMaster() const{
SetUserAction(new RunAction);
}

void ActionInitialization::Build() const{
  SetUserAction(new PrimaryGeneratorAction);
  SetUserAction(new RunAction);
  SetUserAction(new EventAction());
  SetUserAction(new SteppingAction());
}

PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
  SetDefaultCutValue(1.0*micrometer);
  SetVerboseLevel(1);
  RegisterPhysics(new G4EmDNAPhysics());
}


PhysicsList::~PhysicsList(){
}
