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

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4Material.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "PDBbarycenter.hh"
#include "dft.hh"

class DetectorMessenger;
class G4LogicalVolume;
class G4VPhysicalVolume;

class DetectorConstruction : public G4VUserDetectorConstruction{
public:
  DetectorConstruction();
  virtual ~DetectorConstruction();

  virtual G4VPhysicalVolume* Construct();

private:

  G4String fPdbFileName;

  unsigned short int fChosenOption;
  unsigned short int fPdbFileStatus;

  PDBlib fPDBlib;
  Molecule *fpMoleculeList;
  Barycenter *fpBarycenterList;
  G4Material *fpDefaultMaterial;
  G4Material *fpWaterMaterial;
  Molecule Mol;

  void   ConstructMaterials();
  void   CheckMaterials();
  G4VPhysicalVolume* ConstructWorld();

  G4VPhysicalVolume* DefineVolumes(Molecule&
    fpMoleculeList, unsigned short int option, int& typ);
  Molecule Charge_molecule(G4String filename, int& isM);
  void AtomisticView(G4LogicalVolume*,Molecule , double atomSizeFactor);
  void BarycenterView(G4LogicalVolume* ,Barycenter *);
  void ResiduesView(G4LogicalVolume* ,Barycenter *);
  void DrawBoundingVolume(G4LogicalVolume* ,Molecule );

  DetectorMessenger* fpMessenger; 
  G4bool  fCheckOverlaps; 

public:

  Barycenter *GetBarycenterList();
  PDBlib GetPDBlib();
  Molecule *GetMoleculeList();

  void LoadPDBfile(G4String fileName);
  void DrawAtoms_();
  void DrawNucleotides_();
  void DrawResidues_();
  void BuildBoundingVolume();
  void DrawAtomsWithBounding_();
  void DrawNucleotidesWithBounding_();
  void DrawResiduesWithBounding_();
  G4double Molecule_energy_Analysis(Molecule&  fpMol );


};

#endif
