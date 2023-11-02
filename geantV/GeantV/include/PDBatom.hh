

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


#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <G4Pow.hh>
#include <vector>
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "Geant/Config.h"
#include <Geant/VectorTypes.h>
#include "Geant/math_wrappers.h"
#include "Geant/TaskData.h"
#include "Geant/Typedefs.h"


using namespace vecCore::math; 
using vecCore::Get;
using vecCore::Set;
using vecCore::Mask_v;
//using vecCore::AssignMaskLane;
using vecCore::MaskEmpty;
using vecCore::MaskFull;

using Double_v = geant::Double_v;
using Int_v = geant::Double_v;

struct  EML {
G4int Nomega;
std::vector<std::vector<G4double>> leb;
std::vector<std::vector<Double_v>> lebV;
};


class Atom{
public:

  Atom(int serial,
       const std::string& name,
       const std::string& resName,
       int numInRes,
       int resSeq,
       double xInit,
       double yInit,
       double zInit,
       double occupancy,
       double tempFactor,
       const std::string& Chain,
       const std::string& e);

  ~Atom(){};

  Atom *GetNext();

  double GetX();
  double GetY();
  double GetZ();

  int GetID();

  const std::string& GetName();
  const std::string& GetElementName();

  double GetVanDerWaalsRadius();

  void SetNext(Atom *);
  G4double atom_distance(G4double x, G4double y, G4double z);
  void Set_data();
  void Scale_units();
  G4int Get_atomic_number();



  void Gaussian_sto3g();

  int fSerial;       
  int fNumInRes;   
  std::string fName;   
  std::string fResName;  
  int fResSeq;       
  G4double fX;      
  G4double fY;  
  G4double fZ;   
  double fVdwRadius; 
  double fOccupancy; 
  std::string fElement;
  G4int fZat;
  double fTempFactor; 
 
  G4int fId;
  G4String fChainame;
  G4int sizeZ;
  G4int sizeL;
  G4double Bohr_radius;
  G4double Bragg_radius;

  struct Norm_alpha_coff {
  G4double C;
  G4double N;
  G4double A;
  };

  struct orb {
  std::vector<G4int> l;
  Norm_alpha_coff NAC[3];
  Double_v orbitalsAlpha;
  Double_v orbitalsNormals;
  Double_v orbitalsConst;

  std::vector<Double_v> EvalV;
  std::vector<Double_v> X_G;
  std::vector<Double_v> Y_G;
  std::vector<Double_v> Z_G;

  std::vector<Double_v> DXX;
  std::vector<Double_v> DYY;
  std::vector<Double_v> DZZ;

  std::vector<Double_v> DXY;
  std::vector<Double_v> DXZ;
  std::vector<Double_v> DYX;
  std::vector<Double_v> DYZ;
  std::vector<Double_v> DZX;
  std::vector<Double_v> DZY;
  };

 std::vector<orb> orbitals;
 std::vector<EML> atomic_grid;

 std::vector<orb>::iterator shell1;
 std::vector<orb>::iterator shell2;
 std::vector<orb>::iterator shell3;
 std::vector<orb>::iterator shell4;

 std::vector<EML>::iterator grid1;


private:
Atom * fpNext;      
};
#endif

