
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


struct  EML {
G4int Nomega;
std::vector<std::vector<G4double>> leb;
};


class Atom{
public:

  Atom(G4int id,
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
  G4int fNumInRes;    
  std::string fName;      
  std::string fResName;  
  int fResSeq;       
  G4double fX, fY, fZ, fVdwRadius, fOccupancy; 
  std::string fElement;  
  G4int fZat;
  double fTempFactor; 

  G4int fId;
  G4String fChainame;
  G4int sizeZ;
  //G4int sizeBasis; 
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

  std::vector<G4double> Func;
  std::vector<G4double> x_grad;
  std::vector<G4double> y_grad;
  std::vector<G4double> z_grad;
  std::vector<G4double> dxx;
  std::vector<G4double> dyy;
  std::vector<G4double> dzz;
  std::vector<G4double> dxy;
  std::vector<G4double> dxz;
  std::vector<G4double> dyx;
  std::vector<G4double> dyz;
  std::vector<G4double> dzx;
  std::vector<G4double> dzy;
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

