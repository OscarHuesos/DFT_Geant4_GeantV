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
// J. Comput. Phys. 274 (2014) 841-882
// The Geant4-DNA web site is available at http://geant4-dna.org


#include "PDBatom.hh"
#include <new>
#include <strstream>
#include <fstream>
#include <vector>


#ifndef RESIDUE_H
#define RESIDUE_H


class Residue
{
public:

  Residue();
  Residue(const std::string& resName,int resSeq);
  ~Residue() {};

  Residue *GetNext();
  Atom *GetFirst();

  int GetID();
  void SetNext(Residue *);
  void SetFirst(Atom *);

  std::string fResName;       
  int fResSeq;           
  bool fVisible;        
  bool fSelected;     
  int fCenterX;
  int fCenterY;
  int fCenterZ;
  int fChainN;
  G4String fChainame;
  int fNbAtom;       
  std::vector<Atom *> Lista_atoms;


private:
  Residue *fpNext;
  Atom *fpFirst;
};

#endif


#ifndef MOLECULE_H
#define MOLECULE_H


class Molecule{
public:

  Molecule();

  Molecule(const std::string& resName,int mNum);

  ~Molecule() {};


  Molecule *GetNext();
  Residue *GetFirst();

  int GetID();

  void SetNext(Molecule *);
  void SetFirst(Residue *);

  void isbool(bool t);
  void Setvorbits(G4int orb);
  void check_parity();

  bool retstatus();

void Load( const std::string& filename, int& isProtein); 

  std::string fMolName;
  int fMolNum; 
  bool fExiste;
  double fMinGlobZ;   
  double fMaxGlobZ;
  double fMinGlobX;   
  double fMaxGlobX;
  double fMinGlobY;   
  double fMaxGlobY;

  int fCenterX;      
  int fCenterY;     
  int fCenterZ; 
  int fDistCenterMax;
  int fNbResidue;      
  int No_atomos_en_prot;
  int Nresidue=0;
  int cadena;

  G4int no_orbits=0;
  G4int no_electrons=0;
  G4int ocupados=0;
  G4int abiertos=0;
  G4int No_points=0;

std::vector<int>  N_atoms_per_chain;
std::vector<int>  N_residues_per_chain;
std::vector<Residue *> Residuos_cadena;

private:
  Molecule *fpNext;
  Residue *fpFirst;
};
#endif

