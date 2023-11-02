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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// Delage et al. PDB4DNA: implementation of DNA geometry from the Protein Data
//                  Bank (PDB) description for Geant4-DNA Monte-Carlo
//                  simulations (submitted to Comput. Phys. Commun.)
// The Geant4-DNA web site is available at http://geant4-dna.org


#ifndef BARY_H
#define BARY_H

class Barycenter{
public:

  Barycenter();

  Barycenter(int bNum,double x,double y, double z,
      double Bx,double By, double Bz,
      double Sx,double Sy, double Sz, 
      double Px,double Py, double Pz);

  ~Barycenter() {};

  Barycenter *GetNext();
  int GetID();

  void SetNext(Barycenter *);

  void SetDistance(int i, double);

  double GetDistance(int i);

  void SetRadius(double );

  double GetRadius();
  int fBaryNum;
  double fDistanceTab[33];
  double fRadius;
  double fCenterX;           
  double fCenterY;            
  double fCenterZ;           

  double fCenterBaseX;        
  double fCenterBaseY;      
  double fCenterBaseZ;       

  double fCenterSugarX;  
  double fCenterSugarY;   
  double fCenterSugarZ;  

  double fCenterPhosphateX;   
  double fCenterPhosphateY; 
  double fCenterPhosphateZ;  

private:
Barycenter *fpNext;   
};
#endif

#ifndef pdblibread_h
#define pdblibread 1

#include "PDBmolecule.hh"

class PDBlib{
public:

PDBlib();

~PDBlib() {};

Barycenter* ComputeNucleotideBarycenters(Molecule moleculeListTemp);
Barycenter* ComputeProteinBarycenters(Molecule moleculeListTemp);

  void ComputeBoundingVolumeParams(Molecule moleculeListTemp,
      double &dX,double &dY,double &dZ, double &tX,double &tY,double &tZ);     


void ComputeNbNucleotidsPerStrand(Molecule  moleculeListTemp);
void ComputeProteinPerStrand(Molecule  moleculeListTemp);

  unsigned short int ComputeMatchEdepDNA(Barycenter *,Molecule *,
      double x, double y,double z, int &numStrand, int &numNucleotid, int &codeResidue);

private:

  double DistanceTwo3Dpoints(double xA,double xB, double yA,double yB,
      double zA,double zB);

  int fNbNucleotidsPerStrand;
};

#endif
