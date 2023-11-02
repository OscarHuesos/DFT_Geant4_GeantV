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


#include "PDBbarycenter.hh"

#define GEANT4

#ifdef GEANT4

#include "globals.hh"
#else
#define G4cout std::cout
#define G4cerr std::cerr
#define G4endl std::endl
#define G4String std::string 
#include <cfloat>
#endif
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>
#include <vector>
#include <sstream>
#include <stdlib.h>
using namespace std;

PDBlib::PDBlib() :  fNbNucleotidsPerStrand(0){
}

Barycenter * PDBlib::ComputeNucleotideBarycenters(Molecule  moleculeListTemp){

  Barycenter * BarycenterFirst = NULL;
  Barycenter * BarycenterOld = NULL;
  Barycenter * BarycenterNext = NULL;

  Residue *residueListTemp;
  Atom *AtomTemp;
  int k = 0;
  int j_old = 0;
  printf("antes while ComputeNucleotideBarycenters \n");
  while(moleculeListTemp.retstatus())  {
    residueListTemp = moleculeListTemp.GetFirst();
    k++;
    int j = 0;

    int correctNumerotationNumber = 0;
    if(k == 2 && residueListTemp->fResSeq > 1){
      correctNumerotationNumber = residueListTemp->fResSeq;
    }

    while(residueListTemp){
      AtomTemp = residueListTemp->GetFirst();
      j++;

      if(correctNumerotationNumber != 0)  {
        residueListTemp->fResSeq = residueListTemp->fResSeq
            - correctNumerotationNumber+ 1;
      }

      double baryX = 0., baryY = 0., baryZ = 0.;
      double baryBaseX = 0., baryBaseY = 0., baryBaseZ = 0.;
      double barySugX = 0., barySugY = 0., barySugZ = 0.;
      double baryPhosX = 0., baryPhosY = 0., baryPhosZ = 0.;
      unsigned short int nbAtomInBase = 0;

      for(int i = 0; i < residueListTemp->fNbAtom; i++){
 
        baryX += AtomTemp->fX;
        baryY += AtomTemp->fY;
        baryZ += AtomTemp->fZ;
   
        if(residueListTemp->fResSeq == 1){
          if(i == 0){
            baryPhosX += AtomTemp->fX;
            baryPhosY += AtomTemp->fY;
            baryPhosZ += AtomTemp->fZ;
          }  else if(i < 8){
            barySugX += AtomTemp->fX;
            barySugY += AtomTemp->fY;
            barySugZ += AtomTemp->fZ;
          }else{

            if(AtomTemp->fElement != "H ")  {
              baryBaseX += AtomTemp->fX;
              baryBaseY += AtomTemp->fY;
              baryBaseZ += AtomTemp->fZ;
              nbAtomInBase++;
            }
          }
        }
        else{
          if(i < 4){
            baryPhosX += AtomTemp->fX;
            baryPhosY += AtomTemp->fY;
            baryPhosZ += AtomTemp->fZ;
          }  else if(i < 11)  {
            barySugX += AtomTemp->fX;
            barySugY += AtomTemp->fY;
            barySugZ += AtomTemp->fZ;
          }  else{
            
            if(AtomTemp->fElement != "H "){ 
              baryBaseX += AtomTemp->fX;
              baryBaseY += AtomTemp->fY;
              baryBaseZ += AtomTemp->fZ;
              nbAtomInBase++;
            }
          }
        }
        AtomTemp = AtomTemp->GetNext();
      } 

      baryX = baryX / (double) residueListTemp->fNbAtom;
      baryY = baryY / (double) residueListTemp->fNbAtom;
      baryZ = baryZ / (double) residueListTemp->fNbAtom;

      if(residueListTemp->fResSeq != 1)   {
        baryPhosX = baryPhosX / 4.;
        baryPhosY = baryPhosY / 4.;
        baryPhosZ = baryPhosZ / 4.;
      }
      barySugX = barySugX / 7.;
      barySugY = barySugY / 7.;
      barySugZ = barySugZ / 7.;
      baryBaseX = baryBaseX / (double) nbAtomInBase;
      baryBaseY = baryBaseY / (double) nbAtomInBase;
      baryBaseZ = baryBaseZ / (double) nbAtomInBase;

      if(BarycenterOld == NULL){
        BarycenterOld = new Barycenter(j + j_old, baryX, baryY, baryZ, 
                                       baryBaseX,
                                       baryBaseY,
                                       baryBaseZ,
                                       barySugX,
                                       barySugY,
                                       barySugZ,
                                       baryPhosX,
                                       baryPhosY,
                                       baryPhosZ);
        BarycenterFirst = BarycenterOld;
      }
      else{  BarycenterNext = new Barycenter(j + j_old,baryX,baryY, baryZ,
                                        baryBaseX,
                                        baryBaseY,
                                        baryBaseZ,
                                        barySugX,
                                        barySugY,
                                        barySugZ,
                                        baryPhosX,
                                        baryPhosY,
                                        baryPhosZ);
        BarycenterOld->SetNext(BarycenterNext);
        BarycenterOld = BarycenterNext;
      }

      AtomTemp = residueListTemp->GetFirst();
      double dT3Dp;
      double max = 0.;
      for(int ii = 0; ii < residueListTemp->fNbAtom; ii++)  {
        dT3Dp = DistanceTwo3Dpoints(AtomTemp->fX,
                                    BarycenterOld->fCenterX,
                                    AtomTemp->fY,
                                    BarycenterOld->fCenterY,
                                    AtomTemp->fZ,
                                    BarycenterOld->fCenterZ);
        BarycenterOld->SetDistance(ii, dT3Dp);
        if(dT3Dp > max) max = dT3Dp;
        AtomTemp = AtomTemp->GetNext();
      } 
      BarycenterOld->SetRadius(max + 1.8);
      residueListTemp = residueListTemp->GetNext();
    } 
    j_old += j;
   
  } 

  if(BarycenterNext != NULL){
    BarycenterNext->SetNext(NULL);
  }
  return BarycenterFirst;
}

 Barycenter * PDBlib::ComputeProteinBarycenters(Molecule  moleculeListTemp){

   Barycenter * BarycenterFirst = NULL;
   Barycenter * BarycenterOld = NULL;
   Barycenter * BarycenterNext = NULL;

   std::vector<Atom >::iterator atm;
   std::vector<Atom *>::iterator atr;
   std::vector<Residue *>::iterator red;

   int k = 0;
   int i=0;
   int j_old = 0;
   k++;
  int j = 0;
  double baryX = 0., baryY = 0., baryZ = 0.;
  double baryBaseX = 0., baryBaseY = 0., baryBaseZ = 0.;
  double barySugX = 0., barySugY = 0., barySugZ = 0.;
  double baryPhosX = 0., baryPhosY = 0., baryPhosZ = 0.;
  unsigned short int nbAtomInBase = 0;

 for(red = moleculeListTemp.Residuos_cadena.begin();
   red != moleculeListTemp.Residuos_cadena.end() ; ++red){
     j++;
      baryX = 0., baryY = 0., baryZ = 0.;
      baryBaseX = 0., baryBaseY = 0., baryBaseZ = 0.;
      barySugX = 0., barySugY = 0., barySugZ = 0.;
    baryPhosX = 0., baryPhosY = 0., baryPhosZ = 0.;
     nbAtomInBase = 0;

     printf("molecule list %d \n",(*red)->fResSeq);

for(atr = (*red)->Lista_atoms.begin();
   atr != (*red)->Lista_atoms.end(); ++ atr){

  baryX += (*atr)->fX;
  baryY += (*atr)->fY;
  baryZ += (*atr)->fZ;

  if((*red)->fResSeq == 1){

    if(i == 0){
 
      baryPhosX += (*atr)->fX;
      baryPhosY += (*atr)->fY;
      baryPhosZ += (*atr)->fZ;
    }  else if(i < 8){
   
      barySugX += (*atr)->fX;
      barySugY += (*atr)->fY;
      barySugZ += (*atr)->fZ;
    }else{

      if((*atr)->fElement != "H ")  {

        baryBaseX = baryBaseX+(*atr)->fX;
        baryBaseY = baryBaseY+(*atr)->fY;
        baryBaseZ = baryBaseZ+(*atr)->fZ;
        nbAtomInBase++;
      }
    }
  } else{
  
    if(i < 4){

      baryPhosX = baryPhosX+(*atr)->fX;
      baryPhosY = baryPhosY+(*atr)->fY;
      baryPhosZ = baryPhosZ + (*atr)->fZ;
 
    }  else if(i < 11)  {
      barySugX = barySugX+(*atr)->fX;
      barySugY = barySugY+(*atr)->fY;
      barySugZ = barySugZ+(*atr)->fZ;
    }  else{

      if((*atr)->fElement != "H "){ 

        baryBaseX = baryBaseX+(*atr)->fX;
        baryBaseY = baryBaseY+(*atr)->fY;
        baryBaseZ = baryBaseZ+(*atr)->fZ;
        nbAtomInBase++;
      }
    }
  }

i++;

   baryX = baryX / (double) (*red)->fNbAtom;
   baryY = baryY / (double) (*red)->fNbAtom;
   baryZ = baryZ / (double) (*red)->fNbAtom;

   if((*red)->fResSeq != 1) 
   {
     baryPhosX = baryPhosX / 4.;
     baryPhosY = baryPhosY / 4.;
     baryPhosZ = baryPhosZ / 4.;
   }

   barySugX = barySugX / 7.;
   barySugY = barySugY / 7.;
   barySugZ = barySugZ / 7.;
   baryBaseX = baryBaseX / (double) nbAtomInBase;
   baryBaseY = baryBaseY / (double) nbAtomInBase;
   baryBaseZ = baryBaseZ / (double) nbAtomInBase;

   if(BarycenterOld == NULL){
     BarycenterOld = new Barycenter(j + j_old, baryX, baryY, baryZ, 
                                    baryBaseX,
                                    baryBaseY,
                                    baryBaseZ,
                                    barySugX,
                                    barySugY,
                                    barySugZ,
                                    baryPhosX,
                                    baryPhosY,
                                    baryPhosZ);
     BarycenterFirst = BarycenterOld;
   } else{  BarycenterNext = new Barycenter(j + j_old,baryX,baryY, baryZ,
                                     baryBaseX,
                                     baryBaseY,
                                     baryBaseZ,
                                     barySugX,
                                     barySugY,
                                     barySugZ,
                                     baryPhosX,
                                     baryPhosY,
                                     baryPhosZ);
     BarycenterOld->SetNext(BarycenterNext);
     BarycenterOld = BarycenterNext;
   }

   double dT3Dp;
   double max = 0.0;
   double xx;
   double yy;
   double zz;

   for(int ii = 0; ii < (*red)->fNbAtom; ii++)  {

     xx=(*atr)->fX;
     yy=(*atr)->fY;
     zz=(*atr)->fZ;

     dT3Dp = DistanceTwo3Dpoints(xx,
                                 BarycenterOld->fCenterX,
                                 yy,
                                 BarycenterOld->fCenterY,
                                 zz,
                                 BarycenterOld->fCenterZ);

     BarycenterOld->SetDistance(ii, dT3Dp);
     if(dT3Dp > max) max = dT3Dp;

   }
   BarycenterOld->SetRadius(max + 1.8);


 }
}
  j_old += j;

  if(BarycenterNext != NULL){
    BarycenterNext->SetNext(NULL);
    }
   return BarycenterFirst;
   }


void PDBlib::ComputeBoundingVolumeParams(Molecule moleculeListTemp,
                                         double &dX,
                                         double &dY,
                                         double &dZ, 
                                         double &tX,
                                         double &tY,
                                         double &tZ) 
{
  double minminX, minminY, minminZ; 
  double maxmaxX, maxmaxY, maxmaxZ; 

  minminX = DBL_MAX;
  minminY = DBL_MAX;
  minminZ = DBL_MAX;
  maxmaxX = -DBL_MAX;
  maxmaxY = -DBL_MAX;
  maxmaxZ = -DBL_MAX;

printf("antes while ComputeBoundingVolumeParams \n");

    if(minminX > moleculeListTemp.fMaxGlobX){
      minminX = moleculeListTemp.fMaxGlobX;
    }
    if(minminY > moleculeListTemp.fMaxGlobY){
      minminY = moleculeListTemp.fMaxGlobY;
    }
    if(minminZ > moleculeListTemp.fMaxGlobZ){
      minminZ = moleculeListTemp.fMaxGlobZ;
    }
    if(maxmaxX < moleculeListTemp.fMinGlobX){
      maxmaxX = moleculeListTemp.fMinGlobX;
    }
    if(maxmaxY < moleculeListTemp.fMinGlobY){
      maxmaxY = moleculeListTemp.fMinGlobY;
    }
    if(maxmaxZ < moleculeListTemp.fMinGlobZ){
      maxmaxZ = moleculeListTemp.fMinGlobZ;}



  dX = (maxmaxX - minminX) / 2. + 1.8; 
  dY = (maxmaxY - minminY) / 2. + 1.8;
  dZ = (maxmaxZ - minminZ) / 2. + 1.8;
  tX = minminX + (maxmaxX - minminX) / 2.;
  tY = minminY + (maxmaxY - minminY) / 2.;
  tZ = minminZ + (maxmaxZ - minminZ) / 2.;
}


void PDBlib::ComputeNbNucleotidsPerStrand(Molecule  moleculeListTemp){
  Residue *residueListTemp;
  int k = 0;
  int j_old = 0;
  printf("partes de nucleotidos \n" );
  while(moleculeListTemp.retstatus())  {
    residueListTemp = moleculeListTemp.GetFirst();
    k++;
    int j = 0;

    while(residueListTemp){
      j++;
      residueListTemp = residueListTemp->GetNext();
    } 
    j_old += j;

  }

fNbNucleotidsPerStrand = j_old / 2;
}


void PDBlib::ComputeProteinPerStrand(Molecule moleculeListTemp){
  Residue *residueListTemp;
  int k = 0;
  int j_old = 0;
  while(moleculeListTemp.retstatus())  {
    residueListTemp = moleculeListTemp.GetFirst();
    k++;
    int j = 0;

    while(residueListTemp){
      j++;
      residueListTemp = residueListTemp->GetNext();
    } 
    j_old += j;

  }

fNbNucleotidsPerStrand = j_old / 2;
}

unsigned short int PDBlib::ComputeMatchEdepDNA(Barycenter *BarycenterList,
                                               Molecule *moleculeListTemp,
                                               double x,double y,double z,
                                               int &numStrand,int &numNucleotid,
                                               int &codeResidue){

  unsigned short int matchFound = 0;
  Molecule *mLTsavedPointer = moleculeListTemp;
  Barycenter *BLsavedPointer = BarycenterList;
  short int strandNum = 0; 
  int residueNum = 1; 
  G4String baseName; 
  unsigned short int BSP = 2; 
  double smallestDist;
  double distEdepDNA;
  double distEdepAtom;
  Residue *residueListTemp;
  Atom *AtomTemp;
  int k = 0; 
  moleculeListTemp = mLTsavedPointer;
  BarycenterList = BLsavedPointer;
  smallestDist = 33.0; 
  while(moleculeListTemp){
    k++;
    residueListTemp = moleculeListTemp->GetFirst();
    int j = 0; 
    int j_save = INT_MAX; 

    while(residueListTemp){
      j++;
      if(j - j_save > 2) break;

      distEdepDNA = DistanceTwo3Dpoints(x,  BarycenterList->fCenterX,
                                        y,BarycenterList->fCenterY,
                                        z,  BarycenterList->fCenterZ);

      if(distEdepDNA < BarycenterList->GetRadius())  {

        AtomTemp = residueListTemp->GetFirst();
        for(int iii = 0; iii < residueListTemp->fNbAtom; iii++){

          distEdepAtom = DistanceTwo3Dpoints(x,AtomTemp->GetX(),
                                             y, AtomTemp->GetY(),
                                             z, AtomTemp->GetZ());

          if((distEdepAtom < AtomTemp->GetVanDerWaalsRadius()) && (smallestDist
              > distEdepAtom))  {
            strandNum = k;
            if(k == 2)  {
              residueNum = fNbNucleotidsPerStrand + 1
                  - residueListTemp->fResSeq;
            }
            else  {
              residueNum = residueListTemp->fResSeq;
            }
            baseName = (residueListTemp->fResName)[2];
            if(residueListTemp->fResSeq == 1)  {
              if(iii == 0) BSP = 0; 
              else if(iii < 8) BSP = 1; 
              else BSP = 2; 
            }
            else{
              if(iii < 4) BSP = 0;
              else if(iii < 11) BSP = 1; 
              else BSP = 2; 
            }

            smallestDist = distEdepAtom;
            int j_max_value = INT_MAX;
            if(j_save == j_max_value) j_save = j;
            matchFound = 1;
          }
          AtomTemp = AtomTemp->GetNext();
        } 
      } 
      BarycenterList = BarycenterList->GetNext();
      residueListTemp = residueListTemp->GetNext();
    } 
    moleculeListTemp = moleculeListTemp->GetNext();
  } 

  numStrand = strandNum;
  numNucleotid = residueNum;
  codeResidue = BSP;
  return matchFound;
}

double PDBlib::DistanceTwo3Dpoints(double xA,double xB, double yA,double yB,
                                   double zA,double zB){

  return sqrt((xA - xB) * (xA - xB) + (yA - yB) * (yA - yB)+ (zA - zB) * (zA - zB));
}


Barycenter::Barycenter():fBaryNum(0),fRadius(0),
fCenterX(0),fCenterY(0),fCenterZ(0),
fCenterBaseX(0),fCenterBaseY(0),fCenterBaseZ(0),
fCenterSugarX(0),fCenterSugarY(0),fCenterSugarZ(0),
fCenterPhosphateX(0),fCenterPhosphateY(0),fCenterPhosphateZ(0),
fpNext(0){
  for ( int i = 0; i < 33; ++i ){
    fDistanceTab[i] = 0.; 
  }
}

Barycenter::Barycenter(int bNum,double x,double y,double z,
    double Bx,double By, double Bz, double Sx,double Sy, double Sz, 
    double Px,double Py, double Pz) {
  fBaryNum=bNum;
  fRadius=0;
  fCenterX=x;
  fCenterY=y;
  fCenterZ=z;
  for ( int i = 0; i < 33; ++i ){
  fDistanceTab[i] = 0.; 
  }
  fCenterBaseX=Bx;           
  fCenterBaseY=By;          
  fCenterBaseZ=Bz;         
  fCenterSugarX=Sx;          
  fCenterSugarY=Sy;          
  fCenterSugarZ=Sz; 
  fCenterPhosphateX=Px;     
  fCenterPhosphateY=Py;   
  fCenterPhosphateZ=Pz; 
  fpNext=0;
}

Barycenter *Barycenter::GetNext(){
return fpNext;
}

int Barycenter::GetID(){
return fBaryNum;}

void Barycenter::SetNext(Barycenter *barycenterNext){
fpNext=barycenterNext;
}

void Barycenter::SetDistance(int i, double dist){
fDistanceTab[i]=dist;
}

double Barycenter::GetDistance(int i){
return fDistanceTab[i];
}


void Barycenter::SetRadius(double rds){
fRadius=rds;
}

double Barycenter::GetRadius(){
return fRadius;
}
