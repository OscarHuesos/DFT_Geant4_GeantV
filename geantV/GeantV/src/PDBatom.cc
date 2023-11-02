
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
#include <fstream>

  std::vector<G4int> s{0,0,0};
  std::vector<G4int> px{1,0,0};
  std::vector<G4int> py{0,1,0};
  std::vector<G4int> pz{0,0,1};
  std::vector<G4int> dxx{2,0,0};
  std::vector<G4int> dxy{1,1,0};
  std::vector<G4int> dxz{1,0,1};
  std::vector<G4int> dyy{0,2,0};
  std::vector<G4int> dyz{0,1,1};
  std::vector<G4int> dzz{0,0,2};
  std::vector<G4int> fxxx{3,0,0};
  std::vector<G4int> fxxy{2,1,0};
  std::vector<G4int> fxyz{1,1,1};

Atom::Atom(int id, const std::string& n,
           const std::string& rN,int numInRes,int rS,
           double xInit, double yInit, double zInit,
           double o, double tF, const std::string& Chain,
           const std::string& e){
  fId=id;
  fName=n;
  fResName=rN;
  fNumInRes=numInRes;
  fResSeq=rS;
  fChainame=Chain;
  fX=xInit;
  fY=yInit;
  fZ=zInit;
  fOccupancy=o;
  fTempFactor=tF;
  fElement=e;
  fpNext=0;
}


Atom *Atom::GetNext(){
  return fpNext;
}

double Atom::GetX(){
  return fX;
}

double Atom::GetY(){
  return fY;
}


double Atom::GetZ(){
  return fZ;
}


int Atom::GetID(){
return fSerial;
}

G4int Atom::Get_atomic_number(){
return  fZat;
}


const std::string& Atom::GetName(){
return fName;
}

const std::string& Atom::GetElementName()
{
  return fElement;
}

double Atom::GetVanDerWaalsRadius(){
return fVdwRadius;
}


void Atom::SetNext(Atom *AtomNext)
{
  fpNext=AtomNext;
}

void Atom::Scale_units(){
// # Conversion of length from bohr to angstrom
fX=fX/0.529177249 ;
fY=fY/0.529177249 ;
fZ=fZ/0.529177249 ;
printf("scaled x %f y %f z %f \n", fX, fY, fZ);
return;
}


void Atom::Set_data(){


if(fElement == "H "){

orb H1;
H1.l=s;

Set(H1.orbitalsAlpha, 0, 3.425250914);
Set(H1.orbitalsAlpha, 1, 0.6239137298);
Set(H1.orbitalsAlpha, 2, 0.168855404);
Set(H1.orbitalsAlpha, 3, 0.1);

Set(H1.orbitalsConst, 0, 0.1543289673);
Set(H1.orbitalsConst, 1, 0.5353281423);
Set(H1.orbitalsConst, 2, 0.4446345422);
Set(H1.orbitalsConst, 3, 0.);

  fVdwRadius = 1.2;
  fZat=1;

orbitals.push_back(H1);
Bohr_radius = 1.0;
Bragg_radius = 0.25;

}else if(fElement == "He"){

  orb He1;
  He1.l=s;

Set(He1.orbitalsAlpha, 0, 6.3624214);
Set(He1.orbitalsAlpha, 1, 1.159);
Set(He1.orbitalsAlpha, 2, 0.31365);
//   Set(He1.orbitalsAlpha, 3, 0);
Set(He1.orbitalsAlpha, 3, 0.1);

Set(He1.orbitalsConst, 0, 0.1543289673);
Set(He1.orbitalsConst, 1, 0.5353281423);
Set(He1.orbitalsConst, 2, 0.4446345422);
Set(He1.orbitalsConst, 3, 0.);

  fVdwRadius = 1.7;
  fZat=2;
  orbitals.push_back(He1);
  Bohr_radius = 4.0/7.0;


  }else if(fElement == "C "){

    orb C1;
    orb C2;
    orb C3;
    orb C4;
    orb C5;

    C1.l=s;

    Set(C1.orbitalsAlpha, 0, 71.616837);
    Set(C1.orbitalsAlpha, 1, 13.045096);
    Set(C1.orbitalsAlpha, 2, 3.530512);
  //  Set(C1.orbitalsAlpha, 3, 0);
    Set(C1.orbitalsAlpha, 3, 0.1);

    Set(C1.orbitalsConst, 0, 0.1543289673);
    Set(C1.orbitalsConst, 1, 0.5353281423);
    Set(C1.orbitalsConst, 2, 0.4446345422);
    Set(C1.orbitalsConst, 3, 0.);
  //  Set(C1.orbitalsConst, 3, 0.001);
    C2.l=s;

    Set(C2.orbitalsAlpha, 0, 2.94125);
    Set(C2.orbitalsAlpha, 1, 0.683483);
    Set(C2.orbitalsAlpha, 2, 0.22229);
   // Set(C2.orbitalsAlpha, 3, 0);  
    Set(C2.orbitalsAlpha, 3, 0.1);

    Set(C2.orbitalsConst, 0, -0.0999);
    Set(C2.orbitalsConst, 1, 0.3995);
    Set(C2.orbitalsConst, 2, 0.70011);
    Set(C2.orbitalsConst, 3, 0.);

    C3.l=px;


    Set(C3.orbitalsAlpha, 0, 2.94125);
    Set(C3.orbitalsAlpha, 1, 0.683483);
    Set(C3.orbitalsAlpha, 2, 0.22229);
   // Set(C3.orbitalsAlpha, 3, 0);
    Set(C3.orbitalsAlpha, 3, 0.1);

    Set(C3.orbitalsConst, 0, 0.1559);
    Set(C3.orbitalsConst, 1, 0.60768);
    Set(C3.orbitalsConst, 2, 0.39195);
    Set(C3.orbitalsConst, 3, 0.);

    C4.l=py;



    Set(C4.orbitalsAlpha, 0, 2.94125);
    Set(C4.orbitalsAlpha, 1, 0.683483);
    Set(C4.orbitalsAlpha, 2, 0.22229);
  //  Set(C4.orbitalsAlpha, 3, 0);
    Set(C4.orbitalsAlpha, 3, 0.1);

    Set(C4.orbitalsConst, 0,  0.1559);
    Set(C4.orbitalsConst, 1, 0.60768);
    Set(C4.orbitalsConst, 2, 0.39195);
    Set(C4.orbitalsConst, 3, 0.);


    C5.l=pz;

    Set(C5.orbitalsAlpha, 0, 2.94125);
    Set(C5.orbitalsAlpha, 1, 0.683483);
    Set(C5.orbitalsAlpha, 2, 0.22229);
  //  Set(C5.orbitalsAlpha, 3, 0);
    Set(C5.orbitalsAlpha, 3, 0.1);

    Set(C5.orbitalsConst, 0,  0.1559);
    Set(C5.orbitalsConst, 1, 0.60768);
    Set(C5.orbitalsConst, 2, 0.39195);
    Set(C5.orbitalsConst, 3, 0.);

    fVdwRadius = 1.7;
    fZat=6;

    orbitals.push_back(C1);
    orbitals.push_back(C2);
    orbitals.push_back(C3);
    orbitals.push_back(C4);
    orbitals.push_back(C5);


Bohr_radius = 9.0/7.0;
Bragg_radius = 0.7;

  } else if(  fElement == "N ")  {

  orb N1;
  orb N2;
  orb N3;
  orb N4;
  orb N5;

  N1.l=s;

  Set(N1.orbitalsAlpha, 0, 99.10617);
  Set(N1.orbitalsAlpha, 1, 18.052312);
  Set(N1.orbitalsAlpha, 2, 4.88566);
  //Set(N1.orbitalsAlpha, 3, 0);
  Set(N1.orbitalsAlpha, 3, 0.1);

  Set(N1.orbitalsConst, 0, 0.1543289673);
  Set(N1.orbitalsConst, 1, 0.5353281423);
  Set(N1.orbitalsConst, 2, 0.4446345422);
  Set(N1.orbitalsConst, 3, 0.);

  N2.l=s;

  Set(N2.orbitalsAlpha, 0, 3.780456);
  Set(N2.orbitalsAlpha, 1, 0.8785);
  Set(N2.orbitalsAlpha, 2, 0.2857144);
 //   Set(N2.orbitalsAlpha, 3, 0);
  Set(N2.orbitalsAlpha, 3, 0.1);

  Set(N2.orbitalsConst, 0, -0.1);
  Set(N2.orbitalsConst, 1, 0.4);
  Set(N2.orbitalsConst, 2, 0.70011);
  Set(N2.orbitalsConst, 3, 0);

  N3.l=px;

  Set(N3.orbitalsAlpha, 0, 3.780456);
  Set(N3.orbitalsAlpha, 1, 0.8785);
  Set(N3.orbitalsAlpha, 2, 0.285714);
 // Set(N3.orbitalsAlpha, 3, 0);
  Set(N3.orbitalsAlpha, 3, 0.1);

  Set(N3.orbitalsConst, 0, 0.156);
  Set(N3.orbitalsConst, 1, 0.607684);
  Set(N3.orbitalsConst, 2, 0.392);
  Set(N3.orbitalsConst, 3, 0.);

  N4.l=py;

  Set(N4.orbitalsAlpha, 0, 3.780456);
  Set(N4.orbitalsAlpha, 1, 0.8785);
  Set(N4.orbitalsAlpha, 2, 0.285714);
//  Set(N4.orbitalsAlpha, 3, 0);
  Set(N4.orbitalsAlpha, 3, 0.1);


  Set(N4.orbitalsConst, 0, 0.156);
  Set(N4.orbitalsConst, 1, 0.607684);
  Set(N4.orbitalsConst, 2, 0.392);
  Set(N4.orbitalsConst, 3, 0.);

  N5.l=pz;
  
  Set(N5.orbitalsAlpha, 0, 3.780456);
  Set(N5.orbitalsAlpha, 1, 0.8785);
  Set(N5.orbitalsAlpha, 2, 0.285714);
  Set(N5.orbitalsAlpha, 3, 0.1);
//  Set(N5.orbitalsAlpha, 3, 0);

  Set(N5.orbitalsConst, 0, 0.156);
  Set(N5.orbitalsConst, 1, 0.607684);
  Set(N5.orbitalsConst, 2, 0.392);
  Set(N5.orbitalsConst, 3, 0.);

  fVdwRadius = 1.55;
  fZat=7;

  orbitals.push_back(N1);
  orbitals.push_back(N2);
  orbitals.push_back(N3);
  orbitals.push_back(N4);
  orbitals.push_back(N5);

Bohr_radius = 1.0563;
Bragg_radius = 0.65;


} else if(fElement == "N+"){

 orb N_m1;
 orb N_m2;
 orb N_m3;
 orb N_m4;
 orb N_m5;

 N_m1.l=s;

  Set(N_m1.orbitalsAlpha, 0, 99.10617);
  Set(N_m1.orbitalsAlpha, 1, 18.052312);
  Set(N_m1.orbitalsAlpha, 2, 4.88566);
  //Set(N1.orbitalsAlpha, 3, 0);
  Set(N_m1.orbitalsAlpha, 3, 0.1);

  Set(N_m1.orbitalsConst, 0, 0.1543289673);
  Set(N_m1.orbitalsConst, 1, 0.535328);
  Set(N_m1.orbitalsConst, 2, 0.444634542);
  Set(N_m1.orbitalsConst, 3, 0.);

  N_m2.l=s;
  Set(N_m2.orbitalsAlpha, 0, 3.780456);
  Set(N_m2.orbitalsAlpha, 1, 0.878497);
  Set(N_m2.orbitalsAlpha, 2, 0.285714);
  //Set(N1.orbitalsAlpha, 3, 0);
  Set(N_m2.orbitalsAlpha, 3, 0.1);

  Set(N_m2.orbitalsConst, 0, -0.099);
  Set(N_m2.orbitalsConst, 1,  0.444634542);
  Set(N_m2.orbitalsConst, 2,  0.700115);
  Set(N_m2.orbitalsConst, 3,  0.);


  N_m3.l=px;

  Set(N_m3.orbitalsAlpha, 0, 3.780456);
  Set(N_m3.orbitalsAlpha, 1, 0.8785);
  Set(N_m3.orbitalsAlpha, 2, 0.285714);
  //Set(N1.orbitalsAlpha, 3, 0);
  Set(N_m3.orbitalsAlpha, 3, 0.1);

  Set(N_m3.orbitalsConst, 0,  0.155916);
  Set(N_m3.orbitalsConst, 1,  0.607684);
  Set(N_m3.orbitalsConst, 2,  0.392);
  Set(N_m3.orbitalsConst, 3,  0.);

 N_m4.l=py;

  Set(N_m4.orbitalsAlpha, 0, 3.780456);
  Set(N_m4.orbitalsAlpha, 1, 0.8785);
  Set(N_m4.orbitalsAlpha, 2, 0.285714);
  //Set(N1.orbitalsAlpha, 3, 0);
  Set(N_m4.orbitalsAlpha, 3, 0.1);

  Set(N_m4.orbitalsConst, 0,  0.155916);
  Set(N_m4.orbitalsConst, 1,  0.607684);
  Set(N_m4.orbitalsConst, 2,  0.392);
  Set(N_m4.orbitalsConst, 3,  0.);

  N_m5.l=pz;

  Set(N_m5.orbitalsAlpha, 0, 3.780456);
  Set(N_m5.orbitalsAlpha, 1, 0.8785);
  Set(N_m5.orbitalsAlpha, 2, 0.285714);
  //Set(N1.orbitalsAlpha, 3, 0);
  Set(N_m5.orbitalsAlpha, 3, 0.1);

  Set(N_m5.orbitalsConst, 0,  0.155916);
  Set(N_m5.orbitalsConst, 1,  0.607684);
  Set(N_m5.orbitalsConst, 2,  0.392);
  Set(N_m5.orbitalsConst, 3,  0.);

  fVdwRadius = 1.55;
  fZat=7;

  orbitals.push_back(N_m1);
  orbitals.push_back(N_m2);
  orbitals.push_back(N_m3);
  orbitals.push_back(N_m4);
  orbitals.push_back(N_m5);

Bohr_radius = 0.87;
Bragg_radius = 0.58;


} else if(fElement == "O "){

  orb O1;
  orb O2;
  orb O3;
  orb O4;
  orb O5;

  O1.l=s;

  Set(O1.orbitalsAlpha, 0, 130.71);
  Set(O1.orbitalsAlpha, 1, 23.8089);
  Set(O1.orbitalsAlpha, 2, 6.443608);
 // Set(O1.orbitalsAlpha, 3, 0);
  Set(O1.orbitalsAlpha, 3, 0.1);

  Set(O1.orbitalsConst, 0, 0.1543289673);
  Set(O1.orbitalsConst, 1, 0.535328);
  Set(O1.orbitalsConst, 2, 0.4446345422);
  Set(O1.orbitalsConst, 3, 0.);

  O2.l=s;
  Set(O2.orbitalsAlpha, 0, 5.033151);
  Set(O2.orbitalsAlpha, 1, 1.1696);
  Set(O2.orbitalsAlpha, 2, 0.38039);
 // Set(O2.orbitalsAlpha, 3, 0);
  Set(O2.orbitalsAlpha, 3, 0.1);

  Set(O2.orbitalsConst, 0, -0.1);
  Set(O2.orbitalsConst, 1, 0.3995);
  Set(O2.orbitalsConst, 2, 0.70011);
  Set(O2.orbitalsConst, 3, 0.);

  O3.l=px;

  Set(O3.orbitalsAlpha, 0, 5.033151);
  Set(O3.orbitalsAlpha, 1, 1.1696);
  Set(O3.orbitalsAlpha, 2, 0.38039);
  Set(O3.orbitalsAlpha, 3, 0.1);
// Set(O3.orbitalsAlpha, 3, 0);

  Set(O3.orbitalsConst, 0, 0.156);
  Set(O3.orbitalsConst, 1, 0.607684);
  Set(O3.orbitalsConst, 2, 0.392);
  Set(O3.orbitalsConst, 3, 0.);

  O4.l=py;

  Set(O4.orbitalsAlpha, 0, 5.033151);
  Set(O4.orbitalsAlpha, 1, 1.1696);
  Set(O4.orbitalsAlpha, 2, 0.38039);
  //Set(O4.orbitalsAlpha, 3, 0);
  Set(O4.orbitalsAlpha, 3, 0.1);

  Set(O4.orbitalsConst, 0, 0.156);
  Set(O4.orbitalsConst, 1, 0.607684);
  Set(O4.orbitalsConst, 2, 0.392);
  Set(O4.orbitalsConst, 3, 0.);

  O5.l=pz;
  Set(O5.orbitalsAlpha, 0, 5.033151);
  Set(O5.orbitalsAlpha, 1, 1.1696);
  Set(O5.orbitalsAlpha, 2, 0.38039);
  //Set(O5.orbitalsAlpha, 3, 0);
  Set(O5.orbitalsAlpha, 3, 0.1);

  Set(O5.orbitalsConst, 0, 0.156);
  Set(O5.orbitalsConst, 1, 0.607684);
  Set(O5.orbitalsConst, 2, 0.392);
  Set(O5.orbitalsConst, 3, 0.);

  fVdwRadius = 1.52;
  fZat=8;

  orbitals.push_back(O1);
  orbitals.push_back(O2);
  orbitals.push_back(O3);
  orbitals.push_back(O4);
  orbitals.push_back(O5);

Bohr_radius = 0.9;
Bragg_radius = 0.6;

  }
  else if(fElement == "P "){

  orb P1;
  orb P2;
  orb P3;
  orb P4;
  orb P5;
  orb P6;
  orb P7;
  orb P8;
  orb P9;

  P1.l=s;

  Set(P1.orbitalsAlpha, 0, 468.36564);
  Set(P1.orbitalsAlpha, 1,  85.3134);
  Set(P1.orbitalsAlpha, 2,  23.09);
  Set(P1.orbitalsAlpha, 3, 0.1);
//   Set(P1.orbitalsAlpha, 3, 0);

  Set(P1.orbitalsConst, 0, 0.15433);
  Set(P1.orbitalsConst, 1, 0.53533);
  Set(P1.orbitalsConst, 2, 0.444635);
  Set(P1.orbitalsConst, 3, 0.);

  P2.l=s;

  Set(P2.orbitalsAlpha, 0,  28.03264);
  Set(P2.orbitalsAlpha, 1,   6.5142);
  Set(P2.orbitalsAlpha, 2,   2.118614);
  //Set(P2.orbitalsAlpha, 3, 0);
  Set(P2.orbitalsAlpha, 3, 0.1);

  Set(P2.orbitalsConst, 0, -0.1);
  Set(P2.orbitalsConst, 1, 0.4);
  Set(P2.orbitalsConst, 2, 0.7);
  Set(P2.orbitalsConst, 3, 0.);

  Set(P3.orbitalsAlpha, 0, 28.03264);
  Set(P3.orbitalsAlpha, 1,  6.514183);
  Set(P3.orbitalsAlpha, 2,  2.118614);
 // Set(P3.orbitalsAlpha, 3, 0);
  Set(P3.orbitalsAlpha, 3, 0.1);

  Set(P3.orbitalsConst, 0, 0.156);
  Set(P3.orbitalsConst, 1, 0.6077);
  Set(P3.orbitalsConst, 2, 0.392);
  Set(P3.orbitalsConst, 3, 0.);

  P3.l=px;
  Set(P4.orbitalsAlpha, 0, 28.03264);
  Set(P4.orbitalsAlpha, 1,  6.514183);
  Set(P4.orbitalsAlpha, 2,  2.118614);
  Set(P4.orbitalsAlpha, 3, 0.1);
  //  Set(P4.orbitalsAlpha, 3, 0);

  Set(P4.orbitalsConst, 0, 0.156);
  Set(P4.orbitalsConst, 1, 0.6077);
  Set(P4.orbitalsConst, 2, 0.392);
  Set(P4.orbitalsConst, 3, 0.);

  P4.l=py;

  Set(P5.orbitalsAlpha, 0, 28.03264);
  Set(P5.orbitalsAlpha, 1,  6.514183);
  Set(P5.orbitalsAlpha, 2,  2.118614);
 // Set(P5.orbitalsAlpha, 3, 0);
  Set(P5.orbitalsAlpha, 3, 0.1);


  Set(P5.orbitalsConst, 0, 0.156);
  Set(P5.orbitalsConst, 1, 0.6077);
  Set(P5.orbitalsConst, 2, 0.392);
  Set(P5.orbitalsConst, 3, 0.);

  P5.l=pz;
  P6.l=s;

  Set(P6.orbitalsAlpha, 0, 1.7431);
  Set(P6.orbitalsAlpha, 1, 0.48632);
  Set(P6.orbitalsAlpha, 2, 0.190343);
  Set(P6.orbitalsAlpha, 3, 0.1);
 //   Set(P6.orbitalsAlpha, 3, 0);

  Set(P6.orbitalsConst, 0, -0.22);
  Set(P6.orbitalsConst, 1,  0.2256);
  Set(P6.orbitalsConst, 2,  0.9);
  Set(P6.orbitalsConst, 3, 0.);

  P7.l=px;

  Set(P7.orbitalsAlpha, 0, 1.7431);
  Set(P7.orbitalsAlpha, 1, 0.48632);
  Set(P7.orbitalsAlpha, 2, 0.1903);
//  Set(P7.orbitalsAlpha, 3, 0);
  Set(P7.orbitalsAlpha, 3, 0.1);

  Set(P7.orbitalsConst, 0,  0.0106);
  Set(P7.orbitalsConst, 1,  0.595);
  Set(P7.orbitalsConst, 2,  0.462);
  Set(P7.orbitalsConst, 3,  0.);

  P8.l=py;

  Set(P8.orbitalsAlpha, 0, 1.7431);
  Set(P8.orbitalsAlpha, 1, 0.48632);
  Set(P8.orbitalsAlpha, 2, 0.1903);
  Set(P8.orbitalsAlpha, 3, 0.1);
  // Set(P8.orbitalsAlpha, 3, 0);

  Set(P8.orbitalsConst, 0,  0.0106);
  Set(P8.orbitalsConst, 1,  0.595);
  Set(P8.orbitalsConst, 2,  0.462);
  Set(P8.orbitalsConst, 3,  0.);

  P9.l=pz;

  Set(P9.orbitalsAlpha, 0, 1.7431);
  Set(P9.orbitalsAlpha, 1, 0.48632);
  Set(P9.orbitalsAlpha, 2, 0.1903);
  Set(P9.orbitalsAlpha, 3, 0.1);
 //  Set(P9.orbitalsAlpha, 3, 0);

  Set(P9.orbitalsConst, 0,  0.0106);
  Set(P9.orbitalsConst, 1,  0.595);
  Set(P9.orbitalsConst, 2,  0.462);
  Set(P9.orbitalsConst, 3,  0.);

  orbitals.push_back(P1);
  orbitals.push_back(P2);
  orbitals.push_back(P3);
  orbitals.push_back(P4);
  orbitals.push_back(P5);
  orbitals.push_back(P6);
  orbitals.push_back(P7);
  orbitals.push_back(P8);
  orbitals.push_back(P9);

  Bohr_radius = 1.7862;
  Bragg_radius = 1.;
  fVdwRadius = 1.8;
  fZat=15;

} else if( fElement == "S "){

  orb S1;
  orb S2;
  orb S3;
  orb S4;
  orb S5;
  orb S6;
  orb S7;
  orb S8;
  orb S9;

  S1.l=s;
  Set(S1.orbitalsAlpha, 0, 533.12574);
  Set(S1.orbitalsAlpha, 1,  97.11);
  Set(S1.orbitalsAlpha, 2,  26.2816);
  Set(S1.orbitalsAlpha, 3, 0.1);
 // Set(S1.orbitalsAlpha, 3, 0);

  Set(S1.orbitalsConst, 0, 0.15433);
  Set(S1.orbitalsConst, 1, 0.53533);
  Set(S1.orbitalsConst, 2, 0.444635);
  Set(S1.orbitalsConst, 3, 0.);

  S2.l=s;
  Set(S2.orbitalsAlpha, 0,  33.33);
  Set(S2.orbitalsAlpha, 1,  7.74512);
  Set(S2.orbitalsAlpha, 2,  2.519);
 // Set(S2.orbitalsAlpha, 3, 0);
  Set(S2.orbitalsAlpha, 3, 0.1);

  Set(S2.orbitalsConst, 0, -0.1);
  Set(S2.orbitalsConst, 1,  0.4);
  Set(S2.orbitalsConst, 2,  2.0292);
  Set(S2.orbitalsConst, 3, 0.);

  S3.l=px;

  Set(S3.orbitalsAlpha, 0,  33.33);
  Set(S3.orbitalsAlpha, 1,  7.74512);
  Set(S3.orbitalsAlpha, 2,  2.519);
 //   Set(S3.orbitalsAlpha, 3, 0);
  Set(S3.orbitalsAlpha, 3, 0.1);

  Set(S3.orbitalsConst, 0, 0.1559);
  Set(S3.orbitalsConst, 1, 0.607684);
  Set(S3.orbitalsConst, 2, 0.392);
  Set(S3.orbitalsConst, 3, 0.);

  S4.l=py;
  Set(S4.orbitalsAlpha, 0,  33.33);
  Set(S4.orbitalsAlpha, 1,  7.74512);
  Set(S4.orbitalsAlpha, 2,  2.519);
 // Set(S4.orbitalsAlpha, 3, 0);
  Set(S4.orbitalsAlpha, 3, 0.1);

  Set(S4.orbitalsConst, 0, 0.1559);
  Set(S4.orbitalsConst, 1, 0.607684);
  Set(S4.orbitalsConst, 2, 0.392);
  Set(S4.orbitalsConst, 3, 0.);

  S5.l=pz;

  Set(S5.orbitalsAlpha, 0,  33.33);
  Set(S5.orbitalsAlpha, 1,  7.74512);
  Set(S5.orbitalsAlpha, 2,  2.519);
  //Set(S5.orbitalsAlpha, 3, 0);
  Set(S5.orbitalsAlpha, 3, 0.1);

  Set(S5.orbitalsConst, 0, 0.1559);
  Set(S5.orbitalsConst, 1, 0.607684);
  Set(S5.orbitalsConst, 2, 0.392);
  Set(S5.orbitalsConst, 3, 0.);


  S6.l=s;

  Set(S6.orbitalsAlpha, 0,  2.03);
  Set(S6.orbitalsAlpha, 1,  0.566);
  Set(S6.orbitalsAlpha, 2,  0.221583);
  Set(S6.orbitalsAlpha, 3, 0.1);


  Set(S6.orbitalsConst, 0, -0.2196);
  Set(S6.orbitalsConst, 1,  0.2256);
  Set(S6.orbitalsConst, 2,  0.9);
  Set(S6.orbitalsConst, 3, 0.);

  S7.l=px;
  Set(S7.orbitalsAlpha, 0,  2.029);
  Set(S7.orbitalsAlpha, 1,  0.56614);
  Set(S7.orbitalsAlpha, 2,  0.221583);
  //  Set(S7.orbitalsAlpha, 3,  0);
  Set(S7.orbitalsAlpha, 3,  0.1);

  Set(S7.orbitalsConst, 0, 0.010588);
  Set(S7.orbitalsConst, 1, 0.59516);
  Set(S7.orbitalsConst, 2, 0.462);
  Set(S7.orbitalsConst, 3, 0.);

  S8.l=py;

  Set(S8.orbitalsAlpha, 0,  2.029);
  Set(S8.orbitalsAlpha, 1,  0.56614);
  Set(S8.orbitalsAlpha, 2,  0.221583);
  //Set(S8.orbitalsAlpha, 3,  0);
  Set(S8.orbitalsAlpha, 3,  0.1);

  Set(S8.orbitalsConst, 0, 0.010588);
  Set(S8.orbitalsConst, 1, 0.59516);
  Set(S8.orbitalsConst, 2, 0.462);
  Set(S8.orbitalsConst, 3, 0.);

  S9.l=pz;

  Set(S9.orbitalsAlpha, 0,  2.029);
  Set(S9.orbitalsAlpha, 1,  0.56614);
  Set(S9.orbitalsAlpha, 2,  0.221583);
  //  Set(S9.orbitalsAlpha, 3,  0);
  Set(S9.orbitalsAlpha, 3,  0.1);

  Set(S9.orbitalsConst, 0, 0.010588);
  Set(S9.orbitalsConst, 1, 0.59516);
  Set(S9.orbitalsConst, 2, 0.462);
  Set(S9.orbitalsConst, 3, 0.);

  orbitals.push_back(S1);
  orbitals.push_back(S2);
  orbitals.push_back(S3);
  orbitals.push_back(S4);
  orbitals.push_back(S5);
  orbitals.push_back(S6);
  orbitals.push_back(S7);
  orbitals.push_back(S8);
  orbitals.push_back(S9);

fVdwRadius = 1.8;
fZat=16;
Bohr_radius = 1.6489;
Bragg_radius = 1.;


}else{
#ifndef GEANT4
    G4cerr << "Element not recognized : " << fElement << G4endl;
    G4cerr << "Stop now" << G4endl;
    exit(1);
#else
    G4ExceptionDescription errMsg;
    errMsg << "Element not recognized : " << fElement << G4endl;
    G4Exception("PDBlib::Load","ELEM_NOT_RECOGNIZED",
                FatalException,errMsg);
#endif
  }

}







