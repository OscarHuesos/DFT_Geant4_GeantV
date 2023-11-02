
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


Atom::Atom(G4int id, const std::string& n,
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
  return fId;
}

G4int Atom::Get_atomic_number(){
return  fZat;
}


const std::string& Atom::GetName()
{
  return fName;
}


const std::string& Atom::GetElementName()
{
  return fElement;
}


double Atom::GetVanDerWaalsRadius(){
return fVdwRadius;
}


void Atom::SetNext(Atom *AtomNext){
fpNext=AtomNext;
}

void Atom::Scale_units(){
// # Conversion of length from bohr to angstrom
fX=fX/0.529177249 ;
fY=fY/0.529177249 ;
fZ=fZ/0.529177249 ;

return;
}


void Atom::Set_data(){

printf("f element: %s \n", fElement.c_str());
if(fElement == "H "){

  orb H1;
  Norm_alpha_coff s11_H1;
  s11_H1.A=3.425250914;
  s11_H1.C=0.1543289673;
  s11_H1.N=0;

  Norm_alpha_coff s11_H2;
  s11_H2.A=0.6239137298;
  s11_H2.C=0.5353281423;
  s11_H2.N=0;

  Norm_alpha_coff s11_H3;
  s11_H3.A=0.168855404;
  s11_H3.C=0.4446345422;
  s11_H3.N=0;

  H1.l=s;
  H1.NAC[0]=s11_H1;
  H1.NAC[1]=s11_H2;
  H1.NAC[2]=s11_H3;

  fVdwRadius = 1.2;
  fZat=1;

orbitals.push_back(H1);
Bohr_radius = 1.0;
Bragg_radius = 0.25;
//No_electrons = 1;


}else if(fElement == "He"){

  orb He1;

  Norm_alpha_coff s12_He1;
  s12_He1.A=6.3624214;
  s12_He1.C=0.1543289673;
  s12_He1.N=0;

  Norm_alpha_coff s12_He2;
  s12_He2.A=1.159;
  s12_He2.C=0.5353281423;
  s12_He2.N=0;

  Norm_alpha_coff s12_He3;
  s12_He3.A=0.31365;
  s12_He3.C=0.4446345422;
  s12_He3.N=0;

  He1.l=s;
  He1.NAC[0]=s12_He1;
  He1.NAC[1]=s12_He2;
  He1.NAC[2]=s12_He3;

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

    Norm_alpha_coff s12_C1;
    s12_C1.A=71.616837;
    s12_C1.C=0.1543289673;
    s12_C1.N=0;

    Norm_alpha_coff s12_C2;
    s12_C2.A=13.045096;
    s12_C2.C=0.5353281423;
    s12_C2.N=0;

    Norm_alpha_coff s12_C3;
    s12_C3.A=3.530512;
    s12_C3.C=0.4446345422;
    s12_C3.N=0;

    C1.l=s;
    C1.NAC[0]=s12_C1;
    C1.NAC[1]=s12_C2;
    C1.NAC[2]=s12_C3;

    Norm_alpha_coff s22_C1;
    s22_C1.A= 2.94125;
    s22_C1.C= -0.0999;
    s22_C1.N= 0;

    Norm_alpha_coff s22_C2;
    s22_C2.A=0.683483;
    s22_C2.C=0.3995;
    s22_C2.N=0;

    Norm_alpha_coff s22_C3;
    s22_C3.A=0.22229;
    s22_C3.C=0.70011;
    s22_C3.N=0;

    C2.l=s;
    C2.NAC[0]=s22_C1;
    C2.NAC[1]=s22_C2;
    C2.NAC[2]=s22_C3;

    Norm_alpha_coff p22_C1;
    p22_C1.A=2.94125;
    p22_C1.C=0.1559;
    p22_C1.N=0;

    Norm_alpha_coff p22_C2;
    p22_C2.A=0.683483;
    p22_C2.C=0.60768;
    p22_C2.N=0;

    Norm_alpha_coff p22_C3;
    p22_C3.A=0.22229;
    p22_C3.C=0.39195;
    p22_C3.N=0;

    C3.l=px;
    C3.NAC[0]=p22_C1;
    C3.NAC[1]=p22_C2;
    C3.NAC[2]=p22_C3;

    C4.l=py;
    C4.NAC[0]=p22_C1;
    C4.NAC[1]=p22_C2;
    C4.NAC[2]=p22_C3;

    C5.l=pz;
    C5.NAC[0]=p22_C1;
    C5.NAC[1]=p22_C2;
    C5.NAC[2]=p22_C3;

    fVdwRadius = 1.7;
    fZat=6;

    orbitals.push_back(C1);
    orbitals.push_back(C2);
    orbitals.push_back(C3);
    orbitals.push_back(C4);
    orbitals.push_back(C5);

//Bohr_radius.push_back(0.1789);
Bohr_radius = 9.0/7.0;
Bragg_radius = 0.7;


  } else if(  fElement == "N ")  {

  orb N1;
  orb N2;
  orb N3;
  orb N4;
  orb N5;

  Norm_alpha_coff s12_N1;
  s12_N1.A= 99.10617;
  s12_N1.C= 0.1543289673;
  s12_N1.N= 0;

  Norm_alpha_coff s12_N2;
  s12_N2.A= 18.052312;
  s12_N2.C= 0.5353281423;
  s12_N2.N= 0;

  Norm_alpha_coff s12_N3;
  s12_N3.A= 4.88566;
  s12_N3.C= 0.4446345422;
  s12_N3.N= 0;

  N1.l=s;
  N1.NAC[0]=s12_N1;
  N1.NAC[1]=s12_N2;
  N1.NAC[2]=s12_N3;

  Norm_alpha_coff s22_N1;
  s22_N1.A= 3.780456;
  s22_N1.C=-0.1;
  s22_N1.N= 0;

  Norm_alpha_coff s22_N2;
  s22_N2.A= 0.8785;
  s22_N2.C= 0.4;
  s22_N2.N= 0;

  Norm_alpha_coff s22_N3;
  s22_N3.A= 0.2857144;
  s22_N3.C= 0.70011;
  s22_N3.N= 0;

  N2.l=s;
  N2.NAC[0]=s22_N1;
  N2.NAC[1]=s22_N2;
  N2.NAC[2]=s22_N3;

  Norm_alpha_coff p23_N1;
  p23_N1.A= 3.780456;
  p23_N1.C= 0.156;
  p23_N1.N= 0;

  Norm_alpha_coff p23_N2;
  p23_N2.A= 0.8785;
  p23_N2.C= 0.607684;
  p23_N2.N= 0;

  Norm_alpha_coff p23_N3;
  p23_N3.A= 0.285714;
  p23_N3.C= 0.392;
  p23_N3.N= 0;

  N3.l=px;
  N3.NAC[0]=p23_N1;
  N3.NAC[1]=p23_N2;
  N3.NAC[2]=p23_N3;

  N4.l=py;
  N4.NAC[0]=p23_N1;
  N4.NAC[1]=p23_N2;
  N4.NAC[2]=p23_N3;

  N5.l=pz;
  N5.NAC[0]=p23_N1;
  N5.NAC[1]=p23_N2;
  N5.NAC[2]=p23_N3;

  fVdwRadius = 1.55;
  fZat=7;


  orbitals.push_back(N1);
  orbitals.push_back(N2);
  orbitals.push_back(N3);
  orbitals.push_back(N4);
  orbitals.push_back(N5);

Bohr_radius = 1.0563;
Bragg_radius = 0.65;
//Bohr_radius.push_back(0.9863);

} else if(fElement == "N+"){

printf("N+1 dio \n");

 orb N_m1;
 orb N_m2;
 orb N_m3;
 orb N_m4;
 orb N_m5;

 Norm_alpha_coff s12_N_m1;
 s12_N_m1.A= 99.10617;
 s12_N_m1.C= 0.1543289673;
 s12_N_m1.N= 0;

 Norm_alpha_coff s12_N_m2;
 s12_N_m2.A= 18.052312;
 s12_N_m2.C= 0.535328;
 s12_N_m2.N= 0;

 Norm_alpha_coff s12_N_m3;
 s12_N_m3.A= 4.88566;
 s12_N_m3.C= 0.444634542;
 s12_N_m3.N= 0;

 N_m1.l=s;
 N_m1.NAC[0]=s12_N_m1;
 N_m1.NAC[1]=s12_N_m2;
 N_m1.NAC[2]=s12_N_m3;

 Norm_alpha_coff s22_N_m1;
 s22_N_m1.A= 3.780456;
 s22_N_m1.C= -0.099;
 s22_N_m1.N= 0;

 Norm_alpha_coff s22_N_m2;
 s22_N_m2.A= 0.878497;
 s22_N_m2.C= 0.444634542;
 s22_N_m2.N= 0;

 Norm_alpha_coff s22_N_m3;
 s22_N_m3.A= 0.285714;
 s22_N_m3.C= 0.700115;
 s22_N_m3.N= 0;

  N_m2.l=s;
  N_m2.NAC[0]=s22_N_m1;
  N_m2.NAC[1]=s22_N_m2;
  N_m2.NAC[2]=s22_N_m3;

  Norm_alpha_coff p22_N_m1;
  p22_N_m1.A= 3.780456;
  p22_N_m1.C= 0.155916;
  p22_N_m1.N= 0;

  Norm_alpha_coff p22_N_m2;
  p22_N_m2.A= 0.8785;
  p22_N_m2.C= 0.607684;
  p22_N_m2.N= 0;

  Norm_alpha_coff p22_N_m3;
  p22_N_m3.A= 0.285714;
  p22_N_m3.C= 0.392;
  p22_N_m3.N= 0;

  N_m3.l=px;
  N_m3.NAC[0]=p22_N_m1;
  N_m3.NAC[1]=p22_N_m2;
  N_m3.NAC[2]=p22_N_m3;

  N_m4.l=py;
  N_m4.NAC[0]=p22_N_m1;
  N_m4.NAC[1]=p22_N_m2;
  N_m4.NAC[2]=p22_N_m3;

  N_m5.l=pz;
  N_m5.NAC[0]=p22_N_m1;
  N_m5.NAC[1]=p22_N_m2;
  N_m5.NAC[2]=p22_N_m3;

  fVdwRadius = 1.55;
  fZat=7;
  ////////////////////////////////////////
  // int  size_Basis_onda=5;

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

  Norm_alpha_coff s12_O1;
  s12_O1.A= 130.71;
  s12_O1.C= 0.1543289673;
  s12_O1.N= 0.0;

  Norm_alpha_coff s12_O2;
  s12_O2.A= 23.8089;
  s12_O2.C= 0.535328;
  s12_O2.N= 0.0;

  Norm_alpha_coff s12_O3;
  s12_O3.A= 6.443608;
  s12_O3.C= 0.4446345422;
  s12_O3.N= 0.0;

  O1.l=s;
  O1.NAC[0]=s12_O1;
  O1.NAC[1]=s12_O2;
  O1.NAC[2]=s12_O3;

  Norm_alpha_coff s22_O1;
  s22_O1.A=  5.033151;
  s22_O1.C= -0.1;
  s22_O1.N=  0.0;

  Norm_alpha_coff s22_O2;
  s22_O2.A= 1.1696;
  s22_O2.C= 0.3995;
  s22_O2.N= 0.0;

  Norm_alpha_coff s22_O3;
  s22_O3.A= 0.38039;
  s22_O3.C= 0.70011;
  s22_O3.N= 0.0;

  O2.l=s;
  O2.NAC[0]=s22_O1;
  O2.NAC[1]=s22_O2;
  O2.NAC[2]=s22_O3;

  Norm_alpha_coff p24_O1;
  p24_O1.A= 5.033151;
  p24_O1.C= 0.156;
  p24_O1.N= 0.0;

  Norm_alpha_coff p24_O2;
  p24_O2.A= 1.1696;
  p24_O2.C= 0.607684;
  p24_O2.N= 0.0;

  Norm_alpha_coff p24_O3;
  p24_O3.A= 0.38039;
  p24_O3.C= 0.392;
  p24_O3.N= 0.0;

  O3.l=px;
  O3.NAC[0]=p24_O1;
  O3.NAC[1]=p24_O2;
  O3.NAC[2]=p24_O3;

  O4.l=py;
  O4.NAC[0]=p24_O1;
  O4.NAC[1]=p24_O2;
  O4.NAC[2]=p24_O3;

  O5.l=pz;
  O5.NAC[0]=p24_O1;
  O5.NAC[1]=p24_O2;
  O5.NAC[2]=p24_O3;

  fVdwRadius = 1.52;
  //  sizeZ=3*2*2;
  fZat=8;

  orbitals.push_back(O1);
  orbitals.push_back(O2);
  orbitals.push_back(O3);
  orbitals.push_back(O4);
  orbitals.push_back(O5);

//Bohr_radius.push_back(0.1399);
Bohr_radius = 0.9;
Bragg_radius = 0.6;

}else if(fElement == "P "){

  orb P1;
  orb P2;
  orb P3;
  orb P4;
  orb P5;
  orb P6;
  orb P7;
  orb P8;
  orb P9;

  Norm_alpha_coff s12_P1;
  s12_P1.A= 468.36564;
  s12_P1.C= 0.15433;
  s12_P1.N= 0.0;

  Norm_alpha_coff s12_P2;
  s12_P2.A= 85.3134;
  s12_P2.C= 0.53533;
  s12_P2.N= 0.0;

  Norm_alpha_coff s12_P3;
  s12_P3.A= 23.09;
  s12_P3.C= 0.444635;
  s12_P3.N= 0.0;

  P1.l=s;
  P1.NAC[0]=s12_P1;
  P1.NAC[1]=s12_P2;
  P1.NAC[2]=s12_P3;

  Norm_alpha_coff s22_P1;
  s22_P1.A=  28.03264 ;
  s22_P1.C= -0.1;
  s22_P1.N=  0.0;

  Norm_alpha_coff s22_P2;
  s22_P2.A= 6.5142;
  s22_P2.C= 0.4;
  s22_P2.N= 0.0;

  Norm_alpha_coff s22_P3;
  s22_P3.A= 2.118614;
  s22_P3.C= 0.7;
  s22_P3.N= 0.0;

  P2.l=s;
  P2.NAC[0]=s22_P1;
  P2.NAC[1]=s22_P2;
  P2.NAC[2]=s22_P3;


  Norm_alpha_coff p26_P1;
  p26_P1.A= 28.03264 ;
  p26_P1.C= 0.156;
  p26_P1.N= 0.0;

  Norm_alpha_coff p26_P2;
  p26_P2.A= 6.514183;
  p26_P2.C= 0.6077;
  p26_P2.N= 0.0;

  Norm_alpha_coff p26_P3;
  p26_P3.A= 2.118614;
  p26_P3.C= 0.392;
  p26_P3.N= 0.0;

  P3.l=px;
  P3.NAC[0]=p26_P1;
  P3.NAC[1]=p26_P2;
  P3.NAC[2]=p26_P3;

  P4.l=py;
  P4.NAC[0]=p26_P1;
  P4.NAC[1]=p26_P2;
  P4.NAC[2]=p26_P3;

  P5.l=pz;
  P5.NAC[0]=p26_P1;
  P5.NAC[1]=p26_P2;
  P5.NAC[2]=p26_P3;

  Norm_alpha_coff s32_P1;
  s32_P1.A= 1.7431;
  s32_P1.C= -0.22;
  s32_P1.N= 0.0;

  Norm_alpha_coff s32_P2;
  s32_P2.A= 0.48632;
  s32_P2.C= 0.2256;
  s32_P2.N= 0.0;

  Norm_alpha_coff s32_P3;
  s32_P3.A= 0.190343;
  s32_P3.C= 0.9;
  s32_P3.N= 0.0;

  P6.l=s;
  P6.NAC[0]=s32_P1;
  P6.NAC[1]=s32_P2;
  P6.NAC[2]=s32_P3;


 Norm_alpha_coff p33_P1;
  p33_P1.A= 1.7431;
  p33_P1.C= 0.0106;
  p33_P1.N= 0.0;

  Norm_alpha_coff p33_P2;
  p33_P2.A= 0.48632;
  p33_P2.C= 0.595;
  p33_P2.N= 0.0;

  Norm_alpha_coff p33_P3;
  p33_P3.A= 0.1903;
  p33_P3.C= 0.462;
  p33_P3.N= 0.0;

  P7.l=px;
  P7.NAC[0]=p33_P1;
  P7.NAC[1]=p33_P2;
  P7.NAC[2]=p33_P3;

  P8.l=py;
  P8.NAC[0]=p33_P1;
  P8.NAC[1]=p33_P2;
  P8.NAC[2]=p33_P3;

  P9.l=pz;
  P9.NAC[0]=p33_P1;
  P9.NAC[1]=p33_P2;
  P9.NAC[2]=p33_P3;

  orbitals.push_back(P1);
  orbitals.push_back(P2);
  orbitals.push_back(P3);
  orbitals.push_back(P4);
  orbitals.push_back(P5);
  orbitals.push_back(P6);
  orbitals.push_back(P7);
  orbitals.push_back(P8);
  orbitals.push_back(P9);

  fVdwRadius = 1.8;
  fZat=15;

  Bohr_radius = 1.7862;
  Bragg_radius = 1.;
  //Bohr_radius.push_back(1.4311);

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

  Norm_alpha_coff s12_S1;
  s12_S1.A= 533.12574;
  s12_S1.C= 0.15433;
  s12_S1.N= 0.0;

  Norm_alpha_coff s12_S2;
  s12_S2.A= 97.11;
  s12_S2.C= 0.53533;
  s12_S2.N= 0.0;

  Norm_alpha_coff s12_S3;
  s12_S3.A= 26.2816;
  s12_S3.C= 0.444635;
  s12_S3.N= 0.0;

  S1.l=s;
  S1.NAC[0]=s12_S1;
  S1.NAC[1]=s12_S2;
  S1.NAC[2]=s12_S3;

  Norm_alpha_coff s22_S1;
  s22_S1.A=  33.33;
  s22_S1.C=  -0.1;
  s22_S1.N=  0.0;

  Norm_alpha_coff s22_S2;
  s22_S2.A= 7.74512;
  s22_S2.C= 0.4;
  s22_S2.N= 0.0;

  Norm_alpha_coff s22_S3;
  s22_S3.A= 2.519;
  s22_S3.C= 2.0292;
  s22_S3.N= 0.0;

  S2.l=s;
  S2.NAC[0]=s22_S1;
  S2.NAC[1]=s22_S2;
  S2.NAC[2]=s22_S3;

  Norm_alpha_coff p26_S1;
  p26_S1.A= 33.33;
  p26_S1.C= 0.1559;
  p26_S1.N= 0.0;

  Norm_alpha_coff p26_S2;
  p26_S2.A= 7.74512;
  p26_S2.C= 0.607684;
  p26_S2.N= 0.0;

  Norm_alpha_coff p26_S3;
  p26_S3.A= 2.519;
  p26_S3.C= 0.392;
  p26_S3.N= 0.0;

  S3.l=px;
  S3.NAC[0]=p26_S1;
  S3.NAC[1]=p26_S2;
  S3.NAC[2]=p26_S3;

  S4.l=py;
  S4.NAC[0]=p26_S1;
  S4.NAC[1]=p26_S2;
  S4.NAC[2]=p26_S3;

  S5.l=pz;
  S5.NAC[0]=p26_S1;
  S5.NAC[1]=p26_S2;
  S5.NAC[2]=p26_S3;

  Norm_alpha_coff s32_S1;
  s32_S1.A= 2.03;
  s32_S1.C= -0.2196;
  s32_S1.N= 0.0;

  Norm_alpha_coff s32_S2;
  s32_S2.A= 0.566;
  s32_S2.C= 0.2256;
  s32_S2.N= 0.0;

  Norm_alpha_coff s32_S3;
  s32_S3.A= 0.221583;
  s32_S3.C= 0.9;
  s32_S3.N= 0.0;

  S6.l=s;
  S6.NAC[0]=s32_S1;
  S6.NAC[1]=s32_S2;
  S6.NAC[2]=s32_S3;


  Norm_alpha_coff p34_S1;
  p34_S1.A= 2.029;
  p34_S1.C= 0.010588;
  p34_S1.N= 0.0;


  Norm_alpha_coff p34_S2;
  p34_S2.A= 0.56614;
  p34_S2.C= 0.59516;
  p34_S2.N= 0.0;

  Norm_alpha_coff p34_S3;
  p34_S3.A= 0.221583;
  p34_S3.C= 0.462;
  p34_S3.N= 0.0;

  S7.l=px;
  S7.NAC[0]=p34_S1;
  S7.NAC[1]=p34_S2;
  S7.NAC[2]=p34_S3;

  S8.l=py;
  S8.NAC[0]=p34_S1;
  S8.NAC[1]=p34_S2;
  S8.NAC[2]=p34_S3;

  S9.l=pz;
  S9.NAC[0]=p34_S1;
  S9.NAC[1]=p34_S2;
  S9.NAC[2]=p34_S3;

  fVdwRadius = 1.8;
  fZat=16;

  orbitals.push_back(S1);
  orbitals.push_back(S2);
  orbitals.push_back(S3);
  orbitals.push_back(S4);
  orbitals.push_back(S5);
  orbitals.push_back(S6);
  orbitals.push_back(S7);
  orbitals.push_back(S8);
  orbitals.push_back(S9);

Bohr_radius = 1.6489;
Bragg_radius = 1.;
//Bohr_radius.push_back(9.0/8.0);

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





