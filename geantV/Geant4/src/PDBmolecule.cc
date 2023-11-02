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

//MIo tambien

#include "PDBmolecule.hh"

Residue::Residue():
fResName(""),fResSeq(0),fVisible(false),fSelected(false),
fCenterX(0),fCenterY(0),fCenterZ(0),fNbAtom(0),fpNext(0),fpFirst(0)
{
}


Residue::Residue(const std::string& rN,int rS){
  fResName=rN;
  fResSeq=rS; 
  fVisible=false;
  fSelected=false;
  fCenterX=0;
  fCenterY=0;
  fCenterZ=0;
  fNbAtom=0;
  fpNext=0;
  fpFirst=0;
}


Residue *Residue::GetNext(){
return fpNext;
}

Atom *Residue::GetFirst(){
return fpFirst;
}


int Residue::GetID(){
return fResSeq;
}


void Residue::SetNext(Residue *residueNext){
fpNext=residueNext;
}


void Residue::SetFirst(Atom *atomFirst){
fpFirst=atomFirst;
}



Molecule::Molecule():
fMolName(""),fMolNum(0),fMinGlobZ(0),fMaxGlobZ(0),
fMinGlobX(0),fMaxGlobX(0),fMinGlobY(0),fMaxGlobY(0),
fCenterX(0),fCenterY(0),fCenterZ(0),fDistCenterMax(0),fNbResidue(0),
fpNext(0),fpFirst(0)
{
}


Molecule::Molecule(const std::string& mN,int mNum){
  fMolName=mN;  
  fMolNum=mNum; 
  fMinGlobZ=0;
  fMaxGlobZ=0;
  fMinGlobX=0;
  fMaxGlobX=0;
  fMinGlobY=0;
  fMaxGlobY=0;
  fCenterX=0;
  fCenterY=0;
  fCenterZ=0;
  fDistCenterMax=0;
  fNbResidue=0;
  fExiste=true;
  fpNext=0;
  fpFirst=0;
 Nresidue=0;
No_atomos_en_prot=0;
no_orbits=0;
}

Molecule *Molecule::GetNext(){
  return fpNext;
}


Residue *Molecule::GetFirst(){
  return fpFirst;
}


int Molecule::GetID(){
  return fMolNum;
}


void Molecule::SetNext(Molecule *moleculeNext){
fpNext=moleculeNext;
}


void Molecule::SetFirst(Residue *resFirst){
fpFirst=resFirst;
}


void Molecule::isbool(bool t){
fExiste=t;
return;
}

bool Molecule::retstatus(){
return fExiste;
}


void Molecule::Setvorbits(G4int orb){
no_orbits = orb;
return;
}


void Molecule::Load(  const std::string& filename,  int& isProtein ){

    try{

      char fail[5];
      fail[0]='f';
      fail[1]='a';
      fail[2]='i';
      fail[3]='l';
      fail[4]='_';

  puts(filename.c_str());
  std::ifstream Fin;
  Fin.open(filename.c_str());
  if(!Fin)throw strcat(fail,filename.c_str());
  
    G4String buffer;
    G4String ll;
    char Buffer_res[] = "apple";
    char Buffer_cadena[]= "1";
      G4String atomName; 
      G4String chainName=""; 
      G4String element; 
      G4String resName; 
      G4String nameMol;
      G4int res_num;
      unsigned short int numModel = 0;
      unsigned short int model = 0;
      fNbResidue=0;

      int residuo_buffer_per_chain = 0;
    Residue *residuobuffer= new Residue;
    residuobuffer->fNbAtom=0;
    G4String carb("CA");

  //PROTEIN:
  int terMax = 50; 
  isProtein = 1;
  cadena = 0;
  while (getline(Fin, buffer)){

  if((buffer.substr(0, 6)).compare("NUMMDL") == 0){
  std::istringstream((buffer.substr(10, 4))) >> numModel;
  }
  if((numModel > 0) && ((buffer.substr(0, 6)).compare("MODEL ") == 0))  {
  std::istringstream((buffer.substr(10, 4))) >> model;
  //////////////////////////////////////////
  if(model > 1) break;

  }
  
  if((buffer.substr(0, 6)).compare("SEQRES") == 0){

}


if(buffer.substr(0,4)=="ATOM"){
puts(buffer);

Atom* atomo = new Atom(No_atomos_en_prot,
buffer.substr(11,5).c_str(),
buffer.substr(17,3).c_str(),
atoi(buffer.substr(22,6).c_str()),
0, atof(buffer.substr(30,8).c_str()),
atof(buffer.substr(38,8).c_str()),
atof(buffer.substr(46,8).c_str()),
 atoi(buffer.substr(56,1).c_str()),
atof(buffer.substr(60,6).c_str()),
buffer.substr(21,1).c_str(),
buffer.substr(77,2).c_str() );

atomo->Set_data();
no_electrons = no_electrons + atomo->Get_atomic_number();
res_num =  atoi(buffer.substr(22,6).c_str());

No_atomos_en_prot++;
atomo->Scale_units();


if(atomo->fNumInRes!=fNbResidue){
if(fNbResidue==0){
cadena =1;
residuobuffer->fResName=atomo->fResName;
strcpy(Buffer_res, buffer.substr(17,3).c_str());
residuobuffer->fChainame=atomo->fChainame;
residuobuffer->Lista_atoms.push_back(atomo);
residuobuffer->fNbAtom++;
residuobuffer->fChainN=cadena;

}else{

if(strcmp(Buffer_cadena,atomo->fChainame )!=0){
cadena=cadena+1;

N_residues_per_chain.push_back(residuo_buffer_per_chain);
residuo_buffer_per_chain=0;


}
Residuos_cadena.push_back(residuobuffer);
fNbResidue++;
residuo_buffer_per_chain++;

residuobuffer->Lista_atoms.clear();
residuobuffer->fNbAtom=0;
residuobuffer->Lista_atoms.push_back(atomo);
residuobuffer->fNbAtom++;
residuobuffer->fResName=atomo->fResName;
residuobuffer->fChainame=atomo->fChainame;
residuobuffer->fChainN=cadena;
residuobuffer->fResSeq=atomo->fNumInRes;

}
}else{

residuobuffer->Lista_atoms.push_back(atomo);
residuobuffer->fNbAtom++;
residuobuffer->fResName=atomo->fResName;
residuobuffer->fChainame=atomo->fChainame;
residuobuffer->fResSeq=atomo->fNumInRes;
residuobuffer->fChainN=cadena;

if(strcmp(Buffer_cadena,atomo->fChainame )!=0){

cadena=cadena+1;
N_residues_per_chain.push_back(residuo_buffer_per_chain);
residuo_buffer_per_chain=0;

}

}

strcpy(Buffer_res, buffer.substr(17,3).c_str());
strcpy(Buffer_cadena, atomo->fChainame);


}
}

residuobuffer->fResSeq=fNbResidue;
residuobuffer->fResName=Buffer_res;
residuobuffer->fChainame=Buffer_cadena;
residuobuffer->fChainN=cadena;


Residuos_cadena.push_back(residuobuffer);
fNbResidue =  res_num;

printf("Numero de residuos finales: %d y cadenas %d \n",
fNbResidue, cadena);
Fin.close();

}
catch(char* pMsg) { std::cerr << std::endl << "Exception:" << pMsg << std::endl;
}
return;
}


void Molecule::check_parity(){

G4int m =  no_electrons % 2;

if(m == 1){
//impar
abiertos = 1;

}else{

abiertos = 0;
}

G4int op = no_electrons - abiertos; 
ocupados = op/2;

return;  
}











