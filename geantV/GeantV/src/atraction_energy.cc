
#include "atraction_energy.hh"

using namespace std;

Electron_Nuclei_atraction::Electron_Nuclei_atraction(){
}


Double_v Electron_Nuclei_atraction::AtractV(G4double alpha_1, Double_v Alpha2, vector<G4int> fv1,
vector<G4int> fv2, Double_v X1, Double_v  Y1, Double_v Z1, Double_v X2, 
Double_v Y2, Double_v Z2 , Double_v rx, Double_v ry, Double_v rz ){

Operators fOpx;
Double_v N, M, kar;

bool polyn=true;
bool selector = true;

Double_v gam = alpha_1 + Alpha2;
Double_v prd = (X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2);

kar = (2.0*pi/gam)*vecCore::math::Exp (-(alpha_1*Alpha2/gam)*prd);

N=fOpx.Gamma_aproxV(gam, rx, ry, rz, X1, Y1, Z1, X2, Y2, Z2, 
fv1[0], fv1[1], fv1[2], fv2[0],  fv2[1], fv2[2], 
alpha_1, Alpha2, polyn, selector);


M=N*kar;

return M;
}



vector<vector<G4double>> Electron_Nuclei_atraction::NuclearV(Molecule& Mol){

vector<Atom *>::iterator atr_main;
vector<Atom *>::iterator atr_row;
vector<Atom *>::iterator atr_r;
vector<Residue *>::iterator res_main;
vector<Residue *>::iterator res_row;
vector<Residue *>::iterator res_nuc;

G4double cont=0;
G4double cont_nuc=0;
G4double buf_cont=0;
G4int con_aux_main=0;
G4int con_aux_row=0;
G4int atomic_cont=0;
G4int orbital_principal=0;
G4int orbital_sec=0;

Double_v Atrac;
Double_v Vat;

vector<vector<G4double>> Sum;
vector<G4double> buf_Sum;

for(res_main = Mol.Residuos_cadena.begin();
    res_main != Mol.Residuos_cadena.end() ; ++res_main){

for(atr_main = (*res_main)->Lista_atoms.begin();
    atr_main != (*res_main)->Lista_atoms.end(); ++ atr_main){

for((*atr_main)->shell1 = (*atr_main)->orbitals.begin();
    (*atr_main)->shell1 != (*atr_main)->orbitals.end(); ++ (*atr_main)->shell1){


con_aux_main=con_aux_main+ 1;
con_aux_row=0;
orbital_sec = 0;

for(res_row = Mol.Residuos_cadena.begin();
    res_row != Mol.Residuos_cadena.end() ; ++res_row){

for(atr_row = (*res_row)->Lista_atoms.begin();
    atr_row != (*res_row)->Lista_atoms.end(); ++ atr_row){


 for((*atr_row)->shell2 = (*atr_row)->orbitals.begin();
     (*atr_row)->shell2 != (*atr_row)->orbitals.end(); ++ (*atr_row)->shell2){

 con_aux_row=con_aux_row+1;
 atomic_cont=0;

if(   orbital_principal  <= orbital_sec  ){

for (int i = 0; i < 3; i++) {

  atomic_cont = 0;

  for(res_nuc = Mol.Residuos_cadena.begin();
      res_nuc != Mol.Residuos_cadena.end() ; ++res_nuc){

  for(atr_r = (*res_nuc)->Lista_atoms.begin();
      atr_r != (*res_nuc)->Lista_atoms.end(); ++ atr_r){

  atomic_cont = atomic_cont + 1;

Vat = -1.0*((*atr_r)->fZat)*AtractV( Get(  (*(*atr_main)->shell1).orbitalsAlpha  , i ), (*(*atr_row)->shell2).orbitalsAlpha,
(*(*atr_main)->shell1).l, (*(*atr_row)->shell2).l,
(*atr_main)->fX, (*atr_main)->fY, (*atr_main)->fZ,
(*atr_row)->fX,(*atr_row)->fY, (*atr_row)->fZ,
(*atr_r)->fX,(*atr_r)->fY, (*atr_r)->fZ );


Atrac=Vat*( Get( (*(*atr_main)->shell1).orbitalsNormals , i) )*( (*(*atr_row)->shell2).orbitalsNormals );

for (G4int j = 0; j < kVecLenD - 1; ++j) {
cont=cont +  Get(  Atrac , j)* Get(  (*(*atr_row)->shell2).orbitalsConst  , j);

}


cont=cont*( Get( (*(*atr_main)->shell1).orbitalsConst , i ) );

buf_cont = buf_cont+cont;
cont=0;
}
cont_nuc = cont_nuc+ buf_cont ;

buf_cont=0;
}
}


buf_Sum.push_back(cont_nuc);
cont_nuc=0;


}else{


buf_Sum.push_back(  Sum[ orbital_sec ][  orbital_principal  ]   );
cont_nuc=0;

}

orbital_sec++;


}
}
}

Sum.push_back(buf_Sum);
buf_Sum.clear();
orbital_principal++;
}
}

}

return Sum;
}

