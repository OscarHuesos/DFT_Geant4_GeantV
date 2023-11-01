
#include "atraction_energy.hh"

using namespace std;

Electron_Nuclei_atraction::Electron_Nuclei_atraction(){
}


G4double Electron_Nuclei_atraction::Atract(G4double alpha_1, 
G4double alpha_2, vector<G4int> fv1, vector<G4int> fv2, 
G4double x1, G4double  y1, G4double z1,
G4double x2,G4double y2, G4double z2 , G4double rx,
G4double ry, G4double rz ){

Operators fOpx;
bool polyn = true;

// false = old method
bool selector = true;
G4double N, M,  res;
G4double gam = alpha_1 + alpha_2;
G4double prd = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)
+(z1-z2)*(z1-z2);

res=(2.0*pi/gam)*G4Exp(-(alpha_1*alpha_2/gam)*prd);

N=fOpx.Gamma_aprox(gam, rx, ry, rz, x1, y1, z1, x2, y2, z2, 
fv1[0], fv1[1], fv1[2], fv2[0],  fv2[1], fv2[2], alpha_1, alpha_2, polyn, selector);
M=N*res;

return M;
}



vector<vector<G4double>> Electron_Nuclei_atraction::Nuclear(Molecule& Mol){

vector<Atom *>::iterator atr_main;
vector<Atom *>::iterator atr_row;
vector<Atom *>::iterator atr_r;
vector<Residue *>::iterator res_main;
vector<Residue *>::iterator res_row;
vector<Residue *>::iterator res_nuc;

G4double cont=0;
G4double cont_nuc=0;
G4double buf_cont=0;

G4double V=0, Knorms=0;
G4int con_aux_main=0;
G4int con_aux_row=0;
G4int atomic_cont=0;
G4int orbital_principal=0;
G4int orbital_sec=0;

vector<vector<G4double>> Sum;
vector<G4double> buf_Sum;

for(res_main = Mol.Residuos_cadena.begin();
    res_main != Mol.Residuos_cadena.end() ; ++res_main){

for(atr_main = (*res_main)->Lista_atoms.begin();
    atr_main != (*res_main)->Lista_atoms.end(); ++ atr_main){

for((*atr_main)->shell1 = (*atr_main)->orbitals.begin();
    (*atr_main)->shell1 != (*atr_main)->orbitals.end(); ++ (*atr_main)->shell1){

con_aux_main=con_aux_main+1;
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

  atomic_cont=atomic_cont+1;

for (int j = 0; j < 3; j++) {


V= -1.0*((*atr_r)->fZat)*Atract( (*(*atr_main)->shell1).NAC[i].A , (*(*atr_row)->shell2).NAC[j].A ,
(*(*atr_main)->shell1).l, (*(*atr_row)->shell2).l,
(*atr_main)->fX, (*atr_main)->fY, (*atr_main)->fZ,
(*atr_row)->fX,(*atr_row)->fY, (*atr_row)->fZ,
(*atr_r)->fX,(*atr_r)->fY, (*atr_r)->fZ);


Knorms=V*( (*(*atr_main)->shell1).NAC[i].N )*( (*(*atr_row)->shell2).NAC[j].N );
cont=cont+ ( (*(*atr_row)->shell2).NAC[j].C )*Knorms;

}

cont=cont*( (*(*atr_main)->shell1).NAC[i].C );

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

