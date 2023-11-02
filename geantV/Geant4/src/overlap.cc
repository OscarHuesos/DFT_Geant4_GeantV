
#include "overlap.hh"

using namespace std;

Overlap_Int::Overlap_Int(){
}


vector<vector<G4double>>  Overlap_Int::iter(Molecule& Mol){

  Operators fOp;
  vector<Atom *>::iterator atr_main;
  vector<Atom *>::iterator atr_row;
  vector<Residue *>::iterator res_main;
  vector<Residue *>::iterator res_row;
  vector<vector<G4double>> Sum;
  vector<G4double> buf_Sum;

  G4double cont=0;
  G4int con_aux_main=0;
  G4int con_aux_row=0;
  G4double buf_cont=0;
  G4double Ox;
  G4double S;

  G4int orbital_principal=0;
  G4int orbital_sec=0;

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

if(   orbital_principal  <= orbital_sec  ){

cont=0;

 for (int i = 0; i < 3; i++) {

  for (int j = 0; j < 3; j++) {

 Ox=fOp.Oxab( (*(*atr_main)->shell1).NAC[i].A , (*(*atr_row)->shell2).NAC[j].A, 
(*(*atr_main)->shell1).l[0], (*(*atr_main)->shell1).l[1],  (*(*atr_main)->shell1).l[2],
(*(*atr_row)->shell2).l[0], (*(*atr_row)->shell2).l[1], (*(*atr_row)->shell2).l[2], 
(*atr_main)->fX, (*atr_main)->fY, (*atr_main)->fZ, (*atr_row)->fX, (*atr_row)->fY, (*atr_row)->fZ);


S=Ox*( (*(*atr_main)->shell1).NAC[i].N  )*( (*(*atr_row)->shell2).NAC[j].N );
cont= cont +  ((*(*atr_row)->shell2).NAC[j].C)*S;
}


cont=cont*(  (*(*atr_main)->shell1).NAC[i].C );
buf_cont=buf_cont+cont;
cont=0;
}

buf_Sum.push_back(buf_cont);
buf_cont=0;

}else{

buf_Sum.push_back(  Sum[ orbital_sec ][  orbital_principal  ]   );
buf_cont=0;

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

