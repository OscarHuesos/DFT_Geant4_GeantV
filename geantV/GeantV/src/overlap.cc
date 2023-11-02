

#include "overlap.hh"
using namespace std;

Overlap_Int::Overlap_Int(){
}


vector<vector<G4double>>  Overlap_Int::iterVector(Molecule& Mol){

Operators fOp;
vector<Atom *>::iterator atr_main;
vector<Atom *>::iterator atr_row;
vector<Residue *>::iterator res_main;
vector<Residue *>::iterator res_row;
vector<vector<G4double>> SumV;
vector<G4double> buf_SumV;

G4double cont=0;
G4int con_aux_main=0;
G4int con_aux_row=0;
G4double buf_cont=0;
 
Double_v Sv = 0.;
Double_v R= 0.;
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

 for (G4int i = 0; i < 3; i++) {
 
R =fOp.OXABVectortestU( Get( (*(*atr_main)->shell1).orbitalsAlpha , i ) , (*(*atr_row)->shell2).orbitalsAlpha,
(*(*atr_main)->shell1).l[0],   (*(*atr_main)->shell1).l[1],  (*(*atr_main)->shell1).l[2],
(*(*atr_row)->shell2).l[0], (*(*atr_row)->shell2).l[1], (*(*atr_row)->shell2).l[2], 
(*atr_main)->fX,  (*atr_row)->fX, (*atr_main)->fY, (*atr_row)->fY, (*atr_main)->fZ,  (*atr_row)->fZ);


Sv=R*( Get(  (*(*atr_main)->shell1).orbitalsNormals , i ) )*( (*(*atr_row)->shell2).orbitalsNormals );

for (G4int j = 0; j < kVecLenD - 1; ++j) {
cont = cont +  Get(Sv, j)* Get(  (*(*atr_row)->shell2).orbitalsConst  , j);

}



cont= cont*( Get( (*(*atr_main)->shell1).orbitalsConst , i ) );

buf_cont=buf_cont+cont;
cont=0;

}


buf_SumV.push_back(buf_cont);
buf_cont=0;


}else{

buf_SumV.push_back(  SumV[ orbital_sec ][  orbital_principal  ]   );
buf_cont=0;

}

orbital_sec++;


}
}
}

SumV.push_back(buf_SumV);
buf_SumV.clear();
orbital_principal++;

}
}
}


return SumV;
}



