
#include "cinetic_energy.hh"

using namespace std;


Electron_Cinetic::Electron_Cinetic(){
}


Double_v Electron_Cinetic::KinV(G4double Alpha1, Double_v Alpha2,
G4int k1, G4int m1, G4int n1, G4int k2, G4int m2, G4int n2,
Double_v X1, Double_v X2, Double_v Y1, Double_v Y2, 
Double_v Z1, Double_v Z2){

Operators fOpx;

Double_v Q= Double_v(0);
Double_v F= Double_v(0);
Double_v H= Double_v(0);

Double_v LOW3U= Double_v(0);
Double_v LOW3D= Double_v(0);
Double_v LOW3T= Double_v(0);

Q= fOpx.OXABVectortestU(Alpha1, Alpha2, k1+2,m1, n1,
k2, m2, n2, X1, X2, Y1, Y2, Z1, Z2);

F= fOpx.OXABVectortestU(Alpha1, Alpha2,  k1 ,m1+2, n1,
k2, m2, n2, X1 , X2, Y1, Y2, Z1, Z2);

H= fOpx.OXABVectortestU(Alpha1, Alpha2, k1 ,m1, n1+2,
k2, m2, n2, X1 , X2, Y1, Y2, Z1, Z2);

Double_v  L1= 4.0*Alpha1*Alpha1*( Q + F  + H);

Double_v  L2 = -2.0*Alpha1*(2*(k1+ m1 +n1)+3)*(
fOpx.OXABVectortestU(Alpha1, Alpha2, k1 ,m1, n1,
k2, m2, n2, X1 , X2, Y1, Y2, Z1, Z2 )  );

LOW3U = fOpx.OXABVectortestU(Alpha1, Alpha2, k1-2 ,m1, n1,
k2, m2, n2, X1 , X2, Y1, Y2, Z1, Z2);
LOW3D =  fOpx.OXABVectortestU(Alpha1, Alpha2, k1 ,m1-2, n1,
k2, m2, n2, X1 , X2, Y1, Y2, Z1, Z2);
LOW3T = fOpx.OXABVectortestU(Alpha1, Alpha2, k1 ,m1, n1-2,
k2, m2, n2, X1 , X2, Y1, Y2, Z1, Z2);

Double_v  L3 = k1*(k1-1)*LOW3U+m1*(m1-1)*LOW3D + n1*(n1-1)*LOW3T;

return L1+L2+L3;

}    


vector<vector<G4double>> Electron_Cinetic::GradV(Molecule& Mol ){

Operators fOp;
G4double cont=0;
G4double buf_cont=0;
G4int orbital_principal=0;
G4int orbital_sec=0;
Double_v KV = 0.;
Double_v KVN = 0. ;
vector<Atom *>::iterator atr_main;
vector<Atom *>::iterator atr_row;
vector<Residue *>::iterator res_main;
vector<Residue *>::iterator res_row;

vector<vector<G4double>> Sum;
vector<G4double> buf_Sum;

for(res_main = Mol.Residuos_cadena.begin();
    res_main != Mol.Residuos_cadena.end() ; ++res_main){

 for(atr_main = (*res_main)->Lista_atoms.begin();
    atr_main != (*res_main)->Lista_atoms.end(); ++ atr_main){

   for((*atr_main)->shell1 = (*atr_main)->orbitals.begin();
       (*atr_main)->shell1 != (*atr_main)->orbitals.end(); ++ (*atr_main)->shell1){

  orbital_sec = 0;

 for(res_row = Mol.Residuos_cadena.begin();
     res_row != Mol.Residuos_cadena.end() ; ++res_row){

  for(atr_row = (*res_row)->Lista_atoms.begin();
      atr_row != (*res_row)->Lista_atoms.end(); ++ atr_row){

for((*atr_row)->shell2 = (*atr_row)->orbitals.begin();
    (*atr_row)->shell2 != (*atr_row)->orbitals.end(); ++ (*atr_row)->shell2){


if(   orbital_principal  <= orbital_sec  ){

  for (int i = 0; i < 3; i++) {

  KV= KinV(    Get(   (*(*atr_main)->shell1).orbitalsAlpha  , i  )  , (*(*atr_row)->shell2).orbitalsAlpha,
(*(*atr_main)->shell1).l[0], (*(*atr_main)->shell1).l[1],  (*(*atr_main)->shell1).l[2],
(*(*atr_row)->shell2).l[0], (*(*atr_row)->shell2).l[1], (*(*atr_row)->shell2).l[2],
(*atr_main)->fX,  (*atr_row)->fX, (*atr_main)->fY, (*atr_row)->fY, (*atr_main)->fZ,  (*atr_row)->fZ );


KVN=KV*( Get(   (*(*atr_main)->shell1).orbitalsNormals  , i  )  )*( (*(*atr_row)->shell2).orbitalsNormals );


 
for (G4int j = 0; j < kVecLenD - 1; ++j) {
cont = cont +  Get(KVN, j)* Get(  (*(*atr_row)->shell2).orbitalsConst  , j);

}


cont = cont*(Get (  (*(*atr_main)->shell1).orbitalsConst , i) );
buf_cont = buf_cont + cont;
cont=0;
}

buf_Sum.push_back(-0.5*buf_cont);
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







