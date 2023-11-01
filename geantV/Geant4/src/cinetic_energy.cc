
#include "cinetic_energy.hh"

using namespace std;

Electron_Cinetic::Electron_Cinetic(){
}


G4double Electron_Cinetic::Kin(G4double alpha_1, G4double alpha_2, vector<G4int> fv1,
vector<G4int> fv2, G4double x1, G4double  y1, G4double z1,
G4double x2,G4double y2, G4double z2 , bool side){

Operators fOpx;
G4double t1, t2, t3, L1 , L2, L3;

//alpha2 path
if (side == true){

t1=fOpx.Oxab(alpha_1, alpha_2, fv1[0],fv1[1],
fv1[2],fv2[0]+2, fv2[1], fv2[2],x1, y1, z1, x2, y2, z2);
t2=fOpx.Oxab(alpha_1,alpha_2, fv1[0], fv1[1], 
fv1[2],fv2[0], fv2[1]+2, fv2[2],x1, y1,z1, x2, y2, z2);
t3= fOpx.Oxab(alpha_1, alpha_2, fv1[0], fv1[1],
fv1[2],fv2[0], fv2[1], fv2[2]+2,x1, y1, z1, x2, y2, z2);


L1 = 4*alpha_2*alpha_2*(t1 + t2 + t3) ;

L2 = -2.0*alpha_2*( 2.0*(fv2[0]+fv2[1]+fv2[2])+ 3 )*(fOpx.Oxab(alpha_1,
alpha_2, fv1[0], fv1[1], fv1[2], fv2[0], fv2[1], fv2[2], x1, y1, z1, x2, y2, z2))  ;

L3 =   fv2[0]*(fv2[0]-1)*(fOpx.Oxab(alpha_1,
alpha_2, fv1[0] , fv1[1], fv1[2],fv2[0] - 2, fv2[1], fv2[2] ,x1, y1, z1, x2, y2,
z2)) + fv2[1]*(fv2[1]-1)*(fOpx.Oxab(alpha_1,alpha_2, fv1[0],fv1[1], fv1[2],
fv2[0], fv2[1]-2, fv2[2] ,x1, y1, z1, x2, y2, z2)) + fv2[2]*(fv2[2]-1)*(fOpx.Oxab(alpha_1,
alpha_2, fv1[0],fv1[1], fv1[2], fv2[0], fv2[1], fv2[2]- 2 ,x1, y1, z1, x2, y2, z2) )   ;

//alpha1 path
}else{

t1=fOpx.Oxab(alpha_1, alpha_2, fv1[0]+2,fv1[1],
fv1[2],fv2[0], fv2[1], fv2[2],x1, y1, z1, x2, y2, z2);
t2=fOpx.Oxab(alpha_1,alpha_2, fv1[0], fv1[1]+2, 
fv1[2],fv2[0], fv2[1], fv2[2],x1, y1,z1, x2, y2, z2);
t3= fOpx.Oxab(alpha_1, alpha_2, fv1[0], fv1[1],
fv1[2]+2,fv2[0], fv2[1], fv2[2],x1, y1, z1, x2, y2, z2);

L1= 4*alpha_1*alpha_1*(t1+ t2 + t3);

L2=   -2.0*alpha_1*(2*(fv1[0]+fv1[1]+fv1[2])+3)*(fOpx.Oxab(alpha_1,
alpha_2, fv1[0], fv1[1], fv1[2], fv2[0], fv2[1], fv2[2], x1, y1, z1, x2, y2, z2))  ;


L3=  fv1[0]*(fv1[0]-1)*(fOpx.Oxab(alpha_1,
alpha_2, fv1[0]-2,fv1[1], fv1[2],fv2[0], fv2[1], fv2[2] ,x1, y1, z1, x2, y2,
z2))+fv1[1]*(fv1[1]-1)*(fOpx.Oxab(alpha_1,alpha_2, fv1[0],fv1[1]-2, fv1[2],
fv2[0], fv2[1], fv2[2] ,x1, y1, z1, x2, y2, z2))+fv1[2]*(fv1[2]-1)*(fOpx.Oxab(alpha_1,
alpha_2, fv1[0],fv1[1], fv1[2]-2, fv2[0], fv2[1], fv2[2] ,x1, y1, z1, x2, y2, z2) )  ;


}


return L1 + L2 + L3;
}


vector<vector<G4double>> Electron_Cinetic::Grad(Molecule& Mol ){

G4double cont=0;
G4double buf_cont=0;
G4double K=0;
G4double Knorms=0;
G4int orbital_principal=0;
G4int orbital_sec=0;
bool typo = false;
G4int contador = 0;

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

  for (int j = 0; j < 3; j++) {

K=Kin(   (*(*atr_main)->shell1).NAC[i].A ,  (*(*atr_row)->shell2).NAC[j].A , (*(*atr_main)->shell1).l,
(*(*atr_row)->shell2).l, (*atr_main)->fX, (*atr_main)->fY,
(*atr_main)->fZ, (*atr_row)->fX,(*atr_row)->fY, (*atr_row)->fZ, typo);

Knorms=K*( (*(*atr_main)->shell1).NAC[i].N  )*( (*(*atr_row)->shell2).NAC[j].N );

cont=cont +  ( (*(*atr_row)->shell2).NAC[j].C)*Knorms;

}

cont=cont*( (*(*atr_main)->shell1).NAC[i].C );
buf_cont=buf_cont+cont;

cont=0;
}


buf_Sum.push_back(-0.5*buf_cont);
buf_cont=0;

contador++;

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

