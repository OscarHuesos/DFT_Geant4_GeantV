#include "dft.hh"
using namespace std;

DFT::DFT() {
}

Eigen::MatrixXd convert_vvd_to_matrix(vector<vector<G4double> > vvd) {

    size_t n_rows = vvd.size();
    size_t n_cols = vvd.at(0).size();

    Eigen::MatrixXd result(n_rows, n_cols);
    result.row(0) = Eigen::VectorXd::Map(&vvd[0][0], n_cols);

    for (size_t i = 1; i < n_rows; i++) {

        // Make sure that every row of vvd has the same size.
        if (n_cols != vvd.at(i).size()) {
            char buffer[200];
            snprintf(buffer, 200,
                        "vvd[%ld] size (%ld) does not match vvd[0] size (%ld)",
                        i, vvd.at(i).size(), n_cols);
            string err_mesg(buffer);
            throw invalid_argument(err_mesg);
        }

        result.row(i) = Eigen::VectorXd::Map(&vvd[i][0], n_cols);
    }

    return result;
}


void adding( vector<vector<G4double>> V1, vector<vector<G4double>>& V  ){
G4int siz = V1.size();
vector<vector<G4double>>::iterator Col;
vector<G4double>::iterator Row;

for(int i = 0; i < siz ; i++) {
for(int j = 0; j < siz ; j++) {
V[i][j] = V1[i][j] + V[i][j];
}
}

return;
}

void suma( vector<vector<G4double>> V1,
vector<vector<G4double>> V2, vector<vector<G4double>>& Su  ){
G4int siz = V1.size();
vector<vector<G4double>>::iterator Col;
vector<G4double>::iterator Row;

Su.resize(siz);
for(size_t i = 0; i < siz; ++i){
    Su[i].resize(siz);
}

for(int i = 0; i < siz ; i++) {
for(int j = 0; j < siz ; j++) {
Su[i][j] = V1[i][j] + V2[i][j];

}
}
return;
}



G4double suma_completo( vector<vector<G4double>> V  ){

G4int siz = V.size();
G4double sum=0;

for(int i = 0; i < siz ; i++) {
for(int j = 0; j < siz ; j++) {
sum = sum +  V[i][j] ;

}
}

return sum;
}


G4double mult_sum_completo( vector<vector<G4double>> V1, vector<vector<G4double>> V2 ){

G4int size = V1.size();
G4double sum=0;

for(int i = 0; i < size ; i++) {
for(int j = 0; j < size ; j++) {
 sum = V1[i][j]*V2[i][j] + sum;
}
}
return sum;
}


vector<vector<G4double>> DFT::Matriz_D( Eigen::MatrixXd  M,
G4int orb_ocupados,  G4int  orb_abiertos , G4int  orbits ){

G4int it = orb_ocupados + orb_abiertos;
Eigen::MatrixXd  D(orbits, it);
Eigen::MatrixXd  DT;
vector<vector<G4double>> S;
vector<G4double> buffer_D;
G4int cont = 0;
G4double ff =0;

for( int i = 0; i < orbits ; i++) {

 for( int j = 0; j < it ; j++) {

D(i,j) = M(i,j);
}
}

DT = D.transpose();

if(  orb_abiertos == 0    ){

for( int i = 0; i < orbits ; i++) {

 for( int j = 0; j < orbits ; j++) {
ff = 0;
 for( int k = 0; k < orb_ocupados ; k++) {
 
ff = ff + D(i, k)*DT(k,j);
}
 
buffer_D.push_back( 2.0*ff );
}
S.push_back(buffer_D);
buffer_D.clear();
}

}else{

for( int i = 0; i < orbits ; i++) {
cont = 0;
 for( int j = 0; j < orbits ; j++) {
ff = 0;

for( int k = 0; k < it ; k++) {

ff = ff + D(i, k)*DT(k,j);

}


if(cont < orb_ocupados){
ff = 2.0*ff;
}else{

}

cont++;

buffer_D.push_back( ff );

}

S.push_back(buffer_D);
buffer_D.clear();

}
}

return S;
}


vector<vector<G4double>>  DFT::Rep_density (){

vector<vector<vector<vector<G4double>>>>::iterator Colright;
vector<vector<vector<G4double>>>::iterator Rowright;
vector<vector<G4double>>::iterator Colleft;
vector<G4double>::iterator Rowleft;
vector<G4double> iter;
vector<vector<G4double>> DJ;

G4int index2 = 0;
G4int index3 = 0;
G4int index4 = 0;
G4double acum = 0;

for(Colright = J.begin(); Colright != J.end() ; ++Colright){  
index2=0;  
 for(Rowright = Colright->begin(); Rowright != Colright->end(); ++ Rowright){
  index3=0;
  // D:
   for(Colleft = Rowright->begin(); Colleft != Rowright->end() ; ++Colleft){
    index4 = 0;
    for(Rowleft = Colleft->begin(); Rowleft != Colleft->end(); ++ Rowleft){
        acum = acum + ( D[index3][index4] )*(*Rowleft);

index4++;
}
index3++;
}
iter.push_back( acum ) ;
index2++;
acum = 0;
}
DJ.push_back(iter);
iter.clear();
}

return DJ;
}


void DFT::nuclear_nuclear(Molecule& Mol){

vector<Atom *>::iterator atr;
vector<Residue *>::iterator red;

vector<Atom *>::iterator atr_dos;
vector<Residue *>::iterator red_dos;

G4double VV = 0;
G4int index = 0;

for(red = Mol.Residuos_cadena.begin();
    red != Mol.Residuos_cadena.end() ; ++red){

for(atr = (*red)->Lista_atoms.begin();
    atr != (*red)->Lista_atoms.end(); ++ atr){

for(red_dos = Mol.Residuos_cadena.begin();
    red_dos !=  Mol.Residuos_cadena.end() ; ++red_dos ){

for(atr_dos = (*red_dos)->Lista_atoms.begin() + index + 1;
    atr_dos != (*red_dos)->Lista_atoms.end(); ++ atr_dos ){

VV = VV + ((*atr)->fZat*(*atr_dos)->fZat)/fabs(sqrt( ((*atr)->fX - (*atr_dos)->fX )*((*atr)->fX - (*atr_dos)->fX ) +
((*atr)->fY - (*atr_dos)->fY)*((*atr)->fY - (*atr_dos)->fY) + ((*atr)->fZ - (*atr_dos)->fZ)*((*atr)->fZ - (*atr_dos)->fZ )  ));

 }
}
index++;
}
}


Vnn = VV;
return;
}


void DFT::normalizedV(Molecule& Mol){

vector<Atom *>::iterator atr;
vector<Residue *>::iterator red;

G4int orbitales = 0;
vector<G4double> norms;

G4double m;
G4double n1=0;
G4double n2=0;
G4double n3=0;

G4int l1=0;
G4int l2=0;
G4int l3=0;

for(red = Mol.Residuos_cadena.begin();
  red != Mol.Residuos_cadena.end() ; ++red){

 for(atr = (*red)->Lista_atoms.begin();
  atr != (*red)->Lista_atoms.end(); ++ atr){


for((*atr)->shell1 = (*atr)->orbitals.begin();
    (*atr)->shell1 != (*atr)->orbitals.end(); ++ (*atr)->shell1){

  l1= (*(*atr)->shell1).l[0];
  l2= (*(*atr)->shell1).l[1];
  l3= (*(*atr)->shell1).l[2];

orbitales =  orbitales + 1;

  n1=fOperators.doublefactorial( 2*l1-1);
  n2=fOperators.doublefactorial( 2*l2-1);
  n3=fOperators.doublefactorial( 2*l3-1);

  for (int i = 0; i < 3; i++) {

m = pow(2.0*(Get (   (*(*atr)->shell1).orbitalsAlpha , i ) )/pi, 0.75)*
  (pow(4.0*(Get (   (*(*atr)->shell1).orbitalsAlpha , i ) ) , 0.5*(l1+l2+l3))/sqrt(n1*n2*n3) );

Set(   (*(*atr)->shell1).orbitalsNormals , i , m );


}

}

}

}

Mol.Setvorbits(orbitales);
return;
}


void DFT::over(Molecule& Mol ){
S = fOv.iterVector(Mol);
return;
}

void DFT::cin(Molecule& Mol ){
T=fElectron_K.GradV(Mol);
return;
}

void DFT::atrac(Molecule& Mol ){

V=fAtract.NuclearV(Mol);

return;
}

void DFT::rep(Molecule& Mol ){

J=fRepulsion.pairingV(Mol);

return;
}


void DFT::correlation(Molecule& Mol){

// func  0 = LDA, 1 = GGA,
G4int func = 1;
W = fCorrelation.GridV(Mol, func);

return;
}


G4double DFT::exchange(Molecule& Mol ){

// 0 sera lda, 1 gga,
G4int functional = 1;
Vxc = fExchange.XCpotentialV(Mol, D, W, functional);
return fExchange.Energy_xc ;
}


void DFT::SCF(Molecule& Mol,   G4double& dft_energy ){

size_t siz = S.size();
G4int iteraciones = 100;
G4double tol = 0.00001;

vector<vector<G4double>>::iterator Col;
vector<G4double>::iterator Row;
vector<vector<G4double>> one_electron;
vector<vector<G4double>> buffer;

buffer.resize(siz);
for(size_t i = 0; i < siz; ++i){
buffer[i].resize(siz);
}


suma(V , T, one_electron);
G4double exh;
G4double Ein = -0.;

Eigen::MatrixXd B = convert_vvd_to_matrix(one_electron);
Eigen::MatrixXd Os = convert_vvd_to_matrix(S);

Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es( B,
Os, Eigen::ComputeEigenvectors);

D = Matriz_D(es.eigenvectors(), Mol.ocupados, Mol.abiertos, siz);

G4double dif ;
G4int i = 0;

/////////////////////////////////////////////////////////////////////////////////
while(   i <= iteraciones    ){

JD = Rep_density();
exh = exchange(Mol);


buffer = one_electron;

adding(JD ,buffer);
adding(Vxc, buffer);

//energy calc:
energia = Vnn + mult_sum_completo( one_electron, D ) + 
0.5*mult_sum_completo( JD, D ) + exh;

Eigen::MatrixXd F = convert_vvd_to_matrix(buffer);

Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> Eng (F,
Os, Eigen::ComputeEigenvectors);

dif = fabs( energia  - Ein );
if(   tol >   dif){
break;
}


D = Matriz_D(Eng.eigenvectors(), Mol.ocupados, Mol.abiertos, siz);

buffer.clear();
Ein = energia;
i++;
}

dft_energy = energia;
//printf("energia final: %f \n",  energia);
return;
}



