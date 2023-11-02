#include "repulsion_energy.hh"

using namespace std;


Electron_electron_repulsion::Electron_electron_repulsion(){
}


// a0m00 o a0mc0
Double_v Electron_electron_repulsion::repulsion_a0mc0V( G4double alf_1, G4double alf_2, G4double alf_3,
Double_v Alf4, G4int L10, G4int L11, G4int L12, G4int L20, G4int L21, G4int L22,
Double_v X1, Double_v Y1, Double_v Z1,
Double_v X2, Double_v Y2, Double_v Z2,
Double_v X3, Double_v Y3, Double_v Z3,
Double_v X4, Double_v Y4, Double_v Z4,
G4double n1, G4double n2, G4double n3, 
Double_v N4, G4int cc ){

Double_v AC = 0.;
G4int L = L10 + L11 + L12 + L20 + L21 + L22 ;

G4double gam1 = alf_1 + alf_2;

Double_v Gam2 = alf_3 + Alf4;

Double_v PX = (alf_1*X1 + alf_2*X2)/gam1;
Double_v PY = (alf_1*Y1 + alf_2*Y2)/gam1;
Double_v PZ = (alf_1*Z1 + alf_2*Z2)/gam1;

Double_v QX = (alf_3*X3 + Alf4*X4)/Gam2;
Double_v QY = (alf_3*Y3 + Alf4*Y4)/Gam2;
Double_v QZ = (alf_3*Z3 + Alf4*Z4)/Gam2;
Double_v G = gam1 + Gam2;

Operators* fOpx = new Operators();

Double_v WX = (gam1*PX + Gam2*QX)/(G);
Double_v WY = (gam1*PY + Gam2*QY)/(G);
Double_v WZ = (gam1*PZ + Gam2*QZ)/(G);

Double_v Norms = n1*n2*n3*N4;

Double_v D1 = fOpx->NaivePowerVector( X1 - X2, 2) +
fOpx->NaivePowerVector( Y1 - Y2, 2) + fOpx->NaivePowerVector( Z1 - Z2, 2) ;

Double_v D2 =  fOpx->NaivePowerVector( X3 - X4, 2) + 
fOpx->NaivePowerVector( Y3 - Y4, 2) + fOpx->NaivePowerVector( Z3 - Z4, 2) ;


Double_v KP = vecCore::math::Exp( -D1*(alf_1*alf_2)/gam1);
Double_v KQ = vecCore::math::Exp( -D2*(alf_3*Alf4)/Gam2);

Double_v RX = PX - QX;
Double_v RY = PY - QY;
Double_v RZ = PZ - QZ;


Double_v R = RX*RX + RY*RY + RZ*RZ;
Double_v OM = gam1*Gam2/G;


Double_v T = OM*R;
Double_v F = 0. ;

bool polyn = false;

// el viejo es selec = false;
bool selector = true;

// check if 00mc0:
if (  (L10 == 0)  &&  (L11 == 0) &&  (L12 == 0)  ){

// check caso base  00m00:
if ( (L20 == 0) &&  (L21 == 0) &&  (L22 == 0)  )  {


Double_v* OOOMOOO = new Double_v[L+1];
MaskD_v mt = ( T > 100 );

if( MaskFull(mt)  ){

OOOMOOO[0] = vecCore::math::Sqrt( pi/(4*T)  );


if (L > 0){

for (int i = 0; i < L ; i++) {
OOOMOOO[i + 1] = (2.0*i+1)*OOOMOOO[i]/(2.0*T);
}

}

}else{

OOOMOOO[L] = fOpx->Gamma_aproxV( OM , RX, RY, RZ, 1, 1, 1,
1, 1, 1, L10, L11, L12, L20, L21, L22, gam1, Gam2, polyn, selector );

delete fOpx;

if (L > 0){
for(int i = L ; i > 0; i--) { 
OOOMOOO[i - 1] = (2.0*T*OOOMOOO[i] + vecCore::math::Exp(-T) )/( 2.0*(i-1)+1 );
}

}

} //else T > 100

for (int j = 0; j < L + 1 ; j++) {

OOOMOOO[j] = ( OOOMOOO[j]*Norms*pow(2,1)*pow(pi,5.0/2.0)*KP*KQ )/ 
( gam1*Gam2*( vecCore::math::Sqrt(G)  )  );

}

AC = OOOMOOO[0];

delete [] OOOMOOO;

return AC;

}else{ 

F = 0. ;

Double_v ****ARRC;
ARRC = new Double_v***[L20 + 1]; 
for(int i = 0; i < L20 + 1 ; i++){
ARRC[i] = new Double_v**[ L21 + 1 ];
 for(int j=0; j <  L21 + 1 ; j++) {
ARRC[i][j] = new Double_v*[L22 + 1];  
 for(int k=0; k <  L22 + 1 ; k++) {
ARRC[i][j][k] = new Double_v[L + 1];
}
}
}


MaskD_v mt = ( T > 100 );
if( MaskFull(mt)  ){

//F =  vecCore::math::Sqrt( pi/(4*T)  );
ARRC[0][0][0][0] = vecCore::math::Sqrt( pi/(4*T)  );
//f= sqrt(pi/(4*T));

if (L > 0){
for (int i = 0; i < L ; i++) {
ARRC[0][0][0][i + 1] = (2.0*i+1)*ARRC[0][0][0][i]/(2.0*T);
}

}

}else{

ARRC[0][0][0][L] = fOpx->Gamma_aproxV( OM , RX, RY, RZ, 1, 1, 1,
1, 1, 1, L10, L11, L12, L20, L21, L22, gam1, Gam2, polyn, selector );

delete fOpx;

if (L > 0){
for(int i = L ; i > 0; i--) { 

ARRC[0][0][0][i - 1] = (2.0*T*ARRC[0][0][0][i] + vecCore::math::Exp(-T) )/( 2.0*(i-1)+1 );
}

}

} //else T > 100

for (int j = 0; j < L + 1 ; j++) {

ARRC[0][0][0][j] = ( ARRC[0][0][0][j]*Norms*pow(2,1)*pow(pi,5.0/2.0)*KP*KQ )/ 
( gam1*Gam2*( vecCore::math::Sqrt(G)  )  );

}


if (  L20 > 0  ){

for(int i = 1; i < L20 + 1 ; i++){
 for (int l = 0; l < L  ; l++) {


F = (QX -X3)*ARRC[i-1][0][0][l] + (WX - QX)*ARRC[i-1][0][0][l+1] ; 

if( (i-1) > 0  ){
F = F + ( (i-1)/(2.0*Gam2))*( ARRC[i-2][0][0][l] - (OM/Gam2)*ARRC[i-2][0][0][l+1] ) ;

}
ARRC[i][0][0][l] = F;
}
}
}


if (  L21 > 0   ){
for(int j = 1; j < L21 + 1 ; j++){
 for (int l = 0; l < L  ; l++) {

F = (QY - Y3)*ARRC[0][j-1][0][l] + (WY - QY)*ARRC[0][j-1][0][l+1];

if( (j-1) > 0  ){

F = F + ( (j-1)/(2.0*Gam2))*( ARRC[0][j-2][0][l] - (OM/Gam2)*ARRC[0][j-2][0][l+1] ) ;

}

ARRC[0][j][0][l] = F;
}
}
}


if (  L22 > 0   ){
for(int k = 1; k < L22 + 1 ; k++){
 for (int l = 0; l < L  ; l++) {

F = (QZ - Z3)*ARRC[0][0][k-1][l] + (WZ - QZ)*ARRC[0][0][k-1][l+1]  ;

if( (k-1) > 0  ){
F = F + ( (k-1)/(2.0*Gam2))*(ARRC[0][0][k-2][l] - (OM/Gam2)*ARRC[0][0][k-2][l+1] ) ;

}
ARRC[0][0][k][l] = F;

}
}
}

if(   (L20 >= 1) && (  L21 >= 1 )  ) {
F=0;
//000 | nm0
for(int jj = 1; jj < L21 + 1 ; jj++){

  for(int ii = 1; ii < L20 + 1 ; ii++){
for (int l = 0; l < L  ; l++) {

F = (QY - Y3)*ARRC[ii][jj-1][0][l] + (WY - QY)*ARRC[ii][jj-1][0][l+1]  ;

if( (jj-1) > 0  ){

F = F + ( (jj-1)/(2.0*Gam2))*( ARRC[ii][jj-2][0][l] - (OM/Gam2)*ARRC[ii][jj-2][0][l+1] ) ;

}

ARRC[ii][jj][0][l] = F;
}
}
}
}


if(   (L20 >= 1) && (  L22 >= 1 )  ) {
F=0.;
// 000 | n0m
for(int kk = 1; kk < L22 + 1  ; kk++){
  for(int ii = 1; ii < L20 + 1; ii++){
for (int l = 0; l < L  ; l++) {

F = (QZ -Z3)*ARRC[ii][0][kk-1][l] + (WZ - QZ)*ARRC[ii][0][kk-1][l+1]  ;

if( (kk-1) > 0  ){
F = F + ( (kk-1)/(2.0*Gam2))*( ARRC[ii][0][kk-2][l] - (OM/Gam2)*ARRC[ii][0][kk-2][l+1] ) ;

}
ARRC[ii][0][kk][l] = F;
}
}
}
}


if(   (L21 >= 1) && (  L22 >= 1 )  ) {
F=0;
// 000 | 0nm
for(int kk = 1; kk < L22 + 1  ; kk++){
  for(int jj = 1; jj < L21 + 1; jj++){
for (int l = 0; l < L  ; l++) {
F = (QZ - Z3)*ARRC[0][jj][kk-1][l] + (WZ - QZ)*ARRC[0][jj][kk-1][l+1]  ;

if( (kk-1) > 0  ){
F = F + ( (kk-1)/(2.0*Gam2))*( ARRC[0][jj][kk-2][l] - (OM/Gam2)*ARRC[0][jj][kk-2][l+1] ) ;
}

ARRC[0][jj][kk][l] = F;
}
}
}
}


if(   (L20 >= 1) && (  L21 >= 1 )  &&  (L22 >= 1 )  ) {
F=0;

for(int kk = 1; kk < L22 +1 ; kk++){
  for(int jj = 1; jj < L21 +1 ; jj++){
   for(int ii = 1; ii < L20 + 1 ; ii++){
for (int l = 0; l < L ; l++) {

F = (QZ - Z3)*ARRC[ii][jj][kk-1][l] + (WZ - QZ)*ARRC[ii][jj][kk-1][l+1]  ;

if( (kk-1) > 0  ){
F = F +  ( (kk-1)/(2.0*Gam2))*(ARRC[ii][jj][kk-2][l] - (OM/Gam2)*ARRC[ii][jj][kk-2][l+1] ) ;

}
ARRC[ii][jj][kk][l] = F;

}
}
}
}
}

AC = ARRC[L20][L21][L22][0];

for(int i = 0; i < (L20 + 1); i++){ 

  for(int j = 0; j < (L21 + 1); j++){   
    for(int n = 0; n < (L22 +1); n++)  {  

      delete[] ARRC[i][j][n];
        }
        delete [] ARRC[i][j];
      }
     delete [] ARRC[i];

    }

delete[] ARRC;
return AC;
}

// case a0mc0:
}else{

if ( (L20 == 0) &&  (L21 == 0) &&  (L22 == 0)  )  {

Double_v F = 0.;

Double_v ****ARR0;
ARR0 = new Double_v***[L10 + 1]; 
for(int i = 0; i < L10 + 1 ; i++){
ARR0[i] = new Double_v**[ L11 + 1 ];
 for(int j=0; j <  L11 + 1 ; j++) {
ARR0[i][j] = new Double_v*[L12 + 1];  
 for(int k=0; k <  L12 +1 ; k++) {
ARR0[i][j][k] = new Double_v[L+1];
}
}
}


MaskD_v mt = ( T > 100 );
if( MaskFull(mt)  ){

ARR0[0][0][0][0] = vecCore::math::Sqrt( pi/(4*T)  );

if (L > 0){
for (int i = 0; i < L ; i++) {
ARR0[0][0][0][i + 1] = (2.0*i+1)*ARR0[0][0][0][i]/(2.0*T);
}

}

}else{


ARR0[0][0][0][L] = fOpx->Gamma_aproxV( OM , RX, RY, RZ, 1, 1, 1,
1, 1, 1, L10, L11, L12, L20, L21, L22, gam1, Gam2, polyn, selector );

delete fOpx;

if (L > 0){
for(int i = L ; i > 0; i--) { 

ARR0[0][0][0][i - 1] = (2.0*T*ARR0[0][0][0][i] + vecCore::math::Exp(-T) )/( 2.0*(i-1)+1 );
}

}

} //else T > 100

for (int j = 0; j < L + 1 ; j++) {

ARR0[0][0][0][j] = ( ARR0[0][0][0][j]*Norms*pow(2,1)*pow(pi,5.0/2.0)*KP*KQ )/ 
( gam1*Gam2*( vecCore::math::Sqrt(G)  )  );
}

//ARR0[0][0][0] = OOOMOOO;
if (  L10 > 0   ){

for(int k = 1; k < L10 + 1 ; k++){
for (int l = 0; l < L  ; l++) {

F = (PX - X1)*ARR0[k-1][0][0][l] + (WX - PX)*ARR0[k-1][0][0][l+1] ; 

if( (k-1) > 0  ){
F = F + ( (k-1)/(2.0*gam1))*( ARR0[k-2][0][0][l] - (OM/gam1)*ARR0[k-2][0][0][l+1] ) ;
}

ARR0[k][0][0][l] = F;

}
}
}


if (  L11 > 0   ){

for(int k = 1; k < L11 + 1 ; k++){
for (int l = 0; l < L  ; l++) {
F = (PY - Y1)*ARR0[0][k-1][0][l] + (WY - PY)*ARR0[0][k-1][0][l+1]  ;

if( (k-1) > 0  ){
F = F + ( (k-1)/(2.0*gam1))*( ARR0[0][k-2][0][l] - (OM/gam1)*ARR0[0][k-2][0][l+1] ) ;
}

ARR0[0][k][0][l] = F;

}
}
}


if (  L12 > 0   ){

for(int k = 1; k < L12 + 1 ; k++){
for (int l = 0; l < L  ; l++) {

F = (PZ - Z1)*ARR0[0][0][k-1][l] + (WZ - PZ)*ARR0[0][0][k-1][l+1]  ;
if( (k-1) > 0  ){
F = F + ( (k-1)/(2.0*gam1))*( ARR0[0][0][k-2][l] - (OM/gam1)*ARR0[0][0][k-2][l+1] ) ;

}

ARR0[0][0][k][l] = F;

}
}
}

if(   (L10 >= 1) && (  L11  >= 1 )  ) {
F=0;

for(int jj = 1; jj < L11 + 1 ; jj++){
  for(int ii = 1; ii < L10 + 1 ; ii++){
for (int l = 0; l < L  ; l++) {

F = (PY - Y1)*ARR0[ii][jj-1][0][l] + (WY - PY)*ARR0[ii][jj-1][0][l+1]  ;

if( (jj-1) > 0  ){
F = F + ( (jj-1)/(2.0*gam1))*( ARR0[ii][jj-2][0][l] - (OM/gam1)*ARR0[ii][jj-2][0][l+1] ) ;

}

ARR0[ii][jj][0][l] = F;

}
}
}
}

if(   (L10 >= 1) && (  L12 >= 1 )  ) {
F=0;

for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
for (int l = 0; l < L  ; l++) {
//if positive sign f = (Pz - z1)*arr0[ii][0][kk-1][l] - (Wz - Pz)*arr0[ii][0][kk-1][l+1]  ;
F = (PZ - Z1)*ARR0[ii][0][kk-1][l] + (WZ - PZ)*ARR0[ii][0][kk-1][l+1]  ;

if( (kk-1) > 0  ){
F = F + ( (kk-1)/(2.0*gam1))*( 
ARR0[ii][0][kk-2][l] - (OM/gam1)*ARR0[ii][0][kk-2][l+1] ) ;

}
ARR0[ii][0][kk][l] = F;
}
}
}
}


if(   (L11 >= 1) && (  L12 >= 1 )  ) {
F=0;
// 0nm

for(int kk = 1; kk < L12 + 1  ; kk++){
  for(int jj = 1; jj < L11 + 1; jj++){
for (int l = 0; l < L  ; l++) {

F = (PZ - Z1)*ARR0[0][jj][kk-1][l] + (WZ - PZ)*ARR0[0][jj][kk-1][l+1]  ;

if( (kk-1) > 0  ){

F = F + ( (kk-1)/(2.0*gam1))*( ARR0[0][jj][kk-2][l] - (OM/gam1)*ARR0[0][jj][kk-2][l+1] ) ;

}
ARR0[0][jj][kk][l] = F;
}
}
}
}


if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {
F = 0.;

for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
    for (int l = 0; l < L ; l++) {

F = (PZ -Z1)*ARR0[ii][jj][kk-1][l] + (WZ - PZ)*ARR0[ii][jj][kk-1][l+1]  ;

if( (kk-1) > 0  ){
F = F + ( (kk-1)/(2.0*gam1))*(ARR0[ii][jj][kk-2][l] - (OM/gam1)*ARR0[ii][jj][kk-2][l+1] ) ;

}
ARR0[ii][jj][kk][l] = F;

}
}
}
}
}

AC = ARR0[L10][L11][L12][0];

for (int i = 0; i < L10 + 1; i++) {
    for (int j = 0; j < L11 + 1; j++) {
        for (int k = 0; k < L12 + 1; k++) {
            delete[] ARR0[i][j][k];
        }
        delete[] ARR0[i][j];
    }
    delete[] ARR0[i];
}
delete[] ARR0;

return AC;

}else{

if (  (L20 >= 1)  &&  (L21 == 0) &&  (L22 == 0)    ){

F = 0;
// FORMA DE APARTAR   [C1][A1][A2[A3]

Double_v *****ARR4;
ARR4 = new Double_v****[L20 + 1]; 
for(int i = 0; i < L20+1 ; i++){
ARR4[i] = new Double_v***[ L10 + 1];
for(int j = 0; j < L10+1 ; j++){
ARR4[i][j] = new Double_v**[ L11 + 1];
for(int k = 0; k <  L11 + 1 ; k++) {
ARR4[i][j][k] =  new Double_v*[ L12 + 1];
for(int l = 0; l <  L12 + 1 ; l++) {
ARR4[i][j][k][l] =  new Double_v[L + 1];
}
}
}
}

MaskD_v mt = ( T > 100 );
if( MaskFull(mt)  ){
ARR4[0][0][0][0][0] = vecCore::math::Sqrt( pi/(4*T)  );

if (L > 0){
for (int i = 0; i < L ; i++) {
ARR4[0][0][0][0][i + 1] = (2.0*i+1)*ARR4[0][0][0][0][i]/(2.0*T);
}

}

}else{


ARR4[0][0][0][0][L] = fOpx->Gamma_aproxV( OM , RX, RY, RZ, 1, 1, 1,
1, 1, 1, L10, L11, L12, L20, L21, L22, gam1, Gam2, polyn, selector );

delete fOpx;

if (L > 0){
for(int i = L ; i > 0; i--) { 
// ooomooo[i-1] =  (2.0*T*ooomooo[i] + G4Exp(-T) )/(2.0*(i-1)+1);
ARR4[0][0][0][0][i - 1] = (2.0*T*ARR4[0][0][0][0][i] + vecCore::math::Exp(-T) )/( 2.0*(i-1)+1 );
}

}

} //else T > 100

for (int j = 0; j < L + 1 ; j++) {

ARR4[0][0][0][0][j] = ( ARR4[0][0][0][0][j]*Norms*pow(2,1)*pow(pi,5.0/2.0)*KP*KQ )/ 
( gam1*Gam2*( vecCore::math::Sqrt(G)  )  );

}

if (  L10 > 0   ){
for(int i = 1; i < L10 + 1 ; i++){
 for (int l = 0; l < L  ; l++) {
// if positive sign f =  (Px - x1)*arr4[0][i-1][0][0][l] - (Wx - Px)*arr4[0][i-1][0][0][l+1] ; 
F = (PX - X1)*ARR4[0][i-1][0][0][l] + (WX - PX)*ARR4[0][i-1][0][0][l+1] ; 


if( (i-1) > 0  ){
F = F + ( (i-1)/(2.0*gam1))*( ARR4[0][i-2][0][0][l] - (OM/gam1)*ARR4[0][i-2][0][0][l+1] ) ;

}

ARR4[0][i][0][0][l] = F;

}
}
}

if (  L11 > 0   ){
for(int j = 1; j < L11 + 1 ; j++){
 for (int l = 0; l < L  ; l++) {
  //just los Q y W:
  // sign positiv f = (Py - y1)*arr4[0][0][j-1][0][l] - (Wy - Py)*arr4[0][0][j-1][0][l+1];
F = (PY - Y1)*ARR4[0][0][j-1][0][l] + (WY - PY)*ARR4[0][0][j-1][0][l+1];

if( (j-1) > 0  ){
F = F +  ( (j-1)/(2.0*gam1))*( ARR4[0][0][j-2][0][l] - (OM/gam1)*ARR4[0][0][j-2][0][l+1] ) ;

}

ARR4[0][0][j][0][l] = F;

}
}
}

if (  L12 > 0   ){
for(int k = 1; k < L12 + 1 ; k++){
 for (int l = 0; l < L  ; l++) {

F = (PZ - Z1)*ARR4[0][0][0][k-1][l] + (WZ - PZ)*ARR4[0][0][0][k-1][l+1]  ;

if( (k-1) > 0  ){
F = F + ( (k-1)/(2.0*gam1))*(ARR4[0][0][0][k-2][l] - (OM/gam1)*ARR4[0][0][0][k-2][l+1] ) ;

}
ARR4[0][0][0][k][l] = F;
}
}
}

if(   (L10 >= 1) && (  L11 >= 1 )  ) {
F=0;
for(int jj = 1; jj < L11 + 1 ; jj++){
  for(int ii = 1; ii < L10 + 1 ; ii++){
for (int l = 0; l < L  ; l++) {
F = (PY - Y1)*ARR4[0][ii][jj-1][0][l] + (WY - PY)*ARR4[0][ii][jj-1][0][l+1]  ;

if( (jj-1) > 0  ){
F = F + ( (jj-1)/(2.0*gam1))*( ARR4[0][ii][jj-2][0][l] - (OM/gam1)*ARR4[0][ii][jj-2][0][l+1] ) ;

}
ARR4[0][ii][jj][0][l] = F;
}
}
}
}

if(   (L10 >= 1) && (  L12 >= 1 )  ) {
F=0;
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
for (int l = 0; l < L  ; l++) {
F = (PZ -Z1)*ARR4[0][ii][0][kk-1][l] + (WZ - PZ)*ARR4[0][ii][0][kk-1][l+1]  ;

if( (kk-1) > 0  ){
F = F + ( (kk-1)/(2.0*gam1))*( 
ARR4[0][ii][0][kk-2][l] - (OM/gam1)*ARR4[0][ii][0][kk-2][l+1] ) ;

}
ARR4[0][ii][0][kk][l] = F;
}
}
}
}


if(   (L11 >= 1) && (  L12 >= 1 )  ) {
F=0;
// 0nm
for(int kk = 1; kk < L12 + 1  ; kk++){
  for(int jj = 1; jj < L11 + 1; jj++){
for (int l = 0; l < L  ; l++) {

F = (PZ -Z1)*ARR4[0][0][jj][kk-1][l] + (WZ - PZ)*ARR4[0][0][jj][kk-1][l+1]  ;

if( (kk-1) > 0  ){
F = F + ( (kk-1)/(2.0*gam1))*( ARR4[0][0][jj][kk-2][l] - (OM/gam1)*ARR4[0][0][jj][kk-2][l+1] ) ;

}
ARR4[0][0][jj][kk][l] = F;
}
}
}
}


if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {
F=0;

for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
for (int l = 0; l < L ; l++) {

F = (PZ - Z1)*ARR4[0][ii][jj][kk-1][l] + (WZ - PZ)*ARR4[0][ii][jj][kk-1][l+1]  ;
//f = (Pz - z1)*arr4[0][ii][jj][kk-1][l] + (Wz - Pz)*arr4[0][ii][jj][kk-1][l+1]  ;
if( (kk-1) > 0  ){
F = F + ( (kk-1)/(2.0*gam1))*(ARR4[0][ii][jj][kk-2][l] - (OM/gam1)*ARR4[0][ii][jj][kk-2][l+1] ) ;

}
ARR4[0][ii][jj][kk][l] = F;

}
}
}
}
}

// aqui L2[0] > 0
// [C1][A1][A2[A3]
// parte fourth tensor 000 | i00

for(int i = 1; i < L20 + 1 ; i++){
F = 0;
// case   000 | i00
 for (int l = 0; l < L  ; l++) {
F = (QX - X3)*ARR4[i-1][0][0][0][l] + (WX - QX)*ARR4[i-1][0][0][0][l+1] ;

if( (i-1) > 0  ){
F = F + ( (i-1)/(2.0*Gam2))*( ARR4[i-2][0][0][0][l] - (OM/Gam2)*ARR4[i-2][0][0][0][l+1] ) ;

}
ARR4[i][0][0][0][l] = F;
}

if (  L10 > 0   ){
for(int ii = 1; ii < L10 + 1 ; ii++){
 for (int l = 0; l < L  ; l++) {

F = (QX - X3)*ARR4[i-1][ii][0][0][l] + (WX - QX)*ARR4[i-1][ii][0][0][l+1] ;

if( (i-1) > 0  ){
F = F + ( (i-1)/(2.0*Gam2))*( ARR4[i-2][ii][0][0][l] - (OM/Gam2)*ARR4[i-2][ii][0][0][l+1] ) ;
}

F = F + ( (ii)/( 2.0*G ) )*ARR4[i-1][ii-1][0][0][l+1] ;
ARR4[i][ii][0][0][l] = F;
}
}
}

F = 0.;
if (  L11 > 0   ){
for(int jj = 1; jj < L11 + 1 ; jj++){
 for (int l = 0; l < L  ; l++) {

F = (QX - X3)*ARR4[i-1][0][jj][0][l] + (WX - QX)*ARR4[i-1][0][jj][0][l+1] ;

if( (i-1) > 0  ){
F = F + ( (i-1)/(2.0*Gam2))*( ARR4[i-2][0][jj][0][l] - (OM/Gam2)*ARR4[i-2][0][jj][0][l+1] ) ;
}

ARR4[i][0][jj][0][l] = F;

}
}
}

F = 0.;
if (  L12 > 0   ){
for(int kk = 1; kk < L12 + 1 ; kk++){
 for (int l = 0; l < L  ; l++) {

F = (QX - X3)*ARR4[i-1][0][0][kk][l] + (WX - QX)*ARR4[i-1][0][0][kk][l+1] ;

if( (i-1) > 0  ){
F = F + ( (i-1)/(2.0*Gam2))*( ARR4[i-2][0][0][kk][l] - (OM/Gam2)*ARR4[i-2][0][0][kk][l+1] ) ;

}

ARR4[i][0][0][kk][l] = F;

}
}
}


if(   (L10 >= 1) && (  L11 >= 1 )  ) {
F=0;

// i | nm0
for(int jj = 1; jj < L11 + 1 ; jj++){
//jj = 0 -> n00 -> previously, as well as 0m0
  for(int ii = 1; ii < L10 + 1 ; ii++){
    for (int l = 0; l < L  ; l++) { 
F = (QX - X3)*ARR4[i-1][ii][jj][0][l] + (WX - QX)*ARR4[i-1][ii][jj][0][l+1] ;

if( (i-1) > 0  ){
// almenos  i = 2:
F = F + ( (i-1)/(2.0*Gam2))*( ARR4[i-2][ii][jj][0][l] - (OM/Gam2)*ARR4[i-2][ii][jj][0][l+1] ) ;
}
F = F + ( (ii)/( 2.0*G ) )*ARR4[i-1][ii-1][jj][0][l+1] ;
ARR4[i][ii][jj][0][l] = F;
}
}
}
}

if(   (L10 >= 1) && (  L12 >= 1 )  ) {
F=0;
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
for (int l = 0; l < L  ; l++) {
//   n0m | i00
F = (QX - X3)*ARR4[i-1][ii][0][kk][l] + (WX - QX)*ARR4[i-1][ii][0][kk][l+1] ;

if( (i-1) > 0  ){
F = F + ( (i-1)/(2.0*Gam2))*( ARR4[i-2][ii][0][kk][l] - (OM/Gam2)*ARR4[i-2][ii][0][kk][l+1] ) ;
}
F = F +  ( (ii)/( 2.0*G ) )*ARR4[i-1][ii-1][0][kk][l+1] ;
ARR4[i][ii][0][kk][l] = F;
}
}
}
}

if(   (L11 >= 1) && (  L12 >= 1 )  ) {
F=0;
// 0nm | i00
for(int kk = 1; kk < L12 + 1  ; kk++){
  for(int jj = 1; jj < L11 + 1; jj++){
for (int l = 0; l < L  ; l++) {
F = (QX - X3)*ARR4[i-1][0][jj][kk][l] + (WX - QX)*ARR4[i-1][0][jj][kk][l+1]  ;

if( (i-1) > 0  ){
F = F + ( (i-1)/(2.0*Gam2))*( ARR4[i-2][0][jj][kk][l] - (OM/Gam2)*ARR4[i-2][0][jj][kk][l+1] ) ;

}
ARR4[i][0][jj][kk][l] = F;
}
}
}
}

if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {
F=0;
// i | nmk
for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
for (int l = 0; l < (L - i); l++) {
F = (QX - X3)*ARR4[i-1][ii][jj][kk][l] + (WX - QX)*ARR4[i-1][ii][jj][kk][l+1] ;

if( (i-1) > 0  ){
F = F + ( (i-1)/(2.0*Gam2))*(ARR4[i-2][ii][jj][kk][l] - (OM/Gam2)*ARR4[i-2][ii][jj][kk][l+1] ) ;

}
F = F +  ( (ii)/( 2.0*G ) )*ARR4[i-1][ii-1][jj][kk][l+1] ;
ARR4[i][ii][jj][kk][l] = F;

}
}
}
}
}

}// iterate c1 fin

// Check if the pointer is already freed (null)
AC = ARR4[L20][L10][L11][L12][0];

for (int i = 0; i < L20 + 1; i++) {
    for (int j = 0; j < L10 + 1; j++) {
        for (int k = 0; k < L11 + 1; k++) {
            for (int l = 0; l < L12 + 1; l++) {
                delete[] ARR4[i][j][k][l];
            }
            delete[] ARR4[i][j][k];
        }
        delete[] ARR4[i][j];
    }
    delete[] ARR4[i];
}
delete[] ARR4;

return AC;


}else{


if (    (L21 >= 1) &&  (L22 == 0)    ){
Double_v F = 0.;

// FORm   [C2][C1][A1][A2[A3]
Double_v  ******ARR5;

ARR5 = new Double_v*****[L21+1]; 
for(int i = 0; i < L21+1 ; i++){
ARR5[i] = new Double_v****[ L20 +1 ];
for(int j = 0; j < L20+1 ; j++){
ARR5[i][j] = new Double_v***[ L10 +1];
for(int k=0; k <  L10 +1 ; k++) {
ARR5[i][j][k] = new Double_v**[L11 +1 ];  
for(int l = 0; l <  L11 +1 ; l++) {
ARR5[i][j][k][l] = new Double_v*[L12+1];
for(int m=0; m <  L12 +1 ; m++) {
ARR5[i][j][k][l][m] =  new Double_v[L+1];
}
}
}
}
}

MaskD_v mt = ( T > 100 );
if( MaskFull(mt)  ){

ARR5[0][0][0][0][0][0] = vecCore::math::Sqrt( pi/(4*T)  );

if (L > 0){
for (int i = 0; i < L ; i++) {
ARR5[0][0][0][0][0][i + 1] = (2.0*i+1)*ARR5[0][0][0][0][0][i]/(2.0*T);
}

}

}else{

ARR5[0][0][0][0][0][L] = fOpx->Gamma_aproxV( OM , RX, RY, RZ, 1, 1, 1,
1, 1, 1, L10, L11, L12, L20, L21, L22, gam1, Gam2, polyn, selector );

delete fOpx;

if (L > 0){
for(int i = L ; i > 0; i--) { 
ARR5[0][0][0][0][0][i - 1] = (2.0*T*ARR5[0][0][0][0][0][i] + vecCore::math::Exp(-T) )/( 2.0*(i-1)+1 );
}

}

} //else T > 100

for (int j = 0; j < L + 1 ; j++) {

ARR5[0][0][0][0][0][j] = ( ARR5[0][0][0][0][0][j]*Norms*pow(2,1)*pow(pi,5.0/2.0)*KP*KQ )/ 
( gam1*Gam2*( vecCore::math::Sqrt(G)  )  );
}

//ARR5[0][0][0][0][0] = OOOMOOO;

if (  L10 > 0   ){
for(int i = 1; i < L10 + 1 ; i++){
 for (int l = 0; l < L  ; l++) {
F = (PX - X1)*ARR5[0][0][i-1][0][0][l] + (WX - PX)*ARR5[0][0][i-1][0][0][l+1] ;

if( (i-1) > 0  ){
F = F + ( (i-1)/(2.0*gam1) )*( ARR5[0][0][i-2][0][0][l] - (OM/gam1)*ARR5[0][0][i-2][0][0][l+1] ) ;

}
ARR5[0][0][i][0][0][l] = F;
}
}
}


if (  L11 > 0   ){
for(int j = 1; j < L11 + 1 ; j++){
 for (int l = 0; l < L  ; l++) {
F = (PY - Y1)*ARR5[0][0][0][j-1][0][l] + (WY - PY)*ARR5[0][0][0][j-1][0][l+1];

if( (j-1) > 0  ){
F = F + ( (j-1)/(2.0*gam1))*( ARR5[0][0][0][j-2][0][l] - (OM/gam1)*ARR5[0][0][0][j-2][0][l+1] ) ;

}

ARR5[0][0][0][j][0][l] = F;
}
}
}

if (  L12 > 0   ){
for(int k = 1; k < L12 + 1 ; k++){
 for (int l = 0; l < L  ; l++) {
F = (PZ - Z1)*ARR5[0][0][0][0][k-1][l] + (WZ - PZ)*ARR5[0][0][0][0][k-1][l+1]  ;

if( (k-1) > 0  ){
F = F + ( (k-1)/(2.0*gam1))*(ARR5[0][0][0][0][k-2][l] - (OM/gam1)*ARR5[0][0][0][0][k-2][l+1] ) ;

}
ARR5[0][0][0][0][k][l] = F;
}
}
}


if(   (L10 >= 1) && (  L11 >= 1 )  ) {
F=0;
//  ab0 | 000
for(int jj = 1; jj < L11 + 1 ; jj++){
  for(int ii = 1; ii < L10+ 1 ; ii++){
for (int l = 0; l < L  ; l++) {

F = (PY - Y1)*ARR5[0][0][ii][jj-1][0][l] + (WY - PY)*ARR5[0][0][ii][jj-1][0][l+1]  ;

if( (jj-1) > 0  ){
F = F + ( (jj-1)/(2.0*gam1))*( ARR5[0][0][ii][jj-2][0][l] - (OM/gam1)*ARR5[0][0][ii][jj-2][0][l+1] ) ;

}
ARR5[0][0][ii][jj][0][l] = F;
}
}
}
}


if(   (L10 >= 1) && (  L12 >= 1 )  ) {
F=0;
// a0b | 000
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
for (int l = 0; l < L  ; l++) {
F = (PZ - Z1)*ARR5[0][0][ii][0][kk-1][l] + (WZ - PZ)*ARR5[0][0][ii][0][kk-1][l+1]  ;

if( (kk-1) > 0  ){

F = F + ( (kk-1)/(2.0*gam1))*( 
ARR5[0][0][ii][0][kk-2][l] - (OM/gam1)*ARR5[0][0][ii][0][kk-2][l+1] ) ;

}
ARR5[0][0][ii][0][kk][l] = F;
}
}
}
}


if(   (L11 >= 1) && (  L12 >= 1 )  ) {
F=0;
// 0ab | 000
for(int kk = 1; kk < L12 + 1  ; kk++){
  for(int jj = 1; jj < L11 + 1; jj++){
for (int l = 0; l < L  ; l++) {
F = (PZ - Z1)*ARR5[0][0][0][jj][kk-1][l] + (WZ - PZ)*ARR5[0][0][0][jj][kk-1][l+1]  ;

if( (kk-1) > 0  ){

F = F + ( (kk-1)/(2.0*gam1))*( ARR5[0][0][0][jj][kk-2][l] - (OM/gam1)*ARR5[0][0][0][jj][kk-2][l+1] ) ;

}
ARR5[0][0][0][jj][kk][l] = F;
}
}
}
}

if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {

F=0;
for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
for (int l = 0; l < L ; l++) {

F = (PZ - Z1)*ARR5[0][0][ii][jj][kk-1][l] + (WZ - PZ)*ARR5[0][0][ii][jj][kk-1][l+1]  ;

if( (kk-1) > 0  ){
F = F + ( (kk-1)/(2.0*gam1))*( ARR5[0][0][ii][jj][kk-2][l] - (OM/gam1)*ARR5[0][0][ii][jj][kk-2][l+1] ) ;
}
ARR5[0][0][ii][jj][kk][l] = F;

}
}
}
}
}

// aqui C1[0] > 0
// [C2][C1][A1][A2[A3]
//inicial  000 | i00
if ( L20 > 0  ){
F = 0;
for(int i = 1; i < L20 + 1 ; i++){
 for (int l = 0; l < L  ; l++) {

F = (QX - X3)*ARR5[0][i-1][0][0][0][l] + (WX - QX)*ARR5[0][i-1][0][0][0][l+1] ;

if( (i-1) > 0  ){

F = F + ( (i-1)/(2.0*Gam2))*( ARR5[0][i-2][0][0][0][l] - (OM/Gam2)*ARR5[0][i-2][0][0][0][l+1] ) ;

}

ARR5[0][i][0][0][0][l] = F;
}
}

for(int i = 1; i < L20 + 1 ; i++){


if (  L10 > 0   ){
for(int ii = 1; ii < L10 + 1 ; ii++){
 for (int l = 0; l < L  ; l++) {
F = (QX - X3)*ARR5[0][i-1][ii][0][0][l] + (WX - QX)*ARR5[0][i-1][ii][0][0][l+1] ;

if( (i-1) > 0  ){

F = F + ( (i-1)/(2.0*Gam2))*( ARR5[0][i-2][ii][0][0][l] - (OM/Gam2)*ARR5[0][i-2][ii][0][0][l+1] ) ;

}

F = F + ( (ii)/( 2.0*G ) )*ARR5[0][i-1][ii-1][0][0][l+1] ;

ARR5[0][i][ii][0][0][l] = F;
}
}
}

F = 0.;
if (  L11 > 0   ){
for(int jj = 1; jj < L11 + 1 ; jj++){
 for (int l = 0; l < L  ; l++) {
F = (QX - X3)*ARR5[0][i-1][0][jj][0][l] + (WX - QX)*ARR5[0][i-1][0][jj][0][l+1] ;

if( (i-1) > 0  ){
F = F + ( (i-1)/(2.0*Gam2))*( ARR5[0][i-2][0][jj][0][l] - (OM/Gam2)*ARR5[0][i-2][0][jj][0][l+1] ) ;

}
ARR5[0][i][0][jj][0][l] = F;
}
}
}

F = 0.;
if (  L12 > 0   ){
for(int kk = 1; kk < L12 + 1 ; kk++){
 for (int l = 0; l < L  ; l++) {

F = (QX - X3)*ARR5[0][i-1][0][0][kk][l] + (WX - QX)*ARR5[0][i-1][0][0][kk][l+1] ;

if( (i-1) > 0  ){

F = F + ( (i-1)/(2.0*Gam2))*( ARR5[0][i-2][0][0][kk][l] - (OM/Gam2)*ARR5[0][i-2][0][0][kk][l+1] ) ;

}
ARR5[0][i][0][0][kk][l] = F;

}
}
}


if(   (L10 >= 1) && (  L11 >= 1 )  ) {
F=0;
// ab0 | i00
for(int jj = 1; jj < L11 + 1 ; jj++){
  for(int ii = 1; ii < L10 + 1 ; ii++){
    for (int l = 0; l < L  ; l++) { 
F = (QX - X3)*ARR5[0][i-1][ii][jj][0][l] + (WX - QX)*ARR5[0][i-1][ii][jj][0][l+1] ;

if( (i-1) > 0  ){
// almenos  i = 2:
F = F + ( (i-1)/(2.0*Gam2))*( ARR5[0][i-2][ii][jj][0][l] - (OM/Gam2)*ARR5[0][i-2][ii][jj][0][l+1] ) ;

}
F = F + ( (ii)/( 2.0*G ) )*ARR5[0][i-1][ii-1][jj][0][l+1] ;
ARR5[0][i][ii][jj][0][l] = F;
}
}
}
}


if(   (L10 >= 1) && (  L12 >= 1 )  ) {
F=0;  //  a0b | i00
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
for (int l = 0; l < L  ; l++) {
F = (QX - X3)*ARR5[0][i-1][ii][0][kk][l] + (WX - QX)*ARR5[0][i-1][ii][0][kk][l+1]  ;

if( (i-1) > 0  ){
F = F + ( (i-1)/(2.0*Gam2))*( ARR5[0][i-2][ii][0][kk][l] - (OM/Gam2)*ARR5[0][i-2][ii][0][kk][l+1] ) ;

}
F = F + ( (ii)/( 2.0*G ) )*ARR5[0][i-1][ii-1][0][kk][l+1] ;
ARR5[0][i][ii][0][kk][l] = F;
}
}
}
}


if(   (L11 >= 1) && (  L12 >= 1 )  ) {
F=0;  //  0ab | i00 0 nab | i00
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int jj = 1; jj < L11 +1  ; jj++){
for (int l = 0; l < L  ; l++) {
F = (QX - X3)*ARR5[0][i-1][0][jj][kk][l] + (WX - QX)*ARR5[0][i-1][0][jj][kk][l+1]  ;

if( (i-1) > 0  ){
F = F + ( (i-1)/(2.0*Gam2))*( ARR5[0][i-2][0][jj][kk][l] - (OM/Gam2)*ARR5[0][i-2][0][jj][kk][l+1] ) ;

}
ARR5[0][i][0][jj][kk][l] = F;
}
}
}
}

// abc | i00
if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {
F=0;
for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
for (int l = 0; l < L ; l++) {
F = (QX - X3)*ARR5[0][i-1][ii][jj][kk][l] + (WX - QX)*ARR5[0][i-1][ii][jj][kk][l+1] ;

if( (i-1) > 0  ){
F = F + ( (i-1)/(2.0*Gam2))*(ARR5[0][i-2][ii][jj][kk][l] - (OM/Gam2)*ARR5[0][i-2][ii][jj][kk][l+1] ) ;

}
F = F + ( (ii)/( 2.0*G ) )*ARR5[0][i-1][ii-1][jj][kk][l+1] ;

ARR5[0][i][ii][jj][kk][l] = F;
}
}
}
}
}

// ciera c1 mayor a 0
}
// cierre L20 >0
}

// abc | 0j0
// [C2][C1][A1][A2][A3]
for(int j = 1; j < L21 + 1 ; j++){
F = 0;
 for (int l = 0; l < L  ; l++) {
  // el inicial 000 | 0j0
F = (QY - Y3)*ARR5[j-1][0][0][0][0][l] + (WY - QY)*ARR5[j-1][0][0][0][0][l+1] ;

if( (j-1) > 0  ){
F = F + ( (j-1)/(2.0*Gam2))*( ARR5[j-2][0][0][0][0][l] - (OM/Gam2)*ARR5[j-2][0][0][0][0][l+1] ) ;

}
ARR5[j][0][0][0][0][l] = F;
}

// a00 | 0j
if(  L10 > 0  ){  
F = 0;
for(int ii = 1; ii < L10 + 1 ; ii++){
  for (int l = 0; l < L  ; l++) {
F = (QY - Y3)*ARR5[j-1][0][ii][0][0][l] + (WY - QY)*ARR5[j-1][0][ii][0][0][l+1] ;

if( (j-1) > 0  ){
F = F + ( (j-1)/(2.0*Gam2))*( ARR5[j-2][0][ii][0][0][l] - (OM/Gam2)*ARR5[j-2][0][ii][0][0][l+1] ) ;
}
ARR5[j][0][ii][0][0][l] = F;
}
}
}

F = 0;
if (  L11 > 0   ){
for(int jj = 1; jj < L11 + 1 ; jj++){
 for (int l = 0; l < L  ; l++) {
F = (QY - Y3)*ARR5[j-1][0][0][jj][0][l] + (WY - QY)*ARR5[j-1][0][0][jj][0][l+1] ;

if( (j-1) > 0  ){
F = F + ( (j-1)/(2.0*Gam2))*( ARR5[j-2][0][0][jj][0][l] - (OM/Gam2)*ARR5[j-2][0][0][jj][0][l+1] ) ;

}
F = F + ( (jj)/( 2.0*G ) )*ARR5[j-1][0][0][jj-1][0][l+1] ;

ARR5[j][0][0][jj][0][l] = F;
}
}
}

F = 0;
if (  L12 > 0   ){
for(int kk = 1; kk < L12 + 1 ; kk++){
 for (int l = 0; l < L  ; l++) {
F = (QY - Y3)*ARR5[j-1][0][0][0][kk][l] + (WY - QY)*ARR5[j-1][0][0][0][kk][l+1] ;

if( (j-1) > 0  ){
F = F + ( (j-1)/(2.0*Gam2))*( ARR5[j-2][0][0][0][kk][l] - (OM/Gam2)*ARR5[j-2][0][0][0][kk][l+1] ) ;

}
ARR5[j][0][0][0][kk][l] = F;
}
}
}

if(   (L10 >= 1) && (  L11 >= 1 )  ) {
F=0;
// ab0 | 0j0
for(int jj = 1; jj < L11 + 1 ; jj++){
  for(int ii = 1; ii < L10 + 1 ; ii++){ 
   for (int l = 0; l < (L - L10  ); l++) { 
F = (QY - Y3)*ARR5[j-1][0][ii][jj][0][l] + (WY - QY)*ARR5[j-1][0][ii][jj][0][l+1] ;

if( (j-1) > 0  ){
F = F + ( (j-1)/(2.0*Gam2))*( ARR5[j-2][0][ii][jj][0][l] - (OM/Gam2)*ARR5[j-2][0][ii][jj][0][l+1] ) ;

}
F = F + ( (jj)/( 2.0*G) )*ARR5[j-1][0][ii][jj-1][0][l+1] ;
ARR5[j][0][ii][jj][0][l] = F;
}
}
}
}

if(   (L10 >= 1) && (  L12 >= 1 )  ) {
F=0;  
//  a0b | 0j0
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
for (int l = 0; l <  (L - L10  )  ; l++) {
F = (QY - Y3)*ARR5[j-1][0][ii][0][kk][l] + (WY - QY)*ARR5[j-1][0][ii][0][kk][l+1]  ;

if( (j-1) > 0  ){
F = F + ( (j-1)/(2.0*Gam2))*( ARR5[j-2][0][ii][0][kk][l] - (OM/Gam2)*ARR5[j-2][0][ii][0][kk][l+1] ) ;

}
ARR5[j][0][ii][0][kk][l] = F;
}
}
}
}


if(   (L11 >= 1) && (  L12 >= 1 )  ) {
F=0;  //  0ab | 0j0
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int jj = 1; jj < L11 +1  ; jj++){
for (int l = 0; l < (L -L11)  ; l++) {
F = (QY - Y3)*ARR5[j-1][0][0][jj][kk][l] + (WY - QY)*ARR5[j-1][0][0][jj][kk][l+1]  ;

if( (j-1) > 0  ){
F = F + ( (j-1)/(2.0*Gam2))*( ARR5[j-2][0][0][jj][kk][l] - (OM/Gam2)*ARR5[j-2][0][0][jj][kk][l+1] ) ;

}
F = F + ( (jj)/( 2.0*G ) )*ARR5[j-1][0][0][jj-1][kk][l+1] ;
ARR5[j][0][0][jj][kk][l] = F;
}
}
}
}

// abc | 0j0
if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {
F=0;
for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
     for (int l = 0; l < L ; l++) {
      F = (QY - Y3)*ARR5[j-1][0][ii][jj][kk][l] + (WY - QY)*ARR5[j-1][0][ii][jj][kk][l+1] ;

if( (j-1) > 0  ){

F = F + ( (j-1)/(2.0*Gam2))*(ARR5[j-2][0][ii][jj][kk][l] - (OM/Gam2)*ARR5[j-2][0][ii][jj][kk][l+1] ) ;

}
F = F + ( (jj)/( 2.0*G) )*ARR5[j-1][0][ii][jj-1][kk][l+1] ;

ARR5[j][0][ii][jj][kk][l] = F;
}
}
}
}
}


}//

// abc | mn0
// arr5[j][i][n][m][k][l] 
// if  C1 > 0:
if( L20 > 0 )  {
for(int j = 1; j < L21 + 1 ; j++){
for(int i = 1; i < L20 + 1 ; i++){

// inicial 000 | ij0:
for (int l = 0; l < L ; l++) {
F = (QY - Y3)*ARR5[j-1][i][0][0][0][l] + (WY - QY)*ARR5[j-1][i][0][0][0][l+1] ;

if( (j-1) > 0  ){
F = F + ( (j-1)/(2.0*Gam2) )*(ARR5[j-2][i][0][0][0][l] - (OM/Gam2)*ARR5[j-2][i][0][0][0][l+1] ) ;

}
ARR5[j][i][0][0][0][l] = F;
}

// a00 | ij0
if(  L10 > 0  ){  
F= 0;
for(int ii = 1; ii < L10 + 1 ; ii++){
  for (int l = 0; l < L  ; l++) {
F = (QY - Y3)*ARR5[j-1][i][ii][0][0][l] + (WY - QY)*ARR5[j-1][i][ii][0][0][l+1] ;

if( (j-1) > 0  ){
F = F + ( (j-1)/(2.0*Gam2))*( ARR5[j-2][i][ii][0][0][l] - (OM/Gam2)*ARR5[j-2][i][ii][0][0][l+1] ) ;

}
ARR5[j][i][ii][0][0][l] = F;
}
}
}

// 0b0 | ij0
if (  L11 > 0   ){
F = 0;
for(int jj = 1; jj < L11 + 1 ; jj++){
 for (int l = 0; l < L  ; l++) {

F = (QY - Y3)*ARR5[j-1][i][0][jj][0][l] + (WY - QY)*ARR5[j-1][i][0][jj][0][l+1] ;
if( (j-1) > 0  ){
F = F + ( (j-1)/(2.0*Gam2))*( ARR5[j-2][i][0][jj][0][l] - (OM/Gam2)*ARR5[j-2][i][0][jj][0][l+1] ) ;

}
F = F + ( (jj)/( 2.0*G ) )*ARR5[j-1][i][0][jj-1][0][l+1] ;
ARR5[j][i][0][jj][0][l] = F;
}
}
}

// 00c | ij0
if (  L12 > 0   ){
F = 0;
for(int kk = 1; kk < L12 + 1 ; kk++){
 for (int l = 0; l < L  ; l++) {
    F = (QY - Y3)*ARR5[j-1][i][0][0][kk][l] + (WY - QY)*ARR5[j-1][i][0][0][kk][l+1] ;

if( (j-1) > 0  ){
F = F + ( (j-1)/(2.0*Gam2))*( ARR5[j-2][i][0][0][kk][l] - (OM/Gam2)*ARR5[j-2][i][0][0][kk][l+1] ) ;

}
ARR5[j][i][0][0][kk][l] = F;
}
}
}

// ab0 | ij0
if(   (L10 >= 1) && (  L11 >= 1 )  ) {
F=0;
for(int jj = 1; jj < L11 + 1 ; jj++){
  for(int ii = 1; ii < L10 + 1 ; ii++){
    for (int l = 0; l < L  ; l++) { 

F = (QY - Y3)*ARR5[j-1][i][ii][jj][0][l] + (WY - QY)*ARR5[j-1][i][ii][jj][0][l+1] ;
if( (j-1) > 0  ){

F = F + ( (j-1)/(2.0*Gam2))*( ARR5[j-2][i][ii][jj][0][l] - (OM/Gam2)*ARR5[j-2][i][ii][jj][0][l+1] ) ;

}

F = F + ( (jj)/( 2.0*G ) )*ARR5[j-1][i][ii][jj-1][0][l+1] ;
ARR5[j][i][ii][jj][0][l] = F;
}
}
}
}

if(   (L10 >= 1) && (  L12 >= 1 )  ) {
F=0;  
//  a0b | ij0
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
for (int l = 0; l < L  ; l++) {
F = (QY - Y3)*ARR5[j-1][i][ii][0][kk][l] + (WY - QY)*ARR5[j-1][i][ii][0][kk][l+1];

if( (j-1) > 0  ){
F = F + ( (j-1)/(2.0*Gam2))*( ARR5[j-2][i][ii][0][kk][l] - (OM/Gam2)*ARR5[j-2][i][ii][0][kk][l+1] ) ;

}
ARR5[j][i][ii][0][kk][l] = F;
}
}
}
}


if(   (L11 >= 1) && (  L12 >= 1 )  ) {
F=0;  //  0ab | ij0
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int jj = 1; jj < L11 +1  ; jj++){
for (int l = 0; l < L  ; l++) {
F = (QY - Y3)*ARR5[j-1][i][0][jj][kk][l] + (WY - QY)*ARR5[j-1][i][0][jj][kk][l+1]  ;

if( (j-1) > 0  ){
F = F + ( (j-1)/(2.0*Gam2))*( ARR5[j-2][i][0][jj][kk][l] - (OM/Gam2)*ARR5[j-2][i][0][jj][kk][l+1] ) ;

}
F = F + ( (jj)/( 2.0*G ) )*ARR5[j-1][i][0][jj-1][kk][l+1] ;

ARR5[j][i][0][jj][kk][l] = F;
}
}
}
}

// abc | ij0
if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {
F=0;
for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
    for (int l = 0; l < L ; l++) {
     F = (QY - Y3)*ARR5[j-1][i][ii][jj][kk][l] + (WY - QY)*ARR5[j-1][i][ii][jj][kk][l+1] ;

if( (j-1) > 0  ){
F = F + ( (j-1)/(2.0*Gam2))*(ARR5[j-2][i][ii][jj][kk][l] - (OM/Gam2)*ARR5[j-2][i][ii][jj][kk][l+1] ) ;

}
F = F + ( (jj)/( 2.0*G ) )*ARR5[j-1][i][ii][jj-1][kk][l+1] ;

ARR5[j][i][ii][jj][kk][l] = F;
}
}
}
}
}


}// i
}// j

}//fin if  L20 > 0

AC = ARR5[L21][L20][L10][L11][L12][0];  

for (int i = 0; i < L21 + 1; i++) {
    for (int j = 0; j < L20 + 1; j++) {
        for (int k = 0; k < L10 + 1; k++) {
            for (int l = 0; l < L11 + 1; l++) {
                for (int m = 0; m < L12 + 1; m++) {
                    delete[] ARR5[i][j][k][l][m];
                }
                delete[] ARR5[i][j][k][l];
            }
            delete[] ARR5[i][j][k];
        }
        delete[] ARR5[i][j];
    }
    delete[] ARR5[i];
}
delete[] ARR5;

return AC;

}else{

// FORM :  [C3[C2][C1][A1][A2[A3]
F = 0.;

Double_v*******ARR6;
ARR6 = new Double_v******[L22 + 1]; 
for(int i = 0; i < L22 + 1 ; i++){
ARR6[i] = new Double_v*****[ L21 + 1];
for(int j = 0; j < L21 + 1 ; j++){
ARR6[i][j] = new Double_v****[ L20 + 1];
for(int k = 0; k <  L20 + 1 ; k++) {
ARR6[i][j][k] =  new Double_v***[ L10 + 1];
for(int l = 0; l <  L10 + 1 ; l++) {
ARR6[i][j][k][l] = new Double_v**[L11 + 1];
for(int m = 0; m <  L11 + 1 ; m++) {
ARR6[i][j][k][l][m] =  new Double_v*[L12 + 1];
for(int n = 0; n <  L12 + 1 ; n++) {
ARR6[i][j][k][l][m][n] =  new Double_v[L + 1];
}
}
}
}
}
}


MaskD_v mt = ( T > 100 );
if( MaskFull(mt)  ){

ARR6[0][0][0][0][0][0][0] = vecCore::math::Sqrt( pi/(4*T)  );

if (L > 0){
for (int i = 0; i < L ; i++) {
ARR6[0][0][0][0][0][0][i + 1] = (2.0*i+1)*ARR6[0][0][0][0][0][0][i]/(2.0*T);
}

}

}else{

ARR6[0][0][0][0][0][0][L] = fOpx->Gamma_aproxV( OM , RX, RY, RZ, 1, 1, 1,
1, 1, 1, L10, L11, L12, L20, L21, L22, gam1, Gam2, polyn, selector );

delete fOpx;

if (L > 0){
for(int i = L ; i > 0; i--) { 
ARR6[0][0][0][0][0][0][i - 1] = (2.0*T*ARR6[0][0][0][0][0][0][i] + vecCore::math::Exp(-T) )/( 2.0*(i-1)+1 );
}

}

} //else T > 100

for (int j = 0; j < L + 1 ; j++) {
ARR6[0][0][0][0][0][0][j] = ( ARR6[0][0][0][0][0][0][j]*Norms*pow(2,1)*pow(pi,5.0/2.0)*KP*KQ )/ 
( gam1*Gam2*( vecCore::math::Sqrt(G)  )  );
}


if (  L10 > 0   ){
for(int i = 1; i < L10 + 1 ; i++){
 for (int l = 0; l < L  ; l++) {

F = (PX - X1)*ARR6[0][0][0][i-1][0][0][l] + (WX - PX)*ARR6[0][0][0][i-1][0][0][l+1] ; 

if( (i-1) > 0  ){
F = F + ( (i-1)/(2.0*gam1) )*( ARR6[0][0][0][i-2][0][0][l] - (OM/gam1)*ARR6[0][0][0][i-2][0][0][l+1] ) ;

}
ARR6[0][0][0][i][0][0][l] = F;
}
}
}

F = 0.;
if (  L11 > 0   ){
for(int j = 1; j < L11 + 1 ; j++){
 for (int l = 0; l < L  ; l++) {
  F = (PY - Y1)*ARR6[0][0][0][0][j-1][0][l] + (WY - PY)*ARR6[0][0][0][0][j-1][0][l+1];

if( (j-1) > 0  ){
 F = F + ( (j-1)/(2.0*gam1))*( ARR6[0][0][0][0][j-2][0][l] - (OM/gam1)*ARR6[0][0][0][0][j-2][0][l+1] ) ;

}
ARR6[0][0][0][0][j][0][l] = F;
}
}
}


if (  L12 > 0   ){
for(int k = 1; k < L12 + 1 ; k++){
 for (int l = 0; l < L  ; l++) {

F = (PZ - Z1)*ARR6[0][0][0][0][0][k-1][l] + (WZ - PZ)*ARR6[0][0][0][0][0][k-1][l+1]  ;

if( (k-1) > 0  ){
F = F + ( (k-1)/(2.0*gam1))*(ARR6[0][0][0][0][0][k-2][l] - (OM/gam1)*ARR6[0][0][0][0][0][k-2][l+1] ) ;

}
ARR6[0][0][0][0][0][k][l] = F;
}
}
}

if(   (L10 >= 1) && (  L11 >= 1 )  ) {
F=0;
//  ab0 | 000
for(int jj = 1; jj < L11 + 1 ; jj++){
  for(int ii = 1; ii < L10 + 1 ; ii++){
for (int l = 0; l < L  ; l++) {

F = (PY - Y1)*ARR6[0][0][0][ii][jj-1][0][l] + (WY - PY)*ARR6[0][0][0][ii][jj-1][0][l+1]  ;

if( (jj-1) > 0  ){
F = F + ( (jj-1)/(2.0*gam1))*( ARR6[0][0][0][ii][jj-2][0][l] - (OM/gam1)*ARR6[0][0][0][ii][jj-2][0][l+1] ) ;

}
ARR6[0][0][0][ii][jj][0][l] = F;
}
}
}
}


if(   (L10 >= 1) && (  L12 >= 1 )  ) {
F=0;
// a0b | 000
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
    for (int l = 0; l < L  ; l++) {
F = (PZ - Z1)*ARR6[0][0][0][ii][0][kk-1][l] + (WZ - PZ)*ARR6[0][0][0][ii][0][kk-1][l+1]  ;

if( (kk-1) > 0  ){
F = F + ( (kk-1)/(2.0*gam1))*( 
ARR6[0][0][0][ii][0][kk-2][l] - (OM/gam1)*ARR6[0][0][0][ii][0][kk-2][l+1] ) ;

}
ARR6[0][0][0][ii][0][kk][l] = F;
}
}
}
}


if(   (L11 >= 1) && (  L12 >= 1 )  ) {
F=0;
// 0ab | 000
for(int kk = 1; kk < L12 + 1  ; kk++){
  for(int jj = 1; jj < L11 + 1; jj++){
for (int l = 0; l < L  ; l++) {
    F = (PZ - Z1)*ARR6[0][0][0][0][jj][kk-1][l] + (WZ - PZ)*ARR6[0][0][0][0][jj][kk-1][l+1]  ;
if( (kk-1) > 0  ){
F = F + ( (kk-1)/(2.0*gam1))*( ARR6[0][0][0][0][jj][kk-2][l] - (OM/gam1)*ARR6[0][0][0][0][jj][kk-2][l+1] ) ;

}
ARR6[0][0][0][0][jj][kk][l] = F;
}
}
}
}


if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {
F=0;
for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
     for (int l = 0; l < L ; l++) {
        F = (PZ - Z1)*ARR6[0][0][0][ii][jj][kk-1][l] + (WZ - PZ)*ARR6[0][0][0][ii][jj][kk-1][l+1]  ;

if( (kk-1) > 0  ){
    F = F + ( (kk-1)/(2.0*gam1))*( ARR6[0][0][0][ii][jj][kk-2][l] - (OM/gam1)*ARR6[0][0][0][ii][jj][kk-2][l+1] ) ;

}
ARR6[0][0][0][ii][jj][kk][l] = F;
}
}
}
}
}

if( L20 >  0   ){
//  C1[0] > 0
// [C3[C2][C1][A1][A2[A3]
F = 0.;
//  abc | i00
for(int i = 1; i < L20 + 1 ; i++){
F = 0.;
 for (int l = 0; l < (L -L20)  ; l++) {
    F = (QX - X3)*ARR6[0][0][i-1][0][0][0][l] + (WX - QX)*ARR6[0][0][i-1][0][0][0][l+1] ;

if( (i-1) > 0  ){
F = F + ( (i-1)/(2.0*Gam2))*( ARR6[0][0][i-2][0][0][0][l] - (OM/Gam2)*ARR6[0][0][i-2][0][0][0][l+1] ) ;

}
ARR6[0][0][i][0][0][0][l] = F;
}

if (  L10 > 0   ){
for(int ii = 1; ii < L10 + 1 ; ii++){
 for (int l = 0; l < L  ; l++) {
F = (QX - X3)*ARR6[0][0][i-1][ii][0][0][l] + (WX - QX)*ARR6[0][0][i-1][ii][0][0][l+1] ;

if( (i-1) > 0  ){
F = F + ( (i-1)/(2.0*Gam2))*( ARR6[0][0][i-2][ii][0][0][l] - (OM/Gam2)*ARR6[0][0][i-2][ii][0][0][l+1] ) ;

}
F = F + ( (ii)/( 2.0*G ) )*ARR6[0][0][i-1][ii-1][0][0][l+1] ;

ARR6[0][0][i][ii][0][0][l] = F;
}
}
}


if (  L11 > 0   ){
F = 0;
// 0b0 | i00
for(int jj = 1; jj < L11 + 1 ; jj++){
 for (int l = 0; l < L  ; l++) {
    F = (QX - X3)*ARR6[0][0][i-1][0][jj][0][l] + (WX - QX)*ARR6[0][0][i-1][0][jj][0][l+1] ;

if( (i-1) > 0  ){
F = F + ( (i-1)/(2.0*Gam2))*( ARR6[0][0][i-2][0][jj][0][l] - (OM/Gam2)*ARR6[0][0][i-2][0][jj][0][l+1] ) ;

}
ARR6[0][0][i][0][jj][0][l] = F;
}
}
}

// 00c | i00
if (  L12 > 0   ){
for(int kk = 1; kk < L12 + 1 ; kk++){
F = 0;
 for (int l = 0; l < L  ; l++) {

F = (QX - X3)*ARR6[0][0][i-1][0][0][kk][l] + (WX - QX)*ARR6[0][0][i-1][0][0][kk][l+1] ;

if( (i-1) > 0  ){
F = F + ( (i-1)/(2.0*Gam2))*( ARR6[0][0][i-2][0][0][kk][l] - (OM/Gam2)*ARR6[0][0][i-2][0][0][kk][l+1] ) ;

}
ARR6[0][0][i][0][0][kk][l] = F;
}
}
}

if(   (L10 >= 1) && (  L11 >= 1 )  ) {
F=0;
// ab0 | i00
for(int jj = 1; jj < L11 + 1 ; jj++){
  for(int ii = 1; ii < L10 + 1 ; ii++){
    for (int l = 0; l < L  ; l++) { 
     F = (QX - X3)*ARR6[0][0][i-1][ii][jj][0][l] + (WX - QX)*ARR6[0][0][i-1][ii][jj][0][l+1] ;

if( (i-1) > 0  ){
// at least  i = 2:
    F = F + ( (i-1)/(2.0*Gam2))*( ARR6[0][0][i-2][ii][jj][0][l] - (OM/Gam2)*ARR6[0][0][i-2][ii][jj][0][l+1] ) ;

}
F = F + ( (ii)/( 2.0*G ) )*ARR6[0][0][i-1][ii-1][jj][0][l+1] ;
ARR6[0][0][i][ii][jj][0][l] = F;
}
}
}
}

if(   (L10 >= 1) && (  L12 >= 1 )  ) {
F=0;  //  a0b | i00
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
for (int l = 0; l < L  ; l++) {
   
    F = (QX - X3)*ARR6[0][0][i-1][ii][0][kk][l] + (WX - QX)*ARR6[0][0][i-1][ii][0][kk][l+1]  ;

if( (i-1) > 0  ){
F = F + ( (i-1)/(2.0*Gam2))*( ARR6[0][0][i-2][ii][0][kk][l] - (OM/Gam2)*ARR6[0][0][i-2][ii][0][kk][l+1] ) ;

}
F = F + ( (ii)/( 2.0*G ) )*ARR6[0][0][i-1][ii-1][0][kk][l+1] ;

ARR6[0][0][i][ii][0][kk][l] = F;
}
}
}
}


if(   (L11 >= 1) && (  L12 >= 1 )  ) {
F=0;  //  0ab | i00 0 nab | i00
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int jj = 1; jj < L11 +1  ; jj++){
    for (int l = 0; l < L  ; l++) {

        F = (QX - X3)*ARR6[0][0][i-1][0][jj][kk][l] + (WX - QX)*ARR6[0][0][i-1][0][jj][kk][l+1]  ;

if( (i-1) > 0  ){
F = F + ( (i-1)/(2.0*Gam2))*( ARR6[0][0][i-2][0][jj][kk][l] - (OM/Gam2)*ARR6[0][0][i-2][0][jj][kk][l+1] ) ;

}
ARR6[0][0][i][0][jj][kk][l] = F;
}
}
}
}

// abc | i00
if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {
F=0;
for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
    for (int l = 0; l < L ; l++) {
     F = (QX - X3)*ARR6[0][0][i-1][ii][jj][kk][l] + (WX - QX)*ARR6[0][0][i-1][ii][jj][kk][l+1] ;

if( (i-1) > 0  ){
    F = F + ( (i-1)/(2.0*Gam2))*(ARR6[0][0][i-2][ii][jj][kk][l] - (OM/Gam2)*ARR6[0][0][i-2][ii][jj][kk][l+1] ) ;

}
F = F + ( (ii)/( 2.0*G ) )*ARR6[0][0][i-1][ii-1][jj][kk][l+1] ;

ARR6[0][0][i][ii][jj][kk][l] = F;
}
}
}
}
}


} // i

} // L20

if( L21 >  0   ){
F = 0;
// arr6[0][j][0][0][0][0][l] 
for(int j = 1; j < L21 + 1 ; j++){
// init 000 | 0j0
 for (int l = 0; l < (L - L21)  ; l++) {
F = (QY - Y3)*ARR6[0][j-1][0][0][0][0][l] + (WY - QY)*ARR6[0][j-1][0][0][0][0][l+1] ;

if( (j-1) > 0  ){
    F = F + ( (j-1)/(2.0*Gam2))*( ARR6[0][j-2][0][0][0][0][l] - (OM/Gam2)*ARR6[0][j-2][0][0][0][0][l+1] ) ;

}
ARR6[0][j][0][0][0][0][l] = F;
}

// inicial 000 | ij0:
if (  L20 > 0   ){
  for(int i = 1; i < L20 + 1 ; i++){
   for (int l = 0; l < L  ; l++) {
    F = (QY - Y3)*ARR6[0][j-1][i][0][0][0][l] + (WY - QY)*ARR6[0][j-1][i][0][0][0][l+1] ;

    if( (j-1) > 0  ){
     F = F + ( (j-1)/(2.0*Gam2))*(ARR6[0][j-2][i][0][0][0][l] - (OM/Gam2)*ARR6[0][j-2][i][0][0][0][l+1] ) ;

}
ARR6[0][j][i][0][0][0][l] = F;
}
}
}

// abc | 0j0
// [C3][C2][C1][A1][A2][A3]
// a00 | 0j0
if(  L10 > 0  ){  
F = 0.;
for(int ii = 1; ii < L10 + 1 ; ii++){
  for (int l = 0; l < L  ; l++) {
    F = (QY - Y3)*ARR6[0][j-1][0][ii][0][0][l] + (WY - QY)*ARR6[0][j-1][0][ii][0][0][l+1] ;

if( (j-1) > 0  ){
F = F + ( (j-1)/(2.0*Gam2))*( ARR6[0][j-2][0][ii][0][0][l] - (OM/Gam2)*ARR6[0][j-2][0][ii][0][0][l+1] ) ;

}
ARR6[0][j][0][ii][0][0][l] = F;
}
}
}

F = 0;
if (  L11 > 0   ){
for(int jj = 1; jj < L11 + 1 ; jj++){
 for (int l = 0; l < L  ; l++) {
    F = (QY - Y3)*ARR6[0][j-1][0][0][jj][0][l] + (WY - QY)*ARR6[0][j-1][0][0][jj][0][l+1] ;

if( (j-1) > 0  ){
F = F + ( (j-1)/(2.0*Gam2))*( ARR6[0][j-2][0][0][jj][0][l] - (OM/Gam2)*ARR6[0][j-2][0][0][jj][0][l+1] ) ;

}
F = F + ( (jj)/( 2.0*G ) )*ARR6[0][j-1][0][0][jj-1][0][l+1] ;

ARR6[0][j][0][0][jj][0][l] = F;
}
}
}

F = 0.;
// 00c | 0j0
if (  L12 > 0   ){
for(int kk = 1; kk < L12 + 1 ; kk++){
 for (int l = 0; l < L  ; l++) {

    F = (QY - Y3)*ARR6[0][j-1][0][0][0][kk][l] + (WY - QY)*ARR6[0][j-1][0][0][0][kk][l+1] ;

    if( (j-1) > 0  ){
    F = F + ( (j-1)/(2.0*Gam2))*( ARR6[0][j-2][0][0][0][kk][l] - (OM/Gam2)*ARR6[0][j-2][0][0][0][kk][l+1] ) ;

    }
    ARR6[0][j][0][0][0][kk][l] = F;

}
}
}


if(   (L10 >= 1) && (  L11 >= 1 )  ) {
F=0;
// ab0 | 0j0
for(int jj = 1; jj < L11 + 1 ; jj++){
  for(int ii = 1; ii < L10 + 1 ; ii++){
    for (int l = 0; l < L  ; l++) { 
     F = (QY - Y3)*ARR6[0][j-1][0][ii][jj][0][l] + (WY - QY)*ARR6[0][j-1][0][ii][jj][0][l+1] ;

if( (j-1) > 0  ){
    F = F + ( (j-1)/(2.0*Gam2))*( ARR6[0][j-2][0][ii][jj][0][l] - (OM/Gam2)*ARR6[0][j-2][0][ii][jj][0][l+1] ) ;

   }
   F = F + ( (jj)/( 2.0*G ) )*ARR6[0][j-1][0][ii][jj-1][0][l+1] ;
ARR6[0][j][0][ii][jj][0][l] = F;
}
}
}
}


if(   (L10 >= 1) && (  L12 >= 1 )  ) {
F=0;  
//  a0b | 0j0
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
   for (int l = 0; l < L  ; l++) {
   F = (QY - Y3)*ARR6[0][j-1][0][ii][0][kk][l] + (WY - QY)*ARR6[0][j-1][0][ii][0][kk][l+1]  ;

if( (j-1) > 0  ){
F = F + ( (j-1)/(2.0*Gam2))*( ARR6[0][j-2][0][ii][0][kk][l] - (OM/Gam2)*ARR6[0][j-2][0][ii][0][kk][l+1] ) ;
}
ARR6[0][j][0][ii][0][kk][l] = F;
}
}
}
}


if(   (L11 >= 1) && (  L12 >= 1 )  ) {
F=0;  //  0ab | 0j0
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int jj = 1; jj < L11 +1  ; jj++){
    for (int l = 0; l < (L - L21); l++) {
      F = (QY - Y3)*ARR6[0][j-1][0][0][jj][kk][l] + (WY - QY)*ARR6[0][j-1][0][0][jj][kk][l+1]  ;

    if( (j-1) > 0  ){
    F = F + ( (j-1)/(2.0*Gam2))*( ARR6[0][j-2][0][0][jj][kk][l] - (OM/Gam2)*ARR6[0][j-2][0][0][jj][kk][l+1] ) ;

}
F = F + ( (jj)/( 2.0*G ) )*ARR6[0][j-1][0][0][jj-1][kk][l+1] ;

ARR6[0][j][0][0][jj][kk][l] = F;
}
}
}
}

// abc | 0j0
if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {
F=0;
for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
     for (int l = 0; l < L ; l++) {
 
      F = (QY - Y3)*ARR6[0][j-1][0][ii][jj][kk][l] + (WY - QY)*ARR6[0][j-1][0][ii][jj][kk][l+1] ;

    if( (j-1) > 0  ){
      F = F + ( (j-1)/(2.0*Gam2))*(ARR6[0][j-2][0][ii][jj][kk][l] - (OM/Gam2)*ARR6[0][j-2][0][ii][jj][kk][l+1] ) ;
     }
     F = F + ( (jj)/( 2.0*G ) )*ARR6[0][j-1][0][ii][jj-1][kk][l+1] ;

ARR6[0][j][0][ii][jj][kk][l] = F;
}
}
}
}
}

} // Closing j

} // L21 > 0

// abc | mn0
// arr5[j][i][n][m][k][l] 
// if C1 > 0:
if( (L20 > 0) && (  L21 > 0  )   )  {
for(int j = 1; j < L21 + 1 ; j++){
  for(int i = 1; i < L20 + 1 ; i++){

// a00 | ij0
if(  L10 > 0  ){  
F= 0;
for(int ii = 1; ii < L10 + 1 ; ii++){
  for (int l = 0; l < L  ; l++) {
    F = (QY - Y3)*ARR6[0][j-1][i][ii][0][0][l] + (WY - QY)*ARR6[0][j-1][i][ii][0][0][l+1] ;

  if( (j-1) > 0  ){
    F = F + ( (j-1)/(2.0*Gam2))*( ARR6[0][j-2][i][ii][0][0][l] - (OM/Gam2)*ARR6[0][j-2][i][ii][0][0][l+1] ) ;

  }
ARR6[0][j][i][ii][0][0][l] = F;
}
}
}

// 0b0 | ij0
if (  L11 > 0  ){
F = 0;
for(int jj = 1; jj < L11 + 1 ; jj++){
 for (int l = 0; l < L  ; l++) {
  F = (QY - Y3)*ARR6[0][j-1][i][0][jj][0][l] + (WY - QY)*ARR6[0][j-1][i][0][jj][0][l+1] ;

if( (j-1) > 0  ){
  F = F + ( (j-1)/(2.0*Gam2))*( ARR6[0][j-2][i][0][jj][0][l] - (OM/Gam2)*ARR6[0][j-2][i][0][jj][0][l+1] ) ;
 }
  F = F + ( (jj)/( 2.0*G ) )*ARR6[0][j-1][i][0][jj-1][0][l+1] ;

ARR6[0][j][i][0][jj][0][l] = F;
}
}
}

// 00c | ij0
if (  L12 > 0   ){
F = 0;
for(int kk = 1; kk < L12 + 1 ; kk++){
 for (int l = 0; l < L  ; l++) {
  F = (QY - Y3)*ARR6[0][j-1][i][0][0][kk][l] + (WY - QY)*ARR6[0][j-1][i][0][0][kk][l+1] ;

if( (j-1) > 0  ){
  F = F + ( (j-1)/(2.0*Gam2))*( ARR6[0][j-2][i][0][0][kk][l] - (OM/Gam2)*ARR6[0][j-2][i][0][0][kk][l+1] ) ;

}
ARR6[0][j][i][0][0][kk][l] = F;
}
}
}

// ab0 | ij0
if(   (L10 >= 1) && (  L11 >= 1 )  ) {
F=0;
for(int jj = 1; jj < L11 + 1 ; jj++){
  for(int ii = 1; ii < L10 + 1 ; ii++){
    for (int l = 0; l < (L - L10)  ; l++) { 
     F = (QY - Y3)*ARR6[0][j-1][i][ii][jj][0][l] + (WY - QY)*ARR6[0][j-1][i][ii][jj][0][l+1] ;

     if( (j-1) > 0  ){
    F = F + ( (j-1)/(2.0*Gam2))*( ARR6[0][j-2][i][ii][jj][0][l] - (OM/Gam2)*ARR6[0][j-2][i][ii][jj][0][l+1] ) ;
    }
    F = F + ( (jj)/( 2.0*G ) )*ARR6[0][j-1][i][ii][jj-1][0][l+1] ;
ARR6[0][j][i][ii][jj][0][l] = F;
}
}
}
}


if(   (L10 >= 1) && (  L12 >= 1 )  ) {
F=0;  
//  a0b | ij0
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
    for (int l = 0; l < L  ; l++) {
     F = (QY - Y3)*ARR6[0][j-1][i][ii][0][kk][l] + (WY - QY)*ARR6[0][j-1][i][ii][0][kk][l+1];

if( (j-1) > 0  ){
    F = F + ( (j-1)/(2.0*Gam2))*( ARR6[0][j-2][i][ii][0][kk][l] - (OM/Gam2)*ARR6[0][j-2][i][ii][0][kk][l+1] ) ;

}
ARR6[0][j][i][ii][0][kk][l] = F;
}
}
}
}


if(   (L11 >= 1) && (  L12 >= 1 )  ) {
F=0;  //  0ab | ij0
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int jj = 1; jj < L11 +1  ; jj++){
    for (int l = 0; l <  L  ; l++) {
     F = (QY - Y3)*ARR6[0][j-1][i][0][jj][kk][l] + (WY - QY)*ARR6[0][j-1][i][0][jj][kk][l+1]  ;

if( (j-1) > 0  ){
F = F + ( (j-1)/(2.0*Gam2))*( ARR6[0][j-2][i][0][jj][kk][l] - (OM/Gam2)*ARR6[0][j-2][i][0][jj][kk][l+1] ) ;
}
F = F + ( (jj)/( 2.0*G ) )*ARR6[0][j-1][i][0][jj-1][kk][l+1] ;
ARR6[0][j][i][0][jj][kk][l] = F;
}
}
}
}

// abc | ij0
if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {
F=0;
for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
    for (int l = 0; l < L ; l++) {
     F = (QY - Y3)*ARR6[0][j-1][i][ii][jj][kk][l] + (WY - QY)*ARR6[0][j-1][i][ii][jj][kk][l+1] ;

        if( (j-1) > 0  ){
         F = F + ( (j-1)/(2.0*Gam2))*(ARR6[0][j-2][i][ii][jj][kk][l] - (OM/Gam2)*ARR6[0][j-2][i][ii][jj][kk][l+1] ) ;
     }

F = F + ( (jj)/( 2.0*G ) )*ARR6[0][j-1][i][ii][jj-1][kk][l+1] ;
ARR6[0][j][i][ii][jj][kk][l] = F;
}
}
}
}
}


} // i
} // j
} // L20 > 0 && L21 > 0

// FROM   [C3[C2][C1][A1][A2[A3]
//  abc | 00k
for(int k = 1; k < L22 + 1; k++){

for (int l = 0; l < L  ; l++) {

F = (QZ - Z3)*ARR6[k-1][0][0][0][0][0][l] + (WZ - QZ)*ARR6[k-1][0][0][0][0][0][l+1] ;

if( (k-1) > 0  ){
 F = F + ( (k-1)/(2.0*Gam2))*( ARR6[k-2][0][0][0][0][0][l] - (OM/Gam2)*ARR6[k-2][0][0][0][0][0][l+1] ) ;
}
ARR6[k][0][0][0][0][0][l] = F;
}

// a00 | 00k
if(  L10 > 0  ){  
F= 0;
for(int ii = 1; ii < L10 + 1 ; ii++){
  for (int l = 0; l < L  ; l++) {
    F = (QZ - Z3)*ARR6[k-1][0][0][ii][0][0][l] + (WZ - QZ)*ARR6[k-1][0][0][ii][0][0][l+1] ;

if( (k-1) > 0  ){
F = F + ( (k-1)/(2.0*Gam2))*( ARR6[k-2][0][0][ii][0][0][l] - (OM/Gam2)*ARR6[k-2][0][0][ii][0][0][l+1] ) ;
}
ARR6[k][0][0][ii][0][0][l] = F;
}
}
}

// 0b0 | 00k
if (  L11 > 0   ){
F = 0;
for(int jj = 1; jj < L11 + 1 ; jj++){
 for (int l = 0; l < L  ; l++) {
  F = (QZ - Z3)*ARR6[k-1][0][0][0][jj][0][l] + (WZ - QZ)*ARR6[k-1][0][0][0][jj][0][l+1] ;

if( (k-1) > 0  ){
  F = F + ( (k-1)/(2.0*Gam2))*( ARR6[k-2][0][0][0][jj][0][l] - (OM/Gam2)*ARR6[k-2][0][0][0][jj][0][l+1] ) ;
}
ARR6[k][0][0][0][jj][0][l] = F;
}
}
}


// 00c | 00k
if (  L12 > 0   ){
F = 0;
for(int kk = 1; kk < L12 + 1 ; kk++){
 for (int l = 0; l < (L - L12)  ; l++) {

    F = (QZ - Z3)*ARR6[k-1][0][0][0][0][kk][l] + (WZ - QZ)*ARR6[k-1][0][0][0][0][kk][l+1] ;

    if( ( k - 1) > 0  ){
    F = F + ( (k-1)/(2.0*Gam2))*( ARR6[k-2][0][0][0][0][kk][l] - (OM/Gam2)*ARR6[k-2][0][0][0][0][kk][l+1] ) ;

}
F = F + ( (kk)/( 2.0*G ) )*ARR6[k-1][0][0][0][0][kk-1][l+1] ;
ARR6[k][0][0][0][0][kk][l] = F;
}
}
}

// ab0 | 00k
if(   (L10 >= 1) && (  L11 >= 1 )  ) {
F=0;
for(int jj = 1; jj < L11 + 1 ; jj++){
  for(int ii = 1; ii < L10 + 1 ; ii++){
    for (int l = 0; l < L  ; l++) { 
      F = (QZ - Z3)*ARR6[k-1][0][0][ii][jj][0][l] + (WZ - QZ)*ARR6[k-1][0][0][ii][jj][0][l+1] ;

if( (k-1) > 0  ){
 F = F + ( (k-1)/(2.0*Gam2))*( ARR6[k-2][0][0][ii][jj][0][l] - (OM/Gam2)*ARR6[k-2][0][0][ii][jj][0][l+1] ) ;

}
ARR6[k][0][0][ii][jj][0][l] = F;
}
}
}
}

if(   (L10 >= 1) && (  L12 >= 1 )  ) {
F=0;  
//  a0c | 00k
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
    for (int l = 0; l < L  ; l++) {
     F = (QZ - Z3)*ARR6[k-1][0][0][ii][0][kk][l] + (WZ - QZ)*ARR6[k-1][0][0][ii][0][kk][l+1];

if( (k-1) > 0  ){
     F = F + ( (k-1)/(2.0*Gam2))*( ARR6[k-2][0][0][ii][0][kk][l] - (OM/Gam2)*ARR6[k-2][0][0][ii][0][kk][l+1] ) ;

}
F = F + ( (kk)/( 2.0*G ) )*ARR6[k-1][0][0][ii][0][kk-1][l+1] ;
ARR6[k][0][0][ii][0][kk][l] = F;
}
}
}
}


if(   (L11 >= 1) && (  L12 >= 1 )  ) {
F=0;  //  0bc | 00k
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int jj = 1; jj < L11 +1  ; jj++){
    for (int l = 0; l < L  ; l++) {
      F = (QZ - Z3)*ARR6[k-1][0][0][0][jj][kk][l] + (WZ - QZ)*ARR6[k-1][0][0][0][jj][kk][l+1]  ;

    if( (k-1) > 0  ){
     F = F + ( (k-1)/(2.0*Gam2))*( ARR6[k-2][0][0][0][jj][kk][l] - (OM/Gam2)*ARR6[k-2][0][0][0][jj][kk][l+1] ) ;

}
F = F + ( (kk)/( 2.0*G ) )*ARR6[k-1][0][0][0][jj][kk-1][l+1] ;

ARR6[k][0][0][0][jj][kk][l] = F;
}
}
}
}

// abc | 00k
if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {
F=0;
for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
    for (int l = 0; l < L ; l++) {
        F = (QZ - Z3)*ARR6[k-1][0][0][ii][jj][kk][l] + (WZ - QZ)*ARR6[k-1][0][0][ii][jj][kk][l+1] ;

if( (k-1) > 0  ){
  F = F + ( (k-1)/(2.0*Gam2))*(ARR6[k-2][0][0][ii][jj][kk][l] - (OM/Gam2)*ARR6[k-2][0][0][ii][jj][kk][l+1] ) ;

}
F = F + ( (kk)/( 2.0*G ) )*ARR6[k-1][0][0][ii][jj][kk-1][l+1] ;
ARR6[k][0][0][ii][jj][kk][l] = F;
}
}
}
}
}


} // fin k for


//  abc | i0k si i>0
if (L20 > 0 ){  
//abc | i0k
for(int k = 1; k < L22 + 1; k++){
    for(int i = 1; i < L20 + 1 ; i++){
        for (int l = 0; l < L ; l++) {
F=0;
// inicial 000 | i0k

F = (QZ - Z3)*ARR6[k-1][0][i][0][0][0][l] + (WZ - QZ)*ARR6[k-1][0][i][0][0][0][l+1] ;

if( (k-1) > 0  ){
F = F + ( (k-1)/(2.0*Gam2))*(ARR6[k-2][0][i][0][0][0][l] - (OM/Gam2)*ARR6[k-2][0][i][0][0][0][l+1] ) ;

}
ARR6[k][0][i][0][0][0][l] = F;
}


if(  L10 > 0  ){  
F= 0;
// a00 | i0k
for(int ii = 1; ii < L10 + 1 ; ii++){
  for (int l = 0; l < L  ; l++) {
    F = (QZ - Z3)*ARR6[k-1][0][i][ii][0][0][l] + (WZ - QZ)*ARR6[k-1][0][i][ii][0][0][l+1] ;

if( (k-1) > 0  ){
    F = F + ( (k-1)/(2.0*Gam2))*( ARR6[k-2][0][i][ii][0][0][l] - (OM/Gam2)*ARR6[k-2][0][i][ii][0][0][l+1] ) ;

}
ARR6[k][0][i][ii][0][0][l] = F;
}
}
}


// 0b0 | i0k
if (  L11 > 0   ){
F = 0.;
for(int jj = 1; jj < L11 + 1 ; jj++){
 for (int l = 0; l < L  ; l++) {
    F = (QZ - Z3)*ARR6[k-1][0][i][0][jj][0][l] + (WZ - QZ)*ARR6[k-1][0][i][0][jj][0][l+1] ;

if( (k-1) > 0  ){
    F = F + ( (k-1)/(2.0*Gam2))*( ARR6[k-2][0][i][0][jj][0][l] - (OM/Gam2)*ARR6[k-2][0][i][0][jj][0][l+1] ) ;

}
ARR6[k][0][i][0][jj][0][l] = F;
}
}
}


// 00c | i0k
if (  L12 > 0   ){
F = 0;
for(int kk = 1; kk < L12 + 1 ; kk++){
 for (int l = 0; l < L  ; l++) {
    F = (QZ - Z3)*ARR6[k-1][0][i][0][0][kk][l] + (WZ - QZ)*ARR6[k-1][0][i][0][0][kk][l+1] ;

if( (k-1) > 0  ){
F = F + ( (k-1)/(2.0*Gam2))*( ARR6[k-2][0][i][0][0][kk][l] - (OM/Gam2)*ARR6[k-2][0][i][0][0][kk][l+1] ) ;

}
F = F + ( (kk)/( 2.0*G ) )*ARR6[k-1][0][i][0][0][kk-1][l+1] ;
ARR6[k][0][i][0][0][kk][l] = F;
}
}
}

// ab0 | i0k
if(   (L10 >= 1) && (  L11 >= 1 )  ) {
F=0;
for(int jj = 1; jj < L11 + 1 ; jj++){
  for(int ii = 1; ii < L10 + 1 ; ii++){
    for (int l = 0; l < L  ; l++) { 
     F = (QZ - Z3)*ARR6[k-1][0][i][ii][jj][0][l] + (WZ - QZ)*ARR6[k-1][0][i][ii][jj][0][l+1] ;

if( (k-1) > 0  ){
F = F + ( (k-1)/(2.0*Gam2))*( ARR6[k-2][0][i][ii][jj][0][l] - (OM/Gam2)*ARR6[k-2][0][i][ii][jj][0][l+1] ) ;

}
ARR6[k][0][i][ii][jj][0][l] = F;
}
}
}
}


if(   (L10 >= 1) && (  L12 >= 1 )  ) {
F=0;  
//  a0c | i0k
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
for (int l = 0; l < L  ; l++) {
  F = (QZ - Z3)*ARR6[k-1][0][i][ii][0][kk][l] + (WZ - QZ)*ARR6[k-1][0][i][ii][0][kk][l+1];

if( (k-1) > 0  ){
F = F + ( (k-1)/(2.0*Gam2))*( ARR6[k-2][0][i][ii][0][kk][l] - (OM/Gam2)*ARR6[k-2][0][i][ii][0][kk][l+1] ) ;

}
F = F + ( (kk)/( 2.0*G ) )*ARR6[k-1][0][i][ii][0][kk-1][l+1] ;
ARR6[k][0][i][ii][0][kk][l] = F;
}
}
}
}

if(   (L11 >= 1) && (  L12 >= 1 )  ) {
F=0;  //  0bc | i0k
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int jj = 1; jj < L11 +1  ; jj++){
    for (int l = 0; l < L  ; l++) {

F = (QZ - Z3)*ARR6[k-1][0][i][0][jj][kk][l] + (WZ - QZ)*ARR6[k-1][0][i][0][jj][kk][l+1]  ;
if( (k-1) > 0  ){
F = F + ( (k-1)/(2.0*Gam2))*( ARR6[k-2][0][i][0][jj][kk][l] - (OM/Gam2)*ARR6[k-2][0][i][0][jj][kk][l+1] ) ;

}
F = F + ( (kk)/( 2.0*G ) )*ARR6[k-1][0][i][0][jj][kk-1][l+1] ;

ARR6[k][0][i][0][jj][kk][l] = F;
}
}
}
}


//  abc | i0k
if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {
F=0;
for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
    for (int l = 0; l < L ; l++) {
     F = (QZ - Z3)*ARR6[k-1][0][i][ii][jj][kk][l] + (WZ - QZ)*ARR6[k-1][0][i][ii][jj][kk][l+1] ;

if( (k-1) > 0  ){
F = F + ( (k-1)/(2.0*Gam2))*(ARR6[k-2][0][i][ii][jj][kk][l] - (OM/Gam2)*ARR6[k-2][0][i][ii][jj][kk][l+1] ) ;

}
F = F + ( (kk)/( 2.0*G ) )*ARR6[k-1][0][i][ii][jj][kk-1][l+1] ;

ARR6[k][0][i][ii][jj][kk][l] = F;
}
}
}
}
}


// closing k , i
}
}

}// L20 > 0

//  abc | 0jk
if (L21 > 0 ){  
for(int k = 1; k < L22 + 1; k++){
 for(int j = 1; j < L21 + 1 ; j++){
// initial: 000 | 0jk  F=0;
    for (int l = 0; l < L ; l++) {
      F = (QZ - Z3)*ARR6[k-1][j][0][0][0][0][l] + (WZ - QZ)*ARR6[k-1][j][0][0][0][0][l+1] ;

if( (k-1) > 0  ){
    F = F + ( (k-1)/(2.0*Gam2))*(ARR6[k-2][j][0][0][0][0][l] - (OM/Gam2)*ARR6[k-2][j][0][0][0][0][l+1] ) ;
}
ARR6[k][j][0][0][0][0][l] = F;
}

// a00 | 0jk
if(  L10 > 0  ){  
F= 0;
for(int ii = 1; ii < L10 + 1 ; ii++){
  for (int l = 0; l < L  ; l++) {
    F = (QZ - Z3)*ARR6[k-1][j][0][ii][0][0][l] + (WZ - QZ)*ARR6[k-1][j][0][ii][0][0][l+1] ;

if( (k-1) > 0  ){
    F = F + ( (k-1)/(2.0*Gam2))*( ARR6[k-2][j][0][ii][0][0][l] - (OM/Gam2)*ARR6[k-2][j][0][ii][0][0][l+1] ) ;

}
ARR6[k][j][0][ii][0][0][l] = F;
}
}
}

// 0b0 | 0jk
if (  L11 > 0   ){
F = 0.;
for(int jj = 1; jj < L11 + 1 ; jj++){
 for (int l = 0; l < L  ; l++) {
    F = (QZ - Z3)*ARR6[k-1][j][0][0][jj][0][l] + (WZ - QZ)*ARR6[k-1][j][0][0][jj][0][l+1] ;

if( (k-1) > 0  ){
  F = F + ( (k-1)/(2.0*Gam2))*( ARR6[k-2][j][0][0][jj][0][l] - (OM/Gam2)*ARR6[k-2][j][0][0][jj][0][l+1] ) ;
}
ARR6[k][j][0][0][jj][0][l] = F;
}
}
}


// 00c | 0jk
if (  L12 > 0   ){
F = 0;
for(int kk = 1; kk < L12 + 1 ; kk++){
 for (int l = 0; l < L  ; l++) {

    F = (QZ- Z3)*ARR6[k-1][j][0][0][0][kk][l] + (WZ - QZ)*ARR6[k-1][j][0][0][0][kk][l+1] ;

if( (k-1) > 0  ){
    F = F + ( (k-1)/(2.0*Gam2))*( ARR6[k-2][j][0][0][0][kk][l] - (OM/Gam2)*ARR6[k-2][j][0][0][0][kk][l+1] ) ;

}
F = F + ( (kk)/( 2.0*G ) )*ARR6[k-1][j][0][0][0][kk-1][l+1] ;
ARR6[k][j][0][0][0][kk][l] = F;
}
}
}

// ab0 | 0jk
if(   (L10 >= 1) && (  L11 >= 1 )  ) {
F=0;
for(int jj = 1; jj < L11 + 1 ; jj++){
  for(int ii = 1; ii < L10 + 1 ; ii++){
    for (int l = 0; l < L  ; l++) { 

      F = (QZ - Z3)*ARR6[k-1][j][0][ii][jj][0][l] + (WZ - QZ)*ARR6[k-1][j][0][ii][jj][0][l+1] ;

if( (k-1) > 0  ){
      F = F + ( (k-1)/(2.0*Gam2))*( ARR6[k-2][j][0][ii][jj][0][l] - (OM/Gam2)*ARR6[k-2][j][0][ii][jj][0][l+1] ) ;

}
ARR6[k][j][0][ii][jj][0][l] = F;
}
}
}
}

//  a0c | 0jk
if(   (L10 >= 1) && (  L12 >= 1 )  ) {
F=0;  
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
    for (int l = 0; l < L  ; l++) {

    F = (QZ - Z3)*ARR6[k-1][j][0][ii][0][kk][l] + (WZ - QZ)*ARR6[k-1][j][0][ii][0][kk][l+1];

if( (k-1) > 0  ){
    F = F + ( (k-1)/(2.0*Gam2))*( ARR6[k-2][j][0][ii][0][kk][l] - (OM/Gam2)*ARR6[k-2][j][0][ii][0][kk][l+1] ) ;

}
F = F + ( (kk)/( 2.0*G) )*ARR6[k-1][j][0][ii][0][kk-1][l+1] ;

ARR6[k][j][0][ii][0][kk][l] = F;
}
}
}
}

 //  0bc | 0jk
if(   (L11 >= 1) && (  L12 >= 1 )  ) {
F=0; 
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int jj = 1; jj < L11 +1  ; jj++){
    for (int l = 0; l < L  ; l++) {
      F = (QZ - Z3)*ARR6[k-1][j][0][0][jj][kk][l] + (WZ - QZ)*ARR6[k-1][j][0][0][jj][kk][l+1]  ;

if( (k-1) > 0  ){
    F = F + ( (k-1)/(2.0*Gam2))*( ARR6[k-2][j][0][0][jj][kk][l] - (OM/Gam2)*ARR6[k-2][j][0][0][jj][kk][l+1] ) ;

}
F = F + ( (kk)/( 2.0*G) )*ARR6[k-1][j][0][0][jj][kk-1][l+1] ;
ARR6[k][j][0][0][jj][kk][l] = F;
}
}
}
}

//  abc | 0jk
if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {
F=0;
for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
    for (int l = 0; l < L ; l++) {

    F = (QZ - Z3)*ARR6[k-1][j][0][ii][jj][kk][l] + (WZ - QZ)*ARR6[k-1][j][0][ii][jj][kk][l+1] ;

if( (k-1) > 0  ){
F = F + ( (k-1)/(2.0*Gam2))*(ARR6[k-2][j][0][ii][jj][kk][l] - (OM/Gam2)*ARR6[k-2][j][0][ii][jj][kk][l+1] ) ;
}
F = F + ( (kk)/( 2.0*G) )*ARR6[k-1][j][0][ii][jj][kk-1][l+1] ;

ARR6[k][j][0][ii][jj][kk][l] = F;
}
}
}
}
}

// cierre k , j
}
}

}// L21 > 0

// el completo :
//  abc | ijk
//seguro anti memria fallida:
if (  (L21 > 0 )  &&   (L20 > 0 )      ){  

for(int k = 1; k < L22 + 1; k++){
    for(int j = 1; j < L21 + 1 ; j++){
      for(int i = 1; i < L20 + 1 ; i++){
// caso inicial: 000 | ijk
 for (int l = 0; l < L ; l++) {
    F = (QZ - Z3)*ARR6[k-1][j][i][0][0][0][l] + (WZ - QZ)*ARR6[k-1][j][i][0][0][0][l+1] ;

if( (k-1) > 0  ){
    F = F + ( (k-1)/(2.0*Gam2))*(ARR6[k-2][j][i][0][0][0][l] - (OM/Gam2)*ARR6[k-2][j][i][0][0][0][l+1] ) ;

}
ARR6[k][j][i][0][0][0][l] = F;
}


// a00 | ijk
if(  L10 > 0  ){  
F= 0;
for(int ii = 1; ii < L10 + 1 ; ii++){
  for (int l = 0; l < L  ; l++) {
    F = (QZ - Z3)*ARR6[k-1][j][i][ii][0][0][l] + (WZ - QZ)*ARR6[k-1][j][i][ii][0][0][l+1] ;

if( (k-1) > 0  ){
  F = F + ( (k-1)/(2.0*Gam2))*( ARR6[k-2][j][i][ii][0][0][l] - (OM/Gam2)*ARR6[k-2][j][i][ii][0][0][l+1] ) ;

}
ARR6[k][j][i][ii][0][0][l] = F;
}
}
}

// 0b0 | ijk
if (  L11 > 0   ){
F = 0;
for(int jj = 1; jj < L11 + 1 ; jj++){
 for (int l = 0; l < L  ; l++) {

  F = (QZ - Z3)*ARR6[k-1][j][i][0][jj][0][l] + (WZ - QZ)*ARR6[k-1][j][i][0][jj][0][l+1] ;

if( (k-1) > 0  ){
  F = F + ( (k-1)/(2.0*Gam2))*( ARR6[k-2][j][i][0][jj][0][l] - (OM/Gam2)*ARR6[k-2][j][i][0][jj][0][l+1] ) ;

}
ARR6[k][j][i][0][jj][0][l] = F;
}
}
}

// 00c | ijk
if (  L12 > 0   ){
F = 0;
for(int kk = 1; kk < L12 + 1 ; kk++){
 for (int l = 0; l < L  ; l++) {
    F = (QZ- Z3)*ARR6[k-1][j][i][0][0][kk][l] + (WZ - QZ)*ARR6[k-1][j][i][0][0][kk][l+1] ;

if( (k-1) > 0  ){
    F = F + ( (k-1)/(2.0*Gam2))*( ARR6[k-2][j][i][0][0][kk][l] - (OM/Gam2)*ARR6[k-2][j][i][0][0][kk][l+1] ) ;

}
F = F + ( (kk)/( 2.0*G) )*ARR6[k-1][j][i][0][0][kk-1][l+1] ;
ARR6[k][j][i][0][0][kk][l] = F;

}
}
}

// ab0 | ijk
if(   (L10 >= 1) && (  L11 >= 1 )  ) {
F=0;
for(int jj = 1; jj < L11 + 1 ; jj++){
  for(int ii = 1; ii < L10 + 1 ; ii++){
    for (int l = 0; l < L  ; l++) { 
     F = (QZ - Z3)*ARR6[k-1][j][i][ii][jj][0][l] + (WZ - QZ)*ARR6[k-1][j][i][ii][jj][0][l+1] ;

if( (k-1) > 0  ){
    F = F + ( (k-1)/(2.0*Gam2))*( ARR6[k-2][j][i][ii][jj][0][l] - (OM/Gam2)*ARR6[k-2][j][i][ii][jj][0][l+1] ) ;

}
ARR6[k][j][i][ii][jj][0][l] = F;
}
}
}
}

//  a0c | ijk
if(   (L10 >= 1) && (  L12 >= 1 )  ) {
F=0;  
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
    for (int l = 0; l < L  ; l++) {

     F = (QZ - Z3)*ARR6[k-1][j][i][ii][0][kk][l] + (WZ - QZ)*ARR6[k-1][j][i][ii][0][kk][l+1];

if( (k-1) > 0  ){
F = F + ( (k-1)/(2.0*Gam2))*( ARR6[k-2][j][i][ii][0][kk][l] - (OM/Gam2)*ARR6[k-2][j][i][ii][0][kk][l+1] ) ;
}
F = F + ( (kk)/( 2.0*G ) )*ARR6[k-1][j][i][ii][0][kk-1][l+1] ;

ARR6[k][j][i][ii][0][kk][l] = F;
}
}
}
}

 //  0bc | ijk
if(   (L11 >= 1) && (  L12 >= 1 )  ) {
F=0; 
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int jj = 1; jj < L11 +1  ; jj++){
    for (int l = 0; l < L  ; l++) {

    F = (QZ - Z3)*ARR6[k-1][j][i][0][jj][kk][l] + (WZ - QZ)*ARR6[k-1][j][i][0][jj][kk][l+1]  ;

if( (k-1) > 0  ){
    F = F + ( (k-1)/(2.0*Gam2))*( ARR6[k-2][j][i][0][jj][kk][l] - (OM/Gam2)*ARR6[k-2][j][i][0][jj][kk][l+1] ) ;
}
F = F + ( (kk)/( 2.0*G) )*ARR6[k-1][j][i][0][jj][kk-1][l+1] ;

ARR6[k][j][i][0][jj][kk][l] = F;
}
}
}
}

//  abc | ijk
if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {
F=0;
for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
    for (int l = 0; l < (L - k); l++) {

    F = (QZ - Z3)*ARR6[k-1][j][i][ii][jj][kk][l] + (WZ - QZ)*ARR6[k-1][j][i][ii][jj][kk][l+1] ;

if( (k-1) > 0  ){
F = F + ( (k-1)/(2.0*Gam2))*(ARR6[k-2][j][i][ii][jj][kk][l] - (OM/Gam2)*ARR6[k-2][j][i][ii][jj][kk][l+1] ) ;
}
    F = F + ( (kk)/( 2.0*G) )*ARR6[k-1][j][i][ii][jj][kk-1][l+1] ;

ARR6[k][j][i][ii][jj][kk][l] = F;
}
}
}
}
}


} // i
} // j
} // k
// closing if i,j > 0
}

AC = ARR6[L22][L21][L20][L10][L11][L12][0]; 

for (int i = 0; i < L22 + 1; i++) {
    for (int j = 0; j < L21 + 1; j++) {
        for (int k = 0; k < L20 + 1; k++) {
            for (int l = 0; l < L10 + 1; l++) {
                for (int m = 0; m < L11 + 1; m++) {
                    for (int n = 0; n < L12 + 1; n++) {
                        delete[] ARR6[i][j][k][l][m][n];
                    }
                    delete[] ARR6[i][j][k][l][m];
                }
                delete[] ARR6[i][j][k][l];
            }
            delete[] ARR6[i][j][k];
        }
        delete[] ARR6[i][j];
    }
    delete[] ARR6[i];
}
delete[] ARR6;

return AC;

}

}

}

}

//delete [] OOOMOOO;
return AC;
}


Double_v Electron_electron_repulsion::selectorV( G4double alf_1, G4double alf_2, G4double alf_3,
Double_v Alf4, G4int fv1a, G4int fv1b, G4int fv1c,
G4int fv2a, G4int fv2b, G4int fv2c, G4int fv3a, G4int fv3b,
G4int fv3c, G4int fv4a, G4int fv4b, G4int fv4c, 
Double_v X1, Double_v Y1, Double_v Z1, Double_v X2,
Double_v Y2, Double_v Z2, Double_v X3, Double_v Y3, Double_v Z3,
Double_v X4, Double_v Y4, Double_v Z4, G4double n1, G4double n2,
G4double n3, Double_v N4, G4int cont) {


Double_v JepV = 0;

if(   (fv4a == 0) &&   (fv4b == 0) &&   (fv4c == 0)  ){

//check if a0mc0:
if(   (fv2a == 0) &&   (fv2b == 0) &&   ( fv2c == 0)  ){

JepV = JepV + repulsion_a0mc0V(alf_1, alf_2, alf_3, Alf4, fv1a, fv1b, fv1c, fv3a, fv3b, fv3c, 
X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, n1, n2, n3, N4, cont);

return JepV;

}else{


if (    (fv2a > 0)  &&   (fv2b == 0)     &&   (fv2c == 0)     ){

fv1a= fv1a + 1;
fv2a= fv2a - 1;

JepV = selectorV( alf_1, alf_2, alf_3, Alf4, fv1a, fv1b, fv1c,
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c, 
X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, n1, n2, n3, N4, cont);

//checking:
MaskD_v fx = ( vecCore::math::Abs( X1 - X2 ) > 0.0 );
if( MaskFull(fx)  ){

fv1a= fv1a - 1;

JepV = JepV + ( X1 - X2 )*selectorV(alf_1,  alf_2, alf_3, Alf4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
X1, Y1, Z1, X2, Y2, Z2,  X3, Y3,  Z3, X4, Y4, Z4, n1, n2, n3,  N4 , cont);

fv1a=fv1a + 1 ;

}

return JepV;

}else{

if (    (fv2a == 0)  &&   (fv2b > 0)     &&   (fv2c == 0)     ){

fv1b = fv1b + 1;
fv2b = fv2b - 1;

JepV = selectorV( alf_1, alf_2, alf_3, Alf4, fv1a, fv1b, fv1c,
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c, 
X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, n1, n2, n3, N4, cont);

MaskD_v fy = ( vecCore::math::Abs( Y1 - Y2 ) > 0.0 );
if( MaskFull(fy)  ){

fv1b =fv1b - 1;

JepV = JepV + ( Y1 - Y2 )*selectorV(alf_1,  alf_2, alf_3, Alf4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
X1, Y1, Z1, X2, Y2, Z2,  X3, Y3,  Z3, X4, Y4, Z4, n1, n2, n3,  N4 , cont);

fv1b =fv1b + 1;

}

return JepV;

}else{

// 00v
if (   (fv2a == 0)  &&   (fv2b == 0)     &&   (fv2c > 0)      ){

fv1c =fv1c+1;
fv2c =fv2c-1;

JepV = selectorV( alf_1, alf_2, alf_3, Alf4, fv1a, fv1b, fv1c,
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c, 
X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, n1, n2, n3, N4, cont);

MaskD_v fz = ( vecCore::math::Abs( Z1 - Z2 ) > 0.0 );
if( MaskFull(fz)  ){

fv1c=fv1c-1;

JepV = JepV + ( Z1 - Z2 )*selectorV(alf_1,  alf_2, alf_3, Alf4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
X1, Y1, Z1, X2, Y2, Z2,  X3, Y3,  Z3, X4, Y4, Z4, n1, n2, n3,  N4 , cont);

fv1c = fv1c+1;

}

return JepV;

}else{

if (   (fv2a > 0)  &&   (fv2b > 0)     &&   (fv2c == 0)    ){

fv1a = fv1a + 1;
fv2a = fv2a -1;

JepV = selectorV( alf_1, alf_2, alf_3, Alf4, fv1a, fv1b, fv1c,
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c, 
X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, n1, n2, n3, N4, cont);

MaskD_v fx = ( vecCore::math::Abs( X1 - X2 ) > 0.0 );
if( MaskFull(fx)  ){

fv1a = fv1a-1;

JepV = JepV + ( X1 - X2 )*selectorV(alf_1,  alf_2, alf_3, Alf4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
X1, Y1, Z1, X2, Y2, Z2,  X3, Y3,  Z3, X4, Y4, Z4, n1, n2, n3,  N4 , cont);

fv1a=fv1a+1;

}

}else{

if (   (fv2a == 0)  &&   (fv2b > 0)     &&   (fv2c > 0)       ){

fv1b= fv1b + 1;
fv2b= fv2b - 1;

JepV = selectorV( alf_1, alf_2, alf_3, Alf4, fv1a, fv1b, fv1c,
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c, 
X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, n1, n2, n3, N4, cont);

MaskD_v fy = ( vecCore::math::Abs( Y1 - Y2 ) > 0.0 );
if( MaskFull(fy)  ){

fv1b = fv1b - 1;

JepV = JepV + ( Y1 - Y2 )*selectorV(alf_1,  alf_2, alf_3, Alf4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
X1, Y1, Z1, X2, Y2, Z2,  X3, Y3,  Z3, X4, Y4, Z4, n1, n2, n3,  N4 , cont);

fv1b =fv1b + 1;

}

}else{

if (   (fv2a > 0)  &&   (fv2b == 0)     &&   (fv2c > 0)       ){

fv1a=fv1a+1;
fv2a=fv2a-1;

JepV = JepV + selectorV( alf_1, alf_2, alf_3, Alf4, fv1a, fv1b, fv1c,
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c, 
X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, n1, n2, n3, N4, cont);

MaskD_v fx = ( vecCore::math::Abs( X1 - X2 ) > 0.0 );
if( MaskFull(fx)  ){

fv1a = fv1a-1;

JepV = JepV + ( X1 - X2 )*selectorV(alf_1,  alf_2, alf_3, Alf4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
X1, Y1, Z1, X2, Y2, Z2,  X3, Y3,  Z3, X4, Y4, Z4, n1, n2, n3,  N4 , cont);

fv1a = fv1a+1;

}


}else{

fv1a=fv1a+1;
fv2a=fv2a-1;

JepV = JepV + selectorV(alf_1,  alf_2, alf_3, Alf4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
X1, Y1, Z1, X2, Y2, Z2,  X3, Y3,  Z3, X4, Y4, Z4, n1, n2, n3,  N4 , cont);

MaskD_v fx = ( vecCore::math::Abs( X1 - X2 ) > 0.0 );
if( MaskFull(fx)  ){

fv1a=fv1a-1;

JepV = JepV + ( X1 - X2 )*selectorV(alf_1,  alf_2, alf_3, Alf4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
X1, Y1, Z1, X2, Y2, Z2,  X3, Y3,  Z3, X4, Y4, Z4, n1, n2, n3,  N4 , cont);


fv1a = fv1a + 1;

}
}
}
}
}
}
}
}

// check abmcd:
}else{

if (    (fv4a > 0)  &&   (fv4b == 0)     &&   (fv4c == 0)     ){
// a | b | c |n00
fv3a=fv3a+1;
fv4a=fv4a-1;

JepV = JepV + selectorV( alf_1, alf_2, alf_3, Alf4, fv1a, fv1b, fv1c,
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c, 
X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, n1, n2, n3, N4, cont);

MaskD_v fx2 = ( vecCore::math::Abs( X3 - X4 ) > 0.0 );
if( MaskFull( fx2 )  ){

fv3a = fv3a-1;

JepV = JepV + ( X3 - X4 )*selectorV(alf_1,  alf_2, alf_3, Alf4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
X1, Y1, Z1, X2, Y2, Z2,  X3, Y3,  Z3, X4, Y4, Z4, n1, n2, n3,  N4 , cont);

fv3a = fv3a+1;

}

}else{

if (    (fv4a == 0)  &&   (fv4b > 0)     &&   (fv4c == 0)     ){
// a | b | c |0m0
fv3b = fv3b+1;
fv4b = fv4b-1;

JepV = JepV + selectorV( alf_1, alf_2, alf_3, Alf4, fv1a, fv1b, fv1c,
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c, 
X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, n1, n2, n3, N4, cont);

MaskD_v fy2 = ( vecCore::math::Abs( Y3 - Y4 ) > 0.0 );
if( MaskFull( fy2 )  ){

fv3b=fv3b-1;

JepV = JepV + ( Y3 - Y4 )*selectorV(alf_1,  alf_2, alf_3, Alf4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
X1, Y1, Z1, X2, Y2, Z2,  X3, Y3,  Z3, X4, Y4, Z4, n1, n2, n3,  N4 , cont);

fv3b=fv3b+1;


}

}else{

if (    (fv4a == 0)  &&   (fv4b == 0)     &&   (fv4c > 0)     ){
// a | b | c |00k
fv3c=fv3c+1;
fv4c=fv4c-1;

JepV = JepV + selectorV( alf_1, alf_2, alf_3, Alf4, fv1a, fv1b, fv1c,
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c, 
X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, n1, n2, n3, N4, cont);

MaskD_v fz2 = ( vecCore::math::Abs( Z3 - Z4 ) > 0.0 );
if( MaskFull( fz2 )  ){

fv3c =fv3c-1;

JepV = JepV + ( Z3 - Z4 )*selectorV(alf_1,  alf_2, alf_3, Alf4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
X1, Y1, Z1, X2, Y2, Z2,  X3, Y3,  Z3, X4, Y4, Z4, n1, n2, n3,  N4 , cont);


fv3c =fv3c+1;

}


}else{

if (   (fv4a > 0)  &&  (fv4b > 0)  &&  (fv4c == 0)       ){
// ab |c nh0
fv3a =fv3a +1;
fv3a =fv4a -1;

JepV = JepV + selectorV( alf_1, alf_2, alf_3, Alf4, fv1a, fv1b, fv1c,
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c, 
X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, n1, n2, n3, N4, cont);

MaskD_v fx2 = ( vecCore::math::Abs( X3 - X4 ) > 0.0 );
if( MaskFull( fx2 )  ){

fv3a =fv3a-1;

JepV = JepV + ( X3 - X4 )*selectorV(alf_1,  alf_2, alf_3, Alf4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
X1, Y1, Z1, X2, Y2, Z2,  X3, Y3,  Z3, X4, Y4, Z4, n1, n2, n3,  N4 , cont);


fv3a=fv3a+1;

}

}else{

if (   (fv4a == 0)  &&  (fv4b > 0)  &&  (fv4c > 0)   ){
// abc | 0hm
fv3b=fv3b+1;
fv3b=fv4b-1;

JepV = JepV + selectorV( alf_1, alf_2, alf_3, Alf4, fv1a, fv1b, fv1c,
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c, 
X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, n1, n2, n3, N4, cont);

MaskD_v fy2 = ( vecCore::math::Abs( Y3 - Y4 ) > 0.0 );
if( MaskFull( fy2 )  ){

fv3b=fv3b-1;

JepV = JepV + ( Y3 - Y4 )*selectorV(alf_1,  alf_2, alf_3, Alf4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
X1, Y1, Z1, X2, Y2, Z2,  X3, Y3,  Z3, X4, Y4, Z4, n1, n2, n3,  N4 , cont);

fv3b =fv3b +1;


}


}else{

if (   (fv4a > 0)  &&  (fv4b== 0)  &&  (fv4c > 0)   ){
// abc | n0m
fv3a=fv3a+1;
fv3a=fv4a-1;

MaskD_v fx2 = ( vecCore::math::Abs( X3 - X4 ) > 0.0 );
if( MaskFull( fx2 )  ){

fv3a= fv3a - 1;

JepV = JepV + ( X3 - X4 )*selectorV(alf_1,  alf_2, alf_3, Alf4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
X1, Y1, Z1, X2, Y2, Z2,  X3, Y3,  Z3, X4, Y4, Z4, n1, n2, n3,  N4 , cont);


fv3a =fv3a+1;

}

}else{

fv3a=fv3a+1;
fv3a=fv4a-1;

MaskD_v fx2 = ( vecCore::math::Abs( X3 - X4 ) > 0.0 );
if( MaskFull( fx2 )  ){

fv3a=fv3a-1;

JepV = JepV + ( X3 - X4 )*selectorV(alf_1,  alf_2, alf_3, Alf4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
X1, Y1, Z1, X2, Y2, Z2,  X3, Y3,  Z3, X4, Y4, Z4, n1, n2, n3,  N4 , cont);

fv3a = fv3a+1;


}

}

}

}

}

}

}

}

return JepV;
}



vector<vector<vector<vector<G4double>>>> Electron_electron_repulsion::pairingV(Molecule& Mol ){

vector<Atom *>::iterator atr_left_main;
vector<Atom *>::iterator atr_left_row;
vector<Atom *>::iterator atr_right;
vector<Atom *>::iterator atr_right_row;

vector<Residue *>::iterator res_left_main;
vector<Residue *>::iterator res_left_row;
vector<Residue *>::iterator right;
vector<Residue *>::iterator right_row;

vector<G4double> S (Mol.no_orbits, 0.0);
vector<vector<G4double>> Con ( Mol.no_orbits, 
vector<G4double>( Mol.no_orbits, 0.0   )  );


vector<vector<vector<G4double>>> Ven ( Mol.no_orbits, 
vector<vector<G4double>> ( Mol.no_orbits, vector<G4double>( Mol.no_orbits, 0.0 )  )  );

vector<vector<vector<vector<G4double>>>> Rep ( Mol.no_orbits,
vector<vector<vector<G4double>>> ( Mol.no_orbits, 
vector<vector<G4double>> ( 
Mol.no_orbits, vector<G4double>( Mol.no_orbits, 0.0 )  )  ) );


G4double buf_alp_1=0;
G4double buf_alp_2=0;
G4double buf_alp_3=0;

G4double buf_C1=0;
G4double buf_C2=0;
G4double buf_C3=0;

G4double coffs=0;
G4double contcofs_inter=0;

G4int cont = 0;
G4int orbital_principal=0;
G4int orbital_secundario=0;
G4int orbital_terciario=0;
G4int orbital_cuaternario=0;
Double_v JJ;


for(res_left_main = Mol.Residuos_cadena.begin();
    res_left_main != Mol.Residuos_cadena.end() ; ++res_left_main){

for(atr_left_main = (*res_left_main)->Lista_atoms.begin();
    atr_left_main != (*res_left_main)->Lista_atoms.end(); ++ atr_left_main){

for((*atr_left_main)->shell1 = (*atr_left_main)->orbitals.begin();
    (*atr_left_main)->shell1 != (*atr_left_main)->orbitals.end(); ++ (*atr_left_main)->shell1){

orbital_secundario=0;

for(res_left_row = Mol.Residuos_cadena.begin();
    res_left_row != Mol.Residuos_cadena.end() ; ++res_left_row){

for(atr_left_row  = (*res_left_row)->Lista_atoms.begin();
    atr_left_row != (*res_left_row)->Lista_atoms.end(); ++ atr_left_row){

for((*atr_left_row)->shell2 = (*atr_left_row)->orbitals.begin();
    (*atr_left_row)->shell2 != (*atr_left_row)->orbitals.end(); ++ (*atr_left_row)->shell2){


if(   orbital_principal  <= orbital_secundario  ){

orbital_terciario=0;

for(right = Mol.Residuos_cadena.begin();
    right != Mol.Residuos_cadena.end() ; ++right){

for(atr_right = (*right)->Lista_atoms.begin();
    atr_right != (*right)->Lista_atoms.end(); ++  atr_right){

for((*atr_right)->shell3 = (*atr_right)->orbitals.begin();
    (*atr_right)->shell3 != (*atr_right)->orbitals.end(); ++ (*atr_right)->shell3){

orbital_cuaternario=0;

  for(right_row = Mol.Residuos_cadena.begin();
      right_row != Mol.Residuos_cadena.end() ; ++right_row){

  for(atr_right_row  = (*right_row)->Lista_atoms.begin();
      atr_right_row != (*right_row)->Lista_atoms.end(); ++ atr_right_row){

  for((*atr_right_row)->shell4 = (*atr_right_row)->orbitals.begin();
      (*atr_right_row)->shell4 != (*atr_right_row)->orbitals.end(); ++ (*atr_right_row)->shell4){

if(   orbital_terciario  <= orbital_cuaternario  ){

if ( orbital_principal + orbital_secundario <=   orbital_terciario + orbital_cuaternario ){

  for (int i = 0; i < 3; i++) {

  buf_alp_1 = Get( (*(*atr_left_main)->shell1).orbitalsAlpha , i);
  buf_C1    = Get( (*(*atr_left_main)->shell1).orbitalsConst , i);

  for (int j = 0; j < 3; j++) {

  buf_alp_2 = Get( (*(*atr_left_row)->shell2).orbitalsAlpha , j);
  buf_C2    = Get( (*(*atr_left_row)->shell2).orbitalsConst , j);

  for (int k = 0; k < 3; k++) {

  buf_alp_3 = Get ( (*(*atr_right)->shell3).orbitalsAlpha , k);
  buf_C3    = Get ( (*(*atr_right)->shell3).orbitalsConst , k);

JJ = Electron_electron_repulsion::selectorV(buf_alp_1, buf_alp_2, buf_alp_3, (*(*atr_right_row)->shell4).orbitalsAlpha,
  (*(*atr_left_main)->shell1).l[0] , (*(*atr_left_main)->shell1).l[1] , (*(*atr_left_main)->shell1).l[2] , 
  (*(*atr_left_row)->shell2).l[0], (*(*atr_left_row)->shell2).l[1], (*(*atr_left_row)->shell2).l[2],
  (*(*atr_right)->shell3).l[0],  (*(*atr_right)->shell3).l[1],  (*(*atr_right)->shell3).l[2],  
  (*(*atr_right_row)->shell4).l[0],  (*(*atr_right_row)->shell4).l[1],  (*(*atr_right_row)->shell4).l[2],
  (*atr_left_main)->fX, (*atr_left_main)->fY,(*atr_left_main)->fZ,
  (*atr_left_row)->fX,(*atr_left_row)->fY, (*atr_left_row)->fZ,
  (*atr_right)->fX, (*atr_right)->fY, (*atr_right)->fZ,
  (*atr_right_row)->fX, (*atr_right_row)->fY, (*atr_right_row)->fZ,
Get(  (*(*atr_left_main)->shell1).orbitalsNormals , i), Get( (*(*atr_left_row)->shell2).orbitalsNormals , j ),
Get(  (*(*atr_right)->shell3).orbitalsNormals , k ) ,  (*(*atr_right_row)->shell4).orbitalsNormals , cont);


for (int l = 0; l < kVecLenD - 1; ++l) {

coffs = buf_C1*buf_C2*buf_C3*Get( (*(*atr_right_row)->shell4).orbitalsConst , l );
contcofs_inter = contcofs_inter + Get(JJ, l)*coffs;

}


}
}
}


Rep[orbital_principal][orbital_secundario][orbital_terciario][orbital_cuaternario] = contcofs_inter;
Rep[orbital_terciario][orbital_cuaternario][orbital_principal][orbital_secundario] = contcofs_inter;
Rep[orbital_secundario][orbital_principal][orbital_cuaternario][orbital_terciario] = contcofs_inter;
Rep[orbital_cuaternario][orbital_terciario][orbital_secundario][orbital_principal] = contcofs_inter;
Rep[orbital_secundario][orbital_principal][orbital_terciario][orbital_cuaternario] = contcofs_inter;
Rep[orbital_cuaternario][orbital_terciario][orbital_principal][orbital_secundario] = contcofs_inter;
Rep[orbital_principal][orbital_secundario][orbital_cuaternario][orbital_terciario] = contcofs_inter;
Rep[orbital_terciario][orbital_cuaternario][orbital_secundario][orbital_principal] = contcofs_inter;

contcofs_inter=0;
cont = cont + 1 ;

}else{


}


}else{


}

orbital_cuaternario++;

}
}
}


orbital_terciario++;
}
}
}

}else{

}
orbital_secundario++;

}
}
}


orbital_principal++;
}
}
}

return Rep;
}


