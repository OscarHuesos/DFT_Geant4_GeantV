
#include "repulsion_energy.hh"

using namespace std;


Electron_electron_repulsion::Electron_electron_repulsion(){}

// llama a0m00 o a0mc0
G4double Electron_electron_repulsion::repulsion_a0mc0( G4double alf_1, G4double alf_2, G4double alf_3,
G4double alf_4, G4int L10, G4int L11, G4int L12, G4int L20, G4int L21, G4int L22,
G4double x1, G4double y1, G4double z1,
G4double x2, G4double y2, G4double z2,
G4double x3, G4double y3, G4double z3,
G4double x4, G4double y4, G4double z4,
G4double N1, G4double N2, G4double N3, 
G4double N4, G4int cc ){


G4double Norminter = 0;
G4int L = L10 + L11 + L12 + L20 + L21 + L22 ;

G4double gam1 = alf_1 + alf_2;
G4double gam2 = alf_3 + alf_4;

G4double Px = (alf_1*x1 + alf_2*x2)/gam1;
G4double Py = (alf_1*y1 + alf_2*y2)/gam1;
G4double Pz = (alf_1*z1 + alf_2*z2)/gam1;

G4double Qx = (alf_3*x3 + alf_4*x4)/gam2;
G4double Qy = (alf_3*y3 + alf_4*y4)/gam2;
G4double Qz = (alf_3*z3 + alf_4*z4)/gam2;
G4double G = gam1 + gam2;

G4double Wx = (gam1*Px + gam2*Qx)/(G);
G4double Wy = (gam1*Py + gam2*Qy)/(G);
G4double Wz = (gam1*Pz + gam2*Qz)/(G);

G4double D1 = pow(x1-x2,2) + pow(y1-y2,2) + pow(z1-z2,2);
G4double D2 = pow(x3-x4,2) + pow(y3-y4,2) + pow(z3-z4,2);

G4double Kp = G4Exp(-D1*(alf_1*alf_2)/gam1);
G4double Kq = G4Exp(-D2*(alf_3*alf_4)/gam2);
// 2*pow(pi,5.0/2))
G4double Rx = Px-Qx;
G4double Ry = Py-Qy;
G4double Rz = Pz-Qz;
G4double R  = Rx*Rx+Ry*Ry+Rz*Rz;
G4double om = gam1*gam2/(G);
G4double T = om*R;
G4double norms=N1*N2*N3*N4;
G4double f = 0;

// poly false = rep
bool polyn = false;
// the old is selec = false;
bool selector = true;

Operators* fOpx = new Operators();

// check if 00mc0:
if (  (L10 == 0)  &&  (L11 == 0) &&  (L12 == 0)  ){

// check base case 00m00:
if ( (L20 == 0) &&  (L21 == 0) &&  (L22 == 0)  )  {

G4double* ooomooo = new G4double[L+1];

if(T > 100){
ooomooo[0] = sqrt(pi/(4*T));

if (L > 0){
for (int i = 0; i < L ; i++) {
//f = 
ooomooo[i+1] = (2.0*i+1)*ooomooo[i]/(2.0*T);
}

}

}else{

ooomooo[L]  = fOpx->Gamma_aprox( om, Rx, Ry, Rz, 1, 1, 1,
1, 1, 1, L10, L11, L12, L20, L21, L22,
gam1, gam2, polyn, selector );

delete fOpx;

if (L > 0){
for(int i = L ; i > 0; i--) {

ooomooo[i-1] =  (2.0*T*ooomooo[i] + G4Exp(-T) )/(2.0*(i-1)+1);
}

}

} //else T > 100

for (int j = 0; j < L + 1 ; j++) {


ooomooo[j]  =  (ooomooo[j]*norms*pow(2,1)*pow(pi,5.0/2.0)*Kp*Kq*
pow(gam1,-1)*pow(gam2,-1) )/ ( sqrt(G) );

}


Norminter = ooomooo[0];
// check de moment todo al valor l = 0;

delete [] ooomooo;

return Norminter;

}else{  

f=0;  
G4double ****arrc;
arrc = new G4double***[L20 + 1]; 
for(int i = 0; i < L20 + 1 ; i++){
arrc[i] = new G4double**[ L21 + 1 ];
 for(int j=0; j <  L21 + 1 ; j++) {
arrc[i][j] = new G4double*[L22 + 1];  
 for(int k=0; k <  L22 + 1 ; k++) {
arrc[i][j][k] = new G4double[L + 1];
}
}
}

if(T > 100){

//f= 
arrc[0][0][0][0] = sqrt(pi/(4*T));

if (L > 0){

for (int i = 0; i < L ; i++) {
//f = 
arrc[0][0][0][i+1] = (2.0*i+1)*arrc[0][0][0][i]/(2.0*T);
}

}

}else{

arrc[0][0][0][L] = fOpx->Gamma_aprox( om, Rx, Ry, Rz, 1, 1, 1,
1, 1, 1, L10, L11, L12, L20, L21, L22,
gam1, gam2, polyn, selector );

delete fOpx;

if (L > 0){
for(int i = L ; i > 0; i--) {
arrc[0][0][0][i-1] =  (2.0*T*arrc[0][0][0][i] + G4Exp(-T) )/(2.0*(i-1)+1);
}
}
} 

for (int j = 0; j < L + 1 ; j++) {

arrc[0][0][0][j]  =  (arrc[0][0][0][j]*norms*pow(2,1)*pow(pi,5.0/2.0)*Kp*Kq*
pow(gam1,-1)*pow(gam2,-1) )/ ( sqrt(G) );

}


if (  L20 > 0  ){


for(int i = 1; i < L20 + 1 ; i++){
 for (int l = 0; l < L  ; l++) {

f = (Qx - x3)*arrc[i-1][0][0][l] + (Wx - Qx)*arrc[i-1][0][0][l+1] ; 

if( (i-1) > 0  ){
f = f + ( (i-1)/(2.0*gam2))*( arrc[i-2][0][0][l] - (om/gam2)*arrc[i-2][0][0][l+1] ) ;
}
arrc[i][0][0][l] = f;
}
}
}

if (  L21 > 0   ){
for(int j = 1; j < L21 + 1 ; j++){
 for (int l = 0; l < L  ; l++) {
  //just los Q y W:
f = (Qy - y3)*arrc[0][j-1][0][l] + (Wy - Qy)*arrc[0][j-1][0][l+1];
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*( arrc[0][j-2][0][l] - (om/gam2)*arrc[0][j-2][0][l+1] ) ;
}

arrc[0][j][0][l] = f;
}
}
}

if (  L22 > 0   ){
for(int k = 1; k < L22 + 1 ; k++){
 for (int l = 0; l < L  ; l++) {
  // 000 | 00n
  //positve sign f = (Qz - z3)*arrc[0][0][k-1][l] - (Wz - Qz)*arrc[0][0][k-1][l+1]  ;
f = (Qz - z3)*arrc[0][0][k-1][l] + (Wz - Qz)*arrc[0][0][k-1][l+1]  ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*(arrc[0][0][k-2][l] - (om/gam2)*arrc[0][0][k-2][l+1] ) ;
}
arrc[0][0][k][l] = f;
}
}
}

if(   (L20 >= 1) && (  L21 >= 1 )  ) {
f=0;
//000 | nm0
for(int jj = 1; jj < L21 + 1 ; jj++){

  for(int ii = 1; ii < L20 + 1 ; ii++){
for (int l = 0; l < L  ; l++) {

f = (Qy - y3)*arrc[ii][jj-1][0][l] + (Wy - Qy)*arrc[ii][jj-1][0][l+1]  ;
if( (jj-1) > 0  ){
f = f + ( (jj-1)/(2.0*gam2))*( arrc[ii][jj-2][0][l] - (om/gam2)*arrc[ii][jj-2][0][l+1] ) ;
}

arrc[ii][jj][0][l] = f;
}
}
}
}

if(   (L20 >= 1) && (  L22 >= 1 )  ) {
f=0;
// 000 | n0m
for(int kk = 1; kk < L22 + 1  ; kk++){
  for(int ii = 1; ii < L20 + 1; ii++){
for (int l = 0; l < L  ; l++) {
    //primer with sign positve f = (Qz - z3)*arrc[ii][0][kk-1][l] - (Wz - Qz)*arrc[ii][0][kk-1][l+1]  ;
f = (Qz - z3)*arrc[ii][0][kk-1][l] + (Wz - Qz)*arrc[ii][0][kk-1][l+1]  ;
if( (kk-1) > 0  ){
f = f + ( (kk-1)/(2.0*gam2))*( arrc[ii][0][kk-2][l] - (om/gam2)*arrc[ii][0][kk-2][l+1] ) ;
}
arrc[ii][0][kk][l] = f;
}
}
}
}

if(   (L21 >= 1) && (  L22 >= 1 )  ) {
f=0;
// 000 | 0nm
for(int kk = 1; kk < L22 + 1  ; kk++){
  for(int jj = 1; jj < L21 + 1; jj++){
for (int l = 0; l < L  ; l++) {
f = (Qz - z3)*arrc[0][jj][kk-1][l] + (Wz - Qz)*arrc[0][jj][kk-1][l+1]  ;
if( (kk-1) > 0  ){
f = f + ( (kk-1)/(2.0*gam2))*( arrc[0][jj][kk-2][l] - (om/gam2)*arrc[0][jj][kk-2][l+1] ) ;
}
arrc[0][jj][kk][l] = f;
}
}
}
}

if(   (L20 >= 1) && (  L21 >= 1 )  &&  (L22 >= 1 )  ) {
f=0;
// 000 | nmk
//if ll22 = 2, will be 1 without suma :  1; 1 < 2 one iteration
// if ll22 = 2, will be 1,2 with suma :  1; 1 < 3 two iterations
for(int kk = 1; kk < L22 +1 ; kk++){
  for(int jj = 1; jj < L21 +1 ; jj++){
   for(int ii = 1; ii < L20 + 1 ; ii++){
for (int l = 0; l < L ; l++) {

f = (Qz - z3)*arrc[ii][jj][kk-1][l] + (Wz - Qz)*arrc[ii][jj][kk-1][l+1]  ;
if( (kk-1) > 0  ){
f = f + ( (kk-1)/(2.0*gam2))*(arrc[ii][jj][kk-2][l] - (om/gam2)*arrc[ii][jj][kk-2][l+1] ) ;
}
arrc[ii][jj][kk][l] = f;

}
}
}
}
}

Norminter = arrc[L20][L21][L22][0];

for(int i = 0; i < (L20 + 1); i++){ 

  for(int j = 0; j < (L21 + 1); j++){   
    for(int n = 0; n < (L22 +1); n++)  {  

      delete[] arrc[i][j][n];
        }
        delete [] arrc[i][j];

      }
     delete [] arrc[i];

    }


delete[] arrc;


return Norminter;

}
//  a0mc0:
}else{ 


//check if a0m00:
if ( (L20 == 0) &&  (L21 == 0) &&  (L22 == 0)  )  {

f=0;  
G4double ****arr0;
arr0 = new G4double***[L10 + 1]; 
for(int i = 0; i < L10 + 1 ; i++){
arr0[i] = new G4double**[ L11 + 1 ];
 for(int j=0; j <  L11 + 1 ; j++) {
arr0[i][j] = new G4double*[L12 + 1];  
 for(int k=0; k <  L12 +1 ; k++) {
arr0[i][j][k] = new G4double[L+1];
}
}
}


if(T > 100){
//  printf("T  >  100:  %f \n", T);

arr0[0][0][0][0] =  sqrt(pi/(4*T));
//ooomooo[0]= f;

if (L > 0){

for (int i = 0; i < L ; i++) {

arr0[0][0][0][i+1] = (2.0*i+1)*arr0[0][0][0][i]/(2.0*T);
}

}

}else{

arr0[0][0][0][L] = fOpx->Gamma_aprox( om, Rx, Ry, Rz, 1, 1, 1,
1, 1, 1, L10, L11, L12, L20, L21, L22,
gam1, gam2, polyn, selector );

delete fOpx;

if (L > 0){
for(int i = L ; i > 0; i--) {
arr0[0][0][0][i-1] =  (2.0*T*arr0[0][0][0][i] + G4Exp(-T) )/(2.0*(i-1)+1);
}
}
} //else T > 100

for (int j = 0; j < L + 1 ; j++) {

arr0[0][0][0][j]  =  (arr0[0][0][0][j]*norms*pow(2,1)*pow(pi,5.0/2.0)*Kp*Kq*
pow(gam1,-1)*pow(gam2,-1) )/ ( sqrt(G) );
////////////////////////////////////////////////////////////////////////////
}


if (  L10 > 0   ){

// if L_x(n) = 5 ->  0, 1, 2, 3, 4
// st -> k = 1,    1, 2, 3, 4   k < 5 + 1 = 6 -> 1 < 6 -> 1,2,3,4,5

for(int k = 1; k < L10 + 1 ; k++){
for (int l = 0; l < L  ; l++) {
 
f = (Px - x1)*arr0[k-1][0][0][l] + (Wx - Px)*arr0[k-1][0][0][l+1] ; 

if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam1))*( arr0[k-2][0][0][l] - (om/gam1)*arr0[k-2][0][0][l+1] ) ;
}

arr0[k][0][0][l] = f;

}
}
}

if (  L11 > 0   ){

for(int k = 1; k < L11 + 1 ; k++){
for (int l = 0; l < L  ; l++) {
f = (Py - y1)*arr0[0][k-1][0][l] + (Wy - Py)*arr0[0][k-1][0][l+1]  ;

if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam1))*( arr0[0][k-2][0][l] - (om/gam1)*arr0[0][k-2][0][l+1] ) ;
}

arr0[0][k][0][l] = f;

}
}
}

if (  L12 > 0   ){

for(int k = 1; k < L12 + 1 ; k++){
for (int l = 0; l < L  ; l++) {
f = (Pz - z1)*arr0[0][0][k-1][l] + (Wz - Pz)*arr0[0][0][k-1][l+1]  ;

if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam1))*( arr0[0][0][k-2][l] - (om/gam1)*arr0[0][0][k-2][l+1] ) ;
}

arr0[0][0][k][l] = f;

}
}
}

if(   (L10 >= 1) && (  L11  >= 1 )  ) {
f=0;


for(int jj = 1; jj < L11 + 1 ; jj++){
  for(int ii = 1; ii < L10 + 1 ; ii++){
for (int l = 0; l < L  ; l++) {

f = (Py - y1)*arr0[ii][jj-1][0][l] + (Wy - Py)*arr0[ii][jj-1][0][l+1]  ;
if( (jj-1) > 0  ){
f = f + ( (jj-1)/(2.0*gam1))*( arr0[ii][jj-2][0][l] - (om/gam1)*arr0[ii][jj-2][0][l+1] ) ;
}

arr0[ii][jj][0][l] = f;


}
}
}
}


if(   (L10 >= 1) && (  L12 >= 1 )  ) {
f=0;

for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
for (int l = 0; l < L  ; l++) {

f = (Pz - z1)*arr0[ii][0][kk-1][l] + (Wz - Pz)*arr0[ii][0][kk-1][l+1]  ;
if( (kk-1) > 0  ){
f = f + ( (kk-1)/(2.0*gam1))*( 
arr0[ii][0][kk-2][l] - (om/gam1)*arr0[ii][0][kk-2][l+1] ) ;
}
arr0[ii][0][kk][l] = f;
}
}
}
}


if(   (L11 >= 1) && (  L12 >= 1 )  ) {
f=0;


for(int kk = 1; kk < L12 + 1  ; kk++){
  for(int jj = 1; jj < L11 + 1; jj++){
for (int l = 0; l < L  ; l++) {

f = (Pz - z1)*arr0[0][jj][kk-1][l] + (Wz - Pz)*arr0[0][jj][kk-1][l+1]  ;
if( (kk-1) > 0  ){
f = f + ( (kk-1)/(2.0*gam1))*( arr0[0][jj][kk-2][l] - (om/gam1)*arr0[0][jj][kk-2][l+1] ) ;
}
arr0[0][jj][kk][l] = f;
}
}
}
}

if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {

f=0;

for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
for (int l = 0; l < L ; l++) {

f = (Pz - z1)*arr0[ii][jj][kk-1][l] + (Wz - Pz)*arr0[ii][jj][kk-1][l+1]  ;

if( (kk-1) > 0  ){
f = f + ( (kk-1)/(2.0*gam1))*(arr0[ii][jj][kk-2][l] - (om/gam1)*arr0[ii][jj][kk-2][l+1] ) ;
}
arr0[ii][jj][kk][l] = f;

}
}
}
}

}

Norminter = arr0[L10][L11][L12][0];


for (int i = 0; i < L10 + 1; i++) {
    for (int j = 0; j < L11 + 1; j++) {
        for (int k = 0; k < L12 + 1; k++) {
            delete[] arr0[i][j][k];
        }
        delete[] arr0[i][j];
    }
    delete[] arr0[i];
}
delete[] arr0;

return Norminter;


}else{  


if (  (L20 >= 1)  &&  (L21 == 0) &&  (L22 == 0)    ){
// tensor de 4
f = 0;

G4double *****arr4;
arr4 = new G4double****[L20 + 1]; 
for(int i = 0; i < L20+1 ; i++){
arr4[i] = new G4double***[ L10 + 1];
for(int j = 0; j < L10+1 ; j++){
arr4[i][j] = new G4double**[ L11 + 1];
for(int k = 0; k <  L11 + 1 ; k++) {
arr4[i][j][k] =  new G4double*[ L12 + 1];
for(int l = 0; l <  L12 + 1 ; l++) {
arr4[i][j][k][l] =  new G4double[L + 1];
}
}
}
}

if(T > 100){
arr4[0][0][0][0][0] =  sqrt(pi/(4*T));

if (L > 0){

for (int i = 0; i < L ; i++) {
 
arr4[0][0][0][0][i+1] = (2.0*i+1)*arr4[0][0][0][0][i]/(2.0*T);
}

}

}else{

arr4[0][0][0][0][L] = fOpx->Gamma_aprox( om, Rx, Ry, Rz, 1, 1, 1,
1, 1, 1, L10, L11, L12, L20, L21, L22,
gam1, gam2, polyn, selector );

delete fOpx;

if (L > 0){
for(int i = L ; i > 0; i--) {
arr4[0][0][0][0][i-1] =  (2.0*T*arr4[0][0][0][0][i] + G4Exp(-T) )/(2.0*(i-1)+1);
}
}
} //else T > 100


for (int j = 0; j < L + 1 ; j++) {
arr4[0][0][0][0][j]  =  (arr4[0][0][0][0][j]*norms*pow(2,1)*pow(pi,5.0/2.0)*Kp*Kq*
pow(gam1,-1)*pow(gam2,-1) )/ ( sqrt(G) );
}


// iterar los l1 (A1,A2,A3)
if (  L10 > 0   ){
for(int i = 1; i < L10 + 1 ; i++){
 for (int l = 0; l < L  ; l++) {

f = (Px - x1)*arr4[0][i-1][0][0][l] + (Wx - Px)*arr4[0][i-1][0][0][l+1] ; 

if( (i-1) > 0  ){
f = f + ( (i-1)/(2.0*gam1))*( arr4[0][i-2][0][0][l] - (om/gam1)*arr4[0][i-2][0][0][l+1] ) ;
}

arr4[0][i][0][0][l] = f;

}
}
}

if (  L11 > 0   ){
for(int j = 1; j < L11 + 1 ; j++){
 for (int l = 0; l < L  ; l++) {

f = (Py - y1)*arr4[0][0][j-1][0][l] + (Wy - Py)*arr4[0][0][j-1][0][l+1];

if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam1))*( arr4[0][0][j-2][0][l] - (om/gam1)*arr4[0][0][j-2][0][l+1] ) ;
}

arr4[0][0][j][0][l] = f;

}
}
}

if (  L12 > 0   ){
for(int k = 1; k < L12 + 1 ; k++){
 for (int l = 0; l < L  ; l++) {

f = (Pz - z1)*arr4[0][0][0][k-1][l] + (Wz - Pz)*arr4[0][0][0][k-1][l+1]  ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam1))*(arr4[0][0][0][k-2][l] - (om/gam1)*arr4[0][0][0][k-2][l+1] ) ;
}
arr4[0][0][0][k][l] = f;
}
}
}


if(   (L10 >= 1) && (  L11 >= 1 )  ) {
f=0;
for(int jj = 1; jj < L11 + 1 ; jj++){
  for(int ii = 1; ii < L10 + 1 ; ii++){
for (int l = 0; l < L  ; l++) {
 
f = (Py - y1)*arr4[0][ii][jj-1][0][l] + (Wy - Py)*arr4[0][ii][jj-1][0][l+1]  ;
if( (jj-1) > 0  ){
f = f + ( (jj-1)/(2.0*gam1))*( arr4[0][ii][jj-2][0][l] - (om/gam1)*arr4[0][ii][jj-2][0][l+1] ) ;
}
arr4[0][ii][jj][0][l] = f;
}
}
}
}


if(   (L10 >= 1) && (  L12 >= 1 )  ) {
f=0;
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
for (int l = 0; l < L  ; l++) {

f = (Pz - z1)*arr4[0][ii][0][kk-1][l] + (Wz - Pz)*arr4[0][ii][0][kk-1][l+1]  ;
if( (kk-1) > 0  ){
f = f + ( (kk-1)/(2.0*gam1))*( 
arr4[0][ii][0][kk-2][l] - (om/gam1)*arr4[0][ii][0][kk-2][l+1] ) ;
}
arr4[0][ii][0][kk][l] = f;
}
}
}
}


if(   (L11 >= 1) && (  L12 >= 1 )  ) {
f=0;
// 0nm
for(int kk = 1; kk < L12 + 1  ; kk++){
  for(int jj = 1; jj < L11 + 1; jj++){
for (int l = 0; l < L  ; l++) {

f = (Pz - z1)*arr4[0][0][jj][kk-1][l] + (Wz - Pz)*arr4[0][0][jj][kk-1][l+1]  ;
if( (kk-1) > 0  ){
f = f + ( (kk-1)/(2.0*gam1))*( arr4[0][0][jj][kk-2][l] - (om/gam1)*arr4[0][0][jj][kk-2][l+1] ) ;
}
arr4[0][0][jj][kk][l] = f;
}
}
}
}

if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {

f=0;

for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
for (int l = 0; l < L ; l++) {

f = (Pz - z1)*arr4[0][ii][jj][kk-1][l] + (Wz - Pz)*arr4[0][ii][jj][kk-1][l+1]  ;
if( (kk-1) > 0  ){
f = f + ( (kk-1)/(2.0*gam1))*(arr4[0][ii][jj][kk-2][l] - (om/gam1)*arr4[0][ii][jj][kk-2][l+1] ) ;
}
arr4[0][ii][jj][kk][l] = f;

}
}
}
}

}



for(int i = 1; i < L20 + 1 ; i++){
f = 0;
// caso base  000 | i00
 for (int l = 0; l < L  ; l++) {
f = (Qx - x3)*arr4[i-1][0][0][0][l] + (Wx - Qx)*arr4[i-1][0][0][0][l+1] ;
if( (i-1) > 0  ){
f = f + ( (i-1)/(2.0*gam2))*( arr4[i-2][0][0][0][l] - (om/gam2)*arr4[i-2][0][0][0][l+1] ) ;
}
arr4[i][0][0][0][l] = f;
}


if (  L10 > 0   ){
for(int ii = 1; ii < L10 + 1 ; ii++){
 for (int l = 0; l < L  ; l++) {

f = (Qx - x3)*arr4[i-1][ii][0][0][l] + (Wx - Qx)*arr4[i-1][ii][0][0][l+1] ;
if( (i-1) > 0  ){
f = f + ( (i-1)/(2.0*gam2))*( arr4[i-2][ii][0][0][l] - (om/gam2)*arr4[i-2][ii][0][0][l+1] ) ;
}

f = f + ( (ii)/( 2.0*(gam1 + gam2) ) )*arr4[i-1][ii-1][0][0][l+1] ;

arr4[i][ii][0][0][l] = f;
}
}
}


f = 0;
if (  L11 > 0   ){
for(int jj = 1; jj < L11 + 1 ; jj++){
 for (int l = 0; l < L  ; l++) {

 f = (Qx - x3)*arr4[i-1][0][jj][0][l] + (Wx - Qx)*arr4[i-1][0][jj][0][l+1] ;
if( (i-1) > 0  ){
f = f + ( (i-1)/(2.0*gam2))*( arr4[i-2][0][jj][0][l] - (om/gam2)*arr4[i-2][0][jj][0][l+1] ) ;
}

arr4[i][0][jj][0][l] = f;

}
}
}

f = 0;
if (  L12 > 0   ){
for(int kk = 1; kk < L12 + 1 ; kk++){
 for (int l = 0; l < L  ; l++) {
  f = (Qx - x3)*arr4[i-1][0][0][kk][l] + (Wx - Qx)*arr4[i-1][0][0][kk][l+1] ;

if( (i-1) > 0  ){
f = f + ( (i-1)/(2.0*gam2))*( arr4[i-2][0][0][kk][l] - (om/gam2)*arr4[i-2][0][0][kk][l+1] ) ;
}

arr4[i][0][0][kk][l] = f;

}
}
}

if(   (L10 >= 1) && (  L11 >= 1 )  ) {
f=0;

// i | nm0
for(int jj = 1; jj < L11 + 1 ; jj++){
//jj = 0 -> n00 -> previously, as well as 0m0
  for(int ii = 1; ii < L10 + 1 ; ii++){
    for (int l = 0; l < L  ; l++) { 

f = (Qx - x3)*arr4[i-1][ii][jj][0][l] + (Wx - Qx)*arr4[i-1][ii][jj][0][l+1] ;
if( (i-1) > 0  ){
// at least  i = 2:
f = f + ( (i-1)/(2.0*gam2))*( arr4[i-2][ii][jj][0][l] - (om/gam2)*arr4[i-2][ii][jj][0][l+1] ) ;
}
f = f + ( (ii)/( 2.0*(gam1 + gam2) ) )*arr4[i-1][ii-1][jj][0][l+1] ;
arr4[i][ii][jj][0][l] = f;
}
}
}
}

if(   (L10 >= 1) && (  L12 >= 1 )  ) {
f=0;
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
for (int l = 0; l < L  ; l++) {

f = (Qx - x3)*arr4[i-1][ii][0][kk][l] + (Wx - Qx)*arr4[i-1][ii][0][kk][l+1] ;
if( (i-1) > 0  ){
f = f + ( (i-1)/(2.0*gam2))*( arr4[i-2][ii][0][kk][l] - (om/gam2)*arr4[i-2][ii][0][kk][l+1] ) ;
}
f = f + ( (ii)/( 2.0*(gam1 + gam2) ) )*arr4[i-1][ii-1][0][kk][l+1] ;
arr4[i][ii][0][kk][l] = f;
}
}
}
}


if(   (L11 >= 1) && (  L12 >= 1 )  ) {
f=0;
// 0nm | i00
for(int kk = 1; kk < L12 + 1  ; kk++){
  for(int jj = 1; jj < L11 + 1; jj++){
for (int l = 0; l < L  ; l++) {
// f = (Qx - x3)*arr4[i-1][0][jj][kk][l] - (Wx - Qx)*arr4[i-1][0][jj][kk][l+1]  ;
f = (Qx - x3)*arr4[i-1][0][jj][kk][l] + (Wx - Qx)*arr4[i-1][0][jj][kk][l+1]  ;
if( (i-1) > 0  ){
f = f + ( (i-1)/(2.0*gam2))*( arr4[i-2][0][jj][kk][l] - (om/gam2)*arr4[i-2][0][jj][kk][l+1] ) ;
}

arr4[i][0][jj][kk][l] = f;
}
}
}
}



if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {
f=0;
// i | nmk
for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
for (int l = 0; l < (L - i); l++) {
 // f = (Qx - x3)*arr4[i-1][ii][jj][kk][l] - (Wx - Qx)*arr4[i-1][ii][jj][kk][l+1] ;
f = (Qx - x3)*arr4[i-1][ii][jj][kk][l] + (Wx - Qx)*arr4[i-1][ii][jj][kk][l+1] ;
if( (i-1) > 0  ){
f = f + ( (i-1)/(2.0*gam2))*(arr4[i-2][ii][jj][kk][l] - (om/gam2)*arr4[i-2][ii][jj][kk][l+1] ) ;
}
f = f + ( (ii)/( 2.0*(gam1 + gam2) ) )*arr4[i-1][ii-1][jj][kk][l+1] ;
arr4[i][ii][jj][kk][l] = f;

}
}
}
}
}


}

Norminter = arr4[L20][L10][L11][L12][0];
// check de moment todo al valor l = 0;

for (int i = 0; i < L20 + 1; i++) {
    for (int j = 0; j < L10 + 1; j++) {
        for (int k = 0; k < L11 + 1; k++) {
            for (int l = 0; l < L12 + 1; l++) {
                delete[] arr4[i][j][k][l];
            }
            delete[] arr4[i][j][k];
        }
        delete[] arr4[i][j];
    }
    delete[] arr4[i];
}
delete[] arr4;

return Norminter;


}else{

if (    (L21 >= 1) &&  (L22 == 0)    ){

f = 0;

G4double ******arr5;
arr5 = new G4double*****[L21+1]; 
for(int i = 0; i < L21+1 ; i++){
arr5[i] = new G4double****[ L20 +1 ];
for(int j = 0; j < L20+1 ; j++){
arr5[i][j] = new G4double***[ L10 +1];
for(int k=0; k <  L10 +1 ; k++) {
arr5[i][j][k] = new G4double**[L11 +1 ];  
for(int l = 0; l <  L11 +1 ; l++) {
arr5[i][j][k][l] = new G4double*[L12+1];
for(int m=0; m <  L12 +1 ; m++) {
arr5[i][j][k][l][m] =  new G4double[L+1];
}
}
}
}
}

if(T > 100){

arr5[0][0][0][0][0][0] = sqrt(pi/(4*T));

if (L > 0){
for (int i = 0; i < L ; i++) {
//f = 
arr5[0][0][0][0][0][i+1]= (2.0*i+1)*arr5[0][0][0][0][0][i]/(2.0*T);
}

}

}else{


arr5[0][0][0][0][0][L] = fOpx->Gamma_aprox( om, Rx, Ry, Rz, 1, 1, 1,
1, 1, 1, L10, L11, L12, L20, L21, L22,
gam1, gam2, polyn, selector );

delete fOpx;

if (L > 0){
for(int i = L ; i > 0; i--) {
arr5[0][0][0][0][0][i-1] =  (2.0*T*arr5[0][0][0][0][0][i] + G4Exp(-T) )/(2.0*(i-1)+1);
}
}
} //else T > 100

for (int j = 0; j < L + 1 ; j++) {
arr5[0][0][0][0][0][j]  =  (arr5[0][0][0][0][0][j]*norms*pow(2,1)*pow(pi,5.0/2.0)*Kp*Kq*
pow(gam1,-1)*pow(gam2,-1) )/ ( sqrt(G) );
}


f =0;
if (  L10 > 0   ){
for(int i = 1; i < L10 + 1 ; i++){
 for (int l = 0; l < L  ; l++) {
 // f = (Px - x1)*arr5[0][0][i-1][0][0][l] - (Wx - Px)*arr5[0][0][i-1][0][0][l+1] ; 
f = (Px - x1)*arr5[0][0][i-1][0][0][l] + (Wx - Px)*arr5[0][0][i-1][0][0][l+1] ; 
if( (i-1) > 0  ){
f = f + ( (i-1)/(2.0*gam1) )*( arr5[0][0][i-2][0][0][l] - (om/gam1)*arr5[0][0][i-2][0][0][l+1] ) ;
}
arr5[0][0][i][0][0][l] = f;
}
}
}

if (  L11 > 0   ){
for(int j = 1; j < L11 + 1 ; j++){
 for (int l = 0; l < L  ; l++) {
 // f = (Py - y1)*arr5[0][0][0][j-1][0][l] - (Wy - Py)*arr5[0][0][0][j-1][0][l+1];
f = (Py - y1)*arr5[0][0][0][j-1][0][l] + (Wy - Py)*arr5[0][0][0][j-1][0][l+1];
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam1))*( arr5[0][0][0][j-2][0][l] - (om/gam1)*arr5[0][0][0][j-2][0][l+1] ) ;
}

arr5[0][0][0][j][0][l] = f;

}
}
}

if (  L12 > 0   ){
for(int k = 1; k < L12 + 1 ; k++){
 for (int l = 0; l < L  ; l++) {
//f = (Pz - z1)*arr5[0][0][0][0][k-1][l] - (Wz - Pz)*arr5[0][0][0][0][k-1][l+1]  ;
f = (Pz - z1)*arr5[0][0][0][0][k-1][l] + (Wz - Pz)*arr5[0][0][0][0][k-1][l+1]  ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam1))*(arr5[0][0][0][0][k-2][l] - (om/gam1)*arr5[0][0][0][0][k-2][l+1] ) ;
}
arr5[0][0][0][0][k][l] = f;
}
}
}


if(   (L10 >= 1) && (  L11 >= 1 )  ) {
f=0;
//  ab0 | 000
for(int jj = 1; jj < L11 + 1 ; jj++){
  for(int ii = 1; ii < L10+ 1 ; ii++){
for (int l = 0; l < L  ; l++) {

f = (Py - y1)*arr5[0][0][ii][jj-1][0][l] + (Wy - Py)*arr5[0][0][ii][jj-1][0][l+1]  ;
if( (jj-1) > 0  ){
f = f + ( (jj-1)/(2.0*gam1))*( arr5[0][0][ii][jj-2][0][l] - (om/gam1)*arr5[0][0][ii][jj-2][0][l+1] ) ;
}
arr5[0][0][ii][jj][0][l] = f;
}
}
}
}

if(   (L10 >= 1) && (  L12 >= 1 )  ) {
f=0;
// a0b | 000
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
for (int l = 0; l < L  ; l++) {
//f = (Pz - z1)*arr5[0][0][ii][0][kk-1][l] - (Wz - Pz)*arr5[0][0][ii][0][kk-1][l+1]  ;
f = (Pz - z1)*arr5[0][0][ii][0][kk-1][l] + (Wz - Pz)*arr5[0][0][ii][0][kk-1][l+1]  ;
if( (kk-1) > 0  ){
f = f + ( (kk-1)/(2.0*gam1))*( 
arr5[0][0][ii][0][kk-2][l] - (om/gam1)*arr5[0][0][ii][0][kk-2][l+1] ) ;
}
arr5[0][0][ii][0][kk][l] = f;
}
}
}
}

if(   (L11 >= 1) && (  L12 >= 1 )  ) {
f=0;
// 0ab | 000
for(int kk = 1; kk < L12 + 1  ; kk++){
  for(int jj = 1; jj < L11 + 1; jj++){
for (int l = 0; l < L  ; l++) {
//f = (Pz - z1)*arr5[0][0][0][jj][kk-1][l] - (Wz - Pz)*arr5[0][0][0][jj][kk-1][l+1]  ;
f = (Pz - z1)*arr5[0][0][0][jj][kk-1][l] + (Wz - Pz)*arr5[0][0][0][jj][kk-1][l+1]  ;
if( (kk-1) > 0  ){
f = f + ( (kk-1)/(2.0*gam1))*( arr5[0][0][0][jj][kk-2][l] - (om/gam1)*arr5[0][0][0][jj][kk-2][l+1] ) ;
}

arr5[0][0][0][jj][kk][l] = f;
}
}
}
}


if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {

f=0;
for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
for (int l = 0; l < L ; l++) {

f = (Pz - z1)*arr5[0][0][ii][jj][kk-1][l] + (Wz - Pz)*arr5[0][0][ii][jj][kk-1][l+1]  ;

if( (kk-1) > 0  ){
f = f + ( (kk-1)/(2.0*gam1))*( arr5[0][0][ii][jj][kk-2][l] - (om/gam1)*arr5[0][0][ii][jj][kk-2][l+1] ) ;
}
arr5[0][0][ii][jj][kk][l] = f;

}
}
}
}
}


if ( L20 > 0  ){
f = 0;
for(int i = 1; i < L20 + 1 ; i++){
 for (int l = 0; l < L  ; l++) {

f = (Qx - x3)*arr5[0][i-1][0][0][0][l] + (Wx - Qx)*arr5[0][i-1][0][0][0][l+1] ;
if( (i-1) > 0  ){
f = f + ( (i-1)/(2.0*gam2))*( arr5[0][i-2][0][0][0][l] - (om/gam2)*arr5[0][i-2][0][0][0][l+1] ) ;
}

arr5[0][i][0][0][0][l] = f;
}
}


for(int i = 1; i < L20 + 1 ; i++){


if (  L10 > 0   ){
for(int ii = 1; ii < L10 + 1 ; ii++){
 for (int l = 0; l < L  ; l++) {

f = (Qx - x3)*arr5[0][i-1][ii][0][0][l] + (Wx - Qx)*arr5[0][i-1][ii][0][0][l+1] ;
if( (i-1) > 0  ){
f = f + ( (i-1)/(2.0*gam2))*( arr5[0][i-2][ii][0][0][l] - (om/gam2)*arr5[0][i-2][ii][0][0][l+1] ) ;
}

f = f + ( (ii)/( 2.0*(gam1 + gam2) ) )*arr5[0][i-1][ii-1][0][0][l+1] ;
arr5[0][i][ii][0][0][l] = f;
}
}
}


f = 0;
if (  L11 > 0   ){
for(int jj = 1; jj < L11 + 1 ; jj++){
 for (int l = 0; l < L  ; l++) {

f = (Qx - x3)*arr5[0][i-1][0][jj][0][l] + (Wx - Qx)*arr5[0][i-1][0][jj][0][l+1] ;
if( (i-1) > 0  ){
f = f + ( (i-1)/(2.0*gam2))*( arr5[0][i-2][0][jj][0][l] - (om/gam2)*arr5[0][i-2][0][jj][0][l+1] ) ;
}
arr5[0][i][0][jj][0][l] = f;
}
}
}


f = 0;
if (  L12 > 0   ){
for(int kk = 1; kk < L12 + 1 ; kk++){
 for (int l = 0; l < L  ; l++) {

f = (Qx - x3)*arr5[0][i-1][0][0][kk][l] + (Wx - Qx)*arr5[0][i-1][0][0][kk][l+1] ;
if( (i-1) > 0  ){
f = f + ( (i-1)/(2.0*gam2))*( arr5[0][i-2][0][0][kk][l] - (om/gam2)*arr5[0][i-2][0][0][kk][l+1] ) ;
}
arr5[0][i][0][0][kk][l] = f;

}
}
}


if(   (L10 >= 1) && (  L11 >= 1 )  ) {
f=0;
// ab0 | i00
for(int jj = 1; jj < L11 + 1 ; jj++){
  for(int ii = 1; ii < L10 + 1 ; ii++){
    for (int l = 0; l < L  ; l++) { 

f = (Qx - x3)*arr5[0][i-1][ii][jj][0][l] + (Wx - Qx)*arr5[0][i-1][ii][jj][0][l+1] ;
if( (i-1) > 0  ){
// almenos  i = 2:
f = f + ( (i-1)/(2.0*gam2))*( arr5[0][i-2][ii][jj][0][l] - (om/gam2)*arr5[0][i-2][ii][jj][0][l+1] ) ;
}
f = f + ( (ii)/( 2.0*(gam1 + gam2) ) )*arr5[0][i-1][ii-1][jj][0][l+1] ;
arr5[0][i][ii][jj][0][l] = f;
}
}
}
}



if(   (L10 >= 1) && (  L12 >= 1 )  ) {
f=0;  //  a0b | i00
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
for (int l = 0; l < L  ; l++) {
f = (Qx - x3)*arr5[0][i-1][ii][0][kk][l] + (Wx - Qx)*arr5[0][i-1][ii][0][kk][l+1]  ;
if( (i-1) > 0  ){
f = f + ( (i-1)/(2.0*gam2))*( arr5[0][i-2][ii][0][kk][l] - (om/gam2)*arr5[0][i-2][ii][0][kk][l+1] ) ;
}
f = f + ( (ii)/( 2.0*(gam1 + gam2) ) )*arr5[0][i-1][ii-1][0][kk][l+1] ;
arr5[0][i][ii][0][kk][l] = f;
}
}
}
}


if(   (L11 >= 1) && (  L12 >= 1 )  ) {
f=0;  //  0ab | i00 0 nab | i00
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int jj = 1; jj < L11 +1  ; jj++){
for (int l = 0; l < L  ; l++) {

f = (Qx - x3)*arr5[0][i-1][0][jj][kk][l] + (Wx - Qx)*arr5[0][i-1][0][jj][kk][l+1]  ;
if( (i-1) > 0  ){
f = f + ( (i-1)/(2.0*gam2))*( arr5[0][i-2][0][jj][kk][l] - (om/gam2)*arr5[0][i-2][0][jj][kk][l+1] ) ;
}

arr5[0][i][0][jj][kk][l] = f;
}
}
}
}


// abc | i00
if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {
f=0;
for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
for (int l = 0; l < L ; l++) {
f = (Qx - x3)*arr5[0][i-1][ii][jj][kk][l] + (Wx - Qx)*arr5[0][i-1][ii][jj][kk][l+1] ;
if( (i-1) > 0  ){
f = f + ( (i-1)/(2.0*gam2))*(arr5[0][i-2][ii][jj][kk][l] - (om/gam2)*arr5[0][i-2][ii][jj][kk][l+1] ) ;
}
f = f + ( (ii)/( 2.0*(gam1 + gam2) ) )*arr5[0][i-1][ii-1][jj][kk][l+1] ;
arr5[0][i][ii][jj][kk][l] = f;

}
}
}
}
}
// ciera c1 mayor a 0
}

}
// cierre L20 >0

// [C2][C1][A1][A2][A3]
for(int j = 1; j < L21 + 1 ; j++){
  f = 0;
 for (int l = 0; l < L  ; l++) {
  // el inicial 000 | 0j0
f = (Qy - y3)*arr5[j-1][0][0][0][0][l] + (Wy - Qy)*arr5[j-1][0][0][0][0][l+1] ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*( arr5[j-2][0][0][0][0][l] - (om/gam2)*arr5[j-2][0][0][0][0][l+1] ) ;
}
arr5[j][0][0][0][0][l] = f;
}

// a00 | 0j
if(  L10 > 0  ){  
f= 0;
for(int ii = 1; ii < L10 + 1 ; ii++){
  for (int l = 0; l < L  ; l++) {
//f = (Qy - y3)*arr5[j-1][0][ii][0][0][l] - (Wy - Qy)*arr5[j-1][0][ii][0][0][l+1] ;
 f = (Qy - y3)*arr5[j-1][0][ii][0][0][l] + (Wy - Qy)*arr5[j-1][0][ii][0][0][l+1] ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*( arr5[j-2][0][ii][0][0][l] - (om/gam2)*arr5[j-2][0][ii][0][0][l+1] ) ;
}

arr5[j][0][ii][0][0][l] = f;
}
}
}

f = 0;
if (  L11 > 0   ){
for(int jj = 1; jj < L11 + 1 ; jj++){
 for (int l = 0; l < L  ; l++) {
f = (Qy - y3)*arr5[j-1][0][0][jj][0][l] + (Wy - Qy)*arr5[j-1][0][0][jj][0][l+1] ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*( arr5[j-2][0][0][jj][0][l] - (om/gam2)*arr5[j-2][0][0][jj][0][l+1] ) ;
}
f = f + ( (jj)/( 2.0*(gam1 + gam2) ) )*arr5[j-1][0][0][jj-1][0][l+1] ;
arr5[j][0][0][jj][0][l] = f;

}
}
}

f = 0;
if (  L12 > 0   ){
for(int kk = 1; kk < L12 + 1 ; kk++){
 for (int l = 0; l < L  ; l++) {

f = (Qy - y3)*arr5[j-1][0][0][0][kk][l] + (Wy - Qy)*arr5[j-1][0][0][0][kk][l+1] ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*( arr5[j-2][0][0][0][kk][l] - (om/gam2)*arr5[j-2][0][0][0][kk][l+1] ) ;
}
arr5[j][0][0][0][kk][l] = f;

}
}
}


if(   (L10 >= 1) && (  L11 >= 1 )  ) {
f=0;
// ab0 | 0j0
for(int jj = 1; jj < L11 + 1 ; jj++){
  for(int ii = 1; ii < L10 + 1 ; ii++){ 
   for (int l = 0; l < (L - L10  ); l++) {  
f = (Qy - y3)*arr5[j-1][0][ii][jj][0][l] + (Wy - Qy)*arr5[j-1][0][ii][jj][0][l+1] ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*( arr5[j-2][0][ii][jj][0][l] - (om/gam2)*arr5[j-2][0][ii][jj][0][l+1] ) ;
}
f = f + ( (jj)/( 2.0*(gam1 + gam2) ) )*arr5[j-1][0][ii][jj-1][0][l+1] ;
arr5[j][0][ii][jj][0][l] = f;
}
}
}
}


if(   (L10 >= 1) && (  L12 >= 1 )  ) {
f=0;  
//  a0b | 0j0
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
for (int l = 0; l <  (L - L10  )  ; l++) {

f = (Qy - y3)*arr5[j-1][0][ii][0][kk][l] + (Wy - Qy)*arr5[j-1][0][ii][0][kk][l+1]  ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*( arr5[j-2][0][ii][0][kk][l] - (om/gam2)*arr5[j-2][0][ii][0][kk][l+1] ) ;
}
arr5[j][0][ii][0][kk][l] = f;
}
}
}
}


if(   (L11 >= 1) && (  L12 >= 1 )  ) {
f=0;  //  0ab | 0j0
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int jj = 1; jj < L11 +1  ; jj++){
for (int l = 0; l < (L -L11)  ; l++) {
f = (Qy - y3)*arr5[j-1][0][0][jj][kk][l] + (Wy - Qy)*arr5[j-1][0][0][jj][kk][l+1]  ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*( arr5[j-2][0][0][jj][kk][l] - (om/gam2)*arr5[j-2][0][0][jj][kk][l+1] ) ;
}
f = f + ( (jj)/( 2.0*(gam1 + gam2) ) )*arr5[j-1][0][0][jj-1][kk][l+1] ;
arr5[j][0][0][jj][kk][l] = f;
}
}
}
}

// abc | 0j0
if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {
f=0;
for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
for (int l = 0; l < L ; l++) {
f = (Qy - y3)*arr5[j-1][0][ii][jj][kk][l] + (Wy - Qy)*arr5[j-1][0][ii][jj][kk][l+1] ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*(arr5[j-2][0][ii][jj][kk][l] - (om/gam2)*arr5[j-2][0][ii][jj][kk][l+1] ) ;
}

f = f + ( (jj)/( 2.0*(gam1 + gam2) ) )*arr5[j-1][0][ii][jj-1][kk][l+1] ;
arr5[j][0][ii][jj][kk][l] = f;
}
}
}
}
}

}

// abc | mn0
// arr5[j][i][n][m][k][l] 
// if C1 > 0:

if( L20 > 0 )  {
for(int j = 1; j < L21 + 1 ; j++){
for(int i = 1; i < L20 + 1 ; i++){

// inicial 000 | ij0:
    for (int l = 0; l < L ; l++) {
f = (Qy - y3)*arr5[j-1][i][0][0][0][l] + (Wy - Qy)*arr5[j-1][i][0][0][0][l+1] ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*(arr5[j-2][i][0][0][0][l] - (om/gam2)*arr5[j-2][i][0][0][0][l+1] ) ;
}

arr5[j][i][0][0][0][l] = f;
}


// a00 | ij0
if(  L10 > 0  ){  
f= 0;
for(int ii = 1; ii < L10 + 1 ; ii++){
  for (int l = 0; l < L  ; l++) {

f = (Qy - y3)*arr5[j-1][i][ii][0][0][l] + (Wy - Qy)*arr5[j-1][i][ii][0][0][l+1] ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*( arr5[j-2][i][ii][0][0][l] - (om/gam2)*arr5[j-2][i][ii][0][0][l+1] ) ;
}

arr5[j][i][ii][0][0][l] = f;
}
}
}

// 0b0 | ij0
if (  L11 > 0   ){
f = 0;
for(int jj = 1; jj < L11 + 1 ; jj++){
 for (int l = 0; l < L  ; l++) {

f = (Qy - y3)*arr5[j-1][i][0][jj][0][l] + (Wy - Qy)*arr5[j-1][i][0][jj][0][l+1] ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*( arr5[j-2][i][0][jj][0][l] - (om/gam2)*arr5[j-2][i][0][jj][0][l+1] ) ;
}
f = f + ( (jj)/( 2.0*(gam1 + gam2) ) )*arr5[j-1][i][0][jj-1][0][l+1] ;
arr5[j][i][0][jj][0][l] = f;

}
}
}

// 00c | ij0
if (  L12 > 0   ){
f = 0;
for(int kk = 1; kk < L12 + 1 ; kk++){
 for (int l = 0; l < L  ; l++) {
f = (Qy - y3)*arr5[j-1][i][0][0][kk][l] + (Wy - Qy)*arr5[j-1][i][0][0][kk][l+1] ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*( arr5[j-2][i][0][0][kk][l] - (om/gam2)*arr5[j-2][i][0][0][kk][l+1] ) ;
}
arr5[j][i][0][0][kk][l] = f;

}
}
}

// ab0 | ij0
if(   (L10 >= 1) && (  L11 >= 1 )  ) {
f=0;
for(int jj = 1; jj < L11 + 1 ; jj++){
  for(int ii = 1; ii < L10 + 1 ; ii++){
    for (int l = 0; l < L  ; l++) { 

f = (Qy - y3)*arr5[j-1][i][ii][jj][0][l] + (Wy - Qy)*arr5[j-1][i][ii][jj][0][l+1] ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*( arr5[j-2][i][ii][jj][0][l] - (om/gam2)*arr5[j-2][i][ii][jj][0][l+1] ) ;
}
f = f + ( (jj)/( 2.0*(gam1 + gam2) ) )*arr5[j-1][i][ii][jj-1][0][l+1] ;
arr5[j][i][ii][jj][0][l] = f;
}
}
}
}


if(   (L10 >= 1) && (  L12 >= 1 )  ) {
f=0;  
//  a0b | ij0
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
for (int l = 0; l < L  ; l++) {

f = (Qy - y3)*arr5[j-1][i][ii][0][kk][l] + (Wy - Qy)*arr5[j-1][i][ii][0][kk][l+1];
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*( arr5[j-2][i][ii][0][kk][l] - (om/gam2)*arr5[j-2][i][ii][0][kk][l+1] ) ;
}
arr5[j][i][ii][0][kk][l] = f;
}
}
}
}

if(   (L11 >= 1) && (  L12 >= 1 )  ) {
f=0;  //  0ab | ij0
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int jj = 1; jj < L11 +1  ; jj++){
for (int l = 0; l < L  ; l++) {
f = (Qy - y3)*arr5[j-1][i][0][jj][kk][l] + (Wy - Qy)*arr5[j-1][i][0][jj][kk][l+1]  ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*( arr5[j-2][i][0][jj][kk][l] - (om/gam2)*arr5[j-2][i][0][jj][kk][l+1] ) ;
}
f = f + ( (jj)/( 2.0*(gam1 + gam2) ) )*arr5[j-1][i][0][jj-1][kk][l+1] ;
arr5[j][i][0][jj][kk][l] = f;
}
}
}
}

// abc | ij0
if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {
f=0;
for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
for (int l = 0; l < L ; l++) {
f = (Qy - y3)*arr5[j-1][i][ii][jj][kk][l] + (Wy - Qy)*arr5[j-1][i][ii][jj][kk][l+1] ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*(arr5[j-2][i][ii][jj][kk][l] - (om/gam2)*arr5[j-2][i][ii][jj][kk][l+1] ) ;
}

f = f + ( (jj)/( 2.0*(gam1 + gam2) ) )*arr5[j-1][i][ii][jj-1][kk][l+1] ;
arr5[j][i][ii][jj][kk][l] = f;
}
}
}
}
}

}


}
}
//fin if  L20 > 0

Norminter = arr5[L21][L20][L10][L11][L12][0];  
// check de moment todo al valor l = 0;

for (int i = 0; i < L21 + 1; i++) {
    for (int j = 0; j < L20 + 1; j++) {
        for (int k = 0; k < L10 + 1; k++) {
            for (int l = 0; l < L11 + 1; l++) {
                for (int m = 0; m < L12 + 1; m++) {
                    delete[] arr5[i][j][k][l][m];
                }
                delete[] arr5[i][j][k][l];
            }
            delete[] arr5[i][j][k];
        }
        delete[] arr5[i][j];
    }
    delete[] arr5[i];
}
delete[] arr5;

return Norminter;

}else{

f = 0;
G4double *******arr6;
arr6 = new G4double******[L22 + 1]; 
for(int i = 0; i < L22 + 1 ; i++){
arr6[i] = new G4double*****[ L21 + 1];
for(int j = 0; j < L21 + 1 ; j++){
arr6[i][j] = new G4double****[ L20 + 1];
for(int k = 0; k <  L20 + 1 ; k++) {
arr6[i][j][k] =  new G4double***[ L10 + 1];
for(int l = 0; l <  L10 + 1 ; l++) {
arr6[i][j][k][l] = new G4double**[L11 + 1];
for(int m = 0; m <  L11 + 1 ; m++) {
arr6[i][j][k][l][m] =  new G4double*[L12 + 1];
for(int n = 0; n <  L12 + 1 ; n++) {
arr6[i][j][k][l][m][n] =  new G4double[L + 1];
}
}
}
}
}
}


if(T > 100){
arr6[0][0][0][0][0][0][0] = sqrt(pi/(4*T));

if (L > 0){
for (int i = 0; i < L ; i++) {
arr6[0][0][0][0][0][0][i+1]= (2.0*i+1)*arr6[0][0][0][0][0][0][i]/(2.0*T);
}

}

}else{

arr6[0][0][0][0][0][0][L] = fOpx->Gamma_aprox( om, Rx, Ry, Rz, 1, 1, 1,
1, 1, 1, L10, L11, L12, L20, L21, L22,
gam1, gam2, polyn, selector );

delete fOpx;

if (L > 0){
for(int i = L ; i > 0; i--) {
arr6[0][0][0][0][0][0][i-1] =  (2.0*T*arr6[0][0][0][0][0][0][i] + G4Exp(-T) )/(2.0*(i-1)+1);
}
}
} //else T > 100

for (int j = 0; j < L + 1 ; j++) {
arr6[0][0][0][0][0][0][j]  =  (arr6[0][0][0][0][0][0][j]*norms*pow(2,1)*pow(pi,5.0/2.0)*Kp*Kq*
pow(gam1,-1)*pow(gam2,-1) )/ ( sqrt(G) );
}



if (  L10 > 0   ){
for(int i = 1; i < L10 + 1 ; i++){
 for (int l = 0; l < L  ; l++) {

f = (Px - x1)*arr6[0][0][0][i-1][0][0][l] + (Wx - Px)*arr6[0][0][0][i-1][0][0][l+1] ; 
if( (i-1) > 0  ){
f = f + ( (i-1)/(2.0*gam1) )*( arr6[0][0][0][i-2][0][0][l] - (om/gam1)*arr6[0][0][0][i-2][0][0][l+1] ) ;
}
arr6[0][0][0][i][0][0][l] = f;
}
}
}

f =0;
if (  L11 > 0   ){
for(int j = 1; j < L11 + 1 ; j++){
 for (int l = 0; l < L  ; l++) {
    // f = (Py - y1)*arr6[0][0][0][0][j-1][0][l] - (Wy - Py)*arr6[0][0][0][0][j-1][0][l+1];
f = (Py - y1)*arr6[0][0][0][0][j-1][0][l] + (Wy - Py)*arr6[0][0][0][0][j-1][0][l+1];
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam1))*( arr6[0][0][0][0][j-2][0][l] - (om/gam1)*arr6[0][0][0][0][j-2][0][l+1] ) ;
}
arr6[0][0][0][0][j][0][l] = f;
}
}
}


if (  L12 > 0   ){
for(int k = 1; k < L12 + 1 ; k++){
 for (int l = 0; l < L  ; l++) {

f = (Pz - z1)*arr6[0][0][0][0][0][k-1][l] + (Wz - Pz)*arr6[0][0][0][0][0][k-1][l+1]  ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam1))*(arr6[0][0][0][0][0][k-2][l] - (om/gam1)*arr6[0][0][0][0][0][k-2][l+1] ) ;
}
arr6[0][0][0][0][0][k][l] = f;
}
}
}

if(   (L10 >= 1) && (  L11 >= 1 )  ) {
f=0;
//  ab0 | 000
for(int jj = 1; jj < L11 + 1 ; jj++){
  for(int ii = 1; ii < L10 + 1 ; ii++){
for (int l = 0; l < L  ; l++) {
//  f = (Py - y1)*arr6[0][0][0][ii][jj-1][0][l] - (Wy - Py)*arr6[0][0][0][ii][jj-1][0][l+1]  ;
f = (Py - y1)*arr6[0][0][0][ii][jj-1][0][l] + (Wy - Py)*arr6[0][0][0][ii][jj-1][0][l+1]  ;
if( (jj-1) > 0  ){
f = f + ( (jj-1)/(2.0*gam1))*( arr6[0][0][0][ii][jj-2][0][l] - (om/gam1)*arr6[0][0][0][ii][jj-2][0][l+1] ) ;
}
arr6[0][0][0][ii][jj][0][l] = f;
}
}
}
}


if(   (L10 >= 1) && (  L12 >= 1 )  ) {
f=0;
// a0b | 000
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
for (int l = 0; l < L  ; l++) {
 //   f = (Pz - z1)*arr6[0][0][0][ii][0][kk-1][l] - (Wz - Pz)*arr6[0][0][0][ii][0][kk-1][l+1]  ;
f = (Pz - z1)*arr6[0][0][0][ii][0][kk-1][l] + (Wz - Pz)*arr6[0][0][0][ii][0][kk-1][l+1]  ;
if( (kk-1) > 0  ){
f = f + ( (kk-1)/(2.0*gam1))*( 
arr6[0][0][0][ii][0][kk-2][l] - (om/gam1)*arr6[0][0][0][ii][0][kk-2][l+1] ) ;
}
arr6[0][0][0][ii][0][kk][l] = f;
}
}
}
}


if(   (L11 >= 1) && (  L12 >= 1 )  ) {
f=0;
// 0ab | 000
for(int kk = 1; kk < L12 + 1  ; kk++){
  for(int jj = 1; jj < L11 + 1; jj++){
for (int l = 0; l < L  ; l++) {

f = (Pz - z1)*arr6[0][0][0][0][jj][kk-1][l] + (Wz - Pz)*arr6[0][0][0][0][jj][kk-1][l+1]  ;
if( (kk-1) > 0  ){
f = f + ( (kk-1)/(2.0*gam1))*( arr6[0][0][0][0][jj][kk-2][l] - (om/gam1)*arr6[0][0][0][0][jj][kk-2][l+1] ) ;
}

arr6[0][0][0][0][jj][kk][l] = f;
}
}
}
}

if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {

f=0;
for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
for (int l = 0; l < L ; l++) {
f = (Pz - z1)*arr6[0][0][0][ii][jj][kk-1][l] + (Wz - Pz)*arr6[0][0][0][ii][jj][kk-1][l+1]  ;
if( (kk-1) > 0  ){
f = f + ( (kk-1)/(2.0*gam1))*( arr6[0][0][0][ii][jj][kk-2][l] - (om/gam1)*arr6[0][0][0][ii][jj][kk-2][l+1] ) ;
}
arr6[0][0][0][ii][jj][kk][l] = f;

}
}
}
}
}


if( L20 >  0   ){
// here C1[0] > 0
// [C3[C2][C1][A1][A2[A3]
f = 0;
//for(int i = 1; i < L20 + 1 ; i++){
//}
//  abc | i00
for(int i = 1; i < L20 + 1 ; i++){
f = 0;
 for (int l = 0; l < (L -L20)  ; l++) {
f = (Qx - x3)*arr6[0][0][i-1][0][0][0][l] + (Wx - Qx)*arr6[0][0][i-1][0][0][0][l+1] ;
if( (i-1) > 0  ){
f = f + ( (i-1)/(2.0*gam2))*( arr6[0][0][i-2][0][0][0][l] - (om/gam2)*arr6[0][0][i-2][0][0][0][l+1] ) ;
}

arr6[0][0][i][0][0][0][l] = f;
}


if (  L10 > 0   ){
for(int ii = 1; ii < L10 + 1 ; ii++){
 for (int l = 0; l < L  ; l++) {

f = (Qx - x3)*arr6[0][0][i-1][ii][0][0][l] + (Wx - Qx)*arr6[0][0][i-1][ii][0][0][l+1] ;
if( (i-1) > 0  ){
f = f + ( (i-1)/(2.0*gam2))*( arr6[0][0][i-2][ii][0][0][l] - (om/gam2)*arr6[0][0][i-2][ii][0][0][l+1] ) ;
}

f = f + ( (ii)/( 2.0*(gam1 + gam2) ) )*arr6[0][0][i-1][ii-1][0][0][l+1] ;
arr6[0][0][i][ii][0][0][l] = f;
}
}
}


if (  L11 > 0   ){
f = 0;
// 0b0 | i00
for(int jj = 1; jj < L11 + 1 ; jj++){
 for (int l = 0; l < L  ; l++) {
// f = (Qx - x3)*arr6[0][0][i-1][0][jj][0][l] - (Wx - Qx)*arr6[0][0][i-1][0][jj][0][l+1] ;
f = (Qx - x3)*arr6[0][0][i-1][0][jj][0][l] + (Wx - Qx)*arr6[0][0][i-1][0][jj][0][l+1] ;
if( (i-1) > 0  ){
f = f + ( (i-1)/(2.0*gam2))*( arr6[0][0][i-2][0][jj][0][l] - (om/gam2)*arr6[0][0][i-2][0][jj][0][l+1] ) ;
}
arr6[0][0][i][0][jj][0][l] = f;

}
}
}

// 00c | i00
if (  L12 > 0   ){
for(int kk = 1; kk < L12 + 1 ; kk++){
f = 0;
 for (int l = 0; l < L  ; l++) {

f = (Qx - x3)*arr6[0][0][i-1][0][0][kk][l] + (Wx - Qx)*arr6[0][0][i-1][0][0][kk][l+1] ;
if( (i-1) > 0  ){
f = f + ( (i-1)/(2.0*gam2))*( arr6[0][0][i-2][0][0][kk][l] - (om/gam2)*arr6[0][0][i-2][0][0][kk][l+1] ) ;
}
arr6[0][0][i][0][0][kk][l] = f;

}
}
}

if(   (L10 >= 1) && (  L11 >= 1 )  ) {
f=0;
// ab0 | i00
for(int jj = 1; jj < L11 + 1 ; jj++){
  for(int ii = 1; ii < L10 + 1 ; ii++){
    for (int l = 0; l < L  ; l++) { 
f = (Qx - x3)*arr6[0][0][i-1][ii][jj][0][l] + (Wx - Qx)*arr6[0][0][i-1][ii][jj][0][l+1] ;
if( (i-1) > 0  ){
// at least  i = 2:
f = f + ( (i-1)/(2.0*gam2))*( arr6[0][0][i-2][ii][jj][0][l] - (om/gam2)*arr6[0][0][i-2][ii][jj][0][l+1] ) ;
}
f = f + ( (ii)/( 2.0*(gam1 + gam2) ) )*arr6[0][0][i-1][ii-1][jj][0][l+1] ;
arr6[0][0][i][ii][jj][0][l] = f;
}
}
}
}

if(   (L10 >= 1) && (  L12 >= 1 )  ) {
f=0;  //  a0b | i00
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
for (int l = 0; l < L  ; l++) {

f = (Qx - x3)*arr6[0][0][i-1][ii][0][kk][l] + (Wx - Qx)*arr6[0][0][i-1][ii][0][kk][l+1]  ;
if( (i-1) > 0  ){
f = f + ( (i-1)/(2.0*gam2))*( arr6[0][0][i-2][ii][0][kk][l] - (om/gam2)*arr6[0][0][i-2][ii][0][kk][l+1] ) ;
}
f = f + ( (ii)/( 2.0*(gam1 + gam2) ) )*arr6[0][0][i-1][ii-1][0][kk][l+1] ;
arr6[0][0][i][ii][0][kk][l] = f;
}
}
}
}


if(   (L11 >= 1) && (  L12 >= 1 )  ) {
f=0;  //  0ab | i00 0 nab | i00
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int jj = 1; jj < L11 +1  ; jj++){
for (int l = 0; l < L  ; l++) {

f = (Qx - x3)*arr6[0][0][i-1][0][jj][kk][l] + (Wx - Qx)*arr6[0][0][i-1][0][jj][kk][l+1]  ;
if( (i-1) > 0  ){
f = f + ( (i-1)/(2.0*gam2))*( arr6[0][0][i-2][0][jj][kk][l] - (om/gam2)*arr6[0][0][i-2][0][jj][kk][l+1] ) ;
}

arr6[0][0][i][0][jj][kk][l] = f;
}
}
}
}

// abc | i00
if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {
f=0;
for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
for (int l = 0; l < L ; l++) {
f = (Qx - x3)*arr6[0][0][i-1][ii][jj][kk][l] + (Wx - Qx)*arr6[0][0][i-1][ii][jj][kk][l+1] ;
if( (i-1) > 0  ){
f = f + ( (i-1)/(2.0*gam2))*(arr6[0][0][i-2][ii][jj][kk][l] - (om/gam2)*arr6[0][0][i-2][ii][jj][kk][l+1] ) ;
}
f = f + ( (ii)/( 2.0*(gam1 + gam2) ) )*arr6[0][0][i-1][ii-1][jj][kk][l+1] ;
arr6[0][0][i][ii][jj][kk][l] = f;

}
}
}
}
}


}

}
//closing if L20 > 0

if( L21 >  0   ){
f = 0;
// arr6[0][j][0][0][0][0][l] 
for(int j = 1; j < L21 + 1 ; j++){
 
// inicial 000 | 0j0
 for (int l = 0; l < (L - L21)  ; l++) {
f = (Qy - y3)*arr6[0][j-1][0][0][0][0][l] + (Wy - Qy)*arr6[0][j-1][0][0][0][0][l+1] ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*( arr6[0][j-2][0][0][0][0][l] - (om/gam2)*arr6[0][j-2][0][0][0][0][l+1] ) ;
}
arr6[0][j][0][0][0][0][l] = f;
}

// inicial 000 | ij0:
//for(int j = 1; j < L21 + 1 ; j++){
if (  L20 > 0   ){
  for(int i = 1; i < L20 + 1 ; i++){
   for (int l = 0; l < L  ; l++) {
f = (Qy - y3)*arr6[0][j-1][i][0][0][0][l] + (Wy - Qy)*arr6[0][j-1][i][0][0][0][l+1] ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*(arr6[0][j-2][i][0][0][0][l] - (om/gam2)*arr6[0][j-2][i][0][0][0][l+1] ) ;
}

arr6[0][j][i][0][0][0][l] = f;
}
}
}

// abc | 0j0
// [C3][C2][C1][A1][A2][A3]
// a00 | 0j0
if(  L10 > 0  ){  
f= 0;
for(int ii = 1; ii < L10 + 1 ; ii++){
  for (int l = 0; l < L  ; l++) {
  //  f = (Qy - y3)*arr6[0][j-1][0][ii][0][0][l] - (Wy - Qy)*arr6[0][j-1][0][ii][0][0][l+1] ;
f = (Qy - y3)*arr6[0][j-1][0][ii][0][0][l] + (Wy - Qy)*arr6[0][j-1][0][ii][0][0][l+1] ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*( arr6[0][j-2][0][ii][0][0][l] - (om/gam2)*arr6[0][j-2][0][ii][0][0][l+1] ) ;
}

arr6[0][j][0][ii][0][0][l] = f;
}
}
}

f = 0;
if (  L11 > 0   ){
for(int jj = 1; jj < L11 + 1 ; jj++){
 for (int l = 0; l < L  ; l++) {
  //  f = (Qy - y3)*arr6[0][j-1][0][0][jj][0][l] - (Wy - Qy)*arr6[0][j-1][0][0][jj][0][l+1] ;
f = (Qy - y3)*arr6[0][j-1][0][0][jj][0][l] + (Wy - Qy)*arr6[0][j-1][0][0][jj][0][l+1] ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*( arr6[0][j-2][0][0][jj][0][l] - (om/gam2)*arr6[0][j-2][0][0][jj][0][l+1] ) ;
}
f = f + ( (jj)/( 2.0*(gam1 + gam2) ) )*arr6[0][j-1][0][0][jj-1][0][l+1] ;
arr6[0][j][0][0][jj][0][l] = f;

}
}
}

f = 0;
// 00c | 0j0
if (  L12 > 0   ){
for(int kk = 1; kk < L12 + 1 ; kk++){
 for (int l = 0; l < L  ; l++) {

f = (Qy - y3)*arr6[0][j-1][0][0][0][kk][l] + (Wy - Qy)*arr6[0][j-1][0][0][0][kk][l+1] ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*( arr6[0][j-2][0][0][0][kk][l] - (om/gam2)*arr6[0][j-2][0][0][0][kk][l+1] ) ;
}
arr6[0][j][0][0][0][kk][l] = f;

}
}
}


if(   (L10 >= 1) && (  L11 >= 1 )  ) {
f=0;
// ab0 | 0j0
for(int jj = 1; jj < L11 + 1 ; jj++){
  for(int ii = 1; ii < L10 + 1 ; ii++){
    for (int l = 0; l < L  ; l++) { 
f = (Qy - y3)*arr6[0][j-1][0][ii][jj][0][l] + (Wy - Qy)*arr6[0][j-1][0][ii][jj][0][l+1] ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*( arr6[0][j-2][0][ii][jj][0][l] - (om/gam2)*arr6[0][j-2][0][ii][jj][0][l+1] ) ;
}
f = f + ( (jj)/( 2.0*(gam1 + gam2) ) )*arr6[0][j-1][0][ii][jj-1][0][l+1] ;
arr6[0][j][0][ii][jj][0][l] = f;
}
}
}
}


if(   (L10 >= 1) && (  L12 >= 1 )  ) {
f=0;  
//  a0b | 0j0
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
for (int l = 0; l < L  ; l++) {

f = (Qy - y3)*arr6[0][j-1][0][ii][0][kk][l] + (Wy - Qy)*arr6[0][j-1][0][ii][0][kk][l+1]  ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*( arr6[0][j-2][0][ii][0][kk][l] - (om/gam2)*arr6[0][j-2][0][ii][0][kk][l+1] ) ;
}
arr6[0][j][0][ii][0][kk][l] = f;
}
}
}
}


if(   (L11 >= 1) && (  L12 >= 1 )  ) {
f=0;  //  0ab | 0j0
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int jj = 1; jj < L11 +1  ; jj++){
for (int l = 0; l < (L - L21); l++) {
f = (Qy - y3)*arr6[0][j-1][0][0][jj][kk][l] + (Wy - Qy)*arr6[0][j-1][0][0][jj][kk][l+1]  ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*( arr6[0][j-2][0][0][jj][kk][l] - (om/gam2)*arr6[0][j-2][0][0][jj][kk][l+1] ) ;
}
f = f + ( (jj)/( 2.0*(gam1 + gam2) ) )*arr6[0][j-1][0][0][jj-1][kk][l+1] ;
arr6[0][j][0][0][jj][kk][l] = f;
}
}
}
}

// abc | 0j0
if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {
f=0;
for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
for (int l = 0; l < L ; l++) {
    // f = (Qy - y3)*arr6[0][j-1][0][ii][jj][kk][l] - (Wy - Qy)*arr6[0][j-1][0][ii][jj][kk][l+1] ;
f = (Qy - y3)*arr6[0][j-1][0][ii][jj][kk][l] + (Wy - Qy)*arr6[0][j-1][0][ii][jj][kk][l+1] ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*(arr6[0][j-2][0][ii][jj][kk][l] - (om/gam2)*arr6[0][j-2][0][ii][jj][kk][l+1] ) ;
}

f = f + ( (jj)/( 2.0*(gam1 + gam2) ) )*arr6[0][j-1][0][ii][jj-1][kk][l+1] ;
arr6[0][j][0][ii][jj][kk][l] = f;
}
}
}
}
}

}
//}
}


// abc | mn0
// arr5[j][i][n][m][k][l] 
// se ejecuta si C1 > 0:
if( (L20 > 0) && (  L21 > 0  )   )  {

for(int j = 1; j < L21 + 1 ; j++){
  for(int i = 1; i < L20 + 1 ; i++){

// a00 | ij0
if(  L10 > 0  ){  
f= 0;
for(int ii = 1; ii < L10 + 1 ; ii++){
  for (int l = 0; l < L  ; l++) {
f = (Qy - y3)*arr6[0][j-1][i][ii][0][0][l] + (Wy - Qy)*arr6[0][j-1][i][ii][0][0][l+1] ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*( arr6[0][j-2][i][ii][0][0][l] - (om/gam2)*arr6[0][j-2][i][ii][0][0][l+1] ) ;
}

arr6[0][j][i][ii][0][0][l] = f;
}
}
}


// 0b0 | ij0
if (  L11 > 0  ){
f = 0;
for(int jj = 1; jj < L11 + 1 ; jj++){
 for (int l = 0; l < L  ; l++) {

f = (Qy - y3)*arr6[0][j-1][i][0][jj][0][l] + (Wy - Qy)*arr6[0][j-1][i][0][jj][0][l+1] ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*( arr6[0][j-2][i][0][jj][0][l] - (om/gam2)*arr6[0][j-2][i][0][jj][0][l+1] ) ;
}
f = f + ( (jj)/( 2.0*(gam1 + gam2) ) )*arr6[0][j-1][i][0][jj-1][0][l+1] ;
arr6[0][j][i][0][jj][0][l] = f;

}
}
}

// 00c | ij0
if (  L12 > 0   ){
f = 0;
for(int kk = 1; kk < L12 + 1 ; kk++){
 for (int l = 0; l < L  ; l++) {
f = (Qy - y3)*arr6[0][j-1][i][0][0][kk][l] + (Wy - Qy)*arr6[0][j-1][i][0][0][kk][l+1] ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*( arr6[0][j-2][i][0][0][kk][l] - (om/gam2)*arr6[0][j-2][i][0][0][kk][l+1] ) ;
}
arr6[0][j][i][0][0][kk][l] = f;

}
}
}

// ab0 | ij0
if(   (L10 >= 1) && (  L11 >= 1 )  ) {
f=0;
for(int jj = 1; jj < L11 + 1 ; jj++){
  for(int ii = 1; ii < L10 + 1 ; ii++){
    for (int l = 0; l < (L - L10); l++) { 
f = (Qy - y3)*arr6[0][j-1][i][ii][jj][0][l] + (Wy - Qy)*arr6[0][j-1][i][ii][jj][0][l+1] ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*( arr6[0][j-2][i][ii][jj][0][l] - (om/gam2)*arr6[0][j-2][i][ii][jj][0][l+1] ) ;
}
f = f + ( (jj)/( 2.0*(gam1 + gam2) ) )*arr6[0][j-1][i][ii][jj-1][0][l+1] ;
arr6[0][j][i][ii][jj][0][l] = f;
}
}
}
}


if(   (L10 >= 1) && (  L12 >= 1 )  ) {
f=0;  
//  a0b | ij0
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
for (int l = 0; l < L  ; l++) {

f = (Qy - y3)*arr6[0][j-1][i][ii][0][kk][l] + (Wy - Qy)*arr6[0][j-1][i][ii][0][kk][l+1];
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*( arr6[0][j-2][i][ii][0][kk][l] - (om/gam2)*arr6[0][j-2][i][ii][0][kk][l+1] ) ;
}
arr6[0][j][i][ii][0][kk][l] = f;
}
}
}
}

if(   (L11 >= 1) && (  L12 >= 1 )  ) {
f=0;  //  0ab | ij0
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int jj = 1; jj < L11 +1  ; jj++){
for (int l = 0; l < L  ; l++) {
f = (Qy - y3)*arr6[0][j-1][i][0][jj][kk][l] + (Wy - Qy)*arr6[0][j-1][i][0][jj][kk][l+1]  ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*( arr6[0][j-2][i][0][jj][kk][l] - (om/gam2)*arr6[0][j-2][i][0][jj][kk][l+1] ) ;
}
f = f + ( (jj)/( 2.0*(gam1 + gam2) ) )*arr6[0][j-1][i][0][jj-1][kk][l+1] ;
arr6[0][j][i][0][jj][kk][l] = f;
}
}
}
}

// abc | ij0
if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {
f=0;
for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
for (int l = 0; l < L ; l++) {
f = (Qy - y3)*arr6[0][j-1][i][ii][jj][kk][l] + (Wy - Qy)*arr6[0][j-1][i][ii][jj][kk][l+1] ;
if( (j-1) > 0  ){
f = f + ( (j-1)/(2.0*gam2))*(arr6[0][j-2][i][ii][jj][kk][l] - (om/gam2)*arr6[0][j-2][i][ii][jj][kk][l+1] ) ;
}

f = f + ( (jj)/( 2.0*(gam1 + gam2) ) )*arr6[0][j-1][i][ii][jj-1][kk][l+1] ;
arr6[0][j][i][ii][jj][kk][l] = f;
}
}
}
}
}



}
}
}


// From   [C3[C2][C1][A1][A2[A3]
//  abc | 00k
for(int k = 1; k < L22 + 1; k++){
// caso base: 000 | 00k
for (int l = 0; l < L  ; l++) {
  //  f = (Qz - z3)*arr6[k-1][0][0][0][0][0][l] - (Wz - Qz)*arr6[k-1][0][0][0][0][0][l+1] ;
f = (Qz - z3)*arr6[k-1][0][0][0][0][0][l] + (Wz - Qz)*arr6[k-1][0][0][0][0][0][l+1] ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*( arr6[k-2][0][0][0][0][0][l] - (om/gam2)*arr6[k-2][0][0][0][0][0][l+1] ) ;
}
arr6[k][0][0][0][0][0][l] = f;
}

// a00 | 00k 
if(  L10 > 0  ){  
f= 0;
for(int ii = 1; ii < L10 + 1 ; ii++){
  for (int l = 0; l < L  ; l++) {
f = (Qz - z3)*arr6[k-1][0][0][ii][0][0][l] + (Wz - Qz)*arr6[k-1][0][0][ii][0][0][l+1] ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*( arr6[k-2][0][0][ii][0][0][l] - (om/gam2)*arr6[k-2][0][0][ii][0][0][l+1] ) ;
}
arr6[k][0][0][ii][0][0][l] = f;
}
}
}

// 0b0 | 00k
if (  L11 > 0   ){
f = 0;
for(int jj = 1; jj < L11 + 1 ; jj++){
 for (int l = 0; l < L  ; l++) {
// f = (Qz - z3)*arr6[k-1][0][0][0][jj][0][l] - (Wz - Qz)*arr6[k-1][0][0][0][jj][0][l+1] ;
f = (Qz - z3)*arr6[k-1][0][0][0][jj][0][l] + (Wz - Qz)*arr6[k-1][0][0][0][jj][0][l+1] ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*( arr6[k-2][0][0][0][jj][0][l] - (om/gam2)*arr6[k-2][0][0][0][jj][0][l+1] ) ;
}
arr6[k][0][0][0][jj][0][l] = f;
}
}
}


// 00c | 00k
if (  L12 > 0   ){
f = 0;
for(int kk = 1; kk < L12 + 1 ; kk++){
 for (int l = 0; l < (L - L12)  ; l++) {

f = (Qz - z3)*arr6[k-1][0][0][0][0][kk][l] + (Wz - Qz)*arr6[k-1][0][0][0][0][kk][l+1] ;
if( ( k - 1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*( arr6[k-2][0][0][0][0][kk][l] - (om/gam2)*arr6[k-2][0][0][0][0][kk][l+1] ) ;
}
f = f + ( (kk)/( 2.0*(gam1 + gam2) ) )*arr6[k-1][0][0][0][0][kk-1][l+1] ;
arr6[k][0][0][0][0][kk][l] = f;
}
}
}

// ab0 | 00k
if(   (L10 >= 1) && (  L11 >= 1 )  ) {
f=0;
for(int jj = 1; jj < L11 + 1 ; jj++){
  for(int ii = 1; ii < L10 + 1 ; ii++){
    for (int l = 0; l < L  ; l++) { 

f = (Qz - z3)*arr6[k-1][0][0][ii][jj][0][l] + (Wz - Qz)*arr6[k-1][0][0][ii][jj][0][l+1] ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*( arr6[k-2][0][0][ii][jj][0][l] - (om/gam2)*arr6[k-2][0][0][ii][jj][0][l+1] ) ;
}
arr6[k][0][0][ii][jj][0][l] = f;
}
}
}
}


if(   (L10 >= 1) && (  L12 >= 1 )  ) {
f=0;  
//  a0c | 00k
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
for (int l = 0; l < L  ; l++) {

f = (Qz - z3)*arr6[k-1][0][0][ii][0][kk][l] + (Wz - Qz)*arr6[k-1][0][0][ii][0][kk][l+1];
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*( arr6[k-2][0][0][ii][0][kk][l] - (om/gam2)*arr6[k-2][0][0][ii][0][kk][l+1] ) ;
}
f = f + ( (kk)/( 2.0*(gam1 + gam2) ) )*arr6[k-1][0][0][ii][0][kk-1][l+1] ;
arr6[k][0][0][ii][0][kk][l] = f;
}
}
}
}

if(   (L11 >= 1) && (  L12 >= 1 )  ) {
f=0;  //  0bc | 00k
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int jj = 1; jj < L11 +1  ; jj++){
for (int l = 0; l < L  ; l++) {

f = (Qz - z3)*arr6[k-1][0][0][0][jj][kk][l] + (Wz - Qz)*arr6[k-1][0][0][0][jj][kk][l+1]  ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*( arr6[k-2][0][0][0][jj][kk][l] - (om/gam2)*arr6[k-2][0][0][0][jj][kk][l+1] ) ;
}
f = f + ( (kk)/( 2.0*(gam1 + gam2) ) )*arr6[k-1][0][0][0][jj][kk-1][l+1] ;
arr6[k][0][0][0][jj][kk][l] = f;
}
}
}
}

// abc | 00k
if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {
f=0;
for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
for (int l = 0; l < L ; l++) {

f = (Qz - z3)*arr6[k-1][0][0][ii][jj][kk][l] + (Wz - Qz)*arr6[k-1][0][0][ii][jj][kk][l+1] ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*(arr6[k-2][0][0][ii][jj][kk][l] - (om/gam2)*arr6[k-2][0][0][ii][jj][kk][l+1] ) ;
}
f = f + ( (kk)/( 2.0*(gam1 + gam2) ) )*arr6[k-1][0][0][ii][jj][kk-1][l+1] ;
arr6[k][0][0][ii][jj][kk][l] = f;
}
}
}
}
}

// cierre de recorrido k
}

//  abc | i0k si i>0
if (L20 > 0 ){  
//abc | i0k
for(int k = 1; k < L22 + 1; k++){
    for(int i = 1; i < L20 + 1 ; i++){

for (int l = 0; l < L ; l++) {
f=0;
// inicial 000 | i0k
f = (Qz - z3)*arr6[k-1][0][i][0][0][0][l] + (Wz - Qz)*arr6[k-1][0][i][0][0][0][l+1] ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*(arr6[k-2][0][i][0][0][0][l] - (om/gam2)*arr6[k-2][0][i][0][0][0][l+1] ) ;
}
arr6[k][0][i][0][0][0][l] = f;
}


if(  L10 > 0  ){  
f= 0;
// a00 | i0k
for(int ii = 1; ii < L10 + 1 ; ii++){
  for (int l = 0; l < L  ; l++) {
  // f = (Qz - z3)*arr6[k-1][0][i][ii][0][0][l] - (Wz - Qz)*arr6[k-1][0][i][ii][0][0][l+1] ;
f = (Qz - z3)*arr6[k-1][0][i][ii][0][0][l] + (Wz - Qz)*arr6[k-1][0][i][ii][0][0][l+1] ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*( arr6[k-2][0][i][ii][0][0][l] - (om/gam2)*arr6[k-2][0][i][ii][0][0][l+1] ) ;
}
arr6[k][0][i][ii][0][0][l] = f;
}
}
}

// 0b0 | i0k
if (  L11 > 0   ){
f = 0;
for(int jj = 1; jj < L11 + 1 ; jj++){
 for (int l = 0; l < L  ; l++) {

f = (Qz - z3)*arr6[k-1][0][i][0][jj][0][l] + (Wz - Qz)*arr6[k-1][0][i][0][jj][0][l+1] ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*( arr6[k-2][0][i][0][jj][0][l] - (om/gam2)*arr6[k-2][0][i][0][jj][0][l+1] ) ;
}
arr6[k][0][i][0][jj][0][l] = f;
}
}
}

// 00c | i0k
if (  L12 > 0   ){
f = 0;
for(int kk = 1; kk < L12 + 1 ; kk++){
 for (int l = 0; l < L  ; l++) {
// f = (Qz- z3)*arr6[k-1][0][i][0][0][kk][l] - (Wz - Qz)*arr6[k-1][0][i][0][0][kk][l+1] ;
f = (Qz- z3)*arr6[k-1][0][i][0][0][kk][l] + (Wz - Qz)*arr6[k-1][0][i][0][0][kk][l+1] ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*( arr6[k-2][0][i][0][0][kk][l] - (om/gam2)*arr6[k-2][0][i][0][0][kk][l+1] ) ;
}
f = f + ( (kk)/( 2.0*(gam1 + gam2) ) )*arr6[k-1][0][i][0][0][kk-1][l+1] ;
arr6[k][0][i][0][0][kk][l] = f;

}
}
}

// ab0 | i0k
if(   (L10 >= 1) && (  L11 >= 1 )  ) {
f=0;
for(int jj = 1; jj < L11 + 1 ; jj++){
  for(int ii = 1; ii < L10 + 1 ; ii++){
    for (int l = 0; l < L  ; l++) { 

f = (Qz - z3)*arr6[k-1][0][i][ii][jj][0][l] + (Wz - Qz)*arr6[k-1][0][i][ii][jj][0][l+1] ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*( arr6[k-2][0][i][ii][jj][0][l] - (om/gam2)*arr6[k-2][0][i][ii][jj][0][l+1] ) ;
}
arr6[k][0][i][ii][jj][0][l] = f;
}
}
}
}


if(   (L10 >= 1) && (  L12 >= 1 )  ) {
f=0;  
//  a0c | i0k
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
for (int l = 0; l < L  ; l++) {
// f = (Qz - z3)*arr6[k-1][0][i][ii][0][kk][l] - (Wz - Qz)*arr6[k-1][0][i][ii][0][kk][l+1];
f = (Qz - z3)*arr6[k-1][0][i][ii][0][kk][l] + (Wz - Qz)*arr6[k-1][0][i][ii][0][kk][l+1];
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*( arr6[k-2][0][i][ii][0][kk][l] - (om/gam2)*arr6[k-2][0][i][ii][0][kk][l+1] ) ;
}
f = f + ( (kk)/( 2.0*(gam1 + gam2) ) )*arr6[k-1][0][i][ii][0][kk-1][l+1] ;
arr6[k][0][i][ii][0][kk][l] = f;
}
}
}
}

if(   (L11 >= 1) && (  L12 >= 1 )  ) {
f=0;  //  0bc | i0k
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int jj = 1; jj < L11 +1  ; jj++){
for (int l = 0; l < L  ; l++) {

f = (Qz - z3)*arr6[k-1][0][i][0][jj][kk][l] + (Wz - Qz)*arr6[k-1][0][i][0][jj][kk][l+1]  ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*( arr6[k-2][0][i][0][jj][kk][l] - (om/gam2)*arr6[k-2][0][i][0][jj][kk][l+1] ) ;
}
f = f + ( (kk)/( 2.0*(gam1 + gam2) ) )*arr6[k-1][0][i][0][jj][kk-1][l+1] ;
arr6[k][0][i][0][jj][kk][l] = f;
}
}
}
}

//  abc | i0k
if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {
f=0;
for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
for (int l = 0; l < L ; l++) {
f = (Qz - z3)*arr6[k-1][0][i][ii][jj][kk][l] + (Wz - Qz)*arr6[k-1][0][i][ii][jj][kk][l+1] ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*(arr6[k-2][0][i][ii][jj][kk][l] - (om/gam2)*arr6[k-2][0][i][ii][jj][kk][l+1] ) ;
}
f = f + ( (kk)/( 2.0*(gam1 + gam2) ) )*arr6[k-1][0][i][ii][jj][kk-1][l+1] ;
arr6[k][0][i][ii][jj][kk][l] = f;
}
}
}
}
}

// closing k , i
}
}

}

//  abc | 0jk
if (L21 > 0 ){  
for(int k = 1; k < L22 + 1; k++){
    for(int j = 1; j < L21 + 1 ; j++){
// caso inicial: 000 | 0jk
for (int l = 0; l < L ; l++) {
f=0;
// f = (Qz - z3)*arr6[k-1][j][0][0][0][0][l] - (Wz - Qz)*arr6[k-1][j][0][0][0][0][l+1] ;
f = (Qz - z3)*arr6[k-1][j][0][0][0][0][l] + (Wz - Qz)*arr6[k-1][j][0][0][0][0][l+1] ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*(arr6[k-2][j][0][0][0][0][l] - (om/gam2)*arr6[k-2][j][0][0][0][0][l+1] ) ;
}
arr6[k][j][0][0][0][0][l] = f;
}

// a00 | 0jk
if(  L10 > 0  ){  
f= 0;
for(int ii = 1; ii < L10 + 1 ; ii++){
  for (int l = 0; l < L  ; l++) {
f = (Qz - z3)*arr6[k-1][j][0][ii][0][0][l] + (Wz - Qz)*arr6[k-1][j][0][ii][0][0][l+1] ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*( arr6[k-2][j][0][ii][0][0][l] - (om/gam2)*arr6[k-2][j][0][ii][0][0][l+1] ) ;
}
arr6[k][j][0][ii][0][0][l] = f;
}
}
}

// 0b0 | 0jk
if (  L11 > 0   ){
f = 0;
for(int jj = 1; jj < L11 + 1 ; jj++){
 for (int l = 0; l < L  ; l++) {
//  f = (Qz - z3)*arr6[k-1][j][0][0][jj][0][l] - (Wz - Qz)*arr6[k-1][j][0][0][jj][0][l+1] ;
f = (Qz - z3)*arr6[k-1][j][0][0][jj][0][l] + (Wz - Qz)*arr6[k-1][j][0][0][jj][0][l+1] ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*( arr6[k-2][j][0][0][jj][0][l] - (om/gam2)*arr6[k-2][j][0][0][jj][0][l+1] ) ;
}
arr6[k][j][0][0][jj][0][l] = f;
}
}
}

// 00c | 0jk
if (  L12 > 0   ){
f = 0;
for(int kk = 1; kk < L12 + 1 ; kk++){
 for (int l = 0; l < L  ; l++) {
// f = (Qz- z3)*arr6[k-1][j][0][0][0][kk][l] - (Wz - Qz)*arr6[k-1][j][0][0][0][kk][l+1] ;
f = (Qz- z3)*arr6[k-1][j][0][0][0][kk][l] + (Wz - Qz)*arr6[k-1][j][0][0][0][kk][l+1] ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*( arr6[k-2][j][0][0][0][kk][l] - (om/gam2)*arr6[k-2][j][0][0][0][kk][l+1] ) ;
}
f = f + ( (kk)/( 2.0*(gam1 + gam2) ) )*arr6[k-1][j][0][0][0][kk-1][l+1] ;
arr6[k][j][0][0][0][kk][l] = f;

}
}
}

// ab0 | 0jk
if(   (L10 >= 1) && (  L11 >= 1 )  ) {
f=0;
for(int jj = 1; jj < L11 + 1 ; jj++){
  for(int ii = 1; ii < L10 + 1 ; ii++){
    for (int l = 0; l < L  ; l++) { 
//  f = (Qz - z3)*arr6[k-1][j][0][ii][jj][0][l] - (Wz - Qz)*arr6[k-1][j][0][ii][jj][0][l+1] ;
f = (Qz - z3)*arr6[k-1][j][0][ii][jj][0][l] + (Wz - Qz)*arr6[k-1][j][0][ii][jj][0][l+1] ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*( arr6[k-2][j][0][ii][jj][0][l] - (om/gam2)*arr6[k-2][j][0][ii][jj][0][l+1] ) ;
}
arr6[k][j][0][ii][jj][0][l] = f;
}
}
}
}

//  a0c | 0jk
if(   (L10 >= 1) && (  L12 >= 1 )  ) {
f=0;  
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
for (int l = 0; l < L  ; l++) {
// f = (Qz - z3)*arr6[k-1][j][0][ii][0][kk][l] - (Wz - Qz)*arr6[k-1][j][0][ii][0][kk][l+1];
f = (Qz - z3)*arr6[k-1][j][0][ii][0][kk][l] + (Wz - Qz)*arr6[k-1][j][0][ii][0][kk][l+1];
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*( arr6[k-2][j][0][ii][0][kk][l] - (om/gam2)*arr6[k-2][j][0][ii][0][kk][l+1] ) ;
}
f = f + ( (kk)/( 2.0*(gam1 + gam2) ) )*arr6[k-1][j][0][ii][0][kk-1][l+1] ;
arr6[k][j][0][ii][0][kk][l] = f;
}
}
}
}

 //  0bc | 0jk
if(   (L11 >= 1) && (  L12 >= 1 )  ) {
f=0; 
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int jj = 1; jj < L11 +1  ; jj++){
for (int l = 0; l < L  ; l++) {
f = (Qz - z3)*arr6[k-1][j][0][0][jj][kk][l] + (Wz - Qz)*arr6[k-1][j][0][0][jj][kk][l+1]  ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*( arr6[k-2][j][0][0][jj][kk][l] - (om/gam2)*arr6[k-2][j][0][0][jj][kk][l+1] ) ;
}
f = f + ( (kk)/( 2.0*(gam1 + gam2) ) )*arr6[k-1][j][0][0][jj][kk-1][l+1] ;
arr6[k][j][0][0][jj][kk][l] = f;
}
}
}
}


//  abc | 0jk
if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {
f=0;
for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
for (int l = 0; l < L ; l++) {

f = (Qz - z3)*arr6[k-1][j][0][ii][jj][kk][l] + (Wz - Qz)*arr6[k-1][j][0][ii][jj][kk][l+1] ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*(arr6[k-2][j][0][ii][jj][kk][l] - (om/gam2)*arr6[k-2][j][0][ii][jj][kk][l+1] ) ;
}
f = f + ( (kk)/( 2.0*(gam1 + gam2) ) )*arr6[k-1][j][0][ii][jj][kk-1][l+1] ;
arr6[k][j][0][ii][jj][kk][l] = f;
}
}
}
}
}


}
}
//closing if j >0
}

if (  (L21 > 0 )  &&   (L20 > 0 )      ){  

for(int k = 1; k < L22 + 1; k++){
    for(int j = 1; j < L21 + 1 ; j++){
      for(int i = 1; i < L20 + 1 ; i++){

// caso inicial: 000 | ijk
for (int l = 0; l < L ; l++) {
f=0;
// f = (Qz - z3)*arr6[k-1][j][i][0][0][0][l] - (Wz - Qz)*arr6[k-1][j][i][0][0][0][l+1] ;
f = (Qz - z3)*arr6[k-1][j][i][0][0][0][l] + (Wz - Qz)*arr6[k-1][j][i][0][0][0][l+1] ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*(arr6[k-2][j][i][0][0][0][l] - (om/gam2)*arr6[k-2][j][i][0][0][0][l+1] ) ;
}
arr6[k][j][i][0][0][0][l] = f;
}

// a00 | ijk
if(  L10 > 0  ){  
f= 0;
for(int ii = 1; ii < L10 + 1 ; ii++){
  for (int l = 0; l < L  ; l++) {

f = (Qz - z3)*arr6[k-1][j][i][ii][0][0][l] + (Wz - Qz)*arr6[k-1][j][i][ii][0][0][l+1] ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*( arr6[k-2][j][i][ii][0][0][l] - (om/gam2)*arr6[k-2][j][i][ii][0][0][l+1] ) ;
}
arr6[k][j][i][ii][0][0][l] = f;
}
}
}

// 0b0 | ijk
if (  L11 > 0   ){
f = 0;
for(int jj = 1; jj < L11 + 1 ; jj++){
 for (int l = 0; l < L  ; l++) {

f = (Qz - z3)*arr6[k-1][j][i][0][jj][0][l] + (Wz - Qz)*arr6[k-1][j][i][0][jj][0][l+1] ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*( arr6[k-2][j][i][0][jj][0][l] - (om/gam2)*arr6[k-2][j][i][0][jj][0][l+1] ) ;
}
arr6[k][j][i][0][jj][0][l] = f;
}
}
}

// 00c | ijk
if (  L12 > 0   ){
f = 0;
for(int kk = 1; kk < L12 + 1 ; kk++){
 for (int l = 0; l < L  ; l++) {
f = (Qz- z3)*arr6[k-1][j][i][0][0][kk][l] + (Wz - Qz)*arr6[k-1][j][i][0][0][kk][l+1] ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*( arr6[k-2][j][i][0][0][kk][l] - (om/gam2)*arr6[k-2][j][i][0][0][kk][l+1] ) ;
}
f = f + ( (kk)/( 2.0*(gam1 + gam2) ) )*arr6[k-1][j][i][0][0][kk-1][l+1] ;
arr6[k][j][i][0][0][kk][l] = f;

}
}
}

// ab0 | ijk
if(   (L10 >= 1) && (  L11 >= 1 )  ) {
f=0;
for(int jj = 1; jj < L11 + 1 ; jj++){
  for(int ii = 1; ii < L10 + 1 ; ii++){
    for (int l = 0; l < L  ; l++) { 
//  f = (Qz - z3)*arr6[k-1][j][i][ii][jj][0][l] - (Wz - Qz)*arr6[k-1][j][i][ii][jj][0][l+1] ;
f = (Qz - z3)*arr6[k-1][j][i][ii][jj][0][l] + (Wz - Qz)*arr6[k-1][j][i][ii][jj][0][l+1] ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*( arr6[k-2][j][i][ii][jj][0][l] - (om/gam2)*arr6[k-2][j][i][ii][jj][0][l+1] ) ;
}
arr6[k][j][i][ii][jj][0][l] = f;
}
}
}
}

//  a0c | ijk
if(   (L10 >= 1) && (  L12 >= 1 )  ) {
f=0;  
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int ii = 1; ii < L10 +1  ; ii++){
for (int l = 0; l < L  ; l++) {
//  f = (Qz - z3)*arr6[k-1][j][i][ii][0][kk][l] - (Wz - Qz)*arr6[k-1][j][i][ii][0][kk][l+1];
f = (Qz - z3)*arr6[k-1][j][i][ii][0][kk][l] + (Wz - Qz)*arr6[k-1][j][i][ii][0][kk][l+1];
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*( arr6[k-2][j][i][ii][0][kk][l] - (om/gam2)*arr6[k-2][j][i][ii][0][kk][l+1] ) ;
}
f = f + ( (kk)/( 2.0*(gam1 + gam2) ) )*arr6[k-1][j][i][ii][0][kk-1][l+1] ;
arr6[k][j][i][ii][0][kk][l] = f;
}
}
}
}

 //  0bc | ijk
if(   (L11 >= 1) && (  L12 >= 1 )  ) {
f=0; 
for(int kk = 1; kk < L12 +1 ; kk++){
  for(int jj = 1; jj < L11 +1  ; jj++){
for (int l = 0; l < L  ; l++) {

f = (Qz - z3)*arr6[k-1][j][i][0][jj][kk][l] + (Wz - Qz)*arr6[k-1][j][i][0][jj][kk][l+1]  ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*( arr6[k-2][j][i][0][jj][kk][l] - (om/gam2)*arr6[k-2][j][i][0][jj][kk][l+1] ) ;
}
f = f + ( (kk)/( 2.0*(gam1 + gam2) ) )*arr6[k-1][j][i][0][jj][kk-1][l+1] ;
arr6[k][j][i][0][jj][kk][l] = f;
}
}
}
}

//  abc | ijk
if(   (L10 >= 1) && (  L11 >= 1 )  &&  (L12 >= 1 )  ) {
f=0;
for(int kk = 1; kk < L12 + 1 ; kk++){
  for(int jj = 1; jj < L11 + 1 ; jj++){
   for(int ii = 1; ii < L10 + 1 ; ii++){
for (int l = 0; l < (L - k); l++) {
// f = (Qz - z3)*arr6[k-1][j][i][ii][jj][kk][l] + (Wz - Qz)*arr6[k-1][j][i][ii][jj][kk][l+1] ;
f = (Qz - z3)*arr6[k-1][j][i][ii][jj][kk][l] + (Wz - Qz)*arr6[k-1][j][i][ii][jj][kk][l+1] ;
if( (k-1) > 0  ){
f = f + ( (k-1)/(2.0*gam2))*(arr6[k-2][j][i][ii][jj][kk][l] - (om/gam2)*arr6[k-2][j][i][ii][jj][kk][l+1] ) ;
}
f = f + ( (kk)/( 2.0*(gam1 + gam2) ) )*arr6[k-1][j][i][ii][jj][kk-1][l+1] ;
arr6[k][j][i][ii][jj][kk][l] = f;
}
}
}
}
}


}
}
}
// closing if i,j > 0
}


Norminter = arr6[L22][L21][L20][L10][L11][L12][0];  

for (int i = 0; i < L22 + 1; i++) {
    for (int j = 0; j < L21 + 1; j++) {
        for (int k = 0; k < L20 + 1; k++) {
            for (int l = 0; l < L10 + 1; l++) {
                for (int m = 0; m < L11 + 1; m++) {
                    for (int n = 0; n < L12 + 1; n++) {
                        delete[] arr6[i][j][k][l][m][n];
                    }
                    delete[] arr6[i][j][k][l][m];
                }
                delete[] arr6[i][j][k][l];
            }
            delete[] arr6[i][j][k];
        }
        delete[] arr6[i][j];
    }
    delete[] arr6[i];
}
delete[] arr6;


return Norminter;

// cierre de 5 y 6
}

// cierre de 4
}

// cierre arr0
}

// cierre L2x > 0
}

return Norminter;
}




G4double Electron_electron_repulsion::selector( G4double alf_1, G4double alf_2, G4double alf_3,
G4double alf_4, G4int fv1a, G4int fv1b, G4int fv1c,
G4int fv2a, G4int fv2b, G4int fv2c, G4int fv3a, G4int fv3b,
G4int fv3c, G4int fv4a, G4int fv4b, G4int fv4c, 
G4double x1, G4double y1, G4double z1, G4double x2,
G4double y2, G4double z2, G4double x3, G4double y3, G4double z3,
G4double x4, G4double y4, G4double z4, G4double N1, G4double N2,
G4double N3, G4double N4, G4int cont) {


G4double Jep = 0;

//check if abmc0:
if(   (fv4a == 0) &&   (fv4b == 0) &&   (fv4c == 0)  ){

//check if a0mc0:
if(   (fv2a == 0) &&   (fv2b == 0) &&   ( fv2c == 0)  ){

Jep = Jep + repulsion_a0mc0(alf_1, alf_2, alf_3, alf_4, fv1a, fv1b, fv1c, fv3a, fv3b, fv3c, 
x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, N1, N2, N3, N4, cont);
return Jep;

// so abmc0:
}else{

if (    (fv2a > 0)  &&   (fv2b == 0)     &&   (fv2c == 0)     ){


fv1a=fv1a+1;
fv2a=fv2a-1;

Jep = selector(alf_1,  alf_2, alf_3, alf_4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
x1, y1, z1, x2, y2, z2,  x3, y3,  z3,
x4, y4, z4, N1,  N2, N3, N4, cont );

if ( fabs(x1 - x2) > 0 ){
//  al*xyz,000m;
fv1a=fv1a - 1;

Jep = Jep + ( x1 - x2 )*selector(alf_1,  alf_2, alf_3, alf_4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
x1, y1, z1, x2, y2, z2,  x3, y3,  z3,
x4, y4, z4, N1,  N2, N3, N4 , cont);

fv1a=fv1a + 1 ;

}

return Jep;

}else{
//0h0
if (    (fv2a == 0)  &&   (fv2b > 0)     &&   (fv2c == 0)     ){

fv1b = fv1b + 1;
fv2b = fv2b - 1;

Jep = selector(alf_1,  alf_2, alf_3, alf_4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
x1, y1, z1, x2, y2, z2,  x3, y3,  z3,
x4, y4, z4, N1,  N2, N3, N4 , cont);

if ( fabs(y1 - y2) > 0 ){

fv1b =fv1b -1;

Jep = Jep + ( y1 - y2 )*selector(alf_1,  alf_2, alf_3, alf_4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
x1, y1, z1, x2, y2, z2,  x3, y3,  z3,
x4, y4, z4, N1,  N2, N3, N4, cont );

fv1b =fv1b+1;

}

return Jep;


}else{
// 00v
if (   (fv2a == 0)  &&   (fv2b == 0)     &&   (fv2c > 0)      ){

fv1c =fv1c+1;
fv2c =fv2c-1;

Jep = selector(alf_1,  alf_2, alf_3, alf_4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
x1, y1, z1, x2, y2, z2,  x3, y3,  z3,
x4, y4, z4, N1,  N2, N3, N4 , cont);


if ( fabs(z1 - z2) > 0 ){

fv1c=fv1c-1;

Jep = Jep + ( z1 - z2 )*selector(alf_1,  alf_2, alf_3, alf_4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
x1, y1, z1, x2, y2, z2,  x3, y3,  z3,
x4, y4, z4, N1,  N2, N3, N4, cont );

fv1c = fv1c+1;

}

return Jep;


}else{

if (   (fv2a > 0)  &&   (fv2b > 0)     &&   (fv2c == 0)    ){

fv1a = fv1a + 1;
fv2a = fv2a -1;

Jep = selector(alf_1,  alf_2, alf_3, alf_4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
x1, y1, z1, x2, y2, z2,  x3, y3,  z3,
x4, y4, z4, N1,  N2, N3, N4, cont );

if ( fabs(x1 - x2) > 0 ){

fv1a = fv1a-1;

Jep = Jep + ( x1 - x2 )*selector(alf_1,  alf_2, alf_3, alf_4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
x1, y1, z1, x2, y2, z2,  x3, y3,  z3,
x4, y4, z4, N1,  N2, N3, N4, cont );

fv1a=fv1a+1;

}

}else{

if (   (fv2a == 0)  &&   (fv2b > 0)     &&   (fv2c > 0)       ){

fv1b= fv1b + 1;
fv2b= fv2b - 1;


Jep = selector(alf_1,  alf_2, alf_3, alf_4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
x1, y1, z1, x2, y2, z2,  x3, y3,  z3,
x4, y4, z4, N1,  N2, N3, N4 , cont );

if ( fabs(y1 - y2) > 0 ){

fv1b = fv1b - 1;

Jep = Jep + ( y1 - y2 )*selector(alf_1,  alf_2, alf_3, alf_4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
x1, y1, z1, x2, y2, z2,  x3, y3,  z3,
x4, y4, z4, N1,  N2, N3, N4, cont );

fv1b =fv1b + 1;

}


}else{

if (   (fv2a > 0)  &&   (fv2b == 0)     &&   (fv2c > 0)       ){

fv1a=fv1a+1;
fv2a=fv2a-1;


Jep = Jep + selector(alf_1,  alf_2, alf_3, alf_4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
x1, y1, z1, x2, y2, z2,  x3, y3,  z3,
x4, y4, z4, N1,  N2, N3, N4, cont );

if ( fabs(x1 - x2) > 0 ){
//  al*xyz,000m;
fv1a = fv1a-1;

Jep = Jep + ( x1 - x2 )*selector(alf_1,  alf_2, alf_3, alf_4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
x1, y1, z1, x2, y2, z2,  x3, y3,  z3,
x4, y4, z4, N1,  N2, N3, N4, cont );

fv1a = fv1a+1;

}

}else{

 // if    (fv2[0] > 0)  &&   (fv2[1] > 0)     &&   (fv2[2] > 0)    

fv1a=fv1a+1;
fv2a=fv2a-1;

Jep = Jep +  selector(alf_1,  alf_2, alf_3, alf_4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
x1, y1, z1, x2, y2, z2,  x3, y3,  z3,
x4, y4, z4, N1,  N2, N3, N4, cont );

if ( fabs(x1 - x2) > 0 ){
//  al*xyz,000m;
fv1a=fv1a-1;

Jep = Jep + ( x1 - x2 )*selector(alf_1,  alf_2, alf_3, alf_4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
x1, y1, z1, x2, y2, z2,  x3, y3,  z3,
x4, y4, z4, N1,  N2, N3, N4, cont );

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

//  
if (    (fv4a > 0) && (fv4b == 0) && (fv4c == 0)     ){
// a | b | c |n00
fv3a=fv3a+1;
fv4a=fv4a-1;

Jep = Jep + selector(alf_1,  alf_2, alf_3, alf_4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
x1, y1, z1, x2, y2, z2,  x3, y3,  z3,
x4, y4, z4, N1,  N2, N3, N4, cont );

if ( fabs(x3 - x4) > 0 ){

fv3a = fv3a-1;

Jep = Jep + ( x3 - x4 )*selector(alf_1,  alf_2, alf_3, alf_4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
x1, y1, z1, x2, y2, z2,  x3, y3,  z3,
x4, y4, z4, N1,  N2, N3, N4, cont );

fv3a = fv3a+1;

}

}else{

if (    (fv4a == 0)  &&   (fv4b > 0)     &&   (fv4c == 0)     ){
// a | b | c |0m0
fv3b = fv3b+1;
fv4b = fv4b-1;

Jep = Jep + selector(alf_1,  alf_2, alf_3, alf_4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
x1, y1, z1, x2, y2, z2,  x3, y3,  z3,
x4, y4, z4, N1,  N2, N3, N4 ,cont );

if ( fabs(y3 - y4) > 0 ){

fv3b=fv3b-1;

Jep = Jep + (y3 - y4)*selector(alf_1,  alf_2, alf_3, alf_4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
x1, y1, z1, x2, y2, z2,  x3, y3,  z3,
x4, y4, z4, N1,  N2, N3, N4,cont  );

fv3b=fv3b+1;

}

}else{

if (    (fv4a == 0)  &&   (fv4b == 0)     &&   (fv4c > 0)     ){
// a | b | c |00k
fv3c=fv3c+1;
fv4c=fv4c-1;

Jep = Jep + selector(alf_1,  alf_2, alf_3, alf_4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
x1, y1, z1, x2, y2, z2,  x3, y3,  z3,
x4, y4, z4, N1,  N2, N3, N4, cont );

if ( fabs(z3 - z4) > 0 ){

fv3c =fv3c-1;

Jep = Jep + (z3 - z4)*selector(alf_1,  alf_2, alf_3, alf_4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
x1, y1, z1, x2, y2, z2,  x3, y3,  z3,
x4, y4, z4, N1,  N2, N3, N4, cont );

fv3c =fv3c+1;

}


}else{


if (   (fv4a > 0)  &&  (fv4b > 0)  &&  (fv4c == 0)       ){

// ab |c nh0
fv3a =fv3a +1;
fv3a =fv4a -1;

Jep = Jep + selector(alf_1,  alf_2, alf_3, alf_4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
x1, y1, z1, x2, y2, z2,  x3, y3,  z3,
x4, y4, z4, N1,  N2, N3, N4, cont );

if ( fabs(x3 - x4) > 0 ){

fv3a =fv3a-1;

Jep = Jep + (x3 - x4)*selector(alf_1,  alf_2, alf_3, alf_4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
x1, y1, z1, x2, y2, z2,  x3, y3,  z3,
x4, y4, z4, N1,  N2, N3, N4, cont );

fv3a=fv3a+1;

}

}else{

if (   (fv4a == 0)  &&  (fv4b > 0)  &&  (fv4c > 0)   ){
// abc | 0hm
fv3b=fv3b+1;
fv3b=fv4b-1;

Jep = Jep + selector(alf_1,  alf_2, alf_3, alf_4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
x1, y1, z1, x2, y2, z2,  x3, y3,  z3,
x4, y4, z4, N1,  N2, N3, N4,  cont );

if ( fabs(y3 - y4) > 0 ){

fv3b=fv3b-1;

Jep = Jep + (y3 - y4)*selector(alf_1,  alf_2, alf_3, alf_4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
x1, y1, z1, x2, y2, z2,  x3, y3,  z3,
x4, y4, z4, N1,  N2, N3, N4 ,cont );

fv3b =fv3b +1;

}

}else{

if (   (fv4a > 0)  &&  (fv4b== 0)  &&  (fv4c > 0)   ){
// abc | n0m
fv3a=fv3a+1;
fv3a=fv4a-1;

if ( fabs(x3 - x4) > 0 ){

fv3a= fv3a - 1;

Jep = Jep + (x3 - x4)*selector(alf_1,  alf_2, alf_3, alf_4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
x1, y1, z1, x2, y2, z2,  x3, y3,  z3,
x4, y4, z4, N1,  N2, N3, N4, cont );

fv3a =fv3a+1;

}

}else{
// abc | nhm

fv3a=fv3a+1;
fv3a=fv4a-1;

if ( fabs(x3 - x4) > 0 ){

fv3a=fv3a-1;

Jep = Jep + (x3 - x4)*selector(alf_1,  alf_2, alf_3, alf_4, fv1a, fv1b, fv1c, 
fv2a, fv2b, fv2c, fv3a, fv3b, fv3c, fv4a, fv4b, fv4c,
x1, y1, z1, x2, y2, z2,  x3, y3,  z3,
x4, y4, z4, N1,  N2, N3, N4, cont );

fv3a = fv3a+1;

}

}

}

}

}

}

}

}

return Jep;
}



vector<vector<vector<vector<G4double>>>> Electron_electron_repulsion::pairing(Molecule& Mol ){

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
G4double buf_alp_4=0;
G4double coffs=0;
G4double contcofs_inter=0;
G4double J=0;
G4int cont = 0;
G4int orbital_principal=0;
G4int orbital_secundario=0;
G4int orbital_terciario=0;
G4int orbital_cuaternario=0;

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

  buf_alp_1 = (*(*atr_left_main)->shell1).NAC[i].A;

  for (int j = 0; j < 3; j++) {

  buf_alp_2 = (*(*atr_left_row)->shell2).NAC[j].A;

  for (int k = 0; k < 3; k++) {

  buf_alp_3 = (*(*atr_right)->shell3).NAC[k].A;

  for (int l = 0; l < 3; l++) {

  buf_alp_4 = (*(*atr_right_row)->shell4).NAC[l].A;

  coffs = ((*(*atr_left_main)->shell1).NAC[i].C)*((*(*atr_left_row)->shell2).NAC[j].C)*
  ((*(*atr_right)->shell3).NAC[k].C)*((*(*atr_right_row)->shell4).NAC[l].C);

 J = Electron_electron_repulsion::selector(buf_alp_1, buf_alp_2, buf_alp_3, buf_alp_4,
  (*(*atr_left_main)->shell1).l[0] , (*(*atr_left_main)->shell1).l[1] , (*(*atr_left_main)->shell1).l[2] , 
  (*(*atr_left_row)->shell2).l[0], (*(*atr_left_row)->shell2).l[1], (*(*atr_left_row)->shell2).l[2],
  (*(*atr_right)->shell3).l[0],  (*(*atr_right)->shell3).l[1],  (*(*atr_right)->shell3).l[2],  
  (*(*atr_right_row)->shell4).l[0],  (*(*atr_right_row)->shell4).l[1],  (*(*atr_right_row)->shell4).l[2],
  (*atr_left_main)->fX, (*atr_left_main)->fY,(*atr_left_main)->fZ,
  (*atr_left_row)->fX,(*atr_left_row)->fY, (*atr_left_row)->fZ,
  (*atr_right)->fX, (*atr_right)->fY, (*atr_right)->fZ,
  (*atr_right_row)->fX, (*atr_right_row)->fY, (*atr_right_row)->fZ,
  (*(*atr_left_main)->shell1).NAC[i].N, (*(*atr_left_row)->shell2).NAC[j].N,
  (*(*atr_right)->shell3).NAC[k].N, (*(*atr_right_row)->shell4).NAC[l].N , cont); 


contcofs_inter=J*coffs+contcofs_inter;


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


contcofs_inter = 0.;
cont = cont + 1 ;

}else{


}


}else{

//

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

