
#include "operadores.hh"
#include <G4Pow.hh>


using Double_v = geant::Double_v;
using Int_v = geant::Double_v;
using geant::IndexD_v;
using MaskD_v= geant::MaskD_v;
using MaskI_v= geant::MaskI_v;

Operators::Operators(){
}


VECCORE_ATT_HOST_DEVICE
Double_v Operators::NaivePowerVector(Double_v base, G4int ex){

MaskI_v p0 (ex == 0 );
if(MaskFull(p0)) {
return 1.0;
}else{
MaskI_v p1 (ex == 1 );

if(MaskFull(p1)) {
return base;
}else{

MaskI_v p2 (ex == 2 );
if(MaskFull(p2)) {
return base*base;
}else{

MaskI_v p3 (ex == 3 );
if(MaskFull(p3)) {
return base*base*base;
}else{
MaskI_v p4 (ex == 4 );
if(MaskFull(p4)) {
return base*base*base*base;
}else{

MaskI_v p5 (ex == 5 );
if(MaskFull(p5)) {
return base*base*base*base*base;
}else{

MaskI_v p6 (ex == 6 );
if(MaskFull(p6)) {
return base*base*base*base*base*base;
}else{

MaskI_v p7 (ex == 7 );
if(MaskFull(p7)) {
return base*base*base*base*base*base*base;
}else{

MaskI_v p8 (ex == 8 );
if(MaskFull(p8)) {
return base*base*base*base*base*base*base*base;


}else{
MaskI_v p9 (ex == 9 );
if(MaskFull(p9 )) {
return base*base*base*base*base*base*base*base*base;

}else{
MaskI_v p10 (ex == 10 );
if(MaskFull(p10 )) {
return base*base*base*base*base*base*base*base*base*base;


}else{

MaskI_v paro (ex > 10 );
if(MaskFull(p10 )) {
printf("naive power out of range\n");
return 1.0;
}

}

}

}
  
}

}

}

}

}

}  

}

}


}

G4double Operators::doublefactorial(G4int n){

if (n == 0 || n==1 || n==-1){
return 1.0;
}else{
return n*doublefactorial(n- 2);
}

}


Double_v Operators::doublefactorialVector(G4int n){

if (n == 0 || n==1 || n==-1){
return 1.0;
}else{
return n*doublefactorialVector(n-2);
}
}

VECCORE_ATT_HOST_DEVICE
Int_v factorialVector(Int_v n){
MaskI_v fac (n < 1);
if(MaskFull(fac)) {
return 1;
}else{
return factorialVector(n-1)*n;
}


}  

G4double Operators::binomial_coefficient(G4int n, G4int k) {
G4double r=1.0;

try{
if(n==0 && k>0){
throw printf("Error, n=0 \n");
}

if(n==k || k==0){return r;}

 G4Pow* g4pow = G4Pow::GetInstance();

r=(g4pow->factorial(n))/( (g4pow->factorial(k))*(g4pow->factorial(n-k))  ) ;

return r;
}
catch(char* pMsg) { std::cerr << std::endl << "Exception:" << pMsg << std::endl;

}

return r;
}



Double_v Operators::NfuncVector( G4int A, G4int B, Double_v RNN, Double_v RR,
Double_v PA, Double_v PB, Double_v GD, Double_v XX,  G4int bander){

Double_v Eta = 1.0;
Double_v tt = (XX+1.0)/2.0;

MaskI_v ab0 (A==0 && B==0);
MaskI_v a1b0 (A==1 && B==0);
MaskI_v b0 (B==0);

if (MaskFull(ab0)){
return Double_v(1.0);
}else{

if (MaskFull(a1b0)){
Eta = -(PA-RNN + tt*tt*RR);

if (bander == 1){
printf("eta con a=1 y b=0  \n");
}

return Eta;

}else{
if(MaskFull(b0)){

Eta = -(PA-RNN + tt*tt*RR)*NfuncVector(A-1,0, RNN,RR,PA,PB,GD, XX, 0)+
((A-1.0)/(2.0*GD))*(1.0-tt*tt)*NfuncVector(A-2,0,RNN,RR,PA,PB,GD,XX,0);



return Eta;

}else{

Eta = NfuncVector(A+1,B-1, RNN,RR,PA,PB,GD, XX, 0) +
(PA-PB)*NfuncVector(A,B-1,RNN,RR,PA,PB,GD,XX,0);

return Eta;

}
}
}

return Eta;
}



Double_v Operators::fauxV(G4int k, G4int k1, G4int k2, Double_v PA, Double_v PB){

G4double f1,f2;
Double_v P = 0;


for(G4int i = 0; i < k1+1 ; ++i ){
  for(G4int j = 0; j< k2+1 ; ++j ){

    if( i+ j ==k){

f1=binomial_coefficient(k1,i);
f2=binomial_coefficient(k2,j);


P=P+f1*f2*NaivePowerVector(PA, k1-i)* NaivePowerVector(PB, k2-j);


 }
 }
 }     

 return P;
  
}


Double_v Operators::OXABVectortestU(G4double Alpha1, Double_v Alpha2,
G4int k1, G4int m1, G4int n1, G4int k2, G4int m2, G4int n2,
Double_v X1, Double_v X2, Double_v Y1, Double_v Y2, Double_v Z1, Double_v Z2){

Double_v gam = Alpha1 + Alpha2;
Double_v prd = (X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2);

Double_v  G1 = (Alpha1*X1+ Alpha2*X2)/gam;
Double_v  G2 = (Alpha1*Y1+ Alpha2*Y2)/gam;
Double_v  G3 = (Alpha1*Z1+ Alpha2*Z2)/gam;


Double_v  gV = (pi/gam)*(vecCore::math::Sqrt(pi/gam))*(vecCore::math::Exp(-Alpha1*Alpha2*prd/gam));
Double_v X = 0;
Double_v PAX= G1-X1;
Double_v PBX= G1-X2;

G4int tx = vecCore::math::Floor((k1+k2)/2.0);
G4int ty = vecCore::math::Floor((m1+m2)/2.0);
G4int tz = vecCore::math::Floor((n1+n2)/2.0);

if(tx == 0) {

X= fauxV(0, k1, k2, PAX, PBX)*doublefactorialVector(-1)/(1);

}else{

if(tx < 10) {

for(G4int i=0; i < 1 + tx; ++i){
X  = X+(fauxV(2*i, k1, k2, PAX, PBX))*(doublefactorialVector(2*i-1))/(NaivePowerVector(2.0*gam,i));

}


}else{

}
}


Double_v Y = 0;
Double_v PAY= G2-Y1;
Double_v PBY= G2-Y2;


if(ty == 0) {
Y = fauxV(0, m1, m2, PAY, PBY)*doublefactorialVector(-1)/(1);
}else{
if(ty < 10) {

for(G4int j=0; j < 1 + ty; ++j){
Y  = Y+(fauxV(2*j, m1, m2, PAY, PBY))*((doublefactorialVector(2*j-1))/(NaivePowerVector(2.0*gam,j )));



}

}else{
//printf("XX \n");
}
}


Double_v Z = 0.;
Double_v PAZ= G3-Z1;
Double_v PBZ= G3-Z2;

if(tz == 0) {
Z = fauxV(0, n1, n2, PAZ, PBZ)*doublefactorialVector(-1)/(1);
}else{
if(tz < 10) {
for(G4int k=0; k < 1+ tz; ++k){
Z  = Z+(fauxV(2*k, n1, n2, PAZ, PBZ))*((doublefactorialVector(2*k-1))/(NaivePowerVector(2.0*gam,k )));

}

}else{
printf("sup Z \n");
}
}

return  gV*X*Y*Z;
}


Double_v Operators::XniVector(Int_v n, Int_v i){
Double_v s = vecCore::math::Sin(i*pi/(n+1));
Double_v c = vecCore::math::Cos(i*pi/(n+1));
Double_v A = (n+ 1.0-2*i)/(n+1);
Double_v P = A + (2.0/pi)*(1+2*s*s/3.0)*c*s;
return P;
}


Double_v Operators::FmVector(Double_v P, Double_v XX, Double_v YY,
Double_v ZZ, Double_v xp){
Double_v W=vecCore::math::Exp(-(P*((xp+1.0)/2.0)*((xp+1.0)/2.0)*(XX*XX+YY*YY+ZZ*ZZ) ) );
return W;
}



Double_v  Operators::Incomplete_gammaVector(Double_v U, Double_v n ){

Double_v R = 0;
Double_v G, AP ,D, S, DEL, C, FD, UL ,B, H, AN, DD, FC;
G4int IT = 100;
G4double EPS = 0.0000003;
G4double FPMIN = 0.0000000000000000000000000000000001;
G4double  b, g, s, c, d, h, del, an, ap;

Double_v Gln = vecCore::math::Log( vecCore::math::TGamma(n)   ); 

MaskD_v tipo_metodo = ( U <  (  n + 1.0  ) ); 
MaskD_v check_u_uno = ( U <= 0.0 ); 
MaskD_v check_u_dos = ( U < 0.0 ); 

if (MaskFull( tipo_metodo  )){

  if (MaskFull( check_u_uno  )){

    if(   MaskFull( check_u_dos )   ){
      G =  0.0;

      }else{

      for (int j = 0; j < kVecLenD - 1; ++j) {

        if( Get( U , j ) < 0.0 ) { 
        Set(  G  ,j, 0.0 );

        }

      }

    }

  }else{

  // select u > 0.
  if( MaskEmpty( check_u_uno ) ) {

  AP = n;
  DEL = S = 1.0/n;


  for (int i = 1; i <=  IT ; i++) {
  AP = AP + 1;
  DEL = DEL*(U/AP);
  S = S + DEL;
  C = (vecCore::math::Abs(S))*EPS;
  FD = vecCore::math::Abs(DEL);
   MaskD_v wer = ( FD < C ); 


  if(MaskFull( wer)   ) {

  G = S*( vecCore::math::Exp( -U + n*( vecCore::math::Log(U) ) - Gln) );
  break; 

  }
}// end for IT

}else{

G4double u;
G4double z;
G4double gl;

for (int j = 0; j < kVecLenD - 1; ++j) {

u = Get( U , j );
z = Get( n , j );
gl = Get(Gln , j); 

ap = z;
del = s = 1.0/z;

    for (int i = 1; i <=  IT ; i++) {
      ap++;
      del = del*(u/ap);
      s = s + del;
      c = fabs(s)*EPS;
        if(  fabs( del )  <  c ){
      //    G4double u = Get(U , j);
          G4double g = s*G4Exp( -1.0*u + z*log(u ) -  gl  );
          Set(G, j , g );
          break;

        }
 
  }// fin for IT

}


}// fin empty check uno

}// fin maskfull check uno

R = G*( vecCore::math::Exp(Gln) );

return R;

}else{

  if( MaskEmpty( tipo_metodo  ) ) {
    B = U + 1.0 - n;
    C = 1.0/FPMIN;
    D = 1.0/B;
    H = D;

    for (int i = 1; i <= IT ; i++) {

      AN = -i*(i-n);
      B = B + 2.0;
      D= AN*D + B;
      FD = vecCore::math::Abs(D);
      MaskD_v Fcond = ( FD < FPMIN  ); 

      if( MaskFull( Fcond )   ) {
       D = FPMIN ;
      }

      C = B + (AN/C);

      FC = vecCore::math::Abs(C);
      MaskD_v Ffc2 = ( FC < FPMIN  ); 

      if( MaskFull( Ffc2 )   ) {
       C = FPMIN ;
      }

      D   = 1.0/D;
      DEL = D*C;
      H   = H*DEL;

      MaskD_v EPCON = ( vecCore::math::Abs(DEL - 1.0)  < EPS ); 

      if( MaskFull( EPCON )   ) {

      G = H*( vecCore::math::Exp(-U + n*( vecCore::math::Log(U) ) - Gln) );
      break;
      }

    }// end for IT

    G = 1.0 - G;

    R = G*( vecCore::math::Exp(Gln) );

    return R;


}else{

G4double u;
G4double z;

  
for (int j = 0; j < kVecLenD - 1; ++j) {

u = Get( U , j );
z = Get( n , j );
G4double gl = Get( Gln , j );

if( u < ( z + 1.0)  ){

if (u <= 0.0) {

if (u < 0.0){
Set(G , j, 0.0);
//g = 0.0;
}

}else{

G4double ap = z;
G4double del = s = 1.0/z;


for (int i = 1; i <=  IT ; i++) {
ap++;
del = del*(u/ap);
s = s + del;
c = fabs(s)*EPS;

if(  fabs( del )  <  c ){
g = s*G4Exp( -u + z*log(u) - gl );
break;
}

}

}


}else{

b = u + 1.0 - z;
c = 1.0/FPMIN;
d = 1.0/b;
h = d;

for (int i = 1; i <= IT ; i++) {
an = -i*(i - z);
b = b + 2.0;
d = an*d + b;


if(  fabs( d )  <  FPMIN ){
d = FPMIN;
}

c = b + (an/c);

if(  fabs( c )  <  FPMIN ){
c = FPMIN;
}

d = 1.0/d;
del = d*c;
h = h*del;

if(  fabs( del - 1.0 )  <  EPS ){

g = h*G4Exp( -u + z*log(u) - gl );
break;
}

}

g = 1 - g;

}


Set (  R , j , g*G4Exp(gl) );

}// closing j vector 

return R;

}

}

}


//Hermitina gauss cofs:
Double_v Operators::Hermitian_recV( G4int a, G4int b, Double_v D, 
G4double alpha1, Double_v Alpha2, Double_v GAM, G4int t){

Double_v M = (alpha1*Alpha2)/GAM;
G4int c = a + b;
Double_v E1 = 0.;
Double_v E2 = 0.;
Double_v E3 = 0.;
Double_v R  = 0.;


if(  (t < 0) || ( t  > c ) ){

return R;

}else{


if(  (t==0) && (a==0) && (b==0)  ){

R = 1.0;
return R;

}else{

if( b == 0){  

E1 = Hermitian_recV(a-1, b, D, alpha1, Alpha2 , GAM , t-1);
E2 = Hermitian_recV(a-1, b, D, alpha1, Alpha2 , GAM , t);
E3 = Hermitian_recV(a-1, b, D, alpha1, Alpha2 , GAM , t+1);


R = (1.0/(2.0*GAM))*E1 - ( (M*D)/alpha1 )*E2 + (t+1)*E3;


return R;

}else{

E1 = Hermitian_recV(a, b-1, D, alpha1, Alpha2 , GAM , t-1);
E2 = Hermitian_recV(a, b-1, D, alpha1, Alpha2 , GAM , t);
E3 = Hermitian_recV(a, b-1, D, alpha1, Alpha2 , GAM , t+1);

R = (1.0/(2.0*GAM))*E1 + ( (M*D)/Alpha2 )*E2 + (t+1)*E3;

return R;

}


}

}

return R;
}



Double_v Operators::gauss_chebyshevVector(Double_v GAM, Double_v rx, Double_v ry, Double_v rz,
Double_v X1, Double_v Y1, Double_v Z1, Double_v X2, Double_v Y2, Double_v Z2, 
G4int fv10, G4int fv11, G4int fv12,
G4int fv20, G4int fv21, G4int fv22,
G4double alpha1, Double_v Alpha2, bool polyn){


G4double   tol =  0.0000000001;
Double_v error =  10;

G4int N = 50000;
Double_v Rx = 0;
Double_v Ry = 0;
Double_v Rz = 0;

Operators fOpx;

Double_v Sal=0;
Double_v VSIN=vecCore::math::Sin(1.0*pi/(2.0+1));
Double_v PW4= VSIN*VSIN*VSIN*VSIN;
Double_v Om=(16.0/(3.0*(2+1)))*PW4;


Double_v  GPC1, GPC2, GPC3, F0, F1, F2, IT, CC, SS, buf, ABV , TM, tV;

Double_v VC0 = vecCore::math::Cos(pi/6.0);
Double_v VS0 = vecCore::math::Sin(pi/6.0);


G4int n=3;

Double_v  Etav=1.0;
Double_v  Vx=1.0;
Double_v  Vy=1.0;
Double_v  Vz=1.0;

Double_v  Etavm=1.0;
Double_v  VXM=1.0;
Double_v  VYM=1.0;
Double_v  VZM=1.0;

Double_v  Eta0 = 1.0;
Double_v  N0X=1.0;
Double_v  N0Y=1.0;
Double_v  N0Z=1.0;

Double_v  Eta1 = 1.0;
Double_v  NX1=1.0;
Double_v  NY1=1.0;
Double_v  NZ1=1.0;

Double_v  Eta2 = 1.0;
Double_v  NX2=1.0;
Double_v  NY2=1.0;
Double_v  NZ2=1.0;



if (polyn==true){


GPC1 = ( alpha1*X1+Alpha2*X2)/GAM;
GPC2 = ( alpha1*Y1+Alpha2*Y2)/GAM;
GPC3 = ( alpha1*Z1+Alpha2*Z2)/GAM;


Rx = GPC1 - rx;
Ry = GPC2 - ry;
Rz = GPC3 - rz; 


N0X = NfuncVector(fv10, fv20, GPC1, Rx, X1, X2, GAM, 0, 0);

N0Y = NfuncVector(fv11, fv21, GPC2, Ry, Y1, Y2, GAM, 0, 0);

N0Z = NfuncVector(fv12, fv22, GPC3, Rz, Z1, Z2, GAM, 0, 0);

Eta0 = N0X*N0Y*N0Z;


NX1 = NfuncVector(fv10, fv20, GPC1, Rx, X1, X2, GAM,  XniVector(2,1)  , 0);

NY1 = NfuncVector(fv11, fv21, GPC2, Ry, Y1, Y2, GAM,  XniVector(2,1)  , 0);

NZ1 = NfuncVector(fv12, fv22, GPC3, Rz, Z1, Z2, GAM,  XniVector(2,1)  , 0);


Eta1 = NX1*NY1*NZ1;

NX2 = NfuncVector(fv10, fv20, GPC1, Rx, X1, X2, GAM,  -1.0*XniVector(2,1)  , 0);

NY2 = NfuncVector(fv11, fv21, GPC2, Ry, Y1, Y2, GAM,  -1.0*XniVector(2,1)  , 0);

NZ2 = NfuncVector(fv12, fv22, GPC3, Rz, Z1, Z2, GAM,  -1.0*XniVector(2,1)  , 0);

Eta2=NX2*NY2*NZ2;

}else{

if (polyn==false){

G4int L= fv10+fv11+fv12+fv20+fv21+fv22;

Double_v  tV = (  XniVector(2,1) + 1)/2.0;

Eta0 = NaivePowerVector(0.5 , 2*L);
Eta1 = NaivePowerVector(tV , 2*L);
Eta2 = NaivePowerVector(-1.0*tV , 2*L);

Rx=rx;
Ry=ry;
Rz=rz;

}

}

F0 =  Eta0*FmVector(GAM, Rx, Ry, Rz, 0)/2.0;
F1 =  Eta1*FmVector(GAM, Rx, Ry, Rz,  XniVector(2,1)  )/2.0;
F2 =  Eta2*FmVector(GAM, Rx, Ry, Rz,  -1.0*XniVector(2,1)  )/2.0;

Double_v VC1=VS0;
Double_v VS1=VC0;
Double_v VFQ = Om*(F1+F2);
Double_v CV = VFQ+F0;


G4int j=0;
IT = (2.0*n*(1-j)+(j*4.0*n/3)-1);

MaskD_v err_cond = (  error >  tol ); 
MaskD_v it_cond = ( IT <= N ); 


j=0;
n=3;


while(  !MaskEmpty( it_cond )  &&  !MaskEmpty(err_cond)   ){

j = 1 - j;

VC1 = j*VC1 + (1 - j)*VC0;
VS1 = j*VS1 + (1 - j)*VS0;
VC0 = j*VC0 + (1 - j)*vecCore::math::Sqrt(  (1.0+VC0)/2.0  );
VS0 = j*VS0 + (1 - j)*VS0/(2.0*VC0 );
CC  = VC0;
SS  = VS0;

for (int i=1; i <= n-1; i+=2) {
buf = 1.0 + (2.0/(3*pi)*SS*CC*(3.0+2.0*SS*SS) ) - (1.0*i/n) ;

if( ceil(i+2*j) > (i+j)  ){
if(polyn==true){

Vx = NfuncVector( fv10, fv20, GPC1, Rx, X1, X2, GAM, buf ,0);
Vy = NfuncVector( fv11, fv21, GPC2, Ry, Y1, Y2, GAM, buf, 0);
Vz = NfuncVector( fv12, fv22, GPC3, Rz, Z1, Z2, GAM, buf, 0);

Etav=Vx*Vy*Vz;

VXM = NfuncVector( fv10, fv20, GPC1, Rx, X1, X2, GAM, -1.0*buf, 0 );
VYM = NfuncVector( fv11, fv21, GPC2, Ry, Y1, Y2, GAM, -1.0*buf, 0);
VZM = NfuncVector( fv12, fv22, GPC3, Rz, Z1, Z2, GAM, -1.0*buf, 0 );

Etavm = VXM*VYM*VZM; 

}else{
if(polyn==false){

G4int L= fv10+fv11+fv12+fv20+fv21+fv22;
tV = (buf + 1.0)/2.0; 
TM = (-1.0*buf + 1.0)/2.0;

Etavm = NaivePowerVector(TM , 2*L);
Etav  = NaivePowerVector(tV , 2*L);

Rx=rx;
Ry=ry;
Rz=rz;

}
}


CV = CV +  ( (Etavm*FmVector(GAM, Rx, Ry, Rz, -1.0*buf )/2.0) + 
(Etav*FmVector( GAM, Rx, Ry, Rz, buf   )/2.0) )*NaivePowerVector(SS  ,4);

}

buf = SS;
SS = SS*VC1 + CC*VS1;
CC = CC*VC1 - buf*VS1;

}


n=(1 + j)*n;

F0 = F0 + (1.0 - j)*(CV - VFQ );
ABV = (1 - j)*(VFQ - 3.0*F0/2.0) + j*(CV - 2.0*VFQ );
error = 16.0*( vecCore::math::Abs(ABV)/(3.0*n) );

VFQ = (1 - j)*VFQ + j*CV;

IT = (2*n*(1-j)+(j*4.0*n/3)-1);

err_cond = ( error >  tol ); 
it_cond  = ( IT <= N );

}


Sal= (16.0*VFQ)/(3*n);

return Sal;
}



Double_v Operators::Rec_nuclearV(Double_v T, G4int i, G4int j,
G4int k, Double_v Rx, Double_v Ry, Double_v Rz, Double_v GAM,
G4double alpha1, Double_v Alpha2, G4int L, bool polyn ){

Double_v F = 0.;
Double_v W = 0.;
Double_v P = 0.;
Double_v TL = 0.;


if(  (i==0) && (j==0) && (k==0)  ){

W = NaivePowerVector(-2.0*GAM, L);

if (  polyn == false     ){


MaskD_v tt = ( T < 0.00000000001  ); 

if(   MaskFull( tt )     ){
T =     0.00000000001;
}else{

if(  MaskEmpty( tt )   ){

}else{

for (int l = 0; l < kVecLenD -1 ; ++l) {

if( Get( T , l ) < 0.00000000001 ) { 
Set(  T, l , 0.00000000001 );
}

}

}
}



TL = NaivePowerVector(T, L);
P = 0.5/( (vecCore::math::Sqrt(T))*TL  );

F = F + W*( P*Incomplete_gammaVector(T , L + 0.5)    );


}else{

F = F + W*gauss_chebyshevVector(GAM, Rx, Ry, Rz, 1,1,1,1,1,1,
L, 0,0,0,0,0, alpha1, Alpha2 , false);

}


}else{

if( (i==0) && (j==0) ){

if( k > 1){
F = F + (k-1)*Rec_nuclearV( T, i,j, k-2, Rx, Ry, Rz, GAM, alpha1, Alpha2, L+1, polyn);

}

F = F + Rz*Rec_nuclearV( T, i,j, k-1, Rx, Ry, Rz, GAM, alpha1, Alpha2, L+1, polyn);

}else{
if( i==0 ){
if( j > 1 ){

F = F + (j-1)*Rec_nuclearV( T, i, j-2, k, Rx, Ry, Rz, GAM, alpha1, Alpha2, L+1, polyn);
}

F = F +Ry*Rec_nuclearV( T, i, j-1, k, Rx, Ry, Rz, GAM, alpha1, Alpha2, L+1, polyn);

}else{
if( i > 1 ){

F = F + (i-1)*Rec_nuclearV( T, i-2, j, k, Rx, Ry, Rz, GAM, alpha1, Alpha2, L+1, polyn);
}

F = F + Rx*Rec_nuclearV( T, i-1, j, k, Rx, Ry, Rz, GAM, alpha1, Alpha2, L+1, polyn);
}

}

}

return F;

}



Double_v Operators::Gamma_aproxV( Double_v GAM,
Double_v rx, Double_v ry, Double_v rz,
Double_v X1, Double_v Y1, Double_v Z1,
Double_v X2, Double_v Y2, Double_v Z2,
G4int fv10, G4int fv11, G4int fv12,
G4int fv20, G4int fv21, G4int fv22,
G4double alpha1, Double_v Alpha2, bool polyn, bool selector ){

Double_v S = 0;
G4int L = fv10 + fv11 + fv12 + fv20 + fv21 + fv22;
Double_v Rx, Ry, Rz;
// old method without new frame
if( selector == true){

// select repulsion or attraction
if( polyn == true){

Double_v Ex = 0.;
Double_v Ey = 0.;
Double_v Ez = 0.;
Double_v R = 0.;

Double_v GPC1 = ( alpha1*X1+Alpha2*X2)/GAM;
Double_v GPC2 = ( alpha1*Y1+Alpha2*Y2)/GAM;
Double_v GPC3 = ( alpha1*Z1+Alpha2*Z2)/GAM;

Rx = GPC1 - rx;
Ry = GPC2 - ry;
Rz = GPC3 - rz; 

Double_v T = GAM*(  Rx*Rx + Ry*Ry + Rz*Rz  );
S = 0;
bool sus_polyn = false;

for (int i = 0; i < ( fv10 + fv20 +1) ; i++) {

  for (int j = 0; j <( fv11 + fv21 +1) ; j++) {

    for (int k = 0; k < ( fv12 + fv22 +1) ; k++) {
    //  v = 0.0;
      Ex = Hermitian_recV( fv10, fv20, X1 - X2, alpha1, Alpha2, GAM, i);
      Ey = Hermitian_recV( fv11, fv21, Y1 - Y2, alpha1, Alpha2, GAM, j); 
      Ez = Hermitian_recV( fv12, fv22, Z1 - Z2, alpha1, Alpha2, GAM, k);

     R =  Rec_nuclearV( T,  i , j, k, Rx, Ry, Rz, GAM, alpha1, Alpha2, 0 , sus_polyn );
     S = S + Ex*Ey*Ez*R;
}
}
}

return S;

}else{

Double_v T = GAM*(  rx*rx + ry*ry + rz*rz  );

// condition to check to use hibrid or pure algortithm:
bool check_gam = true;

if (check_gam == true){ 

MaskD_v tt = ( T < 0.00000000001  ); 

if(   MaskFull( tt )     ){
T =     0.00000000001;
}else{

if (  MaskEmpty( tt )    ){

}else{

for (int k = 0; k < kVecLenD - 1 ; ++k ) {

if (  Get(T, k) < 0.00000000001   ){
Set( T, k,  0.00000000001 );
}

}
}

}


Double_v TL = NaivePowerVector(T, L);
Double_v P = 0.5/( (vecCore::math::Sqrt(T))*TL  );


S = P*Incomplete_gammaVector(T , L + 0.5);


}else{

S = gauss_chebyshevVector(GAM, Rx, Ry, Rz, 1,1,1,1,1,1,
fv10, fv11, fv12, fv20, fv21, fv22,  alpha1, Alpha2, polyn);

}


return S;
}

// old one selector false
}else{


if( polyn == true){

S = gauss_chebyshevVector(GAM, rx, ry, rz, X1, Y1, Z1,
X2, Y2, Z2, fv10, fv11, fv12, fv20, fv21, fv22, alpha1, Alpha2 , polyn);


}else{   

S = gauss_chebyshevVector(GAM, rx, ry, rz, 1, 1, 1,1, 1, 1,
 fv10, fv11, fv12, fv20, fv21, fv22, alpha1, Alpha2 , polyn);

}

return S;

}

  
return S;
}


