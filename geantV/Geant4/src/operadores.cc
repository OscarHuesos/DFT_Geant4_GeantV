
#include "operadores.hh"
#include <G4Pow.hh>


Operators::Operators(){
}

G4double Operators::doublefactorial(int n){

if (n == 0 || n==1 || n==-1){
return 1;
}else{
return n*doublefactorial(n-2);
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


G4double Operators::faux(G4int k, G4int l1,
G4int l2, G4double PA, G4double PB){

G4double f1,f2;
G4double p=0;
for(int i = 0; i < l1+1 ; i++ ){
  for(int j = 0; j< l2+1 ; j++ ){
    if(i+j==k){
f1=binomial_coefficient(l1,i);
f2=binomial_coefficient(l2,j);
p=p+f1*f2*pow(PA,(l1-i))*pow(PB,(l2-j)) ;

}
}
}

return p;
}



G4double Operators::n_func( G4int a, G4int b, G4double RN,  G4double Rr,
G4double A, G4double B, G4double gam, G4double x, G4int bander ){

G4double eta=1.0;
Operators fOpx;
G4double t = (x+1.0)/2.0;

if(a==0 && b==0){
return 1.0;
}else{
if(a==1 && b==0){

eta = -(A-RN + t*t*Rr);

return eta;

} else{
if( b==0){

eta = -(A-RN + t*t*Rr)*n_func(a-1, 0, RN, Rr,A ,B ,gam, x,0)+
((a-1.0)/(2.0*gam))*(1.0-t*t)*n_func(a-2,0, RN, Rr,A ,B ,gam, x, 0);


return eta;

}else{

eta = n_func(a+1,b-1, RN, Rr,A ,B ,gam, x, 0)+
(A-B)*n_func(a,b-1, RN, Rr,A ,B ,gam, x, 0);

return eta;
}
}
}
return eta;
}



G4double Operators::Oxab(G4double alpha_1, G4double alpha_2,
G4int k1, G4int m1, G4int n1,
G4int k2, G4int m2, G4int n2,
G4double x1, G4double y1, G4double z1,
G4double x2, G4double y2, G4double z2){


G4double gam=alpha_1+alpha_2;
G4double prd=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)
+(z1-z2)*(z1-z2);
G4double gpc1=(alpha_1*x1+alpha_2*x2)/gam;
G4double gpc2=(alpha_1*y1+alpha_2*y2)/gam;
G4double gpc3=(alpha_1*z1+alpha_2*z2)/gam;
G4double res=pow(pi/gam,1.5)*G4Exp(-(alpha_1*alpha_2*prd)/gam);


G4double X=0;
G4double PAX=gpc1-x1;
G4double PBX=gpc1-x2;
G4int rangx= floor((k1+k2)/2.0);

for(int x = 0; x < 1 + rangx ; x++ ){

X=X+(faux(2*x, k1, k2, PAX, PBX))*
((doublefactorial(2*x-1))/(pow(2.0*gam,x)));

}


G4double Y=0;
G4double PAY=gpc2-y1;
G4double PBY=gpc2-y2;
G4int rangy= floor((m1+m2)/2.0);

for(int y = 0; y < 1+rangy ; y++ ){
//for(br = 0; br<m2+1; br++ ){
Y=Y+(faux(2*y,m1, m2, PAY, PBY))*
((doublefactorial(2*y-1))/(pow(2*gam,y)));

}


G4double Z=0;
G4double PAZ=gpc3-z1;
G4double PBZ=gpc3-z2;
G4int rangz= floor((n1+n2)/2.0);

for(int z=0; z< 1+rangz; z++ ){
Z=Z+(faux(2*z,n1, n2, PAZ, PBZ))*
((doublefactorial(2*z-1))/(pow(2*gam,z) ) );

}


return res*X*Y*Z;
}


G4double Operators::xni(G4int n, G4int i){

G4double s1=  sin(i*pi/(n+1));
G4double c1= cos(i*pi/(n+1));
G4double av= (n+1.-2*i)/(n+1);
G4double x = av+(2.0/pi)*(1+2*s1*s1/3.0)*c1*s1;

return x;
}

G4double Operators::Fm(G4double p, G4double X,
G4double Y, G4double Z, G4double x){

return  G4Exp(-( p*((x+1.0)/2.0)*((x+1.0)/2.0)*(X*X+Y*Y+Z*Z) ) );
}


G4double Operators::Incomplete_gamma(G4double u, G4double n){

G4int IT = 100;
G4double EPS = 0.0000003;
G4double FPMIN = 0.0000000000000000000000000000000001;
G4double  b, g, s, c, d, h, del, an;
G4double gln =  lgamma(n); 

if( u < ( n + 1.0)  ){

if (u <= 0.0) {

if (u < 0.0){
g = 0.0;
}

}else{

G4double ap = n;
del = s = 1.0/n;
for (int i = 1; i <=  IT ; i++) {
ap++;
del = del*(u/ap);
s = s + del;
c = fabs(s)*EPS;

if(  fabs( del )  <  c ){
g = s*G4Exp( -u + n*log(u) - gln );
break;
}

}


} 

}else{  

b = u + 1.0 - n;
c = 1.0/FPMIN;
d = 1.0/b;
h = d;

for (int i = 1; i <= IT ; i++) {
an = -i*(i - n);
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
g = h*G4Exp( -u + n*log(u) - gln );
break;
}

}

g = 1 - g;

}

return g*G4Exp(gln);
}


G4double Operators::Hermitian_rec( G4int a, G4int b, 
G4double diff, G4double alpha1, G4double alpha2, G4double gam,
G4int t){


G4double m = (alpha1*alpha2)/gam;
G4int c = a + b;
G4double r = 0;
G4double E1 = 0;
G4double E2 = 0;
G4double E3 = 0;

if(  (t < 0) || ( t  > c ) ){
return r;
}else{

if(  (t==0) && (a==0) && (b==0)  ){

r = 1.0;
return r;

}else{

if (b == 0){

E1 = Hermitian_rec(a-1, b, diff, alpha1, alpha2 , gam , t-1);
E2 = Hermitian_rec(a-1, b, diff, alpha1, alpha2 , gam , t);
E3 = Hermitian_rec(a-1, b, diff, alpha1, alpha2 , gam , t+1);

r = (1.0/(2.0*gam))*E1 - ( (m*diff)/alpha1 )*E2 + (t+1)*E3;
return r;

}else{

E1 = Hermitian_rec(a, b-1, diff, alpha1, alpha2 , gam , t-1);
E2 = Hermitian_rec(a, b-1, diff, alpha1, alpha2 , gam , t);
E3 = Hermitian_rec(a, b-1, diff, alpha1, alpha2 , gam , t+1);


r = (1.0/(2.0*gam))*E1 + ( (m*diff)/alpha2 )*E2 + (t+1)*E3;
return r;

}

}

}

return r;
}


//Coulomb auxiliar:
G4double Operators::Rec_nucl(G4double T, G4int i, G4int j,
G4int k, G4double Rx, G4double Ry, G4double Rz, G4double gam,
G4double alpha_1, G4double alpha_2, 
G4int L, bool polyn ){

G4double f = 0;

if(  (i==0) && (j==0) && (k==0)  ){

if (  (polyn == false)   && (L >= 0)  ){

if( T < 0.00000000001    ){
T =     0.00000000001;
}


f = f + pow(-2.0*gam, L)*(  0.5*pow( T , -L - 0.5 )*(  Incomplete_gamma(T , L + 0.5)    )  );

}else{

f = f + pow(-2.0*gam, L)*gauss_chebyshev(gam,  Rx, Ry, Rz, 1, 1, 1, 1, 1, 1, 
L, 0, 0, 0, 0, 0,  alpha_1, alpha_2, false);

}

}else{

if( (i==0) && (j==0) ){
if( k > 1){
f = f + (k-1)*Rec_nucl( T,  i, j, k-2, Rx, Ry, Rz, gam, alpha_1, alpha_2, L+1 , polyn );
}

f = f + Rz*Rec_nucl( T,  i, j, k-1, Rx, Ry, Rz, gam, alpha_1, alpha_2, L+1 , polyn );

}else{
if( i==0 ){
if( j > 1 ){
f = f + (j-1)*Rec_nucl( T,  i, j-2, k, Rx, Ry, Rz, gam, alpha_1, alpha_2, L+1 , polyn );
}

f = f + Ry*Rec_nucl( T,  i, j-1, k, Rx, Ry, Rz, gam, alpha_1, alpha_2, L+1 , polyn );


}else{
if( i > 1 ){
f = f + (i-1)*Rec_nucl( T,  i-2, j, k, Rx, Ry, Rz, gam, alpha_1, alpha_2, L+1 , polyn );
}

f = f + Rx*Rec_nucl( T,  i-1 , j, k, Rx, Ry, Rz, gam, alpha_1, alpha_2, L+1 , polyn );

}

}

}

return f;
}


//this is the master controller:
G4double Operators::Gamma_aprox( G4double gam,
G4double rx, G4double ry, G4double rz,
G4double x1, G4double y1, G4double z1,
G4double x2, G4double y2, G4double z2,
G4int fv10, G4int fv11, G4int fv12,
G4int fv20, G4int fv21, G4int fv22,
G4double alpha_1, G4double alpha_2, bool polyn, bool selector ){

G4double sal = 0;
G4int L = fv10 + fv11 + fv12 + fv20 + fv21 + fv22;

// select direct old method without new frame   
if( selector == true){

// select repulsion or attraction
if( polyn == true){

G4double gpc1 = (alpha_1*x1 + alpha_2*x2)/gam;
G4double gpc2 = (alpha_1*y1 + alpha_2*y2)/gam;
G4double gpc3 = (alpha_1*z1 + alpha_2*z2)/gam;

G4double Rx = gpc1 - rx;
G4double Ry = gpc2 - ry;
G4double Rz = gpc3 - rz;

G4double T = gam*(Rx*Rx + Ry*Ry + Rz*Rz);
G4double ex = 0.;
G4double ey = 0.;
G4double ez = 0.;
G4double r = 0.;

//selection method base on T ():
sal = 0;
// old gamma is true:
bool sus_polyn = false;

for (int i = 0; i < ( fv10 + fv20 +1) ; i++) {

  for (int j = 0; j <( fv11 + fv21 +1) ; j++) {

    for (int k = 0; k < ( fv12 + fv22 +1) ; k++) {

      ex = Hermitian_rec( fv10, fv20, x1 - x2, alpha_1, alpha_2, gam, i);
      ey = Hermitian_rec( fv11, fv21, y1 - y2, alpha_1, alpha_2, gam, j); 
      ez = Hermitian_rec( fv12, fv22, z1 - z2, alpha_1, alpha_2, gam, k); 

 r =  Rec_nucl( T,  i , j, k, Rx, Ry, Rz, gam, alpha_1, alpha_2, 0 , sus_polyn );
sal = sal + ex*ey*ez*r;
}
}
}

return sal;

}else{
//  direct incomplete gamma 
G4double T = gam*( rx*rx + ry*ry + rz*rz );

bool check_gam = true;

if (check_gam == true){  

if( T < 0.00000000001    ){
T =     0.00000000001;
}

sal = 0.5*pow(T , -L - 0.5 )*( Incomplete_gamma(T, L + 0.5) );

}else{

sal = gauss_chebyshev(gam, rx, ry, rz, 1, 1, 1, 1, 1, 1, 
fv10, fv11, fv12, fv20, fv21, fv22,  alpha_1, alpha_2, polyn);

}

return sal;
}

// old one selector false
}else{


if( polyn == true){

sal = gauss_chebyshev(gam, rx, ry, rz, x1, y1, z1, x2, y2, z2,  
fv10, fv11, fv12, fv20, fv21, fv22,  alpha_1,  alpha_2,  polyn );

}else{   

sal = gauss_chebyshev(gam,  rx, ry, rz, 1, 1, 1, 1, 1, 1, 
fv10, fv11, fv12, fv20, fv21, fv22,  alpha_1, alpha_2, polyn);

}

return sal;
}
  
return sal;
}



G4double Operators::gauss_chebyshev( G4double p,
G4double rx, G4double ry, G4double rz,
G4double x1, G4double y1, G4double z1,
G4double x2, G4double y2, G4double z2,
G4int fv10, G4int fv11, G4int fv12,
G4int fv20, G4int fv21, G4int fv22,
G4double alpha_1, G4double alpha_2, bool polyn){

G4double tol=0.0000000001;
G4double error=10;
G4int N=50000;
G4int L;
G4double Rx=0;
G4double Ry=0;
G4double Rz=0;

Operators fOpx;

G4double sal=0;
G4double v,c, s, gpc1, gpc2, gpc3, ab, it, t, tm;
G4double om=(16.0/(3.0*(2+1)))*(pow(sin(1.0*pi/(2+1)),4) );

G4double c0=cos(pi/6.0);
G4double s0=sin(pi/6.0);
G4int n=3;

G4double etav=1.0;
G4double vx=1;
G4double vy=1;
G4double vz=1;

G4double etavm=1.0;
G4double vxm=1;
G4double vym=1;
G4double vzm=1;

G4double eta0=1.0;
G4double n0x=1;
G4double n0y=1;
G4double n0z=1;

G4double eta1=1.0;
G4double nx1=1;
G4double ny1=1;
G4double nz1=1;

G4double eta2=1.0;
G4double nx2=1;
G4double ny2=1;
G4double nz2=1;


if (polyn==true){

gpc1=(alpha_1*x1+alpha_2*x2)/p;
gpc2=(alpha_1*y1+alpha_2*y2)/p;
gpc3=(alpha_1*z1+alpha_2*z2)/p;

Rx=gpc1-rx;
Ry=gpc2-ry;
Rz=gpc3-rz;

  n0x=fOpx.n_func(fv10, fv20,
  gpc1, Rx, x1, x2, p, 0, 0);
  n0y=fOpx.n_func(fv11, fv21,
  gpc2, Ry, y1, y2, p, 0, 0);
  n0z=fOpx.n_func(fv12, fv22,
  gpc3, Rz, z1, z2, p, 0, 0);

  eta0=n0x*n0y*n0z;

  nx1=fOpx.n_func(fv10, fv20,
  gpc1, Rx, x1, x2, p, xni(2,1), 0);
  ny1=fOpx.n_func(fv11, fv21,
  gpc2, Ry, y1, y2, p, xni(2,1), 0);
  nz1=fOpx.n_func(fv12, fv22,
  gpc3, Rz, z1, z2, p, xni(2,1), 0);

  eta1=nx1*ny1*nz1;

  nx2=fOpx.n_func(fv10, fv20,
  gpc1, Rx, x1, x2, p, -1.0*xni(2,1), 1);
  ny2=fOpx.n_func(fv11, fv21,
  gpc2, Ry, y1, y2, p, -1.0*xni(2,1), 0);
  nz2=fOpx.n_func(fv12, fv22,
  gpc3, Rz, z1, z2, p, -1.0*xni(2,1), 1);

  eta2=nx2*ny2*nz2;


}else{

if (polyn==false){

L= fv10+fv11+fv12+fv20+fv21+fv22;
t=(fOpx.xni(2,1)+1)/2.0;

eta0= pow(0.5,2*L);
eta1= pow(t,2*L);
eta2= pow(-1.0*t,2*L);

Rx=rx;
Ry=ry;
Rz=rz;

}
}

G4double F0=eta0*Fm(p, Rx, Ry, Rz, 0)/2.0;
G4double F1=eta1*Fm(p, Rx, Ry, Rz, xni(2,1))/2.0   ;
G4double F2=eta2*Fm(p, Rx, Ry, Rz, -1.0*xni(2,1) )/2.0;

G4double c1=s0;
G4double s1=c0;
G4double Fq = om*(F1+F2);
G4double C=Fq+F0;
int j=0;

it = (2.0*n*(1-j)+(j*4.0*n/3)-1);

while(  (error > tol) && (  it <= N  ) ){

j=1-j;
c1=j*c1+(1-j)*c0;
s1=j*s1+(1-j)*s0;
c0=j*c0+(1-j)*sqrt( (1.0+c0)/2.0);
s0=j*s0+(1-j)*s0/(2.0*c0);
c=c0;
s=s0;
for (int i=1; i <= n-1; i+=2) {

v=1.0+(2.0/(3*pi)*s*c*(3+2.0*s*s))-(1.0*i/n);

if( ceil(i+2*j) > (i+j)  ){

if(polyn==true){
  
vx=fOpx.n_func(fv10, fv20,
  gpc1, Rx, x1, x2, p, v, 0);
vy=fOpx.n_func(fv11, fv21,
  gpc2, Ry, y1, y2, p, v, 0);
vz=fOpx.n_func(fv12, fv22,
  gpc3, Rz, z1, z2, p, v, 0);

etav=vx*vy*vz;

vxm=fOpx.n_func(fv10, fv20,
  gpc1, Rx, x1, x2, p, -1.0*v, 0);
vym=fOpx.n_func(fv11, fv21,
  gpc2, Ry, y1, y2, p, -1.0*v, 0);
vzm=fOpx.n_func(fv12, fv22,
  gpc3, Rz, z1, z2, p, -1.0*v, 0);

etavm=vxm*vym*vzm;
}else{
if(polyn==false){

L= fv10+fv11+fv12+fv20+fv21+fv22;
t=(v+1)/2.0;
tm=(-1.0*v+1)/2.0;

etavm=pow(tm,2*L);
etav=pow(t,2*L);
Rx=rx;
Ry=ry;
Rz=rz;

}
}

C=C+((etavm*Fm(p, Rx, Ry, Rz,-1.0*v)/2.0)+
(etav*Fm(p, Rx, Ry, Rz,v)/2.0))*pow(s,4);

}
v=s;
s=s*c1+c*s1;
c=c*c1-v*s1;
}

n=(1+j)*n;
F0=F0+(1.0-j)*(C-Fq);
ab= (1-j)*(Fq-3.0*F0/2.0)+j*(C-2.0*Fq);
error=16.0*fabs(ab)/(3.0*n);
Fq=(1-j)*Fq+j*C;
it = (2*n*(1-j)+(j*4.0*n/3)-1);

}
sal= (16.0*Fq)/(3*n);

return sal;

}

