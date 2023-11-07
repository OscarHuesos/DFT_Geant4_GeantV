
#include "Funcionals.hh"

using namespace std;

Hybrid_Func::Hybrid_Func()
{
}


void B3LYP(  vector<G4double> n,
vector<G4double>& Pot, vector<G4double>& Ec ){

return;
}


GGA::GGA()
{
}

G4double GGA::Becke_correctionV( vector<Double_v> P, vector<Double_v>& PVXC, vector<Double_v>& EXCV , 
vector<Double_v> GRADV, vector<Double_v> LAP, vector<Double_v> GVG, vector<Double_v> W ){

G4double Acumv = 0;
G4double ec = 0;
G4double ex = 0;
G4int cont = 0;
G4int contv = 0;
G4int sz = P.size();

//VECTOR VARIBLES
Double_v R ,K, G, FF, NU, GDV, XRHO, CUBR, SINHR, VV, BECK, XX, DERSIN;
Double_v AT, LG, ER ,TAT, MA, TT, UU, UDOS, CHU, POL, H , OSO, ZZ, DLN;
Double_v EPXLDA, EPXGGA, EPCLDA, EPCGGA, EC, EX,  EXSC, DXDGDG, DCDGDG;
Double_v DEPXLDA, DVDX, DBECKX, DEPBE, DEPXGGAR, DLNX, DEPP, DEPCLDA ;
Double_v DUDEC, DUDUDEC, DECLDAPOL,  DHECLDA, DTDGDP, DTDP, DUDT, DHDT;
Double_v DEPSGGAR, DEXDG, DGGADG, DVV, DVDXDX, DODX, DSDX, DSIG, DPOLX;
Double_v XPDG, DXDGDP, DJDT, AUX, DLDT, CHE, CTZ, DNDT, DHDGDT, DCDGDP;

//SCALAR VARIABLES
for (int i = 0; i <  sz ; i++) {

MaskD_v mtol ( P[i] < tol );

if( MaskFull( mtol ) ) {

PVXC[i] = Double_v (0.0);
EXCV[i] = Double_v (0.0);

}else{

if( MaskEmpty( mtol ) ) {

R = 2.0*P[i];
GDV = 2.0*GRADV[i];

CUBR = vecCore::math::Cbrt( R );

// con constante Slter Slat
//EPXLDA  = (-1.0)*Slat*CUBR;

//con cnosnte paper becke:
EPXLDA = (-1.0)*StLDA*( CUBR  );

XRHO = GDV/( CUBR*R ) ;
SINHR = vecCore::math::ASinh( XRHO );

VV = 1.0 + 6.0*beta*XRHO*SINHR;
BECK = (beta*XRHO*XRHO)/VV;
XX = XRHO*XRHO + 1.0;
DERSIN = 1.0/( vecCore::math::Sqrt( XX ) ); 

EPXGGA = EPXLDA + CUBR*( -1.0*BECK  ) ;

K = vecCore::math::Cbrt(  3.0/(4.0*pi*R)  );
G = vecCore::math::Sqrt(K);
FF = 2.0*G + B;
NU = K + B*G + C;

AT = vecCore::math::ATan( Q/FF );
LG = vecCore::math::Log( NU ) ;

EPCLDA = A*( (vecCore::math::Log(K)) - LG + Bq*AT -
b1*( 2.0*(vecCore::math::Log(G - x0)) - LG + 2.0*m1*AT ) ) ;
EXSC = vecCore::math::Exp( -1.0*( EPCLDA/gamma ) );

TAT = Beta_Ec_div_gamma/( EXSC - 1.0  );
ER = vecCore::math::Cbrt( 3.0*pi*pi*R ) ;
MA = raizpi/( 4.0*R*( vecCore::math::Sqrt( ER ) ) );
TT = MA*GDV;
UU = TAT*TT*TT;
UDOS = TAT*TAT*TT*TT*TT*TT;
CHU = 1.0 + UU + UDOS;
POL = (1.0 + UU)/(CHU);
OSO = 1.0 + Beta_Ec_div_gamma*TT*TT*POL;
H = gamma*( vecCore::math::Log( OSO ) );
EPCGGA = EPCLDA + H;

//EC = EC + R*EPCGGA*W[i];
//EX = EX + R*EPXGGA*W[i];

EC =  R*EPCGGA*W[i];
EX =  R*EPXGGA*W[i];
EXCV[i] = R*( EPCGGA  +  EPXGGA )*W[i];

////////////////////////////////////////////////////////////s
// con constante Slater Slat:
// DEPXLDA =  (-1.0*Slat)/( 3.0*(CUBR*CUBR )  ) ;


//con cnosnte paper becke:
DEPXLDA =  (-1.0*StLDA)/( 3.0*(CUBR*CUBR )  ) ;


DVDX = 6.0*beta*( XRHO*DERSIN + SINHR );
DBECKX = beta*( 2.0*XRHO*VV - XRHO*XRHO*DVDX )/(VV*VV);
DEPBE =  ( (4.0*XRHO*DBECKX)/( 3.0*CUBR*CUBR ) ) - BECK/( 3.0*CUBR*CUBR ) ; 
DEPXGGAR =  DEPXLDA + DEPBE;

DLNX = FF/NU;
DEPP = A*( (2.0/G) - DLN - (4.0*B/( FF*FF + QQ )) -
b1*( (2.0/(G - x0) ) - DLN - (4.0*mm/( FF*FF + QQ) )  ) );

DEPCLDA = (-1.0)*( G/(6.0*R) )*DEPP;
ZZ = (-1.0)*( EXSC )/gamma;

DUDEC = ( (-1.0)*TT*TT*Beta_Ec_div_gamma*ZZ )/( ( EXSC - 1.0)*( EXSC - 1.0) );
DUDUDEC = (-2.0*Beta_Ec_div_gamma*Beta_Ec_div_gamma*TT*TT*TT*TT*ZZ)/( (EXSC - 1.0)*(EXSC - 1.0)*(EXSC - 1.0)  );
DECLDAPOL =  ( DUDEC*CHU - (1.0 + UU)*( DUDEC + DUDUDEC )  )/(CHU*CHU);
DHECLDA =ll*(TT*TT*DECLDAPOL/OSO);


DTDGDP = -1.0*q*( (q_prima/( R*( vecCore::math::Sqrt( ER ) )*(3.0*pi*pi*R) ) ) + 1.0/( R*R*( vecCore::math::Sqrt( ER ) ) ) );
DTDP = DTDGDP*GDV;
DUDT = ( 2.0*TAT*TT*CHU - (1.0 + UU)*(2.0*TAT*TT*(1.0 + TAT) ) )/( CHU*CHU );
DHDT = ll*(2.0*TT*POL + TT*TT*DUDT )/( OSO );

DEPSGGAR = DEPCLDA + (DHDT*DTDP + DHECLDA*DEPCLDA);

DEXDG = (-1.0/R)*DBECKX;
DGGADG = DEXDG + MA*DHDT;
DVV = 2.0*VV*DVDX;
DVDXDX = 6.0*beta*( (1.0/( vecCore::math::Sqrt( XX*XX*XX ) ) ) + DERSIN );
DODX = XRHO*XRHO*DVDXDX + 2.0*XRHO*DVDX;
DSDX = ( DODX*VV*VV - XRHO*XRHO*DVDX*DVV  )/( VV*VV*VV*VV  );
DSIG = 2.0*( XRHO*DVDX +VV )/( VV*VV );
DPOLX = beta*( DSIG - DSDX);
XPDG = DPOLX*XRHO/(R*CUBR ) + DBECKX/(R*CUBR) ;
DXDGDP = (1.0/( 3.0*CUBR*CUBR ) )*( 4.0*XPDG - (DBECKX/(R*CUBR) ) );
DJDT = 4.0*TAT*TT + 12.0*(UDOS/TT) + 12.0*(UDOS*UU/TT) + 8.0*(UDOS*UDOS/TT) ;
AUX = 1.0 + 2.0*UU + 3.0*UDOS + 2.0*UDOS*UU + UDOS*UDOS;

DLDT = 2.0*TAT*TAT*( ( 5.0*UU*TT*TT - 1.0 -3.0*UU )*AUX - ( UU*TT*TT*TT -TT - UU*TT )*DJDT )/(CHU*CHU*CHU*CHU);  
CHE = 2.0*( POL + TT*DUDT ) + 2.0*TT*DUDT + TT*TT*DLDT ;
CTZ = Beta_Ec_div_gamma*( 2.0*TT*POL + TT*TT*DUDT );
DNDT = (CHE*OSO - ( 2.0*TT*POL + TT*TT*DUDT )*CTZ  )/(OSO*OSO) ;
DHDGDT = DNDT*MA*ll;
DCDGDP = DHDGDT*DTDP + DHDT*DTDGDP ;
DXDGDG = ( -1.0/( R*R*CUBR ) )*DPOLX ;
DCDGDG = DHDGDT*MA;

///////////////////////////////////////////////////////////////////////////////////////////
PVXC[i] = ( EPCGGA + EPXGGA ) + R*( DEPXGGAR  +  DEPSGGAR) +
(R/(4.0*GRADV[i]*GRADV[i] ) )*DGGADG*4.0*GVG[i] - (R/GDV)*( DXDGDG + DCDGDG )*4.0*GVG[i] - 
(R/GDV)*DGGADG*( 2.0*LAP[i] ) - GDV*DGGADG - R*GDV*(  DXDGDP + DCDGDP );
///////////////////////////////////////////////////////////////////////////////////////////

PVXC[i] = PVXC[i]*W[i]; 

for (int j = 0; j < kVecLenD ; ++j) {

Acumv = Acumv  + Get(EXCV[i] , j);
ex = ex + ( Get( EX , j  ) ) ;
ec = ec + ( Get( EC , j  ) ) ;

contv++;
}


cont++;
R = Double_v (0.0);


}else{
// parte scalar:
Double_v buf1 , buf2;
G4double rho, xrho, gam, z, es, g, f, nu, cubr;
G4double sinhmenos, polu, exc , vxc;
G4double rs, buf_dep, t, At , v, u, uu, chu,  m, aux, xx, dvv;
G4double dudec, dvdx, dvdxdx, dodx;
G4double dJdt, dtdgdp; 
G4double H, depsilonBeck;
G4double dsdx, dsigma, dNdt, dlnx,  dtdp, dhdgdt;
G4double dhdt, Beck, dBeckx;
G4double decldapolu, dLdt, dudt;
G4double Dxdgdp, Dcdgdp, depsiloncLDA, Dcdgdg, Dxdgdg;  
G4double dHeclda, dududec, depsilonxlda;
G4double devasinh, epsilonxLDA, epsilonxGGA , epsiloncLDA, epsiloncGGA ;
G4double despilonxggar, ctz, che,  dexdg , dggadg; 
G4double depsggarho, dpoldx, xpoldg;

for (int j = 0; j < kVecLenD ; ++j) {

if( Get( P[i] , j ) < tol) { 

Set(buf1, j,  0.0);
Set(buf2, j , 0.0);

}else{

rho = 2.0*Get( P[i] , j );
g   = 2.0*Get( GRADV[i] , j );

cubr = pow(rho, 1.0/3.0) ;

// con constante Slter Slat
//epsilonxLDA = (-1.0)*Slat*cubr;

// con constante paper becke:
epsilonxLDA = (-1.0)*StLDA*cubr;

xrho = ( g ) / ( pow(rho, 4.0/3.0)  );
sinhmenos = asinh(xrho);
v = 1 + 6*beta*xrho*sinhmenos;
// constante slater absorbida por Beta:
Beck = (beta*xrho*xrho)/(v);
xx = xrho*xrho + 1.0 ;
devasinh = (1.0)/(sqrt( xx ));

epsilonxGGA = epsilonxLDA + cubr*( (-1.0)*Beck  );

rs = pow(  3.0/(4.0*pi*rho), 1.0/3.0);
gam = sqrt(rs);
f = 2.0*gam + B ;
nu = rs + B*gam + C ;

epsiloncLDA = A*(log(rs) - log(nu) + (Bq)*atan(Q/f) -
b1*(2.0*log(gam - x0) - log(nu) + 2.0*m1*atan(Q/f)  )  ); 

At  = Beta_Ec_div_gamma/( exp( (-1.0)*(epsiloncLDA/( gamma ) ) ) - 1.0 );
m = raizpi /( 4.0*rho*(pow( 3.0*pi*pi*rho , 1.0/6.0  )  )  );
t = m*g;
u = At*t*t;
uu = At*At*t*t*t*t ;
chu = 1.0 + u + uu;
polu = (1.0 + u)/chu;
es = 1.0 + Beta_Ec_div_gamma*t*t*polu;
H  = gamma*(log( es ) );

epsiloncGGA = epsiloncLDA + H;

ec = ec + rho*epsiloncGGA*Get( W[i] , j );    
ex = ex + rho*epsilonxGGA*Get( W[i] , j );    
exc = rho*( epsiloncGGA  +  epsilonxGGA   );


exc = exc*Get(  W[i] , j ) ;
Set(buf2, j , exc );
Acumv = Acumv + exc ;
////////////////////////////////////////////////////////////s

// con constante Slater Slat:
//depsilonxlda = (-1.0*Slat)/( 3.0*(pow(rho, 2.0/3.0 ))  ) ;

//con cnosnte paper becke:
depsilonxlda = (-1.0*StLDA)/( 3.0*(pow(rho, 2.0/3.0 ))  ) ;


dvdx = 6.0*beta*( xrho*devasinh + sinhmenos );
dBeckx = beta*( 2.0*xrho*v - xrho*xrho*dvdx )/(v*v);
depsilonBeck = ( (4.0*xrho*dBeckx)/( 3.0*pow(rho, 2.0/3.0 ) )  ) - Beck/( 3.0*pow(rho, 2.0/3.0 ) ) ;          
despilonxggar = depsilonxlda + depsilonBeck ;
dlnx = f/nu;
buf_dep =  A*( (2.0/gam)- dlnx -( (4.0*B)/( f*f + QQ) ) -
b1*( (2.0/(gam-x0)) - dlnx - (4.0*mm)/( f*f + QQ )  ) );

depsiloncLDA = (-1.0)*( gam/(6.0*rho) )*buf_dep;

z = (-1.0)*(exp( -epsiloncLDA/gamma ) )/gamma;
dudec = ( (-1.0)*t*t*Beta_Ec_div_gamma*z) /( (exp( -epsiloncLDA/gamma ) - 1.0)*( exp( -epsiloncLDA/gamma ) - 1.0 )  )  ;
dududec = (-2.0*Beta_Ec_div_gamma*Beta_Ec_div_gamma*t*t*t*t*z)/( pow( exp( -epsiloncLDA/gamma ) - 1.0 , 3.0 ) );
decldapolu = ( dudec*chu - (1.0 + u)*( dudec + dududec ) )/(chu*chu );
dHeclda = ll*( t*t*decldapolu/( es ) );

dtdgdp = -1.0*q*( ( q_prima/( rho*(pow( 3.0*pi*pi*rho, 7.0/6.0 ) ) ) )  + 1.0/( rho*rho*(pow( 3.0*pi*pi*rho, 1.0/6.0 ) ) )  );
//dtdp = (-1.0*q)*( q_prima/( rho*pow(3.0*pi*pi*rho, 7.0/ 6.0)  )  +  1.0/(rho*rho*( pow(3.0*pi*pi*rho, 1.0/6.0 ) )  )  );
dtdp = dtdgdp*g;
dudt = ( 2.0*At*t*chu - (1.0 + u)*(2.0*At*t*(1.0 + At) ) )/( chu*chu );
dhdt = ll*(2.0*t*polu + t*t*dudt )/( es );

depsggarho = depsiloncLDA + (dhdt*dtdp + dHeclda*depsiloncLDA);

// dex^GGA / dg =
dexdg = (-1.0/rho)*dBeckx;
// dec^GGA / dg =
dggadg = dexdg + m*dhdt;

dvv = 2.0*v*dvdx;
dvdxdx = 6.0*beta*( 1.0/( sqrt( pow( xx, 3.0 ) ) ) + devasinh  ) ;

dodx = xrho*xrho*dvdxdx + 2.0*xrho*dvdx;
dsdx = ( dodx*v*v - xrho*xrho*dvdx*dvv )/(v*v*v*v);
dsigma = 2.0*( xrho*dvdx + v )/(v*v);
dpoldx = beta*( dsigma - dsdx);

xpoldg =  dpoldx*xrho/( pow(rho, 4.0/3.0) )  +  dBeckx/( pow(rho, 4.0/3.0) ) ;
Dxdgdp = (1.0/( 3.0*( pow(rho, 2.0/3.0 ) ) ) )*(  4.0*xpoldg - (dBeckx/pow(rho, 4.0/3.0 ))  );

dJdt = 4.0*At*t + 12.0*(uu/t) + 12.0*(uu*u/t) + 8.0*(uu*uu/t) ;
aux = 1.0 + 2.0*u + 3.0*uu + 2.0*uu*u + uu*uu ;

dLdt = 2.0*At*At*( ( 5.0*u*t*t -1.0 -3.0*u )*aux - ( u*t*t*t -t - u*t )*dJdt )/(chu*chu*chu*chu);  

che = 2.0*( polu + t*dudt ) + 2.0*t*dudt + t*t*dLdt ;
ctz = Beta_Ec_div_gamma*( 2.0*t*polu + t*t*dudt );
dNdt = (che*es - ( 2.0*t*polu + t*t*dudt )*ctz  )/(es*es) ;
dhdgdt = dNdt*m*ll;
Dcdgdp = dhdgdt*dtdp + dhdt*dtdgdp ;

Dxdgdg = (-1.0/( rho*pow(rho, 4.0/3.0 ) ) )*dpoldx ;
Dcdgdg = dhdgdt*m;

vxc = ( epsiloncGGA  +  epsilonxGGA ) + rho*(  despilonxggar +  depsggarho  ) + 
(rho/( 4.0*(Get( GRADV[i] , j ))*(Get( GRADV[i] , j ) ) ) )*dggadg*4.0*( Get( GVG[i] ,  j ) ) - 
(rho/g)*(Dxdgdg + Dcdgdg)*4.0*( Get( GVG[i] ,  j ) ) - 
(rho/g)*dggadg*( 2.0*( Get( LAP[i] , j )  ) ) - g*dggadg  - rho*g*( Dxdgdp + Dcdgdp ) ;


vxc = vxc*Get(  W[i] , j ) ;
Set(buf1, j,  vxc  );

rho = 0.;
exc = 0.;
vxc = 0.;

} 

contv++;
}


PVXC[i] = buf1;
EXCV[i] = buf2;
cont++;

}

}

}// end for the points


return Acumv;
}


LDA::LDA()
{
}


G4double LDA::Exc_LDA_Vwn5V(  vector<Double_v> P,
vector<Double_v>& PVXC, vector<Double_v>& EXCV, vector<Double_v> W ){

G4double AcumV = 0.;
G4double EX = 0.;
G4double EC = 0.;

Double_v R ,K, G, AU, TT, EPP, AT ,LG, DLN, DEPP, DMP, LX;
Double_v FF, NU, NN, EN , V;
G4int sz = P.size();

for (int i = 0; i <  sz ; i++) {

MaskD_v mtol ( P[i] < tol );


if( MaskFull( mtol ) ) {

PVXC[i] = Double_v (0.0);
EXCV[i] = Double_v (0.0);

}else{


if( MaskEmpty( mtol ) ) {

R = 2.0*P[i];
K = vecCore::math::Cbrt(  3.0/(4.0*pi*R)  );
G = vecCore::math::Sqrt(K);
FF = 2.0*G + B;
NU = K + B*G + C;

AT = vecCore::math::ATan( Q/FF );
LG = vecCore::math::Log( NU ) ;

TT = A*( (vecCore::math::Log(K)) - LG + Bq*AT - b1*( 2.0*vecCore::math::Log(G - x0) -
LG + 2.0*m1*AT  ) );
EPP = TT*R;

DLN = FF/NU;


DEPP = A*( (2.0/G) - DLN - (4.0*B/( FF*FF + QQ )) -
b1*( (2.0/(G - x0) ) - DLN - (4.0*mm/( FF*FF + QQ) )  ) );
DMP = TT - (G/6.0)*DEPP;


LX = vecCore::math::Cbrt( P[i] );
NN  = Const_Slater*R*LX ;
EN = NN + EPP;
EXCV[i] = EN*W[i];

// Vxc
V = (4.0/3.0)*Const_Slater*LX + DMP ;
PVXC[i] = W[i]*V ;


for (int j = 0; j < kVecLenD ; ++j) {

AcumV = AcumV  + Get(EXCV[i] , j);
EX = EX + ( Get( NN , j  ) )*( Get( W[i] , j ) )  ;
EC = EC + ( Get( EPP , j  ) )*( Get( W[i] , j ) )  ;
//contv++;
}

//cont++;

}else{

// scalar:
//////////////////////////////////////////////////////////////////
Double_v buf1;
Double_v buf2;
G4double f, nu, buf_dep, exc, vxc, rho;
G4double gam, rs, Epp, dEpp, dlnx, tt;

for (int j = 0; j < kVecLenD ; ++j) {

if( Get( P[i] , j ) < tol) { 

Set(buf1, j,  0.0);
Set(buf2, j , 0.0);


}else{

rho = 2.0*Get( P[i] , j );
rs = pow(  3.0/(4.0*pi*rho), 1.0/3.0);
gam = sqrt(rs);
f = 2.0*gam + B ;
nu = rs + B*gam + C ;

tt = A*( log(rs) - log(nu) + (Bq)*atan(Q/f) - b1*( 2.0*log(gam-x0) - 
log(nu) + 2.0*m1*atan(Q/f) )  ); 
Epp = tt*rho;

dlnx = f/nu;
buf_dep =  A*( (2.0/gam)- dlnx -( (4.0*B)/( f*f + QQ) ) -
b1*( (2.0/(gam-x0)) - dlnx - (4.0*mm)/( f*f + QQ )  ) );
dEpp = tt - (gam/6.0)*buf_dep;

exc = Const_Slater*rho*pow( Get( P[i], j) , 1.0/3.0) + Epp;
exc = exc*Get(  W[i] , j ) ;

Set(buf2, j , exc );
AcumV = AcumV + exc ;

EC = EC + Epp*Get(  W[i] , j ) ;
EX = EX + Const_Slater*rho*( pow( Get( P[i], j) , 1.0/3.0) )*Get(  W[i] , j ) ;

// Vxc:
vxc = (4.0/3.0)*Const_Slater*( pow( Get( P[i], j) , 1.0/3.0)  ) + dEpp;
vxc = vxc*Get(  W[i] , j ) ;
Set(buf1, j,  vxc  );

}
//contv++;
}

PVXC[i] = buf1;
EXCV[i] = buf2;
//cont++;

}// if a zero 


}


}//  closing for


return AcumV;
}


Funcionales::Funcionales(){
}

vector<vector<G4double>> Funcionales::MVxcV(Molecule Mol ,  vector<Double_v>  Fvxc){

vector<vector<G4double>> MV;
G4double buf_cont=0;
G4int orbital_principal=0;


vector<Atom *>::iterator atr_main;
vector<Atom *>::iterator atr_row;
vector<Residue *>::iterator res_main;
vector<Residue *>::iterator res_row;

vector<G4double> buffer;
Double_v BU;

for(res_main = Mol.Residuos_cadena.begin();
    res_main != Mol.Residuos_cadena.end() ; ++res_main){

 for(atr_main = (*res_main)->Lista_atoms.begin();
    atr_main != (*res_main)->Lista_atoms.end(); ++ atr_main){

   for((*atr_main)->shell1 = (*atr_main)->orbitals.begin();
       (*atr_main)->shell1 != (*atr_main)->orbitals.end(); ++ (*atr_main)->shell1){
     
 for(res_row = Mol.Residuos_cadena.begin();
     res_row != Mol.Residuos_cadena.end() ; ++res_row){

  for(atr_row = (*res_row)->Lista_atoms.begin();
      atr_row != (*res_row)->Lista_atoms.end(); ++ atr_row){
 
   for((*atr_row)->shell2 = (*atr_row)->orbitals.begin();
       (*atr_row)->shell2 != (*atr_row)->orbitals.end(); ++ (*atr_row)->shell2){
  
G4int vz2 =  (*(*atr_row)->shell2).EvalV.size();
buf_cont =0;
BU = 0;


for (int i = 0; i <  vz2 ; i++) {

BU = ((*(*atr_main)->shell1).EvalV[i])*((*(*atr_row)->shell2).EvalV[i])*Fvxc[i];

for (int j = 0; j < kVecLenD ; ++j) {

buf_cont = buf_cont +  Get(BU, j);

}

}

buffer.push_back(buf_cont);
buf_cont = 0;
//orbital_sec++;
}
}
}

MV.push_back(buffer);
buffer.clear();
orbital_principal++;
}
}
}

return MV;
}



void Funcionales::Elec_dens_grad_unoV ( Molecule& Mol, vector< vector<G4double>> D, 
vector<Double_v>& vec_norm_gradient,  vector<Double_v>& laplaciano,
vector<Double_v>& Gv ){

vector<Atom *>::iterator atr_main;
vector<Atom *>::iterator atr_row;
vector<Residue *>::iterator res_main;
vector<Residue *>::iterator res_row;

G4int iter_uno = 0;
G4int iter_dos =0;
G4double buf_D = 0;
G4double tol = 0.0000000001;

Double_v LAP, VECN, GVG, UV1, UV2, UV3;

Double_v C_BUX = Double_v (0.0);
Double_v C_BUY = Double_v (0.0);
Double_v C_BUZ = Double_v (0.0);

Double_v DER = Double_v (0.0);
Double_v X_IZQ = Double_v (0.0);
Double_v Y_IZQ = Double_v (0.0);
Double_v Z_IZQ = Double_v (0.0);

Double_v CDX = Double_v (0.0);
Double_v CDY = Double_v (0.0);
Double_v CDZ = Double_v (0.0);

Double_v DXDX = Double_v (0.0);
Double_v DYDY = Double_v (0.0);
Double_v DZDZ = Double_v (0.0);

Double_v TDYDX,TDYDZ;
Double_v TDXDZ, TDXDY;
Double_v TDZDX, TDZDY;

Double_v DYDX, DYDZ, DXDZ, DXDY, DZDY, DZDX;

vector<vector<G4double>>::iterator Col;
vector<G4double>::iterator Row;

for (int i = 0; i <  Mol.No_points ; i++) {

iter_uno = 0;
iter_dos = 0;

for(res_main = Mol.Residuos_cadena.begin();
    res_main != Mol.Residuos_cadena.end() ; ++res_main){

 for(atr_main = (*res_main)->Lista_atoms.begin();
    atr_main != (*res_main)->Lista_atoms.end(); ++ atr_main){

   for((*atr_main)->shell1 = (*atr_main)->orbitals.begin();
       (*atr_main)->shell1 != (*atr_main)->orbitals.end(); ++ (*atr_main)->shell1){
   
//   orbital_sec = 0;
     iter_dos = 0;
  
 for(res_row = Mol.Residuos_cadena.begin();
     res_row != Mol.Residuos_cadena.end() ; ++res_row){

  for(atr_row = (*res_row)->Lista_atoms.begin();
      atr_row != (*res_row)->Lista_atoms.end(); ++ atr_row){
   
    for((*atr_row)->shell2 = (*atr_row)->orbitals.begin();
       (*atr_row)->shell2 != (*atr_row)->orbitals.end(); ++ (*atr_row)->shell2){

buf_D =  0.5*D[iter_uno][iter_dos];

X_IZQ = X_IZQ + (  buf_D*(  (*(*atr_row)->shell2).X_G[i]  )  ) ;
Y_IZQ = Y_IZQ + (  buf_D*(  (*(*atr_row)->shell2).Y_G[i]  )  ) ;
Z_IZQ = Z_IZQ + (  buf_D*(  (*(*atr_row)->shell2).Z_G[i]  )  ) ;
DER = DER +( buf_D*(  (*(*atr_row)->shell2).EvalV[i] )     );

CDX = CDX + ( buf_D*((*(*atr_row)->shell2).DXX[i])   );
CDY = CDY + ( buf_D*((*(*atr_row)->shell2).DYY[i])   );
CDZ = CDZ + ( buf_D*((*(*atr_row)->shell2).DZZ[i])   );

TDYDX = TDYDX + ( buf_D*((*(*atr_row)->shell2).DYX[i])  );
TDYDZ = TDYDZ + ( buf_D*((*(*atr_row)->shell2).DYZ[i])  );

TDXDZ = TDXDZ + ( buf_D*((*(*atr_row)->shell2).DXZ[i])  );
TDXDY = TDXDY + ( buf_D*((*(*atr_row)->shell2).DXY[i])  );


TDZDX = TDZDX + ( buf_D*((*(*atr_row)->shell2).DZX[i])  );
TDZDY = TDZDY + ( buf_D*((*(*atr_row)->shell2).DZY[i])  );

iter_dos++;
}
}
}

C_BUX = C_BUX + ( X_IZQ*( (*(*atr_main)->shell1).EvalV[i] ) + DER*( (*(*atr_main)->shell1).X_G[i]   ) );
C_BUY = C_BUY + ( Y_IZQ*( (*(*atr_main)->shell1).EvalV[i] ) + DER*( (*(*atr_main)->shell1).Y_G[i]   ) );
C_BUZ = C_BUZ + ( Z_IZQ*( (*(*atr_main)->shell1).EvalV[i] ) + DER*( (*(*atr_main)->shell1).Z_G[i]   ) );

DXDX = DXDX + ( 2.0*X_IZQ*( (*(*atr_main)->shell1).X_G[i] ) + ( CDX*(*(*atr_main)->shell1).EvalV[i] )  + DER*( (*(*atr_main)->shell1).DXX[i] )  );
DYDY = DYDY + ( 2.0*Y_IZQ*( (*(*atr_main)->shell1).Y_G[i] ) + ( CDY*(*(*atr_main)->shell1).EvalV[i] )  + DER*( (*(*atr_main)->shell1).DYY[i] )  );
DZDZ = DZDZ + ( 2.0*Z_IZQ*( (*(*atr_main)->shell1).Z_G[i] ) + ( CDZ*(*(*atr_main)->shell1).EvalV[i] )  + DER*( (*(*atr_main)->shell1).DZZ[i] )  );

DYDX = DYDX + ( X_IZQ*( (*(*atr_main)->shell1).Y_G[i] ) + ( TDYDX*(*(*atr_main)->shell1).EvalV[i] ) +                                              
( Y_IZQ*(*(*atr_main)->shell1).X_G[i]  ) ) + ( DER*(*(*atr_main)->shell1).DYX[i]  ) ;

DYDZ = DYDZ + ( Z_IZQ*( (*(*atr_main)->shell1).Y_G[i] ) + ( TDYDZ*(*(*atr_main)->shell1).EvalV[i] ) +                                              
( Y_IZQ*(*(*atr_main)->shell1).Z_G[i]  )  +  ( DER*(*(*atr_main)->shell1).DYZ[i]  ) ) ;

DXDZ = DXDZ + ( Z_IZQ*( (*(*atr_main)->shell1).X_G[i] ) + ( TDXDZ*(*(*atr_main)->shell1).EvalV[i] ) +                                              
( X_IZQ*(*(*atr_main)->shell1).Z_G[i]  )  +  ( DER*(*(*atr_main)->shell1).DXZ[i]  ) ) ;

DXDY = DXDY + ( Y_IZQ*( (*(*atr_main)->shell1).X_G[i] ) + ( TDXDY*(*(*atr_main)->shell1).EvalV[i] ) +                                              
( X_IZQ*(*(*atr_main)->shell1).Y_G[i]  )  +  ( DER*(*(*atr_main)->shell1).DXY[i]  ) ) ;


DZDX = DZDX + ( X_IZQ*( (*(*atr_main)->shell1).Z_G[i] ) + ( TDZDX*(*(*atr_main)->shell1).EvalV[i] ) +                                              
( Z_IZQ*(*(*atr_main)->shell1).X_G[i]  )  +  ( DER*(*(*atr_main)->shell1).DZX[i]  ) ) ;

DZDY = DZDY + ( Y_IZQ*( (*(*atr_main)->shell1).Z_G[i] ) + ( TDZDY*(*(*atr_main)->shell1).EvalV[i] ) +                                              
( Z_IZQ*(*(*atr_main)->shell1).Y_G[i]  )  +  ( DER*(*(*atr_main)->shell1).DZY[i]  ) ) ;

X_IZQ = Double_v (0.0);
Y_IZQ = Double_v (0.0);
Z_IZQ = Double_v (0.0);
DER = Double_v (0.0);

CDX = Double_v (0.0);
CDY = Double_v (0.0);
CDZ = Double_v (0.0);

TDYDX = Double_v (0.0);
TDYDZ = Double_v (0.0);
TDXDY = Double_v (0.0);
TDXDZ = Double_v (0.0);
TDZDX = Double_v (0.0);
TDZDY = Double_v (0.0);


iter_uno++;
}
}
}

MaskD_v con_bux ( C_BUX < tol );

if( MaskFull( con_bux ) ) {

C_BUX = tol;

}else{

for (int j = 0; j < kVecLenD ; ++j) {

if( Get( C_BUX , j ) < tol) { 

Set( C_BUX , j , tol );

}
}
}

MaskD_v con_buy ( C_BUY < tol );

if( MaskFull( con_buy ) ) {

C_BUY = tol;

}else{

for (int j = 0; j < kVecLenD ; ++j) {

if( Get( C_BUY , j ) < tol) { 

Set( C_BUY , j , tol );

}
}
}

MaskD_v con_buz ( C_BUZ < tol );

if( MaskFull( con_buz ) ) {

C_BUZ = tol;

}else{

for (int j = 0; j < kVecLenD ; ++j) {

if( Get( C_BUZ , j ) < tol) { 

Set( C_BUZ , j , tol );

}
}
}

VECN =  vecCore::math::Sqrt( C_BUX*C_BUX + C_BUY*C_BUY + C_BUZ*C_BUZ ) ;
LAP = DXDX + DYDY + DZDZ;

UV1 = C_BUX*DXDX + C_BUY*DXDY + C_BUZ*DXDZ;
UV2 = C_BUX*DYDX + C_BUY*DYDY + C_BUZ*DYDZ;
UV3 = C_BUX*DZDX + C_BUY*DZDY + C_BUZ*DZDZ;

MaskD_v norma ( VECN < tol );

if( MaskFull( norma ) ) {

VECN = tol;

}else{

for (int j = 0; j < kVecLenD ; ++j) {

if( Get( VECN , j ) < tol) { 

Set( VECN , j , tol );

}
}
}

GVG = (UV1/VECN)*C_BUX + (UV2/VECN)*C_BUY + (UV3/VECN)*C_BUZ;

vec_norm_gradient.push_back( VECN );
laplaciano.push_back( LAP );
Gv.push_back( GVG );

C_BUX = Double_v (0.0);
C_BUY = Double_v (0.0);
C_BUZ = Double_v (0.0);

DXDX = Double_v (0.0);
DYDY = Double_v (0.0);
DZDZ = Double_v (0.0);

DYDX = Double_v (0.0);
DYDZ = Double_v (0.0);
DXDZ = Double_v (0.0);
DXDY = Double_v (0.0);
DZDX = Double_v (0.0);
DZDY = Double_v (0.0);

UV1 = Double_v (0.0);
UV2 = Double_v (0.0);
UV3 = Double_v (0.0);

VECN = Double_v (0.0);
LAP = Double_v (0.0);
GVG = Double_v (0.0);

}


return;    
}


vector<Double_v> Funcionales::V_densidad ( Molecule& Mol,  vector<vector<G4double>>  D  ){

vector<G4double> P;
vector<Double_v> rhoV;

vector<Atom *>::iterator atr_main;
vector<Atom *>::iterator atr_row;
vector<Residue *>::iterator res_main;
vector<Residue *>::iterator res_row;

G4int iter_uno = 0;
G4int iter_dos =0;

Double_v bufv;
Double_v U;
Double_v CH;
Double_v CHU;

vector<vector<G4double>>::iterator Col;
vector<G4double>::iterator Row;

for (int i = 0; i <  Mol.No_points ; i++) {

iter_uno = 0;
iter_dos = 0;

for(res_main = Mol.Residuos_cadena.begin();
    res_main != Mol.Residuos_cadena.end() ; ++res_main){

 for(atr_main = (*res_main)->Lista_atoms.begin();
    atr_main != (*res_main)->Lista_atoms.end(); ++ atr_main){

   for((*atr_main)->shell1 = (*atr_main)->orbitals.begin();
       (*atr_main)->shell1 != (*atr_main)->orbitals.end(); ++ (*atr_main)->shell1){
   
     iter_dos = 0;

 for(res_row = Mol.Residuos_cadena.begin();
     res_row != Mol.Residuos_cadena.end() ; ++res_row){

  for(atr_row = (*res_row)->Lista_atoms.begin();
      atr_row != (*res_row)->Lista_atoms.end(); ++ atr_row){

    for((*atr_row)->shell2 = (*atr_row)->orbitals.begin();
       (*atr_row)->shell2 != (*atr_row)->orbitals.end(); ++ (*atr_row)->shell2){
  
CH = ( (*(*atr_row)->shell2).EvalV[i] )*0.5*D[iter_uno][iter_dos];

bufv = bufv + CH;
iter_dos++;

}
}
}

CHU = (bufv)*( (*(*atr_main)->shell1).EvalV[i]   );
U = U + CHU;

bufv = Double_v(0.0);
CH = Double_v(0.0);
iter_uno++;

}
}
}

rhoV.push_back(U);
U = Double_v(0.0);
CHU = Double_v(0.0);
}

return rhoV;
}


vector<vector<G4double>> Funcionales::XCpotentialV(Molecule& Mol,  
vector<vector<G4double>>  D,  vector<Double_v> W, G4int functional){


Energy_xc = 0;
vector<Double_v> PV =  V_densidad (Mol, D );
vector<Double_v> Gradiente_norm , Lap, GVG;
//G4double sum_densidad = 0;

G4int vz =  PV.size();

vector<Double_v> PotXCV (vz, 0.0);
vector<Double_v> EXCV   (vz, 0.0);

switch (functional) {

case 0:

Energy_xc = flda.Exc_LDA_Vwn5V( PV, PotXCV, EXCV, W);

break;

case 1:

Elec_dens_grad_unoV ( Mol, D, Gradiente_norm, Lap, GVG );
Energy_xc = fgga.Becke_correctionV(  PV, PotXCV, EXCV , Gradiente_norm, Lap, GVG,  W );

break;

case 2:


break;

}

 PV.clear();
 EXCV.clear();

//printf("energy xc %f  \n", Energy_xc );
return MVxcV(Mol ,  PotXCV);

}









