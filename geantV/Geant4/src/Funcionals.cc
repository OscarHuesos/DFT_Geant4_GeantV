
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


G4double GGA::Becke_correction( vector<G4double> P, vector<G4double>& Potxc, vector<G4double>& Exc,
vector<G4double> gradp, vector<G4double> lap, vector<G4double> G,  vector<G4double> W ){

G4double rho, gam, z, es, g, f, nu;
G4double sinhmenos, polu ,  cubr;
G4double rs, buf_dep, t, At , v, u, uu, chu, m, aux, xx, dvv;
G4double dudec, dvdx, dvdxdx, dodx;
G4double dJdt, dtdgdp; 
G4double H, depsilonBeck;
G4double dsdx, dsigma, dNdt;
G4double dlnx, dtdp, dhdgdt;
G4double dhdt, Beck, dBeckx;
G4double decldapolu, dLdt, dudt;
G4double Dxdgdp, Dcdgdp, depsiloncLDA, Dcdgdg, Dxdgdg;  
G4double dHeclda, dududec, depsilonxlda;
G4double devasinh,  epsilonxLDA;
G4double epsilonxGGA , epsiloncLDA, epsiloncGGA ;
G4double despilonxggar, ctz, che,  dexdg , dggadg; 
G4double depsggarho, dpoldx, xpoldg;

G4int cont = 0;
G4double acum = 0;
G4double xrho = 0;
G4double Ec_acum_peso = 0;
G4double Ex_acum_peso = 0;
//G4double Exc_acum_peso = 0;
G4int Psize =  P.size();

for (int i = 0; i <  Psize ; i++) {

if(  (P[i] < tol )  ) { 

Potxc[i] = 0.;
Exc[i] = 0.;

}else{

rho = 2.0*P[i];
g =   2.0*gradp[i];

// part of Ex^LDA
/////////////////////////////////////////////////////////////////////////////
//Exlda =  Const_Slater*rho*pow(rho, 1.0/3.0) ;

//epsilonxLDA =  Const_Slater*pow(P[i], 1.0/3.0);
cubr = pow(rho, 1.0/3.0) ;

// con constante Slter Slat
epsilonxLDA = (-1.0)*Slat*cubr;

//con cnosnte paper becke:
//epsilonxLDA = (-1.0)*StLDA*cubr;

xrho = ( g ) / ( pow(rho, 4.0/3.0)  );
sinhmenos = asinh(xrho);
v = 1 + 6*beta*xrho*sinhmenos;

Beck = (beta*xrho*xrho)/(v);

xx = xrho*xrho + 1.0 ;
devasinh = (1.0)/(sqrt( xx ));

epsilonxGGA = epsilonxLDA + cubr*( (-1.0)*Beck  );

rs = pow(  3.0/(4.0*pi*rho), 1.0/3.0);
gam = sqrt(rs);
f = 2.0*gam + B;
nu = rs + B*gam + C;
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

//parte LDA:
Ec_acum_peso = Ec_acum_peso + rho*epsiloncGGA*W[i];
Ex_acum_peso = Ex_acum_peso + rho*epsilonxGGA*W[i];

Exc[i] = rho*( epsiloncGGA  +  epsilonxGGA   );

////////////////////////////////////////////////////////////s
// con constante Slater Slat:
depsilonxlda = (-1.0*Slat)/( 3.0*(pow(rho, 2.0/3.0 ))  ) ;

//con cnosnte paper becke:
//depsilonxlda = (-1.0*StLDA)/( 3.0*(pow(rho, 2.0/3.0 ))  ) ;

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

dtdp = dtdgdp*g;
dudt = ( 2.0*At*t*chu - (1.0 + u)*(2.0*At*t*(1.0 + At) ) )/( chu*chu );
dhdt = ll*(2.0*t*polu + t*t*dudt )/( es );

// dec^GGA/ drho
depsggarho = depsiloncLDA + (dhdt*dtdp + dHeclda*depsiloncLDA);

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

// d^2ec^GGA / dpdg 

dJdt = 4.0*At*t + 12.0*(uu/t) + 12.0*(uu*u/t) + 8.0*(uu*uu/t) ;

aux = 1.0 + 2.0*u + 3.0*uu + 2.0*uu*u + uu*uu ;

dLdt = 2.0*At*At*( ( 5.0*u*t*t -1.0 -3.0*u )*aux - ( u*t*t*t -t - u*t )*dJdt )/(chu*chu*chu*chu);  

che = 2.0*( polu + t*dudt ) + 2.0*t*dudt + t*t*dLdt ;

ctz = Beta_Ec_div_gamma*( 2.0*t*polu + t*t*dudt );

dNdt = (che*es - ( 2.0*t*polu + t*t*dudt )*ctz  )/(es*es) ;
dhdgdt = dNdt*m*ll;

Dcdgdp = dhdgdt*dtdp + dhdt*dtdgdp ;

// d^2 ex^GGA / dg^2 
Dxdgdg = (-1.0/( rho*pow(rho, 4.0/3.0 ) ) )*dpoldx ;

Dcdgdg = dhdgdt*m;

//////////////////////////////////////////////////////////////////////////////////////////////
Potxc[i] = ( epsiloncGGA  +  epsilonxGGA ) + rho*(  despilonxggar +  depsggarho  ) +
(rho/( 4.0*gradp[i]*gradp[i] ) )*dggadg*4.0*G[i] - (rho/g)*(Dxdgdg + Dcdgdg)*4.0*G[i] -
(rho/g)*dggadg*( 2.0*lap[i] ) - g*dggadg  - rho*g*( Dxdgdp + Dcdgdp ) ;
///////////////////////////////////////////////////////////////////////////////////////////////

Exc[i] = W[i]*Exc[i];
acum = acum + Exc[i];

// 
Potxc[i] = W[i]*Potxc[i];

cont++;
rho = 0.;
}
}

return acum; 
}



LDA::LDA()
{
}


void LDA::Exc_lda_vwn5_pol( vector<G4double> na , vector<G4double> nb,
vector<G4double>& Pot, vector<G4double>& Ec){


G4double rho, xx ,rs, at, lg, xpc, tt, pn, epc;
G4double feps = 0;
G4double epsil = 0;
//G4double ecvwn = 0;
G4double alp = 0;
G4double efc = 0;
//G4double dlnxepc = 0;
G4double difepc = 0;
//G4int cont = 0;
G4int nsize =  na.size();


for (int i = 0; i < nsize ; i++) {

rho = na[i] + nb[i];
// nup = n

if(rho < 0.0000000001 ) {   

Pot[i] = 0.;
Ec[i] = 0.;

}else{
// 1.25992 = sqrt3(2)
epsil =  (na[i] - nb[i])/rho;
rs = pow(  3.0/(4.0*pi*rho), 1.0/3.0);
xx = sqrt(rs);
lg = log(rs);
pn = 2.0*xx + B;

xpc = rs + B*xx + C;
at  = atan(Q/pn); 

epc = A*(lg -log(xpc) + Bq*at -
b1*(2.0*log(xx - x0) - log(xpc) + (2.0*pn*at )  ) ); 


difepc =   A*( (2.0/xx)- (pn/xpc) -( (4.0*B)/(pn*pn + Q*Q) ) -
b1*( (2.0/(xx - x0)) - (pn/xpc) - 4.0*(2.0*x0 + B)/( pn*pn + Q*Q ) ) );


if( epsil > 0.0000000001  ){

feps = (pow(1 + epsil, 4.0/3.0) + pow(1 - epsil, 4.0/3.0) -2)/(did);

G4double xfc = rs + bfc*xx + cfc;
G4double atfc = atan(Qfc/(2.0*xx + bfc) );

efc = Afc*(lg - log(xfc) + consfc*atfc -
b2*(2.0*log(xx - x0fc) - log(xfc) + (2*(2.0*x0fc + bfc)/Qfc)*atfc  )  ); 

G4double xalp = rs + balp*xx + calp;
G4double atalp = atan(Qalp/(2.0*xx + balp) );

alp = Aalp*( lg - log(xalp) + consalp*atalp -
b3*( 2.0*log(xx - x0alp) - log(xalp) + ( 2*(2.0*x0alp + balp)/Qalp)*atalp )  ); 

}

tt = epc + (alp/conddif)*feps + feps*epsil*epsil*epsil*epsil*(efc - epc - (alp/conddif) );

}
}


return;
}


G4double LDA::Exc_lda_vwn5(   vector<G4double> P, vector<G4double>& Potcx,
vector<G4double>& Ecx , vector<G4double> W ){

G4double rho, gam, rs, buf_dep, f, nu, tt;
G4double Epp, dEpp, dlnx;
G4int cont = 0;
G4double acum = 0;
G4double Ec_acum_peso = 0;
G4double Ex_acum_peso = 0;
G4int Psize = P.size();


for (int i = 0; i <  Psize ; i++) {
// rho = rhoa + rhob
//if a equal nup = ndown, -> Ddensidad = 2D
if(P[i] < tol ) {   

Potcx[i] = 0.;
Ecx[i] = 0.;

}else{


rho = 2*P[i];

// original with rho = 2P:
rs = pow(  3.0/(4.0*pi*rho), 1.0/3.0);
//rs = pow(  3.0/(4.0*pi*P[i]), 1.0/3.0);
gam = sqrt(rs);
f = 2.0*gam + B ;
nu = rs + B*gam + C ;

tt = A*( log(rs) - log(nu) + (Bq)*atan(Q/f) - b1*( 2.0*log(gam-x0) - 
log(nu) + 2.0*m1*atan(Q/f) )  ); 

Epp = tt*rho;

//d(ln(u))/du= du /u;
//d(ln(f(u)))/du = df(u)/du / f(u) 
dlnx = f/nu;
buf_dep =  A*( (2.0/gam)- dlnx -( (4.0*B)/( f*f + QQ) ) -
b1*( (2.0/(gam-x0)) - dlnx - (4.0*mm)/( f*f + QQ )  ) );

dEpp = tt - (gam/6.0)*buf_dep;

Ecx[i] = Const_Slater*rho*pow(P[i], 1.0/3.0) + Epp;

Ec_acum_peso = Ec_acum_peso + Epp*W[i];
Ex_acum_peso = Ex_acum_peso + Const_Slater*rho*pow(P[i], 1.0/3.0)*W[i];

acum = acum + W[i]*Ecx[i];

Potcx[i] = (4.0/3.0)*Const_Slater*pow(P[i], 1.0/3.0) + dEpp;

Potcx[i] = W[i]*Potcx[i];
cont++;
}
}

return acum;
}


Funcionales::Funcionales()
{
}

void Funcionales::elec_dens_grad_uno ( Molecule& Mol, vector<vector<G4double>> D, 
vector<G4double>& vec_norm_gradient,  vector<G4double>& laplaciano, vector<G4double>& Gv ) {  

vector<Atom *>::iterator atr_main;
vector<Atom *>::iterator atr_row;
vector<Residue *>::iterator res_main;
vector<Residue *>::iterator res_row;

G4int cont = 0;
//G4double cont_buf = 0;
G4int iter_uno = 0;
G4int iter_dos =0;

G4double N = 0;
G4double L = 0;
G4double buf_D = 0;
G4double cont_bux = 0;
G4double cont_buy = 0;
G4double cont_buz = 0;

G4double buf_der = 0;
G4double bufx_izq = 0;
G4double bufy_izq = 0;
G4double bufz_izq = 0;

// var double:
G4double  centrodoublex;
G4double  centrodoubley;
G4double  centrodoublez;

G4double dxdx,dydy, dzdz;
 
G4double U1, U2, U3, gVg;
G4double two_dydx, two_dydz, two_dxdy;
G4double two_dxdz, two_dzdx, two_dzdy;

G4double dydx, dydz;
G4double dxdz, dxdy;
G4double dzdx, dzdy;

vector<vector<G4double>>::iterator Col;
vector<G4double>::iterator Row;


for (int i = 0; i <  Mol.No_points ; i++) {

//cont_buf = 0;
//orbital_sec = 0;
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
  

 buf_D =  0.5*D[iter_uno][iter_dos];
 

bufx_izq = bufx_izq + (  buf_D*( (*(*atr_row)->shell2).x_grad[i] ) );

bufy_izq = bufy_izq + ( buf_D*( (*(*atr_row)->shell2).y_grad[i] ) );

bufz_izq = bufz_izq + ( buf_D*( (*(*atr_row)->shell2).z_grad[i] ) );

buf_der = buf_der + ( buf_D*( (*(*atr_row)->shell2).Func[i]) );

centrodoublex = centrodoublex +  ( buf_D*((*(*atr_row)->shell2).dxx[i])   );

centrodoubley = centrodoubley +  ( buf_D*((*(*atr_row)->shell2).dyy[i])   );

centrodoublez = centrodoublez +  ( buf_D*((*(*atr_row)->shell2).dzz[i])   );

two_dydx = two_dydx + ( buf_D*((*(*atr_row)->shell2).dyx[i])  );

two_dydz = two_dydz + ( buf_D*((*(*atr_row)->shell2).dyz[i])  );

two_dxdz = two_dxdz + ( buf_D*((*(*atr_row)->shell2).dxz[i])  );

two_dxdy = two_dxdy + ( buf_D*((*(*atr_row)->shell2).dxy[i])  );

two_dzdx = two_dzdx + ( buf_D*((*(*atr_row)->shell2).dzx[i])  );
two_dzdy = two_dzdy + ( buf_D*((*(*atr_row)->shell2).dzy[i])  );

iter_dos++;
}
}
}


cont_bux = cont_bux + ( bufx_izq*(  (*(*atr_main)->shell1).Func[i]  ) + buf_der*(  (*(*atr_main)->shell1).x_grad[i]  )  );

cont_buy = cont_buy + ( bufy_izq*(  (*(*atr_main)->shell1).Func[i]  ) + buf_der*(  (*(*atr_main)->shell1).y_grad[i]  )  );

cont_buz = cont_buz + ( bufz_izq*(  (*(*atr_main)->shell1).Func[i]  ) + buf_der*(  (*(*atr_main)->shell1).z_grad[i]  )  );

dxdx = dxdx + ( 2.0*bufx_izq*( (*(*atr_main)->shell1).x_grad[i] ) +  ( centrodoublex*(*(*atr_main)->shell1).Func[i] ) + 
( buf_der*(*(*atr_main)->shell1).dxx[i] )  );

dydy = dydy + ( 2.0*bufy_izq*( (*(*atr_main)->shell1).y_grad[i] ) +  ( centrodoubley*(*(*atr_main)->shell1).Func[i] ) + 
( buf_der*(*(*atr_main)->shell1).dyy[i] )  );

dzdz = dzdz + ( 2.0*bufz_izq*( (*(*atr_main)->shell1).z_grad[i] ) +  ( centrodoublez*(*(*atr_main)->shell1).Func[i] ) + 
( buf_der*(*(*atr_main)->shell1).dzz[i] )  );

dydx = dydx + ( bufx_izq*( (*(*atr_main)->shell1).y_grad[i] )  + ( two_dydx*(*(*atr_main)->shell1).Func[i] ) +
 ( bufy_izq*(*(*atr_main)->shell1).x_grad[i] ) +  ( buf_der*(*(*atr_main)->shell1).dyx[i]  )  );

dydz = dydz + ( ( bufz_izq*(*(*atr_main)->shell1).y_grad[i]  ) + ( two_dydz*(*(*atr_main)->shell1).Func[i] ) +
( bufy_izq*(*(*atr_main)->shell1).z_grad[i] ) + (  buf_der*(*(*atr_main)->shell1).dyz[i] ) );

dxdz = dxdz + ( ( bufz_izq*(*(*atr_main)->shell1).x_grad[i] ) + ( two_dxdz*(*(*atr_main)->shell1).Func[i] ) +
( bufx_izq*(*(*atr_main)->shell1).z_grad[i] ) + ( buf_der*(*(*atr_main)->shell1).dxz[i]  ) );

dxdy = dxdy + ( ( bufy_izq*(*(*atr_main)->shell1).x_grad[i] ) + ( two_dxdy*(*(*atr_main)->shell1).Func[i] ) +
(  bufx_izq*(*(*atr_main)->shell1).y_grad[i]  ) + (  buf_der*(*(*atr_main)->shell1).dxy[i]   )  );


dzdx = dzdx + ( ( bufx_izq*(*(*atr_main)->shell1).z_grad[i] ) + ( two_dzdx*(*(*atr_main)->shell1).Func[i] ) +
( bufz_izq*(*(*atr_main)->shell1).x_grad[i] ) + ( buf_der*(*(*atr_main)->shell1).dzx[i]  ) );

dzdy = dzdy + ( ( bufy_izq*(*(*atr_main)->shell1).z_grad[i] ) + ( two_dzdy*(*(*atr_main)->shell1).Func[i]  ) +
( bufz_izq*(*(*atr_main)->shell1).y_grad[i]  ) + (  buf_der*(*(*atr_main)->shell1).dzy[i] ) );


bufx_izq = 0;
bufy_izq = 0;
bufz_izq = 0;
buf_der = 0;

centrodoublex = 0;

centrodoubley = 0; 

centrodoublez = 0;


two_dydx = 0;
two_dydz = 0;
two_dxdy = 0;
two_dxdz = 0;
two_dzdx = 0;
two_dzdy = 0;


iter_uno++;
//cont1_aux_borrar++;
}
}
}

if(cont_bux < 0.0000000001 ) {  
cont_bux = 0.0000000001;
}


if(cont_buy < 0.0000000001 ) {  
cont_buy = 0.0000000001;
}


if(cont_buz < 0.0000000001 ) {  
cont_buz = 0.0000000001;
}

N = sqrt( cont_bux*cont_bux + cont_buy*cont_buy + cont_buz*cont_buz  );
L = dxdx + dydy + dzdz;

U1 = cont_bux*dxdx + cont_buy*dxdy + cont_buz*dxdz;
U2 = cont_bux*dydx + cont_buy*dydy + cont_buz*dydz;
U3 = cont_bux*dzdx + cont_buy*dzdy + cont_buz*dzdz;

if(N < 0.0000000001 ) {  
N = 0.0000000001;
}  

gVg = (U1/N)*cont_bux + (U2/N)*cont_buy + (U3/N)*cont_buz;

vec_norm_gradient.push_back(N);
laplaciano.push_back(L);
Gv.push_back(gVg);


cont_bux = 0;
cont_buy = 0;
cont_buz = 0;

dxdx = 0;
dydy = 0;
dzdz = 0;

dydx = 0;
dydz = 0;
dxdz = 0;
dxdy = 0;
dzdx = 0;
dzdy = 0;

U1= 0;
U2= 0;
U3= 0;
N= 0;
L=0;
gVg = 0;
cont++;
//cont1_aux_borrar++;
}

return;

}


vector<G4double> Funcionales::densidad ( Molecule& Mol,  vector<vector<G4double>>  D ){

vector<G4double> P;

vector<Atom *>::iterator atr_main;
vector<Atom *>::iterator atr_row;
vector<Residue *>::iterator res_main;
vector<Residue *>::iterator res_row;
G4double cont_buf = 0;
G4int iter_uno = 0;
G4int iter_dos =0;
G4double buff = 0;


vector<vector<G4double>>::iterator Col;
vector<G4double>::iterator Row;


for (int i = 0; i <  Mol.No_points ; i++) {

cont_buf = 0;
//orbital_sec = 0;
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
  

buff = buff + ((*(*atr_row)->shell2).Func[i])*0.5*D[iter_uno][iter_dos];

iter_dos++;
}
}
}

cont_buf = cont_buf + buff*(  (*(*atr_main)->shell1).Func[i]  );
buff = 0;
iter_uno++;

}
}
}


P.push_back(cont_buf);
cont_buf = 0;
//cont1_aux_borrar++;
}

return P;
}



vector<vector<G4double>> Funcionales::MVxc(Molecule Mol ,  vector<G4double>  Fvxc){

vector<vector<G4double>> MV;
// G4double cont=0;
G4double buf_cont=0;
//G4int iter = 0;
G4int orbital_principal=0;
//G4int orbital_sec=0;

vector<Atom *>::iterator atr_main;
vector<Atom *>::iterator atr_row;
vector<Residue *>::iterator res_main;
vector<Residue *>::iterator res_row;

vector<G4double> buffer;
G4double ev1, ev2 , f;

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
  
   G4int vz2 =  (*(*atr_row)->shell2).Func.size();


buf_cont =0;
G4double buf_borar = 0;
ev1 = 0;
ev2 = 0;
f = 0;

for (int i = 0; i <  vz2 ; i++) {

ev1 = ((*(*atr_main)->shell1).Func[i]);
ev2 = ((*(*atr_row)->shell2).Func[i]);
f = Fvxc[i];

buf_borar = ev1*ev2*f;
buf_cont = buf_cont + buf_borar;

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


vector<vector<G4double>> Funcionales::XCpotential(Molecule& Mol,  
vector<vector<G4double>>  D,  vector<G4double> W, G4int functional ){


vector<G4double> P =  densidad (Mol, D );
vector<G4double> grd_uno_norm, laplaciano, gvg;

G4int ss = P.size();

vector<G4double> Potcx (ss, 0.0);
vector<G4double> Ecx   (ss, 0.0);

switch (functional) {

case 0:

//LDA flda;
//printf("funcional LDA \n");
Energy_xc = flda.Exc_lda_vwn5( P, Potcx, Ecx, W);

break;

case 1:

//GGA fgga;
//printf("funcional GGA \n");

elec_dens_grad_uno (Mol, D, grd_uno_norm, laplaciano, gvg );

Energy_xc = fgga.Becke_correction( P, Potcx, Ecx, grd_uno_norm, laplaciano, gvg,  W);

break;

case 2:

break;

}

P.clear();
Ecx.clear();   

//printf("energy xc %f \n", Energy_xc);
// potencial matrix construction
return  MVxc(Mol , Potcx);
}

