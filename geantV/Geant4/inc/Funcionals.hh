
#include "PDBmolecule.hh"
#include <G4Pow.hh>
#include <G4Exp.hh>
#include <G4PhysicalConstants.hh>
#include "TStopwatch.h"
#include "VecCore/Timer.h"


class Hybrid_Func{

public:

Hybrid_Func();
~Hybrid_Func() {};

void B3LYP(  std::vector<G4double> n , 
std::vector<G4double>& Pot, std::vector<G4double>& Ec );

private:	
};


class GGA{

public:
	
GGA();
~GGA() {};


G4double Becke_correction( std::vector<G4double> P, std::vector<G4double>& Potcx, std::vector<G4double>& Ecx , 
std::vector<G4double> gradp , std::vector<G4double> lap, std::vector<G4double> G, std::vector<G4double> W );

G4double tol = 0.0000000001;
// LDA part
G4double x0 = -0.105;
G4double A = 0.0310907;
G4double B = 3.72744;
G4double C = 12.9352;
G4double Q = 6.152;
G4double Bq = 2.0*B/Q;
G4double QQ = Q*Q;
G4double mm = 2.0*x0 + B;
G4double m1 = mm/Q;
G4double pol = 12.55484;
G4double b1 = (B*x0)/pol;

// Becke part:
G4double q_prima = 3.0*pi*pi/(6); 
G4double raizpi = sqrt(pi);
G4double q = raizpi/4.0;
G4double ll = 0.066725;

//for funcional:
G4double alpha = 2.0/3.0;
G4double Const_Slater =  (-1.0)*alpha*(9.0/8.0)*(pow(2, 1./3.0))*(pow( 3.0/pi , 1./3.0));
G4double Slat = alpha*(9.0/8.0)*(pow( 3.0/pi , 1./3.0));
G4double StLDA = (3.0/2.0)*pow( ( 3.0/(4.0*pi) ), 1.0/3.0);
//becke:
G4double beta  = 0.0042;
G4double gamma = 0.0310921;
G4double Beta_Ec_div_gamma = ll/gamma;
//G4double omega = 0.046644;

private:
};



class LDA{

public:

LDA();
~LDA() {};

void Exc_lda_vwn5_pol(  std::vector<G4double> na , std::vector<G4double> nb,
std::vector<G4double>& Pot, std::vector<G4double>& Ec );

G4double Exc_lda_vwn5(  std::vector<G4double> P, std::vector<G4double>& Potcx,
std::vector<G4double>& Ecx , std::vector<G4double> W );

G4double tol = 0.0000000001;
G4double x0 = -0.105;
G4double A = 0.0310907;
G4double B = 3.72744;
G4double C = 12.9352;
G4double Q = 6.152;
G4double Bq = 2.0*B/Q;
G4double QQ = Q*Q;
G4double mm = 2.0*x0 + B;
G4double m1 = mm/Q;
G4double pol = 12.55484;
G4double b1 = (B*x0)/pol;


// polarized !¡
G4double Afc =  0.0155453;
G4double x0fc= -0.325;
G4double bfc =  7.06042;
G4double cfc = 18.0578;
G4double Qfc = sqrt(4.0*cfc - bfc*bfc);
G4double polx0fc = 15.86879;
G4double b2 = (bfc*x0fc)/( polx0fc ); 
G4double consfc = (2.0*bfc)/Qfc;

G4double Aalp  = -1.0/( 6.0*pi );
G4double x0alp = -0.0047584;
G4double balp  =  1.13107;
G4double calp  = 13.0045;
G4double Qalp  = sqrt(4.0*calp - balp*balp);
G4double polx0alp = 12.9991;
G4double b3 = (balp*x0alp)/(polx0alp);
G4double consalp = (2.0*balp)/Qalp;
G4double did = 2.0*( 1.25992 - 1);
G4double conddif = 4.0/(9.0*( 1.25992 - 1 ) );

// LDA functional:
G4double alpha = 2.0/3.0;
G4double Const_Slater =  (-1.0)*alpha*(9.0/8.0)*(pow(2, 1./3.0))*(pow( 3.0/pi , 1./3.0));

private:	

};



class Funcionales{

public:

LDA flda;
GGA fgga;

Funcionales();
~Funcionales() {};

std::vector<G4double> densidad(  Molecule& Mol, std::vector<std::vector<G4double>>  D  );

void elec_dens_grad_uno ( Molecule& Mol, std::vector<std::vector<G4double>> D, 
std::vector<G4double>& vec_norm_gradient, std::vector<G4double>& laplaciano, std::vector<G4double>& GV );

std::vector<std::vector<G4double>> XCpotential(Molecule& Mol, 
std::vector<std::vector<G4double>>  D,  std::vector<G4double> W ,G4int functional );

std::vector<std::vector<G4double>> MVxc(Molecule Mol ,  std::vector<G4double>  Fvxc);
G4double Energy_xc = 0;
private:

};






