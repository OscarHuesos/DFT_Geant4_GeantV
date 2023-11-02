
#include "PDBmolecule.hh"
#include <G4Pow.hh>
#include <G4Exp.hh>
#include <G4PhysicalConstants.hh>
#include "TStopwatch.h"
#include "VecCore/Timer.h"
#include "Geant/Config.h"
#include "operadores.hh"
#include <Geant/VectorTypes.h>
#include "Geant/math_wrappers.h"
#include "Geant/TaskData.h"
#include "Geant/Typedefs.h"
#include "Geant/PhysicalConstants.h"

using namespace vecCore::math; 
using vecCore::Get;
using vecCore::Set;
using vecCore::Mask_v;
//using vecCore::AssignMaskLane;
using vecCore::MaskEmpty;
using vecCore::MaskFull;
using geant::kVecLenD;
using MaskD_v= geant::MaskD_v;


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

G4double Becke_correctionV( std::vector<Double_v> P, std::vector<Double_v>& PVXC , std::vector<Double_v>& EXCV , 
std::vector<Double_v> GRADV , std::vector<Double_v> LAP, std::vector<Double_v> GVG, std::vector<Double_v> W );

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

// Becke part:
G4double q_prima = 3.0*pi*pi/(6); 
G4double raizpi = sqrt(pi);
G4double q = raizpi/4.0;
G4double ll = 0.066725;
G4double beta  = 0.0042;
G4double gamma = 0.0310921;
G4double Beta_Ec_div_gamma = ll/gamma;

//para funcional:
G4double alpha = 2.0/3.0;
G4double Const_Slater =  (-1.0)*alpha*(9.0/8.0)*(pow(2.0 , 1.0/3.0))*(pow( 3.0/pi , 1.0/3.0));
G4double Slat = alpha*(9.0/8.0)*(pow( 3.0/pi , 1./3.0));
G4double StLDA = (3.0/2.0)*pow( ( 3.0/(4.0*pi) ), 1.0/3.0);


private:
};


class LDA{

public:

LDA();
~LDA() {};

G4double Exc_LDA_Vwn5V( std::vector<Double_v> P,
std::vector<Double_v>& PVXC, std::vector<Double_v>& EXCV , std::vector<Double_v> W);

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

//para Slater:
G4double alpha = 2.0/3.0;
G4double Const_Slater =  (-1.0)*alpha*(9.0/8.0)*(pow(2, 1./3.0))*(pow( 3.0/pi , 1./3.0));

private:	

};


class Funcionales{

public:

LDA flda;
GGA fgga;
Operators fOpx;

Funcionales();
~Funcionales() {};

std::vector<Double_v> V_densidad(  Molecule& Mol, std::vector<std::vector<G4double>>  D );

std::vector<std::vector<G4double>> XCpotentialV(Molecule& Mol, 
std::vector<std::vector<G4double>>  D,  std::vector<Double_v> W ,  G4int functional);
std::vector<std::vector<G4double>> MVxcV(Molecule Mol ,  std::vector<Double_v>  Fvxc);

void Elec_dens_grad_unoV ( Molecule& Mol, std::vector<std::vector<G4double>> D, 
std::vector<Double_v>& vec_norm_gradient,  std::vector<Double_v>& laplaciano, std::vector<Double_v>& Gv );


G4double Energy_xc = 0;
private:

};

