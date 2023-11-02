

#include "PDBmolecule.hh"
#include <G4Pow.hh>
#include <G4Exp.hh>
#include <G4PhysicalConstants.hh>

using namespace vecCore::math; 
using vecCore::MaskEmpty;
using vecCore::MaskFull;
using vecCore::Set;
using geant::kVecLenD;
using MaskD_v= geant::MaskD_v;
using MaskI_v= geant::MaskI_v;


class Correlation_Interchange{
public:

Correlation_Interchange();
~Correlation_Interchange() {};

VECCORE_ATT_HOST_DEVICE
std::vector<Double_v> GridV(Molecule& Mol ,  G4int funcional );

VECCORE_ATT_HOST_DEVICE
std::vector<std::vector<Double_v>>  lebedev_weightsV ( G4int siz);


VECCORE_ATT_HOST_DEVICE
std::vector<Double_v> becke_eV ( Atom& atomo_puntos ,  Molecule& Mol,  G4int& cont  ,G4int  id , G4int functional);

void Euler_Mac_Lebedev_SG0 ( Atom& atomo);
void Euler_Mac_Lebedev_SG1 ( Atom& atomo);

VECCORE_ATT_HOST_DEVICE
Double_v evaluation_functionV( G4double X_At,
G4double Y_At, G4double Z_At, G4double Alpha, G4double N,
std::vector<G4int> L, Double_v X, Double_v Y, Double_v Z );

VECCORE_ATT_HOST_DEVICE
Double_v evaluation_partial_functionxV( G4double X_At,
G4double Y_At, G4double Z_At, G4double Alpha, G4double N,
std::vector<G4int> L, Double_v X, Double_v Y, Double_v Z ) ;

VECCORE_ATT_HOST_DEVICE
Double_v evaluation_partial_functionyV( G4double X_At,
G4double Y_At, G4double Z_At, G4double Alpha, G4double N,
std::vector<G4int> L, Double_v X, Double_v Y, Double_v Z );

VECCORE_ATT_HOST_DEVICE
Double_v evaluation_partial_functionzV( G4double X_At,
G4double Y_At, G4double Z_At, G4double Alpha, G4double N,
std::vector<G4int> L, Double_v X, Double_v Y, Double_v Z );

VECCORE_ATT_HOST_DEVICE
Double_v evaluation_double_partial_function_xV ( G4double X_At,
G4double Y_At, G4double Z_At, G4double Alpha, G4double N,
std::vector<G4int> L, Double_v X, Double_v Y, Double_v Z );

VECCORE_ATT_HOST_DEVICE
Double_v evaluation_double_partial_function_yV ( G4double X_At,
G4double Y_At, G4double Z_At, G4double Alpha, G4double N,
std::vector<G4int> L, Double_v X, Double_v Y, Double_v Z );

VECCORE_ATT_HOST_DEVICE
Double_v evaluation_double_partial_function_zV ( G4double X_At,
G4double Y_At, G4double Z_At, G4double Alpha, G4double N,
std::vector<G4int> L, Double_v X, Double_v Y, Double_v Z );

VECCORE_ATT_HOST_DEVICE
Double_v evaluation_partial_function_dydxV( G4double X_At,
G4double Y_At, G4double Z_At,  G4double Alpha, G4double N,
std::vector<G4int> L, Double_v X, Double_v Y, Double_v Z );

VECCORE_ATT_HOST_DEVICE
Double_v evaluation_partial_function_dzdxV( G4double X_At,
G4double Y_At, G4double Z_At,  G4double Alpha, G4double N,
std::vector<G4int> L, Double_v X, Double_v Y, Double_v Z );

VECCORE_ATT_HOST_DEVICE
Double_v evaluation_partial_function_dxdyV( G4double X_At,
G4double Y_At, G4double Z_At,  G4double Alpha, G4double N,
std::vector<G4int> L, Double_v X, Double_v Y, Double_v Z );

VECCORE_ATT_HOST_DEVICE
Double_v evaluation_partial_function_dxdzV( G4double X_At,
G4double Y_At, G4double Z_At, G4double Alpha, G4double N,
std::vector<G4int> L, Double_v X, Double_v Y, Double_v Z );

VECCORE_ATT_HOST_DEVICE
Double_v evaluation_partial_function_dzdyV( G4double X_At,
G4double Y_At, G4double Z_At,  G4double Alpha, G4double N,
std::vector<G4int> L, Double_v X, Double_v Y, Double_v Z );

VECCORE_ATT_HOST_DEVICE
Double_v evaluation_partial_function_dydzV( G4double X_At,
G4double Y_At, G4double Z_At, G4double Alpha, G4double N,
std::vector<G4int> L, Double_v X, Double_v Y, Double_v Z );


//void corr_inter();

private:

};

