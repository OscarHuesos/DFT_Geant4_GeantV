

#include "PDBmolecule.hh"
#include <G4Pow.hh>
#include <G4Exp.hh>
#include <G4PhysicalConstants.hh>


class Correlation_Interchange{
public:

Correlation_Interchange();
~Correlation_Interchange() {};


std::vector<std::vector<G4double>> Euler_Mac_Leb (G4double Bohr_rad);
std::vector<G4double> Grid(Molecule& Mol , G4int funcional );
std::vector<std::vector<G4double>> lebedev_weights ( G4int siz);
std::vector<G4double> becke_eval( Atom& atomo_puntos ,  Molecule& Mol,  G4int& cont  ,G4int id , G4int functional );

void Euler_Mac_Lebedev_SG0 ( Atom& atomo);
void Euler_Mac_Lebedev_SG1 ( Atom& atomo);


G4double evaluation_function( G4double X_At,
G4double Y_At, G4double Z_At, G4double alpha, G4double N,
std::vector<G4int> L, G4double x, G4double y, G4double z );

G4double evaluation_partial_function_x ( G4double X_At,
G4double Y_At, G4double Z_At, G4double alpha, G4double N,
std::vector<G4int> L, G4double x, G4double y, G4double z );

G4double evaluation_partial_function_y ( G4double X_At,
G4double Y_At, G4double Z_At, G4double alpha, G4double N,
std::vector<G4int> L, G4double x, G4double y, G4double z );

G4double evaluation_partial_function_z ( G4double X_At,
G4double Y_At, G4double Z_At, G4double alpha, G4double N,
std::vector<G4int> L, G4double x, G4double y, G4double z );

G4double evaluation_double_partial_function_x (  G4double X_At,
G4double Y_At, G4double Z_At, G4double alpha, G4double N,
std::vector<G4int> L, G4double x, G4double y, G4double z );

G4double evaluation_double_partial_function_y (  G4double X_At,
G4double Y_At, G4double Z_At, G4double alpha, G4double N,
std::vector<G4int> L, G4double x, G4double y, G4double z );

G4double evaluation_double_partial_function_z (  G4double X_At,
G4double Y_At, G4double Z_At, G4double alpha, G4double N,
std::vector<G4int> L, G4double x, G4double y, G4double z );

G4double evaluation_partial_function_dxdy (  G4double X_At,
G4double Y_At, G4double Z_At, G4double alpha, G4double N,
std::vector<G4int> L, G4double x, G4double y, G4double z );

G4double evaluation_partial_function_dxdz (  G4double X_At,
G4double Y_At, G4double Z_At, G4double alpha, G4double N,
std::vector<G4int> L, G4double x, G4double y, G4double z );

G4double evaluation_partial_function_dydx (  G4double X_At,
G4double Y_At, G4double Z_At, G4double alpha, G4double N,
std::vector<G4int> L, G4double x, G4double y, G4double z );

G4double evaluation_partial_function_dydz (  G4double X_At,
G4double Y_At, G4double Z_At, G4double alpha, G4double N,
std::vector<G4int> L, G4double x, G4double y, G4double z );

G4double evaluation_partial_function_dzdx (  G4double X_At,
G4double Y_At, G4double Z_At, G4double alpha, G4double N,
std::vector<G4int> L, G4double x, G4double y, G4double z );

G4double evaluation_partial_function_dzdy (  G4double X_At,
G4double Y_At, G4double Z_At, G4double alpha, G4double N,
std::vector<G4int> L, G4double x, G4double y, G4double z );

//void corr_inter();

private:

};



