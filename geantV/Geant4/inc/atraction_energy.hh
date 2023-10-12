
#include "PDBmolecule.hh"
#include <G4Pow.hh>
#include <G4Exp.hh>
#include <G4PhysicalConstants.hh>
#include "operadores.hh"


class Electron_Nuclei_atraction{

public:
Electron_Nuclei_atraction();
~Electron_Nuclei_atraction() {};


std::vector<std::vector<G4double>> Nuclear(Molecule& Mol);

G4double Atract(G4double alpha_1, G4double alpha_2, std::vector<G4int> fv1,
std::vector<G4int> fv2, G4double x1, G4double  y1, G4double z1,
G4double x2,G4double y2, G4double z2 , G4double rx,
G4double ry, G4double rz );

private:

};
