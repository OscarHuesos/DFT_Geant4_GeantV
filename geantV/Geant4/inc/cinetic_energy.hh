
#include "PDBmolecule.hh"
#include "operadores.hh"


class Electron_Cinetic{

public:
Electron_Cinetic();
~Electron_Cinetic() {};

G4double Kin(G4double alpha_1, G4double alpha_2, std::vector<G4int> fv1,
std::vector<G4int> fv2, G4double x1, G4double  y1, G4double z1,
G4double x2,G4double y2, G4double z2, bool side );

std::vector<std::vector<G4double>> Grad(Molecule& Mol);


private:

};
