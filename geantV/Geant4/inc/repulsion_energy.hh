#include <G4Pow.hh>
#include <G4Exp.hh>
#include <G4PhysicalConstants.hh>
#include "PDBmolecule.hh"
#include "operadores.hh"

class Electron_electron_repulsion{

public:
Electron_electron_repulsion();
~Electron_electron_repulsion() {};


G4double repulsion_a0mc0( G4double alf_1, G4double alf_2, G4double alf_3,
G4double alf_4, G4int L10, G4int L11, G4int L12, G4int L20, G4int L21, G4int L22,
G4double x1, G4double y1, G4double z1,
G4double x2, G4double y2, G4double z2,
G4double x3, G4double y3, G4double z3,
G4double x4, G4double y4, G4double z4,
G4double N1, G4double N2, G4double N3, 
G4double N4 , G4int cc);


std::vector<std::vector<std::vector<std::vector<G4double>>>> pairing(Molecule& Mol );

G4double selector(G4double alf_1, G4double alf_2, G4double alf_3, G4double alf_4, 
G4int fv1a, G4int fv1b, G4int fv1c,
G4int fv2a, G4int fv2b, G4int fv2c, 
G4int fv3a, G4int fv3b, G4int fv3c, 
G4int fv4a, G4int fv4b, G4int fv4c, 
G4double x1, G4double y1, G4double z1, G4double x2,
G4double y2, G4double z2, G4double x3, G4double y3, G4double z3,
G4double x4, G4double y4, G4double z4, G4double N1, G4double N2,
G4double N3, G4double N4, G4int cont); 


private:

};



