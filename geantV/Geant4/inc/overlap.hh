
#include "PDBmolecule.hh"
#include <G4Exp.hh>
#include <G4PhysicalConstants.hh>
#include "operadores.hh"

class Overlap_Int{

public:

Overlap_Int();
~Overlap_Int(){};

std::vector<std::vector<G4double>> iter(Molecule& Mol);

private:

};
