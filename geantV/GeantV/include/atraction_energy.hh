
#include "PDBmolecule.hh"
#include <G4Pow.hh>
#include <G4Exp.hh>
#include <G4PhysicalConstants.hh>
#include "operadores.hh"


using Double_v = geant::Double_v;
using Int_v = geant::Double_v;


class Electron_Nuclei_atraction{

public:
Electron_Nuclei_atraction();
~Electron_Nuclei_atraction() {};


Double_v AtractV(G4double alpha_1, Double_v Alpha2, std::vector<G4int> fv1,
std::vector<G4int> fv2, Double_v X1, Double_v  Y1, Double_v Z1, Double_v X2, 
Double_v Y2, Double_v Z2 , Double_v rx, Double_v ry, Double_v rz );


std::vector<std::vector<G4double>> Nuclear(Molecule& Mol);
std::vector<std::vector<G4double>> NuclearV(Molecule& Mol);


//Double_v AtractVector();


private:

};
