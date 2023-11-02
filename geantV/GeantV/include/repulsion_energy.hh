#include <G4Pow.hh>
#include <G4Exp.hh>
#include <G4PhysicalConstants.hh>
#include "PDBmolecule.hh"
#include "operadores.hh"


using Double_v = geant::Double_v;
using Int_v = geant::Double_v;
//using geant::IndexD_v;
using MaskD_v= geant::MaskD_v;
using MaskI_v= geant::MaskI_v;
using vecCore::Mask_v;


class Electron_electron_repulsion{

public:
Electron_electron_repulsion();
~Electron_electron_repulsion() {};


Double_v repulsion_a0mc0V( G4double alf_1, G4double alf_2, G4double alf_3,
Double_v Alf4, G4int L10, G4int L11, G4int L12, G4int L20, G4int L21, G4int L22,
Double_v X1, Double_v Y1, Double_v Z1,
Double_v X2, Double_v Y2, Double_v Z2,
Double_v X3, Double_v Y3, Double_v Z3,
Double_v X4, Double_v Y4, Double_v Z4,
G4double n1, G4double n2, G4double n3, 
Double_v N4, G4int cc);

std::vector<std::vector<std::vector<std::vector<G4double>>>> pairingV(Molecule& Mol);


Double_v selectorV( G4double alf_1, G4double alf_2, G4double alf_3,
Double_v Alf4, G4int fv1a, G4int fv1b, G4int fv1c,
G4int fv2a, G4int fv2b, G4int fv2c, G4int fv3a, G4int fv3b,
G4int fv3c, G4int fv4a, G4int fv4b, G4int fv4c, 
Double_v X1, Double_v Y1, Double_v Z1, Double_v X2,
Double_v Y2, Double_v Z2, Double_v X3, Double_v Y3, Double_v Z3,
Double_v X4, Double_v Y4, Double_v Z4, G4double n1, G4double n2,
G4double n3, Double_v N4, G4int cont);

private:

};

