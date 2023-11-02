
#include "PDBmolecule.hh"
#include "operadores.hh"
#include "Geant/Config.h"
#include <Geant/VectorTypes.h>
#include "Geant/math_wrappers.h"
#include "Geant/TaskData.h"
#include "Geant/Typedefs.h"
#include "Geant/PhysicalConstants.h"



class Electron_Cinetic{

public:
Electron_Cinetic();
~Electron_Cinetic() {};

VECCORE_ATT_HOST_DEVICE
Double_v KinV(G4double Alpha1, Double_v Alpha2,
G4int k1, G4int m1, G4int n1, G4int k2, G4int m2, G4int n2,
Double_v X1, Double_v X2, Double_v Y1, Double_v Y2, 
Double_v Z1, Double_v Z2);


std::vector<std::vector<G4double>> Grad(Molecule& Mol);
std::vector<std::vector<G4double>> GradV(Molecule& Mol);
private:

};
