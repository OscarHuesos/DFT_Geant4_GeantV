

#include "PDBmolecule.hh"
#include <G4Exp.hh>
#include <G4PhysicalConstants.hh>
#include "operadores.hh"
#include "Geant/Config.h"
#include "Geant/math_wrappers.h"
#include "Geant/TaskData.h"
#include "Geant/Typedefs.h"
#include "Geant/PhysicalConstants.h"
#include <Geant/VectorTypes.h>

using namespace vecCore::math; 
using vecCore::Get;
using vecCore::Set;
using vecCore::Mask_v;
//using vecCore::Mask;
using geant::MaskD_v;
//using vecCore::AssignMaskLane;
using geant::MaskI_v;
using vecCore::MaskEmpty;
using vecCore::MaskFull;
using geant::kVecLenD;

using Double_v = geant::Double_v;
using Int_v = geant::Double_v;

//template <typename T>
//using vector_t = std::vector<T>;


class Overlap_Int{

public:

Overlap_Int();
~Overlap_Int(){};


std::vector<std::vector<G4double>> iterVector(Molecule& Mol);

Double_v OverVector(Atom atomo, Molecule& Mol );



private:

};
