
#ifndef OPERATORS_H
#define OPERATORS_H

#include <G4Exp.hh>
#include <G4Pow.hh>
#include <vector>
#include <G4PhysicalConstants.hh>
#include "globals.hh"
#include "Geant/Config.h"
#include <Geant/VectorTypes.h>
#include "Geant/math_wrappers.h"
#include "Geant/TaskData.h"
#include "Geant/Typedefs.h"
#include "Geant/PhysicalConstants.h"


using namespace vecCore::math; 
using vecCore::Get;
using vecCore::Set;
using vecCore::Mask_v;
//using vecCore::AssignMaskLane;
using vecCore::MaskEmpty;
using vecCore::MaskFull;
using geant::kVecLenD;

using Double_v = geant::Double_v;
using Int_v = geant::Double_v;


//template <typename T>
//using vector_t = std::vector<T>;


class Operators{


public:
Operators();
~Operators() {};

G4double doublefactorial(G4int n);

G4double binomial_coefficient(G4int n,  int k);

VECCORE_ATT_HOST_DEVICE
Double_v NaivePowerVector(Double_v base, G4int ex);

VECCORE_ATT_HOST_DEVICE
Double_v doublefactorialVector(G4int n);

VECCORE_ATT_HOST_DEVICE
Double_v FmVector(Double_v P, Double_v XX, Double_v YY,
Double_v ZZ, Double_v xp);  

VECCORE_ATT_HOST_DEVICE
Double_v XniVector(Int_v n, Int_v i);

VECCORE_ATT_HOST_DEVICE
Double_v Gamma_aproxV( Double_v GAM,
Double_v rx, Double_v ry, Double_v rz,
Double_v X1, Double_v Y1, Double_v Z1,
Double_v X2, Double_v Y2, Double_v Z2,
G4int fv10, G4int fv11, G4int fv12,
G4int fv20, G4int fv21, G4int fv22,
G4double alpha1, Double_v Alpha2, bool polyn, bool selector );


VECCORE_ATT_HOST_DEVICE
Double_v gauss_chebyshevVector( Double_v GAM, Double_v rx, Double_v ry, Double_v rz,
Double_v X1, Double_v Y1, Double_v Z1, Double_v X2, Double_v Y2, Double_v Z2, 
G4int fv10, G4int fv11, G4int fv12,
G4int fv20, G4int fv21, G4int fv22, 
G4double alpha1, Double_v Alpha2, bool polyn);

VECCORE_ATT_HOST_DEVICE
Double_v Incomplete_gammaVector(Double_v U, Double_v n );

VECCORE_ATT_HOST_DEVICE
Double_v Hermitian_recV( G4int a, G4int b, 
Double_v D, G4double alpha1, Double_v Alpha2, Double_v GAM, G4int t);

VECCORE_ATT_HOST_DEVICE
Double_v Rec_nuclearV(Double_v T, G4int i, G4int j,
G4int k, Double_v Rx, Double_v Ry, Double_v Rz, Double_v GAM,
G4double alpha1, Double_v Alpha2, G4int L, bool polyn );


VECCORE_ATT_HOST_DEVICE
Double_v NfuncVector(G4int A, G4int B, Double_v RNN, Double_v RR,
Double_v PA, Double_v PB, Double_v GD, Double_v XX,  G4int bander);


VECCORE_ATT_HOST_DEVICE
Double_v OXABVectortestU(G4double Alpha1, Double_v Alpha2,
G4int k1, G4int m1, G4int n1, G4int k2, G4int m2, G4int n2,
Double_v X1, Double_v X2, Double_v Y1, Double_v Y2, Double_v Z1, Double_v Z2);


VECCORE_ATT_HOST_DEVICE
Double_v fauxV(G4int k, G4int k1, G4int k2,  Double_v PA, Double_v PB);

private:

};
#endif

