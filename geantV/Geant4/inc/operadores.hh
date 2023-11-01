
#ifndef OPERATORS_H
#define OPERATORS_H

#include <G4Exp.hh>
#include <G4Pow.hh>
#include <vector>
#include <G4PhysicalConstants.hh>
#include "globals.hh"


class Operators{

public:
Operators();
~Operators() {};

G4double doublefactorial(G4int n);
G4double binomial_coefficient(G4int n,  G4int k);
G4double Fm(G4double p, G4double X, G4double Y, G4double Z, G4double x);
G4double xni(G4int n, G4int i);
G4double Incomplete_gamma(G4double u, G4double n);

G4double Gamma_aprox( G4double gam,
G4double rx, G4double ry, G4double rz,
G4double x1, G4double y1, G4double z1,
G4double x2, G4double y2, G4double z2,
G4int fv10, G4int fv11, G4int fv12,
G4int fv20, G4int fv21, G4int fv22,
G4double alpha_1, G4double alpha_2, bool polyn, bool selector );


G4double Rec_nucl(G4double T, G4int i, G4int j,
G4int k, G4double Rx, G4double Ry, G4double Rz, G4double gam,
G4double alpha_1, G4double alpha_2, 
G4int L, bool polyn );

G4double Hermitian_rec( G4int a, G4int b, 
G4double diff, G4double alpha1, G4double alpha2, G4double gam,
G4int t);


G4double gauss_chebyshev(G4double p,
G4double rx, G4double ry, G4double rz,
G4double x1, G4double y1, G4double z1,
G4double x2, G4double y2, G4double z2,
G4int fv10, G4int fv11, G4int fv12,
G4int fv20, G4int fv21, G4int fv22,
G4double alpha_1, G4double alpha_2, bool polyn);


G4double faux(G4int k, G4int l1,
G4int l2, G4double PA, G4double PB);

G4double n_func( G4int a, G4int b, G4double RN,  G4double Rr,
G4double A, G4double B, G4double gam, G4double x, G4int bander);

G4double Oxab(G4double alpha_1, G4double alpha_2, G4int k1,
G4int m1, G4int n1, G4int k2, G4int m2, G4int n2,
G4double x1, G4double  y1, G4double z1,
G4double x2,G4double y2,G4double z2);

private:

};
#endif

