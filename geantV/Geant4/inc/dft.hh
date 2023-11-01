
#include "PDBmolecule.hh"
#include <G4Pow.hh>
#include <G4Exp.hh>
#include "corr_interchange.hh"
#include "repulsion_energy.hh"
#include "cinetic_energy.hh"
#include "operadores.hh"
#include "atraction_energy.hh"
#include "overlap.hh"
#include "Funcionals.hh"
#include "VecCore/Timer.h"
#include <Eigen/Eigenvalues> 

class DFT{
public:

Overlap_Int fOv;
Electron_Cinetic fElectron_K;
Operators fOperators;
Electron_Nuclei_atraction fAtract;
Electron_electron_repulsion fRepulsion;
Correlation_Interchange fCorrelation;
Funcionales fExchange;


std::vector<std::vector<G4double>> S;
std::vector<std::vector<G4double>> T;
std::vector<std::vector<G4double>> V;
std::vector<std::vector<std::vector<std::vector<G4double>>>> J;
std::vector<G4double> W;
std::vector<std::vector<G4double>> JD;
std::vector<std::vector<G4double>> D;
std::vector<std::vector<G4double>> Vxc;

G4double energia;
G4double Vnn;

DFT();

~DFT() {};


void nuclear_nuclear(Molecule& Mol);
void normalize(Molecule& Mol);
void over(Molecule& Mol);
void cin(Molecule& Mol);
void atrac(Molecule& Mol);
void rep(Molecule& Mol);
void correlation(Molecule& Mol);
G4double exchange(Molecule& Mol);

void SCF(Molecule& Mol,   G4double& dft_energy );
std::vector<std::vector<G4double>> Matriz_D(Eigen::MatrixXd  M,
G4int orb_ocupados, G4int orb_abiertos, G4int  orbits );
std::vector<std::vector<G4double>>  Rep_density ();

private:

};

