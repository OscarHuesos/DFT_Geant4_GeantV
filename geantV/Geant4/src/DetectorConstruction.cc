
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// Delage et al. PDB4DNA: implementation of DNA geometry from the Protein Data
//                  Bank (PDB) description for Geant4-DNA Monte-Carlo
//                  simulations (submitted to Comput. Phys. Commun.)
// The Geant4-DNA web site is available at http://geant4-dna.org


#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Orb.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4PVPlacement.hh"
#include "G4SolidStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include<iterator> 
#include<vector> 


#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fpMessenger(0),fCheckOverlaps(false){

  fPdbFileName=G4String("arginine_neutro.pdb");
  fPdbFileStatus=0;
  fChosenOption=11;
  fpDefaultMaterial=0;
  fpWaterMaterial=0;

  fpMessenger = new DetectorMessenger(this);
}

DetectorConstruction::~DetectorConstruction(){
}


G4VPhysicalVolume* DetectorConstruction::Construct(){
  fChosenOption=11;
  fPdbFileStatus=0;
  G4VPhysicalVolume* worldPV;
  G4double energy = 0.;
  int typ;
  G4double Etime= 0.0;
  Mol=Charge_molecule(fPdbFileName, typ);
  worldPV=DefineVolumes(Mol,fChosenOption,typ);

printf("Entrando DFT \n");

TStopwatch timer;
static Timer<nanoseconds> TTimer;

TTimer.Start();
energy = Molecule_energy_Analysis(Mol);
Etime =TTimer.Elapsed();

printf("Total DFT energy in Geant4: %f with %f ns \n",  energy ,Etime);

return worldPV;
}

void DetectorConstruction::ConstructMaterials(){

  G4NistManager* nistManager = G4NistManager::Instance();
  G4bool fromIsotopes = false;
  nistManager->FindOrBuildMaterial("G4_WATER", fromIsotopes);
  fpWaterMaterial = G4Material::GetMaterial("G4_WATER");

  G4double a;  
  G4double z; 
  G4double density;
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                 kStateGas, 2.73*kelvin, 3.e-18*pascal);
  fpDefaultMaterial = G4Material::GetMaterial("Galactic");
}

void DetectorConstruction::CheckMaterials(){
  if(!fpDefaultMaterial)
    G4Exception("DetectorConstruction::CheckMaterials",
                "DEFAULT_MATERIAL_NOT_INIT_1",
                FatalException,
                "Default material not initialized.");

  if(!fpWaterMaterial)
    G4Exception("DetectorConstruction::CheckMaterials",
                "WATER_MATERIAL_NOT_INIT",
                FatalException,
                "Water material not initialized.");
}

G4VPhysicalVolume* DetectorConstruction::ConstructWorld(){

  G4double worldSize  = 1000*1*angstrom;

  if ( !fpDefaultMaterial ){
     G4Exception("DetectorConstruction::ConstructWorld",
                "DEFAULT_MATERIAL_NOT_INIT_2",
                FatalException,  "Default material not initialized.");
  }

  G4VSolid* worldS= new G4Box("World", worldSize/2, worldSize/2, worldSize/2);

  G4LogicalVolume*  worldLV  = new G4LogicalVolume(  worldS,   
      fpDefaultMaterial,  
      "World");  

  G4VisAttributes *MyVisAtt_ZZ = new G4VisAttributes(
      G4Colour(G4Colour::Gray()));
  MyVisAtt_ZZ ->SetVisibility (false);
  worldLV->SetVisAttributes(MyVisAtt_ZZ);

  G4VPhysicalVolume*  worldPV  = new G4PVPlacement( 0, G4ThreeVector(),  worldLV,  "World",     
      0,   
      false,           
      0,    
      true);            

  return worldPV;
}


Molecule DetectorConstruction::Charge_molecule(G4String filename, int& isM){
Molecule fpMol;
fpMol.Load(filename,isM);
fpMol.check_parity();
return fpMol;
}


G4VPhysicalVolume* DetectorConstruction::DefineVolumes(Molecule&
  fpMoleculeList, unsigned short int option, int& typ){

  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  // Define Materials
  ConstructMaterials();

  G4VPhysicalVolume* worldPV;
  G4LogicalVolume* worldLV;
  worldPV=ConstructWorld();
  worldLV=worldPV->GetLogicalVolume();
  fpBarycenterList=NULL;

printf("N atoms total: %d  \n",fpMoleculeList.No_atomos_en_prot);

if ( typ == 1){
printf("es proteina \n");
//fpBarycenterList=fPDBlib.ComputeProteinBarycenters(fpMoleculeList );
G4cout<<"This PDB file is Protein."<<G4endl;
}

if (typ ==0 ){
fPdbFileStatus=1;
printf("No protein \n");
fpMoleculeList.isbool(false);
}

if (option==1){
AtomisticView( worldLV,fpMoleculeList,1.0);
} else
    if (option==2){
      BarycenterView(worldLV,fpBarycenterList);
    }  else
      if (option==3){
        ResiduesView(worldLV,fpBarycenterList);
      }  else
        if (option==10){
          DrawBoundingVolume( worldLV,fpMoleculeList);
        }  else
          if (option==11){
            printf("esta opcion \n");
            AtomisticView( worldLV,fpMoleculeList,1.0);
            DrawBoundingVolume( worldLV,fpMoleculeList);
          }  else
            if (option==12){
              BarycenterView(worldLV,fpBarycenterList);
              DrawBoundingVolume( worldLV,fpMoleculeList);
            }  else
              if (option==13){
                ResiduesView(worldLV,fpBarycenterList);
                DrawBoundingVolume( worldLV,fpMoleculeList);
              }

  G4cout << "Number of loaded chains = " << fpMoleculeList.cadena <<G4endl;
  return worldPV;
}



G4double DetectorConstruction::Molecule_energy_Analysis(Molecule&  fpMol ){

G4double energy_min;
DFT fDFT;
fDFT.nuclear_nuclear(fpMol);
fDFT.normalize(fpMol);
fDFT.over(fpMol);
fDFT.cin(fpMol);
fDFT.atrac(fpMol);
fDFT.rep(fpMol);
fDFT.correlation(fpMol);
fDFT.SCF(fpMol, energy_min);
return energy_min;
}



PDBlib DetectorConstruction::GetPDBlib(){
return fPDBlib;
}

Barycenter *DetectorConstruction::GetBarycenterList(){
return fpBarycenterList;
}

Molecule *DetectorConstruction::GetMoleculeList(){
return fpMoleculeList;
}

void DetectorConstruction::AtomisticView(G4LogicalVolume* worldLV,
    Molecule moleculeListTemp, double atomSizeFactor){
  CheckMaterials();

  G4double sphereSize  = atomSizeFactor*1*angstrom;
  G4VSolid* atomS_H = new G4Orb("Sphere", sphereSize*1.2);
  G4VSolid* atomS_He = new G4Orb("Sphere", sphereSize*1.4);
  G4VSolid* atomS_C = new G4Orb("Sphere", sphereSize*1.7);
  G4VSolid* atomS_O = new G4Orb("Sphere", sphereSize*1.52);
  G4VSolid* atomS_N = new G4Orb("Sphere", sphereSize*1.55);

  G4VSolid* atomS_Nm = new G4Orb("Sphere", sphereSize*1.55);
  G4VSolid* atomS_S = new G4Orb("Sphere", sphereSize*1.8);
  G4VSolid* atomS_P = new G4Orb("Sphere", sphereSize*1.8);
  G4VSolid* atomS_X = new G4Orb("Sphere", sphereSize); 

  G4LogicalVolume* atomLV_H = new G4LogicalVolume(
      atomS_H, fpWaterMaterial,  "atomLV_H"); 
  G4VisAttributes * MyVisAtt_H = new G4VisAttributes(
      G4Colour(G4Colour::White()));
  MyVisAtt_H->SetForceSolid(true);
  atomLV_H->SetVisAttributes(MyVisAtt_H);


  G4LogicalVolume* atomLV_C = new G4LogicalVolume(
      atomS_C, fpWaterMaterial, "atomLV_C"); 
  G4VisAttributes * MyVisAtt_C = new G4VisAttributes(
      G4Colour(G4Colour::Gray())); 
  MyVisAtt_C->SetForceSolid(true);
  atomLV_C->SetVisAttributes(MyVisAtt_C);

  G4LogicalVolume* atomLV_O = new G4LogicalVolume(
      atomS_O, fpWaterMaterial, "atomLV_O"); // its solid, material, name
  G4VisAttributes * MyVisAtt_O = new G4VisAttributes(
      G4Colour(G4Colour::Red()));
  MyVisAtt_O->SetForceSolid(true);
  atomLV_O->SetVisAttributes(MyVisAtt_O);

  G4LogicalVolume* atomLV_N = new G4LogicalVolume(
      atomS_N, fpWaterMaterial, "atomLV_N"); // its solid, material, name
  G4VisAttributes * MyVisAtt_N = new G4VisAttributes(
      G4Colour(G4Colour(0.,0.,0.5)));//Dark blue
  MyVisAtt_N->SetForceSolid(true);
  atomLV_N->SetVisAttributes(MyVisAtt_N);


G4LogicalVolume* atomLV_Nm = new G4LogicalVolume(
      atomS_N, fpWaterMaterial, "atomLV_N+"); // its solid, material, name
  G4VisAttributes * MyVisAtt_Nm = new G4VisAttributes(
      G4Colour(G4Colour(0.,0.,0.5)));//Dark blue
  MyVisAtt_Nm->SetForceSolid(true);
  atomLV_Nm->SetVisAttributes(MyVisAtt_Nm);

  G4LogicalVolume* atomLV_S = new G4LogicalVolume(
      atomS_S, fpWaterMaterial, "atomLV_S"); // its solid, material, name
  G4VisAttributes * MyVisAtt_S = new G4VisAttributes(G4Colour(
      G4Colour::Yellow()));
  MyVisAtt_S->SetForceSolid(true);
  atomLV_S->SetVisAttributes(MyVisAtt_S);

  G4LogicalVolume* atomLV_P = new G4LogicalVolume(
      atomS_P, fpWaterMaterial, "atomLV_P"); // its solid, material, name
  G4VisAttributes * MyVisAtt_P = new G4VisAttributes(
      G4Colour(G4Colour(1.0,0.5,0.)));//Orange
  MyVisAtt_P->SetForceSolid(true);
  atomLV_P->SetVisAttributes(MyVisAtt_P);

  G4LogicalVolume* atomLV_He = new G4LogicalVolume(
  atomS_He, fpWaterMaterial, "atomLV_He"); // its solid, material, name
  G4VisAttributes * MyVisAtt_He = new G4VisAttributes(
  G4Colour(G4Colour(1.0,0.75,0.8)));//rosa para helio de momento, other elements in CKP
  MyVisAtt_He->SetForceSolid(true);
  atomLV_He->SetVisAttributes(MyVisAtt_He);

  G4LogicalVolume* atomLV_X = new G4LogicalVolume(atomS_X, fpWaterMaterial, "atomLV_X"); // its solid, material, name
  G4VisAttributes * MyVisAtt_X = new G4VisAttributes(
  G4Colour(G4Colour(1.0,0.75,0.8)));//Pink, other elements in CKP
  MyVisAtt_X->SetForceSolid(true);
  atomLV_X->SetVisAttributes(MyVisAtt_X);

  int nbAtomTot=0;  
  int nbAtomH=0, nbAtomC=0, nbAtomO=0, nbAtomN=0, nbAtomS=0, nbAtomP=0, nbAtomHe=0, nbAtomNm=0;
  int nbAtomX=0;

std::vector<Atom *>::iterator atr;
std::vector<Residue *>::iterator red;

for(red = moleculeListTemp.Residuos_cadena.begin();
red != moleculeListTemp.Residuos_cadena.end() ; ++red){

for(atr = (*red)->Lista_atoms.begin();
atr != (*red)->Lista_atoms.end(); ++ atr){

std::cout << (*atr)->fElement << std::endl;

if ((*atr)->fElement.compare("H ") == 0){
      nbAtomH++;
      new G4PVPlacement(0, G4ThreeVector((*atr)->fX*1*angstrom,
                                      (*atr)->fY*1*angstrom,
                                      (*atr)->fZ*1*angstrom),
                                      atomLV_H,
                                      "atomP",
                                      worldLV,
                                      false,
                                      0,fCheckOverlaps);

}else if ((*atr)->fElement.compare("He") == 0)  {
nbAtomHe++;
  new G4PVPlacement(0,  G4ThreeVector((*atr)->fX*1*angstrom,
                                  (*atr)->fY*1*angstrom,
                                  (*atr)->fZ*1*angstrom),
                                  atomLV_He,
                                  "atomP",
                                  worldLV,
                                  false,
                                  0,
                                  fCheckOverlaps);

}else if ((*atr)->fElement.compare("C ") == 0)  {
  nbAtomC++;
  new G4PVPlacement(0,  G4ThreeVector((*atr)->fX*1*angstrom,
                                  (*atr)->fY*1*angstrom,
                                  (*atr)->fZ*1*angstrom),
                                  atomLV_C,
                                  "atomP",
                                  worldLV,
                                  false,
                                  0,
                                  fCheckOverlaps);
 // printf("Paso por C ? \n" );
}else if ((*atr)->fElement.compare("O ") == 0){
    nbAtomO++;
    new G4PVPlacement(0, G4ThreeVector((*atr)->fX*1*angstrom,
                                    (*atr)->fY*1*angstrom,
                                    (*atr)->fZ*1*angstrom),
                                    atomLV_O,
                                    "atomP",
                                    worldLV,
                                    false,
                                    0,
                                    fCheckOverlaps);
  }  else if ((*atr)->fElement.compare("N ") == 0)  {
    nbAtomN++;
    new G4PVPlacement(0,  G4ThreeVector((*atr)->fX*1*angstrom,
                                    (*atr)->fY*1*angstrom,
                                    (*atr)->fZ*1*angstrom),
                                    atomLV_N,
                                    "atomP",
                                    worldLV,
                                    false,
                                    0,
                                    fCheckOverlaps);

  }  else if ((*atr)->fElement.compare("N+") == 0)  {
    nbAtomNm++;
    new G4PVPlacement(0,  G4ThreeVector((*atr)->fX*1*angstrom,
                                    (*atr)->fY*1*angstrom,
                                    (*atr)->fZ*1*angstrom),
                                    atomLV_Nm,
                                    "atomP",
                                    worldLV,
                                    false,
                                    0,
                                    fCheckOverlaps);

  }  else if ((*atr)->fElement.compare("S ") == 0){
    nbAtomS++;
    new G4PVPlacement(0,  G4ThreeVector((*atr)->fX*1*angstrom,
                                    (*atr)->fY*1*angstrom,
                                    (*atr)->fZ*1*angstrom),
                                    atomLV_S,
                                    "atomP",
                                    worldLV,
                                    false,
                                    0,
                                    fCheckOverlaps);
  }
  else if ((*atr)->fElement.compare("P ") == 0){
    nbAtomP++;
    new G4PVPlacement(0,  G4ThreeVector((*atr)->fX*1*angstrom,
                                    (*atr)->fY*1*angstrom,
                                    (*atr)->fZ*1*angstrom),
                                    atomLV_P,
                                    "atomP",
                                    worldLV,
                                    false,
                                    0,
                                    fCheckOverlaps);
  }  else  {
    nbAtomX++;
    new G4PVPlacement(0, G4ThreeVector((*atr)->fX*1*angstrom,
                                    (*atr)->fY*1*angstrom,
                                    (*atr)->fZ*1*angstrom),
                                    atomLV_X, "atomP",
                                    worldLV,
                                    false,0,fCheckOverlaps);
  }
nbAtomTot++;

}
}


G4cout << "**************** atomisticView(...) ****************" <<G4endl;
G4cout << "Number of Atoms      = " << nbAtomTot <<G4endl;
G4cout << "Number of Helios     = " << nbAtomHe <<G4endl;
G4cout << "Number of Hydrogens  = " << nbAtomH <<G4endl;
G4cout << "Number of Carbons    = " << nbAtomC <<G4endl;
G4cout << "Number of Oxygens    = " << nbAtomO <<G4endl;
G4cout << "Number of Nitrogens  = " << nbAtomN <<G4endl;
G4cout << "Number of Nitrogens+ = " << nbAtomNm <<G4endl;
G4cout << "Number of Sulfurs    = " << nbAtomS <<G4endl;
G4cout << "Number of Phosphorus = " << nbAtomP <<G4endl;
G4cout << "Number of undifined atoms =" << nbAtomX <<G4endl<<G4endl;
return;
}


void DetectorConstruction::BarycenterView(G4LogicalVolume* worldLV,
  Barycenter *barycenterListTemp){
  CheckMaterials();
  G4VSolid* atomS_ZZ;
  G4LogicalVolume* atomLV_ZZ;
  G4VisAttributes* MyVisAtt_ZZ;
  int k=0;

  while (barycenterListTemp){
    k++;
    atomS_ZZ = new G4Orb("Sphere", (barycenterListTemp->GetRadius())*angstrom);
    atomLV_ZZ = new G4LogicalVolume(atomS_ZZ,fpWaterMaterial,"atomLV_ZZ");
    MyVisAtt_ZZ = new G4VisAttributes(G4Colour(G4Colour::Magenta()));
    MyVisAtt_ZZ->SetForceSolid(true);
    atomLV_ZZ->SetVisAttributes(MyVisAtt_ZZ);

    new G4PVPlacement(0, G4ThreeVector(
        barycenterListTemp->fCenterX*1*angstrom,
        barycenterListTemp->fCenterY*1*angstrom,
        barycenterListTemp->fCenterZ*1*angstrom),
        atomLV_ZZ,
        "atomZZ",
        worldLV,
        false,
        0,
        fCheckOverlaps);
    barycenterListTemp=barycenterListTemp->GetNext();
  }
  return;
}

///////////////// BEGIN representation for Base Sugar Phosphate barycenters

void DetectorConstruction::ResiduesView(G4LogicalVolume* worldLV,
    Barycenter *barycenterListTemp){
  CheckMaterials();
  G4VisAttributes* MyVisAtt_ZZ;
  G4VSolid* tubS1_ZZ;
  G4LogicalVolume* tubLV1_ZZ;
  G4VSolid* tubS2_ZZ;
  G4LogicalVolume* tubLV2_ZZ;
  G4VSolid* AS_ZZ;
  G4LogicalVolume* ALV_ZZ;
  G4VSolid* BS_ZZ;
  G4LogicalVolume* BLV_ZZ;
  G4VSolid* CS_ZZ;
  G4LogicalVolume* CLV_ZZ;
  int k=0;

  while (barycenterListTemp)  {
    k++;
    AS_ZZ = new G4Orb("Sphere", 1.*angstrom);
    ALV_ZZ = new G4LogicalVolume(AS_ZZ,fpWaterMaterial, "ALV_ZZ");
    MyVisAtt_ZZ = new G4VisAttributes(G4Colour(G4Colour::Blue()));
    MyVisAtt_ZZ->SetForceSolid(true);
    ALV_ZZ->SetVisAttributes(MyVisAtt_ZZ);
    new G4PVPlacement(0,  G4ThreeVector(barycenterListTemp->fCenterBaseX*angstrom,
                                    barycenterListTemp->fCenterBaseY*angstrom,
                                    barycenterListTemp->fCenterBaseZ*angstrom),
                                    ALV_ZZ,
                                    "AZZ",
                                    worldLV,
                                    false,
                                    0,
                                    fCheckOverlaps);

    BS_ZZ = new G4Orb("Sphere", 1.*angstrom);
    BLV_ZZ = new G4LogicalVolume(BS_ZZ,fpWaterMaterial, "BLV_ZZ");
    MyVisAtt_ZZ = new G4VisAttributes(G4Colour(G4Colour::Red()));
    MyVisAtt_ZZ->SetForceSolid(true);
    BLV_ZZ->SetVisAttributes(MyVisAtt_ZZ);
    new G4PVPlacement(0,  G4ThreeVector(
                          (barycenterListTemp->fCenterPhosphateX)*angstrom,
                          (barycenterListTemp->fCenterPhosphateY)*angstrom,
                          (barycenterListTemp->fCenterPhosphateZ)*angstrom),
                          BLV_ZZ,
                          "BZZ",
                          worldLV,
                          false,
                          0,
                          fCheckOverlaps);

    CS_ZZ = new G4Orb("Sphere", 1.*angstrom);
    CLV_ZZ = new G4LogicalVolume(CS_ZZ,fpWaterMaterial, "CLV_ZZ");
    MyVisAtt_ZZ = new G4VisAttributes(G4Colour(G4Colour::Yellow()));
    MyVisAtt_ZZ->SetForceSolid(true);
    CLV_ZZ->SetVisAttributes(MyVisAtt_ZZ);
    new G4PVPlacement(0,  G4ThreeVector(
                          barycenterListTemp->fCenterSugarX*angstrom,
                          barycenterListTemp->fCenterSugarY*angstrom,
                          barycenterListTemp->fCenterSugarZ*angstrom),
                          CLV_ZZ,
                          "CZZ",
                          worldLV,
                          false,
                          0,
                          fCheckOverlaps);

    tubS1_ZZ  = new G4Tubs( "Cylinder",  0.,
                              0.5*angstrom,
                              std::sqrt (
           (barycenterListTemp->fCenterBaseX-barycenterListTemp->fCenterSugarX)
         * (barycenterListTemp->fCenterBaseX-barycenterListTemp->fCenterSugarX)
         + (barycenterListTemp->fCenterBaseY-barycenterListTemp->fCenterSugarY)
         * (barycenterListTemp->fCenterBaseY-barycenterListTemp->fCenterSugarY)
         + (barycenterListTemp->fCenterBaseZ-barycenterListTemp->fCenterSugarZ)
         * (barycenterListTemp->fCenterBaseZ-barycenterListTemp->fCenterSugarZ)
                              ) /2 *angstrom, 0.,  2.*pi );

    tubLV1_ZZ    = new G4LogicalVolume(tubS1_ZZ,fpWaterMaterial, "tubLV_ZZ");
    MyVisAtt_ZZ = new G4VisAttributes(G4Colour(G4Colour::Green()));
    MyVisAtt_ZZ->SetForceSolid(true);
    tubLV1_ZZ->SetVisAttributes(MyVisAtt_ZZ);

    G4double Ux= barycenterListTemp->fCenterBaseX-
        barycenterListTemp->fCenterSugarX;
    G4double Uy= barycenterListTemp->fCenterBaseY-
        barycenterListTemp->fCenterSugarY;
    G4double Uz= barycenterListTemp->fCenterBaseZ-
        barycenterListTemp->fCenterSugarZ;
    G4double llUll=std::sqrt(Ux*Ux+Uy*Uy+Uz*Uz);

    Ux=Ux/llUll;
    Uy=Uy/llUll;
    Uz=Uz/llUll;

    G4ThreeVector direction = G4ThreeVector(Ux,Uy,Uz);
    G4double theta_euler =  direction.theta();
    G4double phi_euler   =  direction.phi();
    G4double psi_euler   = 0;

    //Warning : clhep Euler constructor build inverse matrix !
    G4RotationMatrix rotm1Inv  = G4RotationMatrix(phi_euler+pi/2,
                                                  theta_euler,
                                                  psi_euler);
    G4RotationMatrix rotm1 = rotm1Inv.inverse();
    G4ThreeVector translm1 = G4ThreeVector(
        (barycenterListTemp->fCenterBaseX+barycenterListTemp->fCenterSugarX)/2.
        *angstrom,
        (barycenterListTemp->fCenterBaseY+barycenterListTemp->fCenterSugarY)/2.
        *angstrom,
        (barycenterListTemp->fCenterBaseZ+barycenterListTemp->fCenterSugarZ)/2.
        *angstrom);
    G4Transform3D transform1 = G4Transform3D(rotm1,translm1);
    new G4PVPlacement(transform1,  // rotation translation
                      tubLV1_ZZ,"atomZZ",worldLV,false, 0, fCheckOverlaps);

    tubS2_ZZ    = new G4Tubs( "Cylinder2", 0.,
                              0.5*angstrom,
                              std::sqrt (
      (barycenterListTemp->fCenterSugarX-barycenterListTemp->fCenterPhosphateX)
    * (barycenterListTemp->fCenterSugarX-barycenterListTemp->fCenterPhosphateX)
    + (barycenterListTemp->fCenterSugarY-barycenterListTemp->fCenterPhosphateY)
    * (barycenterListTemp->fCenterSugarY-barycenterListTemp->fCenterPhosphateY)
    + (barycenterListTemp->fCenterSugarZ-barycenterListTemp->fCenterPhosphateZ)
    * (barycenterListTemp->fCenterSugarZ-barycenterListTemp->fCenterPhosphateZ)
                              ) /2 *angstrom, 0.,  2.*pi );

    tubLV2_ZZ    = new G4LogicalVolume(tubS2_ZZ, fpWaterMaterial, "tubLV2_ZZ");
    MyVisAtt_ZZ = new G4VisAttributes(G4Colour(1,0.5,0));
    MyVisAtt_ZZ->SetForceSolid(true);
    tubLV2_ZZ->SetVisAttributes(MyVisAtt_ZZ);

    Ux= barycenterListTemp->fCenterSugarX-
        barycenterListTemp->fCenterPhosphateX;
    Uy= barycenterListTemp->fCenterSugarY-
        barycenterListTemp->fCenterPhosphateY;
    Uz= barycenterListTemp->fCenterSugarZ-
        barycenterListTemp->fCenterPhosphateZ;
    llUll=std::sqrt(Ux*Ux+Uy*Uy+Uz*Uz);

    Ux=Ux/llUll;
    Uy=Uy/llUll;
    Uz=Uz/llUll;

    direction = G4ThreeVector(Ux,Uy,Uz);
    theta_euler = direction.theta();
    phi_euler = direction.phi();
    psi_euler = 0;

    //Warning : clhep Euler constructor build inverse matrix !
    rotm1Inv  = G4RotationMatrix(phi_euler+pi/2,theta_euler,psi_euler);
    rotm1 = rotm1Inv.inverse();
    translm1 = G4ThreeVector(
        (barycenterListTemp->fCenterSugarX+
            barycenterListTemp->fCenterPhosphateX)/2.*angstrom,
            (barycenterListTemp->fCenterSugarY+
                barycenterListTemp->fCenterPhosphateY)/2.*angstrom,
                (barycenterListTemp->fCenterSugarZ+
                    barycenterListTemp->fCenterPhosphateZ)/2.*angstrom);
    transform1 = G4Transform3D(rotm1,translm1);
    new G4PVPlacement(transform1, // rotation translation
                      tubLV2_ZZ,
                      "atomZZ",
                      worldLV,
                      false,
                      0,
                      fCheckOverlaps);

    barycenterListTemp=barycenterListTemp->GetNext();
  }
}


void DetectorConstruction::DrawBoundingVolume(G4LogicalVolume* worldLV,
    Molecule moleculeListTemp){
  CheckMaterials();
  double dX,dY,dZ;
  double tX,tY,tZ;
  printf("drawing vol \n" );
  fPDBlib.ComputeBoundingVolumeParams(moleculeListTemp,dX, dY, dZ,tX, tY, tZ);
  //[Geometry] Build a box :
  G4VSolid* boundingS = new G4Box("Bounding", dX*1*angstrom, dY*1*angstrom, dZ*1*angstrom);
  G4LogicalVolume* boundingLV = new G4LogicalVolume(boundingS, fpWaterMaterial, "BoundingLV");
  G4RotationMatrix *pRot = new G4RotationMatrix();

  G4VisAttributes *MyVisAtt_ZZ = new G4VisAttributes(
      G4Colour(G4Colour::Gray()));
  boundingLV->SetVisAttributes(MyVisAtt_ZZ);

  new G4PVPlacement(pRot,G4ThreeVector(
                        tX*1*angstrom,
                        tY*1*angstrom,
                        tZ*1*angstrom),  
                        boundingLV,
                        "boundingPV",
                        worldLV,false,0,fCheckOverlaps);
return;
}

void DetectorConstruction::LoadPDBfile(G4String fileName){
  G4cout<<"Load PDB file se carga en otro metodo :"  << G4endl<<G4endl;
  fPdbFileName=fileName;
//  int m;
#ifdef G4MULTITHREADED

  G4MTRunManager::GetRunManager()->ReinitializeGeometry();
#else

#endif
}

void DetectorConstruction::BuildBoundingVolume(){
int m;
  if (fPdbFileStatus>0) //a PDB file has been loaded
  {
    G4cout<<"Build only world volume and bounding volume"
        " for computation."<<G4endl<<G4endl;
#ifdef G4MULTITHREADED
    G4MTRunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(Mol,10,m)  );
    G4MTRunManager::GetRunManager()->ReinitializeGeometry();
#else
    G4RunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(Mol,10,m)  );
#endif
  }  else G4cout<<"PDB file not found!"<<G4endl<<G4endl;
}

void DetectorConstruction::DrawAtoms_(){
//int m;
  if (fPdbFileStatus>0)   {
#ifdef G4MULTITHREADED
    G4MTRunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(Mol,10,m)  );
    G4MTRunManager::GetRunManager()->ReinitializeGeometry();
#else
int m;
    G4RunManager::GetRunManager()->DefineWorldVolume(
    DefineVolumes(Mol ,1, m)  );
#endif
  }  else G4cout<<"PDB file not found!"<<G4endl<<G4endl;
}

void DetectorConstruction::DrawNucleotides_(){
  if (fPdbFileStatus>0) //a PDB file has been loaded
  {
    int m;
#ifdef G4MULTITHREADED
    G4MTRunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(Mol,2, m));
    G4MTRunManager::GetRunManager()->ReinitializeGeometry();
#else
    G4RunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(Mol,2, m)  );
#endif
  } else G4cout<<"PDB file not found!"<<G4endl<<G4endl;
}


void DetectorConstruction::DrawResidues_(){
int m;
  if (fPdbFileStatus>0) 
  {
#ifdef G4MULTITHREADED
    G4MTRunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(Mol,3,m));
    G4MTRunManager::GetRunManager()->ReinitializeGeometry();
#else
    G4RunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(Mol,3,m));
#endif
  } else G4cout<<"PDB file not found!"<<G4endl<<G4endl;
}


void DetectorConstruction::DrawAtomsWithBounding_(){
  if (fPdbFileStatus>0) //a PDB file has been loaded
  {
#ifdef G4MULTITHREADED
int m;

    );
    G4MTRunManager::GetRunManager()->ReinitializeGeometry();
#else
// int m;

#endif
  }  else G4cout<<"PDB file not found!"<<G4endl<<G4endl;
}


void DetectorConstruction::DrawNucleotidesWithBounding_(){
int m;
  if (fPdbFileStatus>0) 
  {
#ifdef G4MULTITHREADED
    printf("hay multitread \n");

    G4MTRunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(Mol,12,  m) );
    G4MTRunManager::GetRunManager()->ReinitializeGeometry();
#else

    G4RunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(Mol ,12,  m)  );
#endif
  } else G4cout<<"PDB file not found!"<<G4endl<<G4endl;
}


void DetectorConstruction::DrawResiduesWithBounding_(){
int m;
  if (fPdbFileStatus>0) 
  {
#ifdef G4MULTITHREADED
    G4MTRunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(Mol,13, m));
    G4MTRunManager::GetRunManager()->ReinitializeGeometry();
#else
    G4RunManager::GetRunManager()->DefineWorldVolume(
        DefineVolumes(Mol,13, m)  );
#endif
  }  else G4cout<<"PDB file not found!"<<G4endl<<G4endl;
}


