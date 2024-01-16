# DFT_Geant4_GeantV

Implementation of the DFT-LDA and DFT-GGA in the Geant4 and GeantV frameworks. The vectorized implementations uses VecCore API and the respective dependencies. The basis set functions are STO-3G. The files needed are .pdb extension.


## Requeriments
The implementation requires the same dependencies that the GeantV project. Please, visit :  https://gitlab.cern.ch/GeantV/geant

- [GCC] <= gcc-9, g++-9

- [CMake] >=  3.8.0
  
- [Eigen3] >= 3.3

- [Xerces] >= 3.2.3

## Build and Install
The project can be installed manually. Each module can be compiled independently from  the source.

### UMESIMD
Expression template for vectorization library.

```sh
git clone https://github.com/edanor/umesimd.git
cd umesimd/
mkdir install
mkdir build/ && cd build/
cmake ../ -DCMAKE_INSTALL_PREFIX=$PREFIX \
make
make install
```

### Vc
Library to explicit vectorization in C++.

```sh
git clone https://github.com/VcDevel/Vc.git
cd Vc/
git checkout -b 1.3.3 1.3.3
mkdir install
mkdir build && cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=$PREFIX \
make
make install
```

### Google test
Testing library for the C++.

```sh
git clone https://github.com/google/googletest.git
cd googletest
mkdir install
mkdir build && cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -DCMAKE_BUILD_TYPE=RELEASE
make
sudo make install
```

### Google Benckmark
Library to benchmark codesnippets.

```sh
git clone https://github.com/google/benchmark.git
cd benchmark/
mkdir install
mkdir build && cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -DBENCHMARK_DOWNLOAD_DEPENDENCIES=ON -DCMAKE_BUILD_TYPE=Release 
make
make install
```
### ROOT
Open-source data analysis framework. Tested in the version v6.14.06.

```sh
cd root
mkdir install
mkdir build && cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -DFFTW=0N -DGSL=ON -DOPENGL=ON -DPythia8=ON -Dgdml=ON -Dc++11=ON  -Dccache=ON  -Dmathmore=ON   
make
make install
```

To source, run:
```sh
cd install/bin
source thisroot.sh
```

Extra dependencies may be required. If they are not needed, turn the flags off:

```sh
-Dtmva=ON -DTBB=OFF -Dbuiltin_tbb=OFF -Droofit=ON -Dzlib=ON -DCFITSIO=ON     \
-DOpenSSL=ON -Dimt=OFF -Doracle=OFF -Dasimage=OFF -Dbuiltin_afterimage=OFF   \
-Dvdt=OFF -Ddavix=OFF -Dgviz=ON
```

### HEPMC/3
C++ event record for HEP.

```sh
git clone https://gitlab.cern.ch/hepmc/HepMC3.git
cd HepMC3/
git checkout -b 3.0.0 3.0.0
mkdir install
mkdir build && cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -DROOT_DIR = "path to dir with .../.cmake files of ROOT"
make
make install
```

### VecCore
Abstraction API for vectorization.

```sh
git clone https://gitlab.cern.ch/VecGeom/VecCore.git
cd VecCore/
git tag -l
mkdir install
mkdir build && cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -DROOT=ON -DBACKEND=Vc \
  -DVc_DIR =  "path to dir with /Vc/.../.cmake files of Vc"   \
  -DROOT_DIR = "path to dir with .../.cmake files of ROOT"
make 
make install
```

### VecMath
Library of vectorized math utilities based on VecCore.

```sh
git clone https://github.com/root-project/vecmath.git
cd vecmath/
mkdir install
mkdir build && cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -DVecCore_DIR =  "path to dir with /VecCore/.../.cmake files of VecCore"  \
  -DVc_DIR = "path to dir with /Vc/.../.cmake files of Vc"                  \
  -DROOT_DIR = "path to dir with .../.cmake files of ROOT"
make
make install
```

### VecCoreLib
Additional libraries for VecCore.

```sh
git clone https://gitlab.cern.ch/GeantV/VecCoreLib.git
mkdir install
mkdir build && cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -DROOT=ON -DBACKEND=Vc      \
  -DVecCore_DIR = "path to dir with /VecCore/.../.cmake files of VecCore"  \
  -DVc_DIR = "path to dir with /Vc/.../.cmake files of Vc"                 \
  -DROOT_DIR = "path to dir with .../.cmake files of ROOT"
make
make install
```
### Geant4
Toolkit for HEP. Tested in the version 10.07.04

```sh
cd root
mkdir install
mkdir data &&
mkdir build && cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -DGEANT4_USE_OPENGL_X11=ON -DGEANT4_USE_GDML=ON -DGEANT4_USE_SYSTEM_ZLIB=ON  \  
  -DGEANT4_ROOT=ON -DROOT=ON -DGEANT4_INSTALL_DATA=ON -DGEANT4_USE_QT=ON        \
  -DROOT_DIR =  "path to dir with .../.cmake files of ROOT"  \
  -DGEANT4_INSTALL_DATADIR =  "path to the dir designated for the databases"  \
  -DXERCESC_ROOT_DIR = "path to dir for the Xerces installation"  
make -j4
make install
```

To source, run:
```sh
cd install/bin
. geant4.sh  /  source geant4.sh
```

### VecGeom
Library for vectorized geometry used in HEP. Tested in the version VecGeom 00.05.01

```sh
cd vecgeom
mkdir install
mkdir build && cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -DBACKEND=Vc -DROOT=ON -DVECGEOM_ROOT=ON -DVECGEOM_BACKEND=Vc   \
  -DVc_DIR = "path to dir with /Vc/.../.cmake files of Vc"                 \
  -DROOT_DIR = "path to dir with .../.cmake files of ROOT"
make 
make install
```

### GeantV with the DFT method
Vectorized vesion of Geant4 (under testing). Download  cms2018.gdml and add to the data folder manually. Then:

```sh
cd geant
mkdir install
mkdir build && cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -DVc_DIR = "path to dir with /Vc/.../.cmake files of Vc"        \
  -DROOT_DIR = "path to dir with .../.cmake files of ROOT"        \         
  -DGeant4_DIR =  "path to dir with .../.cmake files of Geant4"        \   
  -DVecCore_DIR = "path to dir with /VecCore/.../.cmake files of VecCore"  \
  -DVecMath_DIR = "path to dir for the VecMath installation"              \
  -DVecCoreLib_DIR =  "path to dir for the VecCoreLib installation"              \
  -DVecGeom_DIR =  "path to dir with .../.cmake files of VecGeom"        \   
  -DHepMC_DIR =  "path to dir with .../.cmake files of HepMC3"        \   
  -DCMAKE_PREFIX_PATH= "path to dir for the Eigen3 installation"  \
  -DWITH_GEANT4=ON -DBUILD_REAL_PHYSICS_TESTS=ON  -DUSE_ROOT=ON -DUSE_NUMA= OFF   \
make
make install
```

If TBB is selected, switch to ON:
```sh
-- USE_TBB =  OFF  
```

To compile GeantV without DFT, download directly the GeantV project from:

```sh
git clone https://gitlab.cern.ch/GeantV/geant.git
```

### Build Options 

|Option|Description=ON|API|
|------|--------------|---|
|OPENGL|Enable OpenGL support|ROOT|
|Pythia8|Build with Pythia8 tool|ROOT|  
|GSL|Enable GSL library|ROOT|
|WITH_GEANT4|Build with Geant4 examples|GeantV|
|BUILD_REAL_PHYSICS_TESTS|Enable tests|GeantV|
|USE_ROOT|Build with ROOT support|GeantV|
|USE_NUMA|Enable NUMA-aware, requires hwloc|GeantV|
|ROOT|Build with ROOT support|VecCore|
|BACKEND=Vc|Build with Vc Backend options|VecCore|
|ROOT|Build with ROOT support|VecCoreLib|
|BACKEND=Vc|Build with Vc Backend options|VecCoreLib|
|GEANT4_USE_OPENGL_X11|OpenGL support|Geant4|
|GEANT4_ROOT|Build with ROOT support|Geant4|
|GEANT4_USE_GDML|Enable GDML|Geant4|
|GEANT4_INSTALL_DATA|Enable databases from Geant4|Geant4|
|GEANT4_USE_QT|Build with QT support|Geant4|
|VECGEOM_BACKEND=Vc|Build with Vc Backend options|VecGeom|
|VECGEOM_ROOT|Build with ROOT support|VecGeom|


### Benchmarking 
Currently, the project is under benchmarking. The molecule is set on the DetectorConstruction.cc file. Either the folder of Geant4 or GeantV contains the available molecules.
To add new molecules, it should be put in the respective folder, an agregate the respective name in the CMakeLists within the Test_SCRIPTS section.
Choose 0 = LDA, 1 = GGA in DFT::correlation and DFT::exchange method from dft.cc to select the appropiate functional.

To launch the test, in "GeantV/build/dft/Geant(4/V)" folder, run:

```sh
source "path to Geant4 dir/install/bin/geant4.sh"
source "path to ROOT dir/install/bin/thisroot.sh"
./dftprotein(4/V)
```
The Benckmarking works in both scalar and vector mode. GUI using .mac files are available.
