# DFT_Geant4_GeantV

Implementation of the DFT-LDA and DFT-GGA in the Geant4 and GeantV frameworks. The vectorized implementations uses VecCore API and the respective dependencies. The basis set functions are STO-3G. The files needed are .pdb extension.

# Requeriments
The implementation requires the same dependencies that the GeantV project. Please, visit :  https://gitlab.cern.ch/GeantV/geant


- [GCC] >= gcc-9 g++-9

- [Clang] >= 10

- [CMake] >=  3.8.0
  
- [Eigen3] >= 3.3

# Build and Install

### Build Options
The table below shows the available CMake options for VecGeom that may be used to customize the build:

|Option|Default|Description|
|------|:-----:|-----------|
|VECGEOM_BACKEND|scalar|Vector backend API to be used|
|VECGEOM_BUILTIN_VECCORE|OFF|Build VecCore and its dependencies from source|
|VECGEOM_CUDA_VOLUME_SPECIALIZATION|OFF|Use specialized volumes for CUDA|
|VECGEOM_DISTANCE_DEBUG|OFF|Enable comparison of calculated distances againt ROOT/Geant4 behind the scenes|

## HEPMC/3:
```sh
git clone https://gitlab.cern.ch/hepmc/HepMC3.git
cd HepMC3/
git checkout -b 3.0.0 3.0.0
mkdir build
mkdir install
```

## Vc

```sh
$ mkdir build && cd build
$ cmake ..
$ cmake --build .
$ cmake --build . --target install
```
## ROOT 6.14

## Geant4



## Google test


## Google Benckmark

# Mode of use



