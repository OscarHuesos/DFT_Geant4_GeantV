# DFT_Geant4_GeantV

Implementation of the DFT-LDA and DFT-GGA in the Geant4 and GeantV frameworks. The vectorized implementations uses VecCore API and the respective dependencies. The basis set functions are STO-3G. The files needed are .pdb extension.

# Requeriments
The implementation requires the same dependencies that the GeantV project. Please, visit :  https://gitlab.cern.ch/GeantV/geant


- [GCC] >= gcc-9 g++-9

- Clang >= 10

- [CMake] >=  3.7.0


# Build and Install

```sh
$ mkdir build && cd build
$ cmake ..
$ cmake --build .
$ cmake --build . --target install
```


## Geant4



# Mode of use

### Build Options
The table below shows the available CMake options for VecGeom that may be used to customize the build:

|Option|Default|Description|
|------|:-----:|-----------|
|VECGEOM_BACKEND|scalar|Vector backend API to be used|
|VECGEOM_BUILTIN_VECCORE|OFF|Build VecCore and its dependencies from source|
|VECGEOM_CUDA_VOLUME_SPECIALIZATION|OFF|Use specialized volumes for CUDA|
|VECGEOM_DISTANCE_DEBUG|OFF|Enable comparison of calculated distances againt ROOT/Geant4 behind the scenes|



