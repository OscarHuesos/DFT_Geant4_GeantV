# DFT_Geant4_GeantV

Implementation of the DFT-LDA and DFT-GGA in the Geant4 and GeantV frameworks. The vectorized implementations uses VecCore API and the respective dependencies. The basis set functions are STO-3G. The files needed are .pdb extension.

C++11

## Requeriments
The implementation requires the same dependencies that the GeantV project. Please, visit :  https://gitlab.cern.ch/GeantV/geant


- [GCC] >= gcc-9 g++-9

- [CMake] >=  3.8.0
  
- [Eigen3] >= 3.3

- [Xerces] >= 3.2.3

## Build and Install

The table below shows the available CMake options for VecGeom that may be used to customize the build:

### UMESIMD
```sh
git clone https://github.com/edanor/umesimd.git
cd umesimd/
mkdir install
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/home/choscar/geantv/umesimd/install \
-DCMAKE_C_COMPILER=gcc-9 -DCMAKE_CXX_COMPILER=g++-9 \
/home/choscar/geantv/umesimd/
make
make install
```

### Vc
```sh
git clone https://github.com/VcDevel/Vc.git
cd Vc/
git checkout -b 1.3.3 1.3.3
mkdir build
mkdir install
cd build
cmake -DCMAKE_INSTALL_PREFIX=/home/choscar/geantv/Vc/install \
-DCMAKE_C_COMPILER=gcc-9 -DCMAKE_CXX_COMPILER=g++-9 CC=gcc-9  -DCMAKE_CXX_FLAGS="-std=c++11"    \
/home/choscar/geantv/Vc/
make
make install
```

### Google test
```sh
git clone https://github.com/google/googletest.git
cd googletest
mkdir build
mkdir install
cd build
cmake -DCMAKE_INSTALL_PREFIX=/home/choscar/geantv/googletest/install \
-DCMAKE_C_COMPILER=gcc-9 -DCMAKE_CXX_COMPILER=g++-9 CC=gcc-9  -DCMAKE_CXX_FLAGS="-std=c++11"  \
-DCMAKE_BUILD_TYPE=RELEASE /home/choscar/geantv/googletest/
make
sudo make install
```

### Google Benckmark
```sh
git clone https://github.com/google/benchmark.git
cd benchmark/
mkdir build
mkdir install
cd build
cmake -DCMAKE_INSTALL_PREFIX=/home/choscar/geantv/benchmark/install \
-DCMAKE_C_COMPILER=gcc-9 -DCMAKE_CXX_COMPILER=g++-9 CC=gcc-9  -DCMAKE_CXX_FLAGS="-std=c++11"  \
-DBENCHMARK_DOWNLOAD_DEPENDENCIES=ON -DCMAKE_BUILD_TYPE=Release \
/home/choscar/geantv/benchmark/
make
make install
```
### ROOT

v6.14.06

```sh
cmake -DCMAKE_INSTALL_PREFIX=/home/choscar/geantv/root/install  \
-DCMAKE_C_COMPILER=gcc-9 -DCMAKE_CXX_COMPILER=g++-9 CC=gcc-9  -DCMAKE_CXX_FLAGS="-std=c++11"    \
-DCFITSIO=ON -DFFTW=0N -DGSL=ON -DOPENGL=ON -DPythia8=ON  \
-Dgdml=ON -Ddavix=OFF -Dc++11=ON  \
-Dtmva=ON -Droofit=ON -Dccache=ON  \
-DTBB=OFF -Dvdt=OFF -Dbuiltin_tbb=OFF  -Dimt=OFF \
-Dgviz=ON -DOpenSSL=ON -Doracle=OFF -Dzlib=ON  -Dmathmore=ON   \
-Dasimage=OFF -Dbuiltin_afterimage=OFF \
/home/choscar/geantv/root/root-6.14.06/
make
make install
source thisroot.sh en install bin
```



### HEPMC/3
```sh
git clone https://gitlab.cern.ch/hepmc/HepMC3.git
cd HepMC3/
git checkout -b 3.0.0 3.0.0
mkdir build
mkdir install
cd build
cmake -DCMAKE_INSTALL_PREFIX=/home/choscar/geantv/HepMC3/install \
-DCMAKE_PREFIX_PATH=/home/choscar/geantv/root/install  \
-DROOT_DIR=/home/choscar/geantv/root/install/cmake \
-DCMAKE_C_COMPILER=gcc-9 -DCMAKE_CXX_COMPILER=g++-9 CC=gcc-9  -DCMAKE_CXX_FLAGS="-std=c++11"    \
/home/choscar/geantv/HepMC3/
make
make install
```



### VecCore
```sh
git clone https://gitlab.cern.ch/VecGeom/VecCore.git
cd VecCore/
git tag -l
mkdir build
mkdir install
cd build
cmake -DCMAKE_INSTALL_PREFIX=/home/choscar/geantv/VecCore/install \
-DVc_DIR=/home/choscar/geantv/Vc/install/lib/cmake/Vc \
-DCMAKE_C_COMPILER=gcc-9 -DCMAKE_CXX_COMPILER=g++-9 CC=gcc-9  -DCMAKE_CXX_FLAGS="-std=c++11"  \
-DCMAKE_PREFIX_PATH=/home/choscar/geantv/root/install  \
-DROOT_DIR=/home/choscar/geantv/root/install/cmake  \
-DROOT=ON -DBACKEND=Vc  /home/choscar/geantv/VecCore
make 
make install
```

### VecMath
```sh
git clone https://github.com/root-project/vecmath.git
cd vecmath/
mkdir build
mkdir install
cd build/
cmake -DCMAKE_INSTALL_PREFIX=/home/choscar/geantv/vecmath/install  \
-DVecCore_DIR=/home/choscar/geantv/VecCore/install/lib/cmake/VecCore  \
-DVc_DIR=/home/choscar/geantv/Vc/install/lib/cmake/Vc \
-DROOT_DIR=/home/choscar/geantv/root/install/cmake  \
-DCMAKE_C_COMPILER=gcc-9 -DCMAKE_CXX_COMPILER=g++-9 CC=gcc-9  -DCMAKE_CXX_FLAGS="-std=c++11"    \
/home/choscar/geantv/vecmath
make
make install
```







### Geant4
```sh

```






### VecGeom

VecGeom 00.05.01


```sh

```

### GeantV
-- USE_TBB OFF  : TBB disabled
git clone https://gitlab.cern.ch/GeantV/geant.git
download: cms2018.gdml y agregar a data



```sh

```


### Build Options 

|Option|Description=ON|API|
|------|--------------|---|
|WITH_GEANT4|Build with Geant4 examples.|GeantV|
|BUILD_REAL_PHYSICS_TESTS|Enable tests.|GeantV|
|USE_ROOT|Build with ROOT support.|GeantV|
|USE_NUMA|Enable NUMA-aware, requires hwloc.|GeantV|
|ROOT|Build with ROOT support.|VecCore|
|DBACKEND=Vc|Build with Vc Backend options.|VecCore|






|VECGEOM_BACKEND|scalar|Vector backend API to be used|
|VECGEOM_BUILTIN_VECCORE|OFF|Build VecCore and its dependencies from source|
|VECGEOM_CUDA_VOLUME_SPECIALIZATION|OFF|Use specialized volumes for CUDA|
|VECGEOM_DISTANCE_DEBUG|OFF|Enable comparison of calculated distances againt ROOT/Geant4 behind the scenes|


### Benckmarking 



