#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

# Eigen:
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

project(dftprotein4)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
#----------------------------------------------------------------------------
# Find Geant4 package, batch mode only executable (i.e. no need ui and vis).
#
#find_package(Geant4 REQUIRED)

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
#include_directories(${PROJECT_SOURCE_DIR}/include)


include_directories(
  ${PROJECT_SOURCE_DIR}/inc
  ${CMAKE_SOURCE_DIR}/core/numa/inc
  ${CMAKE_SOURCE_DIR}/core/concurrency/inc
  ${CMAKE_SOURCE_DIR}/core/interfaces/inc
  ${CMAKE_SOURCE_DIR}/run/userapp/inc
  ${CMAKE_SOURCE_DIR}/magneticfield/inc
  ${CMAKE_SOURCE_DIR}/physics/generators/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/particles/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/material/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/management/inc
  ${CMAKE_SOURCE_DIR}/physics/electromagnetic/models/inc
  ${CMAKE_SOURCE_DIR}/physics/electromagnetic/processes/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/cuts/inc
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm5/GeantV/inc
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm3/GeantV/inc
  ${CMAKE_SOURCE_DIR}/examples/physics/FullLHCb/GeantV/inc
  ${CMAKE_SOURCE_DIR}/examples/physics/FullCMS/GeantV/inc
  )

#set(ExampleScripts naiver.mac README})
#set(ExampleCTestLevel 2)

set(sources
   ${CMAKE_SOURCE_DIR}/core/base/src/Error.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/TaskData.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/Track.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/TrackGeo.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/Basket.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/Handler.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/SimulationStage.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/PreStepHandler.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeomQueryHandler.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/FieldPropagationHandler.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/FieldLookup.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/LinearPropagationHandler.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/SteppingActionsHandler.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/PreStepStage.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeomQueryStage.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/PropagationStage.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/SteppingActionsStage.cxx
   ${CMAKE_SOURCE_DIR}/core/concurrency/src/sync_objects.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/Propagator.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/Event.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/EventSet.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/PhysicsInterface.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/WorkloadManager.cxx
   ${CMAKE_SOURCE_DIR}/run/userapp/src/UserApplication.cxx
   ${CMAKE_SOURCE_DIR}/run/userapp/src/UserDetectorConstruction.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/FactoryStore.cxx
   ${CMAKE_SOURCE_DIR}/run/userapp/src/MyHit.cxx
   ${CMAKE_SOURCE_DIR}/run/userapp/src/StdApplication.cxx
#  ${CMAKE_CURRENT_SOURCE_DIR}/src/TransportManager.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/LocalityManager.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/MCTruthMgr.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/RunManager.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/EventServer.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/TrackStat.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/Region.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/ScalarNavInterfaceVG.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/ScalarNavInterfaceVGM.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/VectorNavInterface.cxx
)

set(headers
   ${CMAKE_SOURCE_DIR}/core/base/inc/Geant/Config.h
   ${CMAKE_SOURCE_DIR}/core/base/inc/Geant/Error.h
   ${CMAKE_SOURCE_DIR}/core/base/inc/Geant/Fwd.h
   ${CMAKE_SOURCE_DIR}/core/concurrency/inc/Geant/dcqueue.h
   ${CMAKE_SOURCE_DIR}/core/concurrency//inc/Geant/mpmc_bounded_queue.h
   ${CMAKE_SOURCE_DIR}/core/concurrency//inc/Geant/priority_queue.h
   ${CMAKE_SOURCE_DIR}/core/concurrency//inc/Geant/condition_locker.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/TaskData.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/RngWrapper.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/Geant/PhysicsInterface.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/Basket.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/Handler.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/Basketizer.h
 #  ${CMAKE_CURRENT_SOURCE_DIR}/inc/StackLikeBuffer.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/SimulationStage.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/PreStepHandler.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/GeomQueryHandler.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/FieldPropagationHandler.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/FieldLookup.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/LinearPropagationHandler.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/SteppingActionsHandler.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/PreStepStage.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/GeomQueryStage.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/PropagationStage.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/SteppingActionsStage.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/Propagator.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/Geant/Track.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/TrackGeo.h
   ${CMAKE_SOURCE_DIR}/run/userapp/inc/Geant/UserApplication.h
   ${CMAKE_SOURCE_DIR}/run/userapp/inc/Geant/UserDetectorConstruction.h
   ${CMAKE_SOURCE_DIR}/run/userapp/inc/Geant/UserFieldConstruction.h
   ${CMAKE_SOURCE_DIR}/run/userapp/inc/Geant/MyHit.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/Geant/Factory.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/Geant/FactoryStore.h
   ${CMAKE_SOURCE_DIR}/run/userapp/inc/Geant/StdApplication.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/Geant/PrimaryGenerator.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/TrackManager.h
#   ${CMAKE_CURRENT_SOURCE_DIR}/inc/TransportManager.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/LocalityManager.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/Geant/MCTruthMgr.h
   ${CMAKE_SOURCE_DIR}/core/concurrency/inc/Geant/cuckoohash_map.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/RunManager.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/Geant/Event.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/Geant/EventSet.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/EventServer.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/GeantConfig.h
#   ${CMAKE_CURRENT_SOURCE_DIR}/inc/TrackStat.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/Geant/Region.h
   ${CMAKE_SOURCE_DIR}/run/userapp/inc/Geant/ExternalFramework.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/Geant/ScalarNavInterfaceVG.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/Geant/ScalarNavInterfaceVGM.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/Geant/VectorNavInterface.h
)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#

#include_directories(${PROJECT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIR})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/inc/*.hh)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project

#include(${Geant4_USE_FILE})
#include_directories(${PROJECT_SOURCE_DIR}/inc)


#esta es la que agregas, es la buena
#include_directories("${PROJECT_SOURCE_DIR}/core/base/inc")

#include_directories( inc
#  ${CMAKE_SOURCE_DIR}/examples/magneticfield/simplifiedCMS/inc
#  ${CMAKE_SOURCE_DIR}/base/inc
#)

#add_library(Vmagfield SHARED ${sources} )
#target_link_libraries(Vmagfield  -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${VECGEOM_LIBRARIES} ${VecCoreLib_LIBRARIES})

#set(OUTPUT bin/examples/pathtodft)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${OUTPUT})


add_executable(dftprotein4 dftprotein4.cc ${sources} ${headers})
target_link_libraries(dftprotein4 ${Geant4_LIBRARIES} Geant_v RealPhysics ${ROOT_LIBRARIES} -lpthread ${Vc_LIBRARIES} ${HEPMC_LIBRARIES} ${NUMA_ALL_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY}  Eigen3::Eigen)

#

set(Test_SCRIPTS
gui.mac
init.mac
init_vis.mac
dftprotein4.in
runInGUI.mac
vis.mac
analysis.C
h2.pdb
h2O.pdb
benzene.pdb
nodos.pdb
glycine.pdb
METH.pdb
#benzene_dos.pdb
lys.pdb
methane.pdb
CO.pdb
alanine.pdb
proline.pdb
ethylene.pdb
arginine_neutro.pdb
asparagine.pdb
aspartic_acid.pdb
cysteine.pdb
ethane.pdb
glu.pdb
glutamine.pdb
histidine_neut.pdb
#histidine_nmas.pdb
isoleucine.pdb
leucine.pdb
methionine.pdb
phenol.pdb
phenylalanine.pdb
threonine.pdb
tryptophano.pdb
tyrosine.pdb
ser.pdb
valine.pdb
pyridine.pdb
no.pdb
h2s.pdb
carb_sol.pdb
ox_sol.pdb
n_sol.pdb
P_sol.pdb
S_sol.pdb
#H_solo.pdb
nitrogen_mol.pdb
oxygen.pdb
#cho_pyridine.pdb
#1plx_model_uno.pdb
#1plx_model_uno_sin_h.pdb
)

foreach(_script ${Test_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()


#include(${GEANTV_USE_FILE})
#include_directories(${PROJECT_SOURCE_DIR}/include)

#set(ExampleScripts test_checker.mac README ${CMAKE_SOURCE_DIR}/data/cms2018.gdml
 #                  ${CMAKE_SOURCE_DIR}/data/cmsmagfield2015.txt )

# ${CMAKE_SOURCE_DIR}/data/cms2018.gdml ${CMAKE_SOURCE_DIR}/data/cmsmagfield2015.txt )

# Set to 2 or more to exclude from the default set
#set(ExampleCTestLevel 2)

## The additional source/header files are aglomerated into a single library in
#     examples/physics/CMakeLists.txt
#  To add files to a GeantV example, currently you must edit that file.
#
# file(GLOB ExampleHeaders inc/*.h)
# file(GLOB ExampleSources src/*.cxx)

# set(ExampleHeaders CMSData.h CMSDetectorConstruction.h CMSFieldConstruction.h CMSFullApp.h CMSParticleGun.h CMSPhysicsList.h CMSmagField.h)
# set(ExampleSources CMSData.cxx CMSDetectorConstruction.cxx CMSFieldConstruction.cxx CMSFullApp.cxx CMSParticleGun.cxx CMSPhysicsList.cxx)

# The following can usually be guessed from the directory name.

#set(ExampleName test_checker)
#set(ExampleMainSource ${ExampleName}.cc)
#set(ExampleOutputDir install)

# Generic Part that can be used for other examples.
# The macro is declared in ./examples/CMakeLists.txt


#set(ExampleName pdbprotein)
#set(ExampleMainSource ${ExampleName}.cc)
#set(ExampleOutputDir bin/dftprotein4)


#add_executable(${ExampleName} ${ExampleMainSource})
#add_dependencies(dftprotein4 Geant_v RealPhysics GeantNuma)


#add_dependencies(${ExampleName} Geant_v RealPhysics GeantExamplesRP GeantNuma)
#add_executable(naiver naiver.cc ${sources} ${headers})
#target_link_libraries(naiver ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS dftprotein4 DESTINATION dft/Geant4/dftprotein4)
install(TARGETS dftprotein4 DESTINATION bin/dftprotein4)
install(FILES DESTINATION bin/dftprotein4)



if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()

