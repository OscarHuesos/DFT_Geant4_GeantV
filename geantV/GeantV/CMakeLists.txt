#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()
#cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

# Eigen:
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

project(dftproteinV)

#TBB
if (USE_TBB)
   include_directories(${TBB_INCLUDE_DIR})
endif()

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})

include_directories(
  ${PROJECT_SOURCE_DIR}/include
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
  ${CMAKE_SOURCE_DIR}/core/base/inc
  ${CMAKE_SOURCE_DIR}/run/scheduler/inc
  ${CMAKE_SOURCE_DIR}/run/userapp/inc
  ${CMAKE_SOURCE_DIR}/core/interfaces/inc
  ${CMAKE_SOURCE_DIR}/core/concurrency/inc
  ${CMAKE_SOURCE_DIR}/core/numa/inc
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
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#

add_executable(dftproteinV dftproteinV.cc ${sources} ${headers})
target_link_libraries(dftproteinV -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Geant_v RealPhysics ${ROOT_LIBRARIES} -lpthread ${Vc_LIBRARIES} ${HEPMC_LIBRARIES} ${NUMA_ALL_LIBRARIES} Eigen3::Eigen)

#add_executable(${ExampleName} ${ExampleMainSource})
#add_executable(dft pdbprotein.cc pdbprotein.cc ${ExampleMainSource} ${sources} ${headers})
#target_link_libraries(dft -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Geant_v RealPhysics ${ROOT_LIBRARIES} -lpthread ${Vc_LIBRARIES} ${HEPMC_LIBRARIES} ${NUMA_ALL_LIBRARIES})


add_dependencies(dftproteinV Geant_v RealPhysics GeantNuma)


#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build pdb4dna. This is so that we can run the executable directly because
# it relies on these scripts being in the current working directory.
#
set(EXAMPLEpdb4dna_SCRIPTS
gui.mac
init.mac
init_vis.mac
dftproteinV.in
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
histidine_nmas.pdb
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
H_solo.pdb
nitrogen_mol.pdb
oxygen.pdb
cho_pyridine.pdb
1plx_model_uno.pdb
1plx_model_uno_sin_h.pdb
)

foreach(_script ${EXAMPLEpdb4dna_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#


install(TARGETS dftproteinV DESTINATION dft/GeantV/dftproteinV)
install(TARGETS dftproteinV DESTINATION bin/dftproteinV)
install(FILES DESTINATION dft/GeantV/dftproteinV)


if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()
