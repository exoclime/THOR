# to build with cmake
# create a build directory and move into it
# $ mkdir build
# $ cd build
# generate the makefile (to do only ones, if we don't add files or change makefiles)
# don't forget the two points at the end of the command '..'.
# It runs cmake in the 'build' directory
# but with the data from the '..' directory.
# the script should autodetect the CUDA architecture, when run alone
# $ cmake ..
# You can specify the architecture with the SM variable:
# $ cmake -DSM=30 ..
# compile (-jX runs X jobs in parallel, useful with multiple cores to speed up)
# $ make -j8
# to debug makefile by showing commands
# $ make VERBOSE=1

# specify the path to physical modules
# $ cmake -DMODULES_SRC="src/physics/managers/multi" ..
# build the Alfrodull module
# $ cmake -DMODULES_SRC=Alfrodull ..
# debug or release
# release build
# $ cmake -DCMAKE_BUILD_TYPE=Release ..
# debug build
# $ cmake -DCMAKE_BUILD_TYPE=Debug ..

#######################################################
# Config for travis-ci

# basic setup for cmake
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON)
set(CMAKE_COLOR_MAKEFILE ON)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_STANDARD 11)
# disable in source builds
# this is only a temporary fix, but for now we need it as cmake will
# otherwise overwrite the existing makefiles
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

#
# PROJECT is QUDA
#
project (THOR CUDA CXX)
  find_package(CUDA REQUIRED)
  set(USING_CUDA_LANG_SUPPORT False)
#######################################################

cmake_minimum_required (VERSION 3.8 FATAL_ERROR)

MESSAGE(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# check command line args, user can specify SM to override manual setting
set(SM "0" CACHE STRING "GPU SM value")
set(MODULES_SRC "0" CACHE STRING "Physics module subdirectory")

if (MODULES_SRC MATCHES "0")
  #set(MODULES_SRC "Alfrodull")
  #set(MODULES_SRC "src/physics/managers/multi/")
  #set(MODULES_SRC "src/physics/managers/empty/")
endif()

# if CUDA wants gcc/g++ 5
# SET(CMAKE_CXX_COMPILER g++-7)
# SET(CMAKE_C_COMPILER gcc-7)

# versions of CMAKE above 3.8
message(STATUS "detected CMake version ${CMAKE_VERSION}")
  
project (THOR CUDA CXX)

if (BUILD_TYPE MATCHES "Debug")
  set(CMAKE_CUDA_FLAGS ${CMAKE_CUDA_FLAGS_DEBUG})
endif()

# some C++ flags
set (CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
string (APPEND CMAKE_CXX_FLAGS " -Wall -Wpedantic ")
string (APPEND CMAKE_CUDA_FLAGS " --compiler-options -Wall ")

include(CheckLanguage)
check_language(CUDA)
  

  # Find CUDA
  #  find_package(CUDA REQUIRED)
  #if (CUDA_FOUND AND SM MATCHES "0")
  #  CUDA_SELECT_NVCC_ARCH_FLAGS(ARCH_FLAGS Auto)
  #else()
  #  message(STATUS "CUDA Architecture manually set to: -arch=sm_${SM}")
  #  set(ARCH_FLAGS "-arch=sm_${SM}")
  #endif()
  

#####################################################################
# output directories
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})


# build flafs
#LIST(APPEND CUDA_NVCC_FLAGS ${ARCH_FLAGS})
#LIST(APPEND CUDA_NVCC_FLAGS ${COMPILE_FLAGS})
#LIST(APPEND CUDA_NVCC_FLAGS "-std c++14")

#LIST(APPEND CUDA_NVCC_FLAGS "-ccbin gcc-5")
#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

string (APPEND CMAKE_CUDA_FLAGS " -cudart shared" )

set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)


######################################################################
# define sources 
# source for esp main program
set(SRC
  src/esp.cu
  src/ESP/grid.cu
  src/ESP/esp_initial.cu
  src/ESP/simulation_setup.cu
  src/ESP/esp_output.cu
  src/ESP/profx_driver.cu
  src/ESP/thor_driver.cu
  src/ESP/insolation.cu
  src/utils/binary_test.cpp
  src/utils/debug_helpers.cu
  src/utils/storage.cpp
  src/utils/config_file.cpp
  src/utils/cmdargs.cpp
  src/utils/directories.cpp
  src/utils/log_writer.cpp
  src/utils/iteration_timer.cpp
  src/ESP/profx_globdiag.cu
  src/ESP/profx_sponge.cu
  src/ESP/reduction_add.cu
  src/ESP/phy_modules_device.cu
  src/ESP/ultrahot_thermo.cu
  src/utils/cuda_device_memory.cu
  src/physics/managers/multi/src/phy_modules.cu
  src/physics/modules/src/radiative_transfer.cu
  src/physics/modules/src/boundary_layer.cu
  src/physics/modules/src/chemistry.cu
  src/utils/diagnostics.cu
  )

# source for test program for storage class
set(SRC_STORAGE
  src/test/storage_test.cpp
  src/utils/storage.cpp
  src/utils/directories.cpp
  src/utils/log_writer.cpp)

set(SRC_GRID
  src/ESP/grid.cu
  src/test/grid_test.cpp
  src/utils/storage.cpp
  src/utils/binary_test.cpp
  src/utils/debug_helpers.cu
  src/utils/directories.cpp
  src/utils/log_writer.cpp)

set(SRC_CONFIG
  src/test/config_test.cpp
  src/utils/config_file.cpp)

set(SRC_CMDARGS
  src/test/cmdargs_test.cpp
  src/utils/cmdargs.cpp)

set(SRC_DIRECTORIES
  src/utils/log_writer.cpp
  src/test/directories_test.cpp
  src/utils/directories.cpp)

set(SRC_GEN_INIT
  src/test/gen_init.cpp
  src/utils/storage.cpp
  src/utils/log_writer.cpp
  src/utils/directories.cpp
  src/ESP/grid.cu
  src/ESP/simulation_setup.cu)

set(SRC_REDUCTION_ADD
  src/ESP/reduction_add.cu
  src/test/reduction_add_test.cu)

######################################################################
# FIND HDF5
find_package(HDF5 COMPONENTS CXX)
if(HDF5_FOUND)
        include_directories(${HDF5_INCLUDE_DIR})
	set(HDF5_LIBS hdf5 hdf5_cpp)
endif()

message(STATUS "HDF5 C++libraries " ${HDF5_LIBRARIES})
include_directories("src/headers")
include_directories("src/physics/modules/inc")
######################################################################
# get git version information
# from https://stackoverflow.com/questions/1435953/how-can-i-pass-git-sha1-to-compiler-as-definition-using-cmake/4318642#4318642

# The 'real' git information file
SET(GITREV_BARE_FILE git-rev.h)
# The temporary git information file
SET(GITREV_BARE_TMP git-rev-tmp.h)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/ )
SET(GITREV_FILE ${CMAKE_CURRENT_BINARY_DIR}/include/${GITREV_BARE_FILE})
SET(GITREV_TMP ${CMAKE_CURRENT_BINARY_DIR}/include/${GITREV_BARE_TMP})

include_directories(${CMAKE_CURRENT_BINARY_DIR}/include/)

set(GIT_BRANCH_COMMAND "git rev-parse --abbrev-ref HEAD | sed \'s/\(.*\)/\"\\1\"/g\'" )


ADD_CUSTOM_COMMAND(
  OUTPUT ${GITREV_TMP}
  COMMAND echo -n \\\#define GIT_BRANCH_RAW\  > ${GITREV_TMP}
  COMMAND echo \\\"\$\$\(git rev-parse --abbrev-ref HEAD\)\\\" >> ${GITREV_TMP} 
  COMMAND echo -n \\\#define GIT_HASH_RAW\  >> ${GITREV_TMP}
  COMMAND echo \\\"\$\$\(git describe --always --dirty --abbrev=40 --match=\"NoTagWithThisName\"\)\\\" >> ${GITREV_TMP} 
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${GITREV_TMP} ${GITREV_FILE}
#  COMMAND ${CMAKE_COMMAND} -E remove ${GITREV_TMP}
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
#  VERBATIM
)
# Finally, the temporary file should be added as a dependency to the target

######################################################################
# define main application to build 
add_executable(esp ${SRC} ${GITREV_TMP})

######################################################################
# physics modules subdirectory

# set(THOR_MODULE "1")

# file(RELATIVE_PATH THOR_ROOT "${PROJECT_SOURCE_DIR}/${MODULES_SRC}" "${PROJECT_SOURCE_DIR}/" )
# message("Relative path to THOR root from physics module: " ${THOR_ROOT})
# set(THOR_ROOT_INCLUDES "${THOR_ROOT}src/headers" )
# set(SHARED_PHY_MODULES_INC_DIR "${THOR_ROOT}src/physics/modules/inc" )
# set(SHARED_PHY_MODULES_SRC_DIR "${THOR_ROOT}src/physics/modules/src" )

# add_subdirectory(${MODULES_SRC})


# set(PHYSICS_LIBRARY_PATH "${MODULES_SRC}/libphy_modules.a")
if (EXISTS "${PROJECT_SOURCE_DIR}/Alfrodull" AND EXISTS "${PROJECT_SOURCE_DIR}/Alfrodull/CMakeLists.txt")
  message("Alfrodull library detected")
  set(THOR_MODULE "1")
  file(RELATIVE_PATH THOR_ROOT "${PROJECT_SOURCE_DIR}/${MODULES_SRC}" "${PROJECT_SOURCE_DIR}/" )
  message("Relative path to THOR root from physics module: " ${THOR_ROOT})
  set(THOR_ROOT_INCLUDES "${THOR_ROOT}src/headers" )
  set(SHARED_PHY_MODULES_INC_DIR "${THOR_ROOT}src/physics/modules/inc" )
  set(SHARED_PHY_MODULES_SRC_DIR "${THOR_ROOT}src/physics/modules/src" )
  add_subdirectory(Alfrodull)
  set(ALFRODULL_INCLUDE_DIRS
    ${PROJECT_SOURCE_DIR}/Alfrodull/src/inc/)
  include_directories(${ALFRODULL_INCLUDE_DIRS})
  target_compile_definitions(esp PRIVATE "HAS_ALFRODULL")
  target_link_libraries (esp alfrodull)
  target_link_directories (esp PRIVATE Alfrodull)
else()
  message("Alfrodull library not detected")
endif()

######################################################################
# build application

target_compile_options(esp  PUBLIC -DBUILD_LEVEL="release" -DDEVICE_SM=${SM})
# options for debugging
#target_compile_options(esp PUBLIC -DBENCHMARKING -DBENCH_NAN_CHECK)
# target_link_libraries (esp ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES} phy_modules)
target_link_libraries (esp ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES})

# build tests
add_executable(test_storage EXCLUDE_FROM_ALL ${SRC_STORAGE})
target_link_libraries (test_storage ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES}  )


add_executable(test_args EXCLUDE_FROM_ALL ${SRC_CMDARGS})
add_executable(test_config EXCLUDE_FROM_ALL ${SRC_CONFIG})
add_executable(test_directories EXCLUDE_FROM_ALL ${SRC_DIRECTORIES})
add_executable(test_reduction_add EXCLUDE_FROM_ALL ${SRC_REDUCTION_ADD})
add_executable(test_gen_init EXCLUDE_FROM_ALL ${SRC_GEN_INIT})
target_link_libraries (test_gen_init ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES}  )

add_executable(test_grid EXCLUDE_FROM_ALL ${SRC_GRID})
target_compile_options(test_grid PUBLIC -DBENCHMARKING -DBENCH_POINT_COMPARE -DBENCH_PRINT_DEBUG)

target_link_libraries (test_grid ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES}  )

set(TEST_NAMES
	test_storage
	test_args
	test_config
	test_directories
	test_gen_init
	test_grid)

add_custom_target(tests EXCLUDE_FROM_ALL DEPENDS ${TEST_NAMES})
