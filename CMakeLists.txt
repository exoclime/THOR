# to build with cmake
# create a build directory and move into it
# $ mkdir build
# $ cd build
# generate the makefile (to do only ones, if we don't add files or change makefiles)
# don't forget the two points at the end of the command '..'.
# It runs cmake in the 'build' directory
# but with the data from the '..' directory
# $ cmake -DCMAKE_CUDA_FLAGS='-arch=sm_30' ..
# compile (-jX runs X jobs in parallel, useful with multiple cores to speed up)
# $ make -j8
# to debug makefile by showing commands
# $ make VERBOSE=1



cmake_minimum_required (VERSION 3.5 FATAL_ERROR)

project (THOR C CXX CUDA)
set (CMAKE_CXX_STANDARD 14)
set ( COMPILE_FLAGS --compiler-options;--Wall) 
find_package(CUDA REQUIRED)

CUDA_SELECT_NVCC_ARCH_FLAGS(ARCH_FLAGS Auto)
LIST(APPEND CUDA_NVCC_FLAGS ${ARCH_FLAGS})
LIST(APPEND CUDA_NVCC_FLAGS ${COMPILE_FLAGS})
#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

string (APPEND CMAKE_CUDA_FLAGS " -cudart shared")


set(SRC
  src/esp.cu
  src/grid/grid.cu
  src/initial/esp_initial.cu
  src/initial/planet.cu
  src/output/esp_output.cu
  src/profx/profx_driver.cu
  src/thor/thor_driver.cu
  src/devel/binary_test.cu
  src/output/storage.cpp)

set(SRC_STORAGE
  src/test/storage_test.cpp
  src/output/storage.cpp)




find_package(HDF5 COMPONENTS CXX)
if(HDF5_FOUND)
        include_directories(${HDF5_INCLUDE_DIR})
	set(HDF5_LIBS hdf5 hdf5_cpp)
endif()

message(STATUS "HDF5 C++libraries " ${HDF5_LIBRARIES})
include_directories("src/headers")
add_executable(esp ${SRC})
target_link_libraries (esp ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES}  )
add_executable(test_storage ${SRC_STORAGE})
target_link_libraries (test_storage ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES}  )

