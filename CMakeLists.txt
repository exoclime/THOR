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
set (CMAKE_CXX_STANDARD 11)
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
  src/test/binary_test.cu)


find_package(HDF5)
if(HDF5_FOUND)
        include_directories(${HDF5_INCLUDE_DIR})
	set(HDF5_LIBS hdf5 hdf5_cpp)
endif()

include_directories("src/headers")
add_executable(esp ${SRC})
target_link_libraries (esp ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES}  )


