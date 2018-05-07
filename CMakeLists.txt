# to build with cmake
# mkdir build
# cd build
# generate the makefile (to do only ones, if we don't add files or change makefiles)
# cmake -DCMAKE_CUDA_FLAGS='-arch=sm_30' ..
# compile
# make -j8
# to debug makefile by showing commands
# make VERBOSE=1



cmake_minimum_required (VERSION 3.9 FATAL_ERROR)

project (THOR C CXX CUDA)


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
  src/thor/thor_driver.cu)


find_package(HDF5)
if(HDF5_FOUND)
        include_directories(${HDF5_INCLUDE_DIR})
	set(HDF5_LIBS hdf5 hdf5_cpp)
endif()
      
add_executable(esp ${SRC})
target_link_libraries (esp ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES}  )


