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



cmake_minimum_required (VERSION 3.5 FATAL_ERROR)

# check command line args, user can specify SM to override manual setting
set(SM "0" CACHE STRING "GPU SM value")


# check what CMAKE version we use, above 3.8, CUDA is recognised as a
# language and fully integrated. UBUNTU 14 used on the clusters still
# uses cmake 3.5
# set project and cuda variables
if(${CMAKE_VERSION} VERSION_LESS 3.8)
  message(STATUS "detected CMake version ${CMAKE_VERSION} below 3.8")



  
  # looks like CUDA wants gcc/g++ 5
  SET(CMAKE_CXX_COMPILER g++-5)
  SET(CMAKE_C_COMPILER gcc-5)

  project (THOR C CXX)
  # some C++ flags
  set (CMAKE_CXX_STANDARD 11)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set ( COMPILE_FLAGS --compiler-options;-Wall)

  # Find CUDA
  find_package(CUDA REQUIRED)


  # check CUDA for CMake 3.5
  # from here:
  # https://wagonhelm.github.io/articles/2018-03/detecting-cuda-capability-with-cmake
  if (CUDA_FOUND AND SM MATCHES "0")
    #Get CUDA compute capability
    set(OUTPUTFILE ${CMAKE_CURRENT_SOURCE_DIR}/build/cuda_script) # No suffix required
    set(CUDAFILE ${CMAKE_CURRENT_SOURCE_DIR}/cmake/check_cuda.cu)
    execute_process(COMMAND nvcc -lcuda ${CUDAFILE} -o ${OUTPUTFILE})
    execute_process(COMMAND ${OUTPUTFILE}
                            RESULT_VARIABLE CUDA_RETURN_CODE
                            OUTPUT_VARIABLE ARCH)

     if(${CUDA_RETURN_CODE} EQUAL 0)
       set(CUDA_SUCCESS "TRUE")
     else()
       set(CUDA_SUCCESS "FALSE")
     endif()

     if (${CUDA_SUCCESS})
       message(STATUS "CUDA Architecture: ${ARCH}")
       message(STATUS "CUDA Version: ${CUDA_VERSION_STRING}")
       message(STATUS "CUDA Path: ${CUDA_TOOLKIT_ROOT_DIR}")
       message(STATUS "CUDA Libraries: ${CUDA_LIBRARIES}")
       message(STATUS "CUDA Performance Primitives: ${CUDA_npp_LIBRARY}")
       
       set(ARCH_FLAGS "${ARCH}")
       add_definitions(-DGPU) #You may not require this

     else()
       message(WARNING ${ARCH})
     endif()
   else()
     message(STATUS "CUDA Architecture manually set to: -arch=sm_${SM}")
     set(ARCH_FLAGS "-arch=sm_${SM}")
   endif()
else()
  # versions of CMAKE above 3.8
  message(STATUS "detected CMake version ${CMAKE_VERSION} above 3.8")
    
 
  # looks like CUDA wants gcc/g++ 5
  SET(CMAKE_CXX_COMPILER g++-5)
  SET(CMAKE_C_COMPILER gcc-5)

  project (THOR C CXX CUDA)

  # some C++ flags
  set (CMAKE_CXX_STANDARD 11)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set ( COMPILE_FLAGS --compiler-options;-Wall)
  

  # Find CUDA
  find_package(CUDA REQUIRED)
  if (CUDA_FOUND AND SM MATCHES "0")
    CUDA_SELECT_NVCC_ARCH_FLAGS(ARCH_FLAGS Auto)
  else()
    message(STATUS "CUDA Architecture manually set to: -arch=sm_${SM}")
    set(ARCH_FLAGS "-arch=sm_${SM}")
  endif()
  
endif()

#####################################################################
# output directories
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})


# build flafs
LIST(APPEND CUDA_NVCC_FLAGS ${ARCH_FLAGS})
LIST(APPEND CUDA_NVCC_FLAGS ${COMPILE_FLAGS})
LIST(APPEND CUDA_NVCC_FLAGS "-std c++11")

#LIST(APPEND CUDA_NVCC_FLAGS "-ccbin gcc-5")
#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

string (APPEND CMAKE_CUDA_FLAGS " -cudart shared" )

######################################################################
# define sources 
# source for esp main program
set(SRC
  src/esp.cu
  src/grid/grid.cu
  src/initial/esp_initial.cu
  src/initial/planet.cu
  src/output/esp_output.cu
  src/profx/profx_driver.cu
  src/thor/thor_driver.cu
  src/devel/binary_test.cpp
  src/output/storage.cpp
  src/input/config_file.cpp
  src/input/cmdargs.cpp
  src/files/directories.cpp)

# source for test program for storage class
set(SRC_STORAGE
  src/test/storage_test.cpp
  src/output/storage.cpp)

set(SRC_GRID
  src/initial/esp_initial.cu
  src/output/esp_output.cu
  src/initial/planet.cu
  src/profx/profx_driver.cu
  src/thor/thor_driver.cu
  src/grid/grid.cu
  src/test/grid_test.cpp
  src/output/storage.cpp
  src/devel/binary_test.cpp)

set(SRC_CONFIG
  src/test/config_test.cpp
  src/input/config_file.cpp)

set(SRC_CMDARGS
  src/test/cmdargs_test.cpp
  src/input/cmdargs.cpp)

set(SRC_DIRECTORIES
  src/test/directories_test.cpp
  src/files/directories.cpp)
######################################################################
# FIND HDF5
find_package(HDF5 COMPONENTS CXX)
if(HDF5_FOUND)
        include_directories(${HDF5_INCLUDE_DIR})
	set(HDF5_LIBS hdf5 hdf5_cpp)
endif()

message(STATUS "HDF5 C++libraries " ${HDF5_LIBRARIES})
include_directories("src/headers")

######################################################################
# build application
cuda_add_executable(esp ${SRC} OPTIONS -DBUILD_LEVEL="release")
target_link_libraries (esp ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES}  )

# build tests
cuda_add_executable(test_storage ${SRC_STORAGE})
target_link_libraries (test_storage ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES}  )


add_executable(test_args ${SRC_CMDARGS})
add_executable(test_config ${SRC_CONFIG})
add_executable(test_directories ${SRC_DIRECTORIES})

cuda_add_executable(test_grid ${SRC_GRID} OPTIONS -DBENCHMARKING -DBENCH_POINT_COMPARE)
target_compile_definitions(test_grid PRIVATE -DBENCHMARKING -DBENCH_POINT_COMPARE)
target_link_libraries (test_grid ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES}  )
