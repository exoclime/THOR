cmake_minimum_required (VERSION 3.8 FATAL_ERROR)

# to build localy with tests
# cmake -DBUILD_TESTS=ON -DTHOR_BASE_FOLDER=../ ..  

# check command line args, user can specify SM to override manual setting
set(SM "0" CACHE STRING "GPU SM value")

# if CUDA wants specific gcc/g++ 
#SET(CMAKE_CXX_COMPILER g++-8)
#SET(CMAKE_C_COMPILER gcc-8)

option(BUILD_TESTS "Build test programs" OFF)

if (DEFINED THOR_MODULE)
  message(STATUS "Building as THOR module")
  
  set(THOR_BASE_FOLDER ${THOR_ROOT})
  message("THOR_BASE_FOLDER: " ${THOR_BASE_FOLDER})
else()
  message(STATUS "Building as python library module")
  project (Alfrodull CUDA CXX)
  # some C++ flags
  set (CMAKE_CXX_STANDARD 14)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set ( COMPILE_FLAGS --compiler-options;-Wall)
  
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

  if (NOT DEFINED THOR_BASE_FOLDER)
    set(THOR_BASE_FOLDER "../")
  endif()
  
  
  include(CheckLanguage)
  check_language(CUDA)
  
  set(CMAKE_CUDA_STANDARD 14)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  
  # Find CUDA
  #find_package(CUDA REQUIRED)
  #if (CUDA_FOUND AND SM MATCHES "0")
  #  CUDA_SELECT_NVCC_ARCH_FLAGS(ARCH_FLAGS Auto)
  #else()
  #  message(STATUS "CUDA Architecture manually set to: -arch=sm_${SM}")
  #  set(ARCH_FLAGS "-arch=sm_${SM}")
  #endif()
  
  #SET(CUDA_SEPARABLE_COMPILATION ON)
  
  # build flags
  #LIST(APPEND CUDA_NVCC_FLAGS ${ARCH_FLAGS})
  #LIST(APPEND CUDA_NVCC_FLAGS ${COMPILE_FLAGS})
  #LIST(APPEND CUDA_NVCC_FLAGS "-std c++11")
  
  #string (APPEND CMAKE_CUDA_FLAGS " -dc -cudart shared -rdc=true -lcudart -lcudadevrt -lcuda" )
  
  find_package(pybind11 REQUIRED)
  #find_package(PythonLibs REQUIRED)
  
  
  ######################################################################
  # FIND HDF5
  find_package(HDF5 REQUIRED COMPONENTS CXX)
  if(HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIR})
    set(HDF5_LIBS hdf5 hdf5_cpp)
  endif()
  
  message(STATUS "HDF5 include dir: " ${HDF5_INCLUDE_DIR})
  message(STATUS "HDF5 C++libraries: " ${HDF5_LIBRARIES})
  ######################################################################
endif()


set(SRC_ALFRODULLIB
	src/alfrodullib.cu
	src/alfrodull_engine.cu
	src/planck_table.cu
	src/gauss_legendre_weights.cpp
	src/opacities/opacities.cu
	src/opacities/cloud_opacities.cu
	src/kernels/calculate_physics.cu
	src/kernels/math_helpers.cu
	src/kernels/interpolate_values.cu
	src/kernels/integrate_flux.cu)


	
include_directories(src/inc)

# THOR utilities
set(SRC_THOR_UTILS
  ${THOR_BASE_FOLDER}/src/utils/storage.cpp
  ${THOR_BASE_FOLDER}/src/utils/cuda_device_memory.cu
  )


######################################################################
# build and link binaries
add_library(alfrodullib STATIC ${SRC_ALFRODULLIB})
set_target_properties(alfrodullib PROPERTIES 
    POSITION_INDEPENDENT_CODE ON
    CUDA_VISIBILITY_PRESET "hidden"
    CUDA_SEPARABLE_COMPILATION ON
)

if (DEFINED THOR_MODULE)
  include_directories(${THOR_ROOT_INCLUDES} ${SHARED_PHY_MODULES_INC_DIR})

  include_directories(${THOR_BASE_FOLDER}/src/physics/modules/inc/)
    
  set(SRC_THOR_CLASS
    src/two_streams_radiative_transfer.cu)

  
  set(LIBRARY_OUTPUT_PATH)

  ######################################################################
  # build and link binaries
  add_library(alfrodull STATIC ${SRC_THOR_CLASS})
  set_target_properties(alfrodull PROPERTIES 
    POSITION_INDEPENDENT_CODE ON
    CUDA_VISIBILITY_PRESET "hidden"
    CUDA_SEPARABLE_COMPILATION ON
	)
      target_link_libraries(alfrodull PRIVATE alfrodullib)
#      export(TARGETS alfrodull FILE libalfrodull.a)
else()
  # specify CGS units
  target_compile_definitions(alfrodullib PRIVATE "CGS_UNITS")
  # build as independent python library
  include_directories(${THOR_BASE_FOLDER}/src/headers/ )
  include_directories(pylfrodull/src/inc)
  set(SRC_PYMODULE
    pylfrodull/src/pylfrodull.cpp)

  # THOR utilities
  set(SRC_THOR_UTILS
    ${THOR_BASE_FOLDER}/src/utils/storage.cpp
    ${THOR_BASE_FOLDER}/src/utils/cuda_device_memory.cu
    )

  add_library(thor_utils STATIC ${SRC_THOR_UTILS})

  set_target_properties(thor_utils PROPERTIES 
    POSITION_INDEPENDENT_CODE ON
    CUDA_VISIBILITY_PRESET "hidden"
    CUDA_SEPARABLE_COMPILATION ON
    )

  
  add_library(pylfrodull SHARED ${SRC_PYMODULE})
  #target_link_libraries (pylfrodull alfrodullib)
  
  set_target_properties(pylfrodull PROPERTIES 
    CXX_VISIBILITY_PRESET "hidden"
    #				         INTERPROCEDURAL_OPTIMIZATION TRUE
    PREFIX "${PYTHON_MODULE_PREFIX}"
    SUFFIX "${PYTHON_MODULE_EXTENSION}"
    )
  
  target_link_libraries(pylfrodull PRIVATE alfrodullib)
  target_link_libraries(pylfrodull PRIVATE ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES})
  target_link_libraries(pylfrodull PRIVATE thor_utils)
  target_link_libraries(pylfrodull PRIVATE pybind11::module )
  #target_link_libraries(pylfrodull PRIVATE ${PYTHON_LIBRARIES} ${CUDA_LIBRARIES})
  # set_target_properties(pylfrodull PROPERTIES
  # 				 POSITION_INDEPENDENT_CODE ON
  # #				 CUDA_VISIBILITY_PRESET "hidden"
  # 				 CXX_VISIBILITY_PRESET "hidden"
  # 				 PREFIX "${PYTHON_MODULE_PREFIX}"
  #                                  SUFFIX "${PYTHON_MODULE_EXTENSION}")
endif()

#######################################################################
# build and execute tests
if (BUILD_TESTS)
  #include_directories(tests/inc)
  include_directories(${THOR_ROOT_INCLUDES})

  include_directories(${THOR_BASE_FOLDER}/src/headers/ )

  
  # build as test
  set(SRC_TEST_INTEGRATE_FLUX
    tests/integrate_flux/flux_int_test.cu
    )


  set(SRC_TEST_TWO_STREAMS_SOLVE
    tests/integrate_flux/two_streams_solve_test.cu
    src/kernels/integrate_flux.cu
    src/kernels/calculate_physics.cu
    src/kernels/math_helpers.cu
    )
  
  
  set(LIBRARY_OUTPUT_PATH)

  add_executable(test_two_streams_solve ${SRC_TEST_TWO_STREAMS_SOLVE} ${SRC_THOR_UTILS})
  set_target_properties(test_two_streams_solve PROPERTIES 
    POSITION_INDEPENDENT_CODE ON
    CUDA_VISIBILITY_PRESET "hidden"
    CUDA_SEPARABLE_COMPILATION ON
    )
  #target_link_libraries(test_integrate_flux PRIVATE alfrodull)
  target_link_libraries(test_two_streams_solve PRIVATE ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES})
  

  # Disable integration comparison test, code to compare to is not
  # supported anymore, so can't be compiled
  # was used in development when porting code from HELIOS to Alf before making it evolve
  # add_executable(test_integrate_flux ${SRC_TEST_INTEGRATE_FLUX} ${SRC_THOR_UTILS})
  # set_target_properties(test_integrate_flux PROPERTIES 
  #   POSITION_INDEPENDENT_CODE ON
  #   CUDA_VISIBILITY_PRESET "hidden"
  #   CUDA_SEPARABLE_COMPILATION ON
  #   )
  # target_link_libraries(test_integrate_flux PRIVATE alfrodull)
  # target_link_libraries(test_integrate_flux PRIVATE ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES})
endif()
