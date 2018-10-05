
project (phy_modules)

set(MODULES_SRCS
  src/phy_modules.cu
  ${SHARED_PHY_MODULES_SRC_DIR}/radiative_transfer.cu
  )


include_directories(${THOR_ROOT_INCLUDES} ${SHARED_PHY_MODULES_INC_DIR})
include_directories(inc)
cuda_add_library(phy_modules STATIC ${MODULES_SRCS})
