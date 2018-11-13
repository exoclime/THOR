#include "dyn/phy_modules_device.h"

#include <cstdio>

const int NUM_PHY_MODULES_DYN_CORE_ARRAYS = 10;


__constant__ device_RK_array dynamical_core_phy_modules_arrays[NUM_PHY_MODULES_DYN_CORE_ARRAYS];
__constant__ int             num_dynamical_arrays;

device_RK_array_manager::device_RK_array_manager() {
}

bool device_RK_array_manager::register_array(double* array_d,
                                             double* arrayk_d,
                                             double* arrayi_d,
                                             int     dimensions) {

    if (data.size() == NUM_PHY_MODULES_DYN_CORE_ARRAYS) {
        printf("Not enough space to allocate array definition for phy_modules\n");
        return false;
    }

    data.emplace_back(array_d, arrayk_d, arrayi_d, dimensions);

    return true;
}

void device_RK_array_manager::allocate_device_array() {
    cudaMemcpyToSymbol(dynamical_core_phy_modules_arrays,
                       data.data(),
                       data.size() * sizeof(device_RK_array));

    int datasize = data.size();

    // maybe this needs a pointer ?
    cudaMemcpyToSymbol(num_dynamical_arrays,
                       &datasize,
                       sizeof(size_t));

    // wait for copy to finish before deallocating local memory
    cudaDeviceSynchronize();
}
