#include "dyn/phy_modules_device.h"

#include <cstdio>

__constant__ device_RK_array dynamical_core_phy_modules_arrays[NUM_PHY_MODULES_DYN_CORE_ARRAYS];
__constant__ int             num_dynamical_arrays[1];

device_RK_array_manager::device_RK_array_manager() {
}

bool device_RK_array_manager::register_array(double* array_d,
                                             double* arrayk_d,
                                             double* arrayi_d,
                                             int     dimensions) {

    if (data.size() == NUM_PHY_MODULES_DYN_CORE_ARRAYS) {
        printf("Not enough space to allocate array definitions for phy_modules\n"
               "  increase NUM_PHY_MODULES_DYN_CORE_ARRAYS\n");
        return false;
    }

    data.emplace_back(array_d, arrayk_d, arrayi_d, dimensions);

    return true;
}

void device_RK_array_manager::allocate_device_array() {
    cudaMemcpyToSymbol(dynamical_core_phy_modules_arrays,
                       data.data(),
                       data.size() * sizeof(device_RK_array));
    {
        cudaError_t err = cudaGetLastError();
        
        // Check device query
        if (err != cudaSuccess) {
            printf("phy: array cuda error: %s\n", cudaGetErrorString(err));
        }
    }
    
    int datasize = data.size();
#ifdef __DEBUG
    printf("Num data: %d\n", datasize);
    
    for (auto & d : data)
        printf("%d %p %p %p\n",
               d.dimensions,
               (void*)d.array_d,
               (void*)d.arrayk_d,
               (void*)d.arrayi_d);
#endif // __DEBUG


    // maybe this needs a pointer ?
    cudaMemcpyToSymbol(num_dynamical_arrays,
                       &datasize,
                       sizeof(int));
    {
        cudaError_t err = cudaGetLastError();
        
        // Check device query
        if (err != cudaSuccess) {
            printf("'phy: num' cuda error: %s\n",  cudaGetErrorString(err));
        }
    }
    
    // wait for copy to finish before deallocating local memory
    cudaDeviceSynchronize();
}
