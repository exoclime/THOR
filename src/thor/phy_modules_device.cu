// ==============================================================================
// This file is part of THOR.
//
//     THOR is free software : you can redistribute it and / or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
//
//     THOR is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
//     GNU General Public License for more details.
//
//     You find a copy of the GNU General Public License in the main
//     THOR directory under <license.txt>.If not, see
//     <http://www.gnu.org/licenses/>.
// ==============================================================================
//
// Manager routines for physics modules that need direct access to dynamical core
//
// Known limitations: - Runs in a single GPU.
//
// Known issues: None
//
//
// If you use this code please cite the following reference:
//
//       [1] Mendonca, J.M., Grimm, S.L., Grosheintz, L., & Heng, K., ApJ, 829, 115, 2016
//
// Current Code Owners: Joao Mendonca (joao.mendonca@space.dtu.dk)
//                      Russell Deitrick (russell.deitrick@csh.unibe.ch)
//                      Urs Schroffenegger (urs.schroffenegger@csh.unibe.ch)
//
// History:
// Version Date       Comment
// ======= ====       =======
// 2.0     30/11/2018 Released version (RD & US)
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////



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
