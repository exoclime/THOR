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
// Description: binary correctness test of output, enabled with compile time switches
//
//
//
// Method: [1] - Dumps output to binary file on a flag
//         [2] - Reads data from binary files on a flag and compare to
//               dumped output
//
// Known limitations: None.
//
//
// Known issues: None.
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
#include <iostream>

#include "binary_test.h"
#include "directories.h"
#include "log_writer.h"

__host__ void getDeviceData(const double *device, double *host, int size) {
    cudaMemcpy(host, device, size, cudaMemcpyDeviceToHost);
}

bool *init_device_mem_check(bool *ptr) {
    cudaMalloc((void **)&ptr, sizeof(bool));

    return ptr;
}

void deinit_device_mem_check(bool *ptr) {
    cudaFree(ptr);
}


// returns true if data contains a NaN
__global__ void isnan_check_device(double *array, int size, bool *check) {
    //
    //  Description: Check for nan in array.

    int idx = threadIdx.x + blockDim.x * blockIdx.x;

    if (idx < size && ::isnan(array[idx])) {
        *check = true;
    }
}


// check data for NaNs
// returns true if data contains a NaN
__host__ bool check_array_for_nan(double *ptr, int size, bool on_device, bool *check_d) {
    // TODO: could probably act on multiple arrays through streams
    if (on_device) {
        bool check_h = false;
        cudaMemcpy(check_d, &check_h, sizeof(bool), cudaMemcpyHostToDevice);
        isnan_check_device<<<size / 256 + 1, 256>>>(ptr, size, check_d);
        cudaMemcpy(&check_h, check_d, sizeof(bool), cudaMemcpyDeviceToHost);

        return check_h;
    }
    else {
        bool isnan = false;
        for (int i = 0; i < size; i++) {
            isnan |= std::isnan(ptr[i]); //
        }
        return isnan;
    }
}

#ifdef BENCHMARKING
void crash_report(const output_def &def, std::string output_dir, const std::string &iteration) {
    int  first = 1;
    path o(output_dir);
    o /= ("crash_report_" + iteration + "_" + def.short_name + ".txt");
    std::fstream crash_report_file;
    crash_report_file.open(o.to_string(), std::ofstream::out);
    if (def.device_ptr) {
        double *ptr_h;
        ptr_h = new double[def.size];
        cudaMemcpy(ptr_h, def.data, def.size * sizeof(double), cudaMemcpyDeviceToHost);
        for (int i = 0; i < def.size; i++) {
            if (std::isnan(ptr_h[i])) {
                crash_report_file << def.fun(i, first) << std::endl;
                first = 0;
            }
        }
    }
    else {
        for (int i = 0; i < def.size; i++) {
            if (std::isnan(def.data[i])) {
                crash_report_file << def.fun(i, first) << std::endl;
                first = 0;
            }
        }
    }
}
#endif

void check_last_cuda_error(std::string ref_name) {
    cudaError_t err = cudaGetLastError();

    // Check device query
    if (err != cudaSuccess) {
        log::printf("'%s' cuda error: %s\n", ref_name.c_str(), cudaGetErrorString(err));
    }
}
