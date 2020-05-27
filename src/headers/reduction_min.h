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
//
//
//
// Description: Computes the reduction minimum value of an array
//
//
// Method: -
//
// Known limitations: None.
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
#pragma once

#include <cmath>
#include <iostream>

#include "log_writer.h"

const double dmax = std::numeric_limits<double>::max();

// GPU reduction sum kernel
template<int BLOCK_SIZE> __global__ void gpu_reduction_min(double *d, double *o, long length) {
    // temporary memory for all tiles in that thread
    __shared__ double ds_in[2 * BLOCK_SIZE];

    // import all the data from global memory
    int mem_offset1 = 2 * (blockDim.x * blockIdx.x + threadIdx.x);

    if (mem_offset1 + 1 < length) {
        *((double2 *)(&(ds_in[2 * threadIdx.x]))) = *((double2 *)(&(d[mem_offset1])));
    }

    else if (mem_offset1 < length) {
        ds_in[2 * threadIdx.x]     = d[mem_offset1];
        ds_in[2 * threadIdx.x + 1] = dmax;
    }
    else {
        ds_in[2 * threadIdx.x]     = dmax;
        ds_in[2 * threadIdx.x + 1] = dmax;
    }


    // loop on stride and add
    for (int stride = blockDim.x; stride > 0; stride /= 2) {
        __syncthreads();
        if (threadIdx.x < stride) {
            ds_in[threadIdx.x] = min(ds_in[threadIdx.x + stride], ds_in[threadIdx.x]);
        }
    }

    __syncthreads();

    // copy to output

    if (threadIdx.x == 0)
        o[blockIdx.x] = ds_in[0];
};

// host function running reduction add on data from device
template<int BLOCK_SIZE> __host__ double gpu_min_on_device(double *in_d, long length) {
    int num_blocks = ceil(double(length) / double(2 * BLOCK_SIZE));

    double *out_h = new double[num_blocks];
    double *out_d;
    // create device temp array
    cudaMalloc((void **)&out_d, num_blocks * sizeof(double));
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
        log::printf("Malloc: %s\n", cudaGetErrorString(err));


    gpu_reduction_min<BLOCK_SIZE><<<num_blocks, BLOCK_SIZE>>>(in_d, out_d, length);
    err = cudaGetLastError();
    if (err != cudaSuccess)
        log::printf("krnl: %s\n", cudaGetErrorString(err));

    cudaMemcpy(out_h, out_d, num_blocks * sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaGetLastError();
    if (err != cudaSuccess)
        log::printf("cpyD2H: %s\n", cudaGetErrorString(err));

    double out = std::numeric_limits<double>::max();
    for (int i = 0; i < num_blocks; i++) {
        out = std::min(out, out_h[i]);
    }

    cudaFree(out_d);
    delete[] out_h;


    return out;
};
