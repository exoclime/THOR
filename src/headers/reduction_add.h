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
// Description: Computes the reduction sum of an array
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

// CPU reduction sum function
double cpu_reduction_sum(double *d, long length);


// GPU reduction sum kernel
template<int BLOCK_SIZE> __global__ void gpu_reduction_sum(double *d, double *o, long length) {
    // temporary memory for all tiles in that thread
    __shared__ double ds_in[2 * BLOCK_SIZE];

    // import all the data from global memory
    int mem_offset1 = 2 * (blockDim.x * blockIdx.x + threadIdx.x);

    if (mem_offset1 + 1 < length) {
        *((double2 *)(&(ds_in[2 * threadIdx.x]))) = *((double2 *)(&(d[mem_offset1])));
    }

    else if (mem_offset1 < length) {
        ds_in[2 * threadIdx.x]     = d[mem_offset1];
        ds_in[2 * threadIdx.x + 1] = 0.0f;
    }
    else {
        ds_in[2 * threadIdx.x]     = 0.0f;
        ds_in[2 * threadIdx.x + 1] = 0.0f;
    }


    // loop on stride and add
    for (int stride = blockDim.x; stride > 0; stride /= 2) {
        __syncthreads();
        if (threadIdx.x < stride)
            ds_in[threadIdx.x] += ds_in[threadIdx.x + stride];
    }

    __syncthreads();

    // copy to output

    if (threadIdx.x == 0)
        o[blockIdx.x] = ds_in[0];
};

// host function running reduction add on data from device
template<int BLOCK_SIZE> __host__ double gpu_sum_on_device(double *in_d, long length) {
    int num_blocks = ceil(double(length) / double(2 * BLOCK_SIZE));

    double *out_h = new double[num_blocks];
    double *out_d;
    // create device temp array
    cudaMalloc((void **)&out_d, num_blocks * sizeof(double));
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
        log::printf("Malloc: %s\n", cudaGetErrorString(err));


    gpu_reduction_sum<BLOCK_SIZE><<<num_blocks, BLOCK_SIZE>>>(in_d, out_d, length);
    err = cudaGetLastError();
    if (err != cudaSuccess)
        log::printf("krnl: %s\n", cudaGetErrorString(err));

    cudaMemcpy(out_h, out_d, num_blocks * sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaGetLastError();
    if (err != cudaSuccess)
        log::printf("cpyD2H: %s\n", cudaGetErrorString(err));

    double out = 0.0;
    for (int i = 0; i < num_blocks; i++) {
        out += out_h[i];
    }

    cudaFree(out_d);
    delete[] out_h;


    return out;
};


// host function running reduction add on data from host
template<int BLOCK_SIZE> __host__ double gpu_sum_from_host(double *d, long length) {
    double *in_d;


    cudaMalloc((void **)&in_d, length * sizeof(double));
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
        log::printf("Malloc: %s\n", cudaGetErrorString(err));


    cudaMemcpy(in_d, d, length * sizeof(double), cudaMemcpyHostToDevice);
    err = cudaGetLastError();
    if (err != cudaSuccess)
        log::printf("cpyH2D: %s\n", cudaGetErrorString(err));

    double out = gpu_sum_on_device<BLOCK_SIZE>(in_d, length);


    cudaFree(in_d);

    return out;
};


template<int BLOCK_SIZE> double cpu_sum(double *d, long length) {
    long num_blocks = ceil(double(length) / double(2 * BLOCK_SIZE));
    long bl         = 2 * BLOCK_SIZE;


    double buf[2 * BLOCK_SIZE];

    double out = 0.0;

    for (long i = 0l; i < num_blocks; i++) {
        for (volatile long j = 0l; j < bl; j++) {
            long idx = i * bl + j;

            if (idx < length)
                buf[j] = d[idx];
            else
                buf[j] = 0.0;
        }


        double o = cpu_reduction_sum(buf, BLOCK_SIZE);
        out += o;
    }

    return out;
};


template<int BLOCK_SIZE> double cpu_linear_sum(double *d, long length) {

    double out = 0.0;
    for (long i = 0l; i < length; i++) {
        out += d[i];
    }


    return out;
};
