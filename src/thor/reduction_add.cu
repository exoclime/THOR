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
// Current Code Owner: Joao Mendonca, EEG. joao.mendonca@csh.unibe.ch,
//                     Russell Deitrick, russell.deitrick@csh.unibe.ch
//                     Urs Schroffenegger, urs.schroffenegger@csh.unibe.ch
//
//
// History:
// Version Date       Comment
// ======= ====       =======
//
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////

#include "reduction_add.h"



const long MAX_BLOCK_SIZE = 4096;


double cpu_reduction_sum(double * d, long length)
{
    
    for (int stride = length; stride > 0; stride /= 2)
    {
        for (int i = 0; i < stride; i++)
            d[i] += d[i+stride];
    }

    return d[0];
}

double buf[MAX_BLOCK_SIZE];

template<int BLOCK_SIZE>
double cpu_sum(double *d, long length)
{
    long num_blocks = ceil(double(length)/double(2*BLOCK_SIZE));
    

    

    double out = 0.0;
    
    for (long i = 0; i < num_blocks; i++)
    {
        for (long j = 0; j < 2*BLOCK_SIZE; j++)
        {
            long idx = i*2*BLOCK_SIZE + j;
            
            if (idx < length)
                buf[j] = d[idx];
            else
                buf[j] = 0.0;
        }
        
            
        double o = cpu_reduction_sum(buf, BLOCK_SIZE);
        //    printf("%d: %g %g\n", i, o, out);
        out += o;
        
        
    }


    

    return out;
}


template<int BLOCK_SIZE>
__global__ void gpu_reduction_sum(double * d,
                       double * o,
                       long length)
{
   // temporary memory for all tiles in that thread
    __shared__ double ds_in[2*BLOCK_SIZE];  
 
    // import all the data from global memory
    int mem_offset1 = 2*(blockDim.x*blockIdx.x + threadIdx.x);

    if (mem_offset1 + 1 < length)
    {
        *((double2*)(&(ds_in[2*threadIdx.x]))) = *((double2*)(&(d[mem_offset1])));
    }

    else if  (mem_offset1 < length)
    {
        ds_in[2*threadIdx.x] = d[mem_offset1];
        ds_in[2*threadIdx.x + 1] = 0.0f;
    }
    else
    {
        ds_in[2*threadIdx.x] = 0.0f;
        ds_in[2*threadIdx.x + 1] = 0.0f;
    }
    
    
    // loop on stride and add
    for (int stride = blockDim.x; stride > 0; stride /= 2)
    {
        __syncthreads();
        if (threadIdx.x < stride)
            ds_in[threadIdx.x] += ds_in[threadIdx.x + stride];
    }
    
    __syncthreads();
    
    // copy to output
    
    if (threadIdx.x == 0)
       o[blockIdx.x] = ds_in[0];
    
}

        
