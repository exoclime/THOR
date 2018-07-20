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
// Current Code Owner: Joao Mendonca, EEG. joao.mendonca@csh.unibe.ch
//
// If you use this code please cite the following reference: 
//
//       [1] Mendonca, J.M., Grimm, S.L., Grosheintz, L., & Heng, K., ApJ, 829, 115, 2016  
//
// History:
// Version Date       Comment
// ======= ====       =======
//
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////
#include <iostream>


__host__ void getDeviceData(const double * device, double * host, int size)
{
    cudaMemcpy(host, device, size, cudaMemcpyDeviceToHost);
}

bool * init_device_mem_check(bool * ptr)
{
    cudaMalloc((void **)&ptr, sizeof (bool));

    return ptr;
    
}

void deinit_device_mem_check(bool *ptr)
{
    cudaFree(ptr);
}



__global__ void isnan_check_device(double *array, int size, bool *check)
{
//
//  Description: Check for nan in array.

  int idx = threadIdx.x+blockDim.x*blockIdx.x;
  
  if (idx < size && ::isnan(array[idx])) {
      *check = true;
  }
 
}



__host__ bool check_array_for_nan(double * ptr, int size, bool on_device, bool * check_d)
{
    // TODO: could probably act on multiple arrays through streams
    if (on_device)
    {
        bool check_h = false;
        cudaMemcpy(check_d, &check_h, sizeof(bool), cudaMemcpyHostToDevice);
        isnan_check_device<<< size/256+1, 256 >>>(ptr, size, check_d);
        cudaMemcpy(&check_h, check_d, sizeof(bool), cudaMemcpyDeviceToHost);

        return check_d;
    }
    else
    {
        bool isnan = false;
        for (int i = 0; i < size; i++)
        {
            isnan |= std::isnan(ptr[i]);
        }
        return isnan;
    }
    
    

}

void check_last_cuda_error(std::string ref_name)
{
    cudaError_t err = cudaGetLastError();
            
    // Check device query
    if (err != cudaSuccess) 
    {
        printf("'%s' cuda error: %s\n", ref_name.c_str(), cudaGetErrorString(err));
    }
}


