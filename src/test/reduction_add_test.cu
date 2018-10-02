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
// Description: test for reduction sum 
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

#include <random>
#include <cmath>
#include <iostream>

#include <chrono>
#include <iomanip>
#include <sstream>

#include "reduction_add.h"

using std::cout;
using std::endl;
using std::abs;




template<int BLOCK_SIZE>
bool cpu_gpu_test(double * s, long size)
{
    bool overall_result = true;
    
    for (    long compute_size = size; compute_size > 0; compute_size /= 2)
    {
        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
        double reduction_sum_CPU = cpu_sum<BLOCK_SIZE>(s, compute_size);
        std::chrono::system_clock::time_point stop = std::chrono::system_clock::now();
        auto duration_cpu = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        
        start = std::chrono::system_clock::now();
        double output_val = gpu_sum_from_host<BLOCK_SIZE>(s, compute_size);
        stop = std::chrono::system_clock::now();
        auto duration_gpu = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);       
        
        
        double output_ref = reduction_sum_CPU;
        
        bool result = output_val == output_ref;
        overall_result &= result;
        
        printf("[%ld] [%s] Computed in: CPU: %ld us, GPU: %ld us, CPU/GPU ratio: %f\n",
               compute_size,
               result?"SUCCESS":"FAIL",
               duration_cpu.count(),
               duration_gpu.count(),
               double(duration_cpu.count())/double(duration_gpu.count())
            );
        
        if (!result)
        {
            
            printf("CPU reduction sum: %32.15f\n", reduction_sum_CPU);
            printf("GPU reduction sum: %32.15f\n", output_val);
        }
    }

    return overall_result;
    
}


int main ()
{
//    long size = 500000000;
    long size = 1000000000;
//    int size = 434567890;
    
    // allocate on heap
    double * s =  new double[size];
    
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);

    
    double lin_sum = 0.0;
    printf("Generating test data\n");
    
    for (int i = 0; i < size; i++)
    {
        s[i] = dis(gen);
        //s[i] = 1.0;
        
        lin_sum += s[i];
    }

    
    printf("Computing sum over %ld elements\n", size);
    
    printf("Linear sum: %32.15f\n", lin_sum);

    bool overall_result = true;
    
    printf("\n");
    printf("Test BLOCK_SIZE = 512\n");
    printf("\n");
    overall_result &= cpu_gpu_test<512>(s, size);

    printf("\n");
    printf("Test BLOCK_SIZE = 1024\n");
    printf("\n");
    overall_result &= cpu_gpu_test<1024>(s, size);

    //bool result = abs(output_val - output_ref) < epsilon;
    
    
    if (overall_result)
        cout << "reduce sum compare SUCCESS" << endl;
    else
    {
        cout << "reduce sum compare FAIL" << endl;
    }
    
    
    delete[] s;

    
    exit(0);
}
