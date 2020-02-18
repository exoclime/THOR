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
// 2.0     30/11/2018 Released version (RD & US)
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <iostream>
#include <random>

#include <chrono>
#include <iomanip>
#include <sstream>

#include "reduction_add.h"

using std::abs;
using std::cout;
using std::endl;


template<int BLOCK_SIZE> bool cpu_gpu_test(double *s, long size) {
    bool overall_result = true;

    for (long compute_size = size; compute_size > 0; compute_size /= 2) {
        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
        double reduction_sum_CPU                    = cpu_sum<BLOCK_SIZE>(s, compute_size);
        std::chrono::system_clock::time_point stop  = std::chrono::system_clock::now();
        auto duration_cpu = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

        start             = std::chrono::system_clock::now();
        double output_val = gpu_sum_from_host<BLOCK_SIZE>(s, compute_size);
        stop              = std::chrono::system_clock::now();
        auto duration_gpu = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);


        double output_ref = reduction_sum_CPU;

        bool result = output_val == output_ref;
        overall_result &= result;

        printf("[%ld] [%s] Computed in: CPU: %ld us, GPU: %ld us, CPU/GPU ratio: %f\n",
               compute_size,
               result ? "SUCCESS" : "FAIL",
               duration_cpu.count(),
               duration_gpu.count(),
               double(duration_cpu.count()) / double(duration_gpu.count()));

        if (!result) {

            printf("CPU reduction sum: %32.15f\n", reduction_sum_CPU);
            printf("GPU reduction sum: %32.15f\n", output_val);
        }
    }

    return overall_result;
}


int main() {
    //    long size = 500000000;
    long size = 1000000;
    //    int size = 434567890;

    // allocate on heap
    double *s = new double[size];

    std::random_device rd;        //Will be used to obtain a seed for the random number engine
    std::mt19937       gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);


    // double lin_sum = 0.0;
    // printf("Generating test data\n");
    //
    // for (int i = 0; i < size; i++) {
    //     s[i] = dis(gen);
    //     //s[i] = 1.0;
    //
    //     lin_sum += s[i];
    // }
    //
    //
    // printf("Computing sum over %ld elements\n", size);
    //
    // printf("Linear sum: %32.15f\n", lin_sum);
    //
    // bool overall_result = true;
    //
    // printf("\n");
    // printf("Test BLOCK_SIZE = 512\n");
    // printf("\n");
    // overall_result &= cpu_gpu_test<512>(s, size);
    //
    // printf("\n");
    // printf("Test BLOCK_SIZE = 1024\n");
    // printf("\n");
    // overall_result &= cpu_gpu_test<1024>(s, size);
    //
    // //bool result = abs(output_val - output_ref) < epsilon;
    //
    //
    // if (overall_result)
    //     cout << "reduce sum compare SUCCESS" << endl;
    // else {
    //     cout << "reduce sum compare FAIL" << endl;
    // }
    //
    //
    // delete[] s;


    // testing for sponge layer optimization
    int     nv = 32, nlat_bins = 20, max_count = 100;
    double *vbar_d, *vbar_h_cpu, *vbar_h_gpu, *vbar_h_gpu2, *utmp, *vtmp, *wtmp, *utmp_h, *vtmp_h,
        *wtmp_h;
    vbar_h_cpu  = (double *)malloc(3 * nv * nlat_bins * sizeof(double));
    vbar_h_gpu  = (double *)malloc(3 * nv * nlat_bins * sizeof(double));
    vbar_h_gpu2 = (double *)malloc(3 * nv * nlat_bins * sizeof(double));
    cudaMalloc((void **)&utmp, nv * nlat_bins * max_count * sizeof(double));
    cudaMalloc((void **)&vtmp, nv * nlat_bins * max_count * sizeof(double));
    cudaMalloc((void **)&wtmp, nv * nlat_bins * max_count * sizeof(double));
    cudaMalloc((void **)&vbar_d, nv * nlat_bins * 3 * sizeof(double));
    utmp_h = (double *)malloc(nv * nlat_bins * max_count * sizeof(double));
    vtmp_h = (double *)malloc(nv * nlat_bins * max_count * sizeof(double));
    wtmp_h = (double *)malloc(nv * nlat_bins * max_count * sizeof(double));

    for (int ilat = 0; ilat < nlat_bins; ilat++) {
        for (int lev = 0; lev < nv; lev++) {
            vbar_h_cpu[ilat * nv * 3 + lev * 3 + 0] = 0.0;
            vbar_h_cpu[ilat * nv * 3 + lev * 3 + 1] = 0.0;
            vbar_h_cpu[ilat * nv * 3 + lev * 3 + 2] = 0.0;

            vbar_h_gpu[ilat * nv * 3 + lev * 3 + 0] = 0.0;
            vbar_h_gpu[ilat * nv * 3 + lev * 3 + 1] = 0.0;
            vbar_h_gpu[ilat * nv * 3 + lev * 3 + 2] = 0.0;

            vbar_h_gpu2[ilat * nv * 3 + lev * 3 + 0] = 0.0;
            vbar_h_gpu2[ilat * nv * 3 + lev * 3 + 1] = 0.0;
            vbar_h_gpu2[ilat * nv * 3 + lev * 3 + 2] = 0.0;

            for (int count = 0; count < max_count; count++) {
                utmp_h[ilat * nv * max_count + lev * max_count + count] = dis(gen);
                vbar_h_cpu[ilat * nv * 3 + lev * 3 + 0] +=
                    utmp_h[ilat * nv * max_count + lev * max_count + count];

                vtmp_h[ilat * nv * max_count + lev * max_count + count] = dis(gen);
                vbar_h_cpu[ilat * nv * 3 + lev * 3 + 1] +=
                    vtmp_h[ilat * nv * max_count + lev * max_count + count];

                wtmp_h[ilat * nv * max_count + lev * max_count + count] = dis(gen);
                vbar_h_cpu[ilat * nv * 3 + lev * 3 + 2] +=
                    wtmp_h[ilat * nv * max_count + lev * max_count + count];
            }
        }
    }

    cudaMemcpy(utmp, utmp_h, nv * nlat_bins * max_count * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(vtmp, vtmp_h, nv * nlat_bins * max_count * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(wtmp, wtmp_h, nv * nlat_bins * max_count * sizeof(double), cudaMemcpyHostToDevice);

    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
    for (int ilat = 0; ilat < nlat_bins; ilat++) {
        for (int lev = 0; lev < nv; lev++) {
            vbar_h_gpu[ilat * nv * 3 + lev * 3 + 0] = gpu_sum_on_device<1024>(
                &(utmp[ilat * nv * max_count + lev * max_count]), max_count);
            vbar_h_gpu[ilat * nv * 3 + lev * 3 + 1] = gpu_sum_on_device<1024>(
                &(vtmp[ilat * nv * max_count + lev * max_count]), max_count);
            vbar_h_gpu[ilat * nv * 3 + lev * 3 + 2] = gpu_sum_on_device<1024>(
                &(wtmp[ilat * nv * max_count + lev * max_count]), max_count);
        }
    }
    cudaMemcpy(vbar_d, vbar_h_gpu, 3 * nlat_bins * nv * sizeof(double), cudaMemcpyHostToDevice);
    std::chrono::system_clock::time_point stop = std::chrono::system_clock::now();
    auto duration_one = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    printf("old method: %ld us\n", duration_one.count());

    //do check for accuracy
    double epsilon = 1e-14;
    bool   sum_ok  = true;
    for (int ilat = 0; ilat < nlat_bins; ilat++) {
        for (int lev = 0; lev < nv; lev++) {
            for (int xyz = 0; xyz < 3; xyz++) {
                if (fabs(vbar_h_cpu[ilat * nv * 3 + lev * 3 + xyz]
                         - vbar_h_gpu[ilat * nv * 3 + lev * 3 + xyz])
                        / vbar_h_cpu[ilat * nv * 3 + lev * 3 + xyz]
                    > epsilon) {
                    sum_ok = false;
                    printf("failure at ilat = %d, lev = %d, xyz = %d\n", ilat, lev, xyz);
                    printf("rel error = %.16e\n",
                           fabs(vbar_h_cpu[ilat * nv * 3 + lev * 3 + xyz]
                                - vbar_h_gpu[ilat * nv * 3 + lev * 3 + xyz])
                               / vbar_h_cpu[ilat * nv * 3 + lev * 3 + xyz]);
                }
            }
            // printf("rel error = %.16e\n",
            //        fabs(vbar_h_cpu[ilat * nv * 3 + lev * 3 + 2]
            //             - vbar_h_gpu[ilat * nv * 3 + lev * 3 + 2])
            //            / vbar_h_cpu[ilat * nv * 3 + lev * 3 + 2]);
        }
    }
    if (sum_ok) {
        printf("old method sum ok\n");
    }
    else {
        printf("old method sum NOT ok\n");
    }

    std::chrono::system_clock::time_point start2 = std::chrono::system_clock::now();
    gpu_sum_on_device_sponge<1024>(utmp, max_count, vbar_h_gpu2, nv, nlat_bins, 0);
    gpu_sum_on_device_sponge<1024>(vtmp, max_count, vbar_h_gpu2, nv, nlat_bins, 1);
    gpu_sum_on_device_sponge<1024>(wtmp, max_count, vbar_h_gpu2, nv, nlat_bins, 2);
    cudaMemcpy(vbar_d, vbar_h_gpu2, 3 * nlat_bins * nv * sizeof(double), cudaMemcpyHostToDevice);
    std::chrono::system_clock::time_point stop2 = std::chrono::system_clock::now();
    auto duration_two = std::chrono::duration_cast<std::chrono::microseconds>(stop2 - start2);

    printf("new method: %ld us\n", duration_two.count());


    sum_ok = true;
    for (int ilat = 0; ilat < nlat_bins; ilat++) {
        for (int lev = 0; lev < nv; lev++) {
            for (int xyz = 0; xyz < 3; xyz++) {
                if (fabs(vbar_h_cpu[ilat * nv * 3 + lev * 3 + xyz]
                         - vbar_h_gpu2[ilat * nv * 3 + lev * 3 + xyz])
                        / vbar_h_cpu[ilat * nv * 3 + lev * 3 + xyz]
                    > epsilon) {
                    sum_ok = false;
                    printf("failure at ilat = %d, lev = %d, xyz = %d\n", ilat, lev, xyz);
                    printf("rel error = %.16e\n",
                           fabs(vbar_h_cpu[ilat * nv * 3 + lev * 3 + xyz]
                                - vbar_h_gpu2[ilat * nv * 3 + lev * 3 + xyz])
                               / vbar_h_cpu[ilat * nv * 3 + lev * 3 + xyz]);
                }
            }
            // printf("rel error = %.16e\n",
            //        fabs(vbar_h_cpu[ilat * nv * 3 + lev * 3 + 2]
            //             - vbar_h_gpu[ilat * nv * 3 + lev * 3 + 2])
            //            / vbar_h_cpu[ilat * nv * 3 + lev * 3 + 2]);
        }
    }
    if (sum_ok) {
        printf("new method sum ok\n");
    }
    else {
        printf("new method sum NOT ok\n");
    }

    exit(0);
}
