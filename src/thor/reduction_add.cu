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



template <int BLOCK_SIZE>
__global__ void reduction_add (double *A_d        ,
                               double * out_d,
                               int len)
{
//    __shared__ double 
/*    
    int x = threadIdx.x;
    int y = threadIdx.y;
    int ib = blockIdx.x;
    int nv = gridDim.y;
    int lev = blockIdx.y;

    int pt1, pt2, pt3, pt4, pt5, pt6;
    int nhl = nl_region + 2;
    int nhl2 = nhl*nhl;

    int ir = (y + 1)*nhl + x + 1;   // Region index
    int iri;

    __shared__ double nflxv_s[3 * NX*NY];
    __shared__ double pressure_s[(NX + 2)*(NY + 2)];

    int pent_ind = 0;
    int ig, ir2, id, twot;
    double vr;
    double alt, rscale;
    double Mx, My, Mz;
    double funcx, funcy, funcz;

    // Load shared memory
    ig = maps_d[ib*nhl2 + ir];
    id = ig;
    if (x == 0 && y == 0) if (maps_d[ib * nhl2] == -1) pent_ind = 1;
    pressure_s[ir] = pressure_d[ig * nv + lev];

    ///////////////////////////////
    //////////// Halo /////////////
    ///////////////////////////////
    if (x == 0) {
        ir2 = (y + 1) * nhl + x;
        ig = maps_d[ib * nhl2 + ir2];
        pressure_s[ir2] = pressure_d[ig * nv + lev];
    }
    if (x == nhl - 3){
        ir2 = (y + 1) * nhl + x + 2;
        ig = maps_d[ib * nhl2 + ir2];
        pressure_s[ir2] = pressure_d[ig * nv + lev];
    }
    if (y == 0){
        twot = 1;
        ir2 = y * nhl + (x + 1);
        if (x == 0) twot = 2;

        for (int k = 0; k < twot; k++){
            if (k == 1) ir2 = y * nhl + x;
            ig = maps_d[ib * nhl2 + ir2];
            if (ig >= 0) pressure_s[ir2] = pressure_d[ig * nv + lev];
            else         pressure_s[ir2] = 0.0;
        }
    }
    if (y == nhl - 3) {
        twot = 1;
        ir2 = (y + 2) * nhl + (x + 1);
        if (x == nhl - 3) twot = 2;
        for (int k = 0; k < twot; k++){
            if (k == 1) ir2 = (y + 2) * nhl + (x + 2);
            ig = maps_d[ib * nhl2 + ir2];
            pressure_s[ir2] = pressure_d[ig * nv + lev];
        }
    }
    __syncthreads();
*/
}

