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
// Description: Helper functions to load variables to memory with halos
//
//
// Method:
//
// Known limitations: None.
//
// Known issues: None.
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

__device__ __forceinline__ bool compute_mem_idx(int*       maps_d,
                                                const int& nhl,
                                                const int& nhl2,
                                                int&       ig,
                                                int&       igh,
                                                int&       ir,
                                                int&       ir2,
                                                bool&      pent_ind) {
    int x  = threadIdx.x;
    int y  = threadIdx.y;
    int ib = blockIdx.x;
    //    int nv  = gridDim.y;
    //    int lev = blockIdx.y;

    bool load_halo = false;

    ir = (y + 1) * nhl + x + 1; // Region index
    ig = maps_d[ib * nhl2 + ir];

    if (x == 0 && y == 0)
        if (maps_d[ib * nhl2] == -1) pent_ind = true;

    ///////////////////////////////
    //////////// Halo /////////////
    ///////////////////////////////

    // x: 0 halo
    //    if (x == 0) {
    //        ir2 = (y + 1) * nhl + x;
    if (y == 0) {
        ir2       = (x + 1) * nhl;
        load_halo = true;
    }
    // x: nhl halo
    //if (x == nhl - 3){
    //        ir2 = (y + 1) * nhl + x + 2;
    else if (y == 3) {
        ir2       = (x + 1) * nhl + nhl - 3 + 2;
        load_halo = true;
    }
    // y: 0 halo
    //if (y == 0){
    //    ir2 = y * nhl + (x + 1);
    else if (y == 7) {
        ir2       = x + 1;
        load_halo = true;
    }


    // x: 0, y: 0 corner point
    //    if (y == 0 && x == 0) {
    //        ir2 = y * nhl + x;
    else if (x == 4 && y == 4) {
        ir2       = 0;
        load_halo = true;
    }

    // y: nhl halo
    //    if (y == nhl - 3) {
    //        ir2 = (y + 2) * nhl + (x + 1);
    else if (y == 11) {

        ir2       = (nhl - 3 + 2) * nhl + (x + 1);
        load_halo = true;
    }
    // x: nhl, y: nhl corner point
    // if (y == nhl - 3 && x == nhl - 3) {
    // ir2 = (y + 2) * nhl + (x + 2);
    else if (y == 4 && x == 5) {

        ir2       = (nhl - 3 + 2) * nhl + (nhl - 3 + 2);
        load_halo = true;
    }

    if (load_halo)
      igh = maps_d[ib * nhl2 + ir2];
    else
        igh = 0;

    return load_halo;
}
