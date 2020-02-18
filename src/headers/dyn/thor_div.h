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
// Description: Computes the 3D divergence damping (see equations 44 to 46 from Mendonca et al. 2016)
//
// Method: The same as for the 4th order hyper-diffusion
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
#include "kernel_halo_helpers.h"


template<int NX, int NY>
__global__ void DivM_Op(double* DivM_d,
                        double* divg_Mh_d,
                        double* Mh_d,
                        double* Wh_d,
                        double* K_d,
                        double* areasTr_d,
                        double* areas_d,
                        double* nvecoa_d,
                        double* nvecti_d,
                        double* nvecte_d,
                        double* func_r_d,
                        double* Altitudeh_d,
                        double* Altitude_d,
                        double  A,
                        int*    maps_d,
                        int     nl_region,
                        bool    laststep,
                        bool    DeepModel) {

    int x = threadIdx.x;
    int y = threadIdx.y;
    //int ib  = blockIdx.x;
    int nv  = gridDim.y;
    int lev = blockIdx.y;

    int nhl  = nl_region + 2;
    int nhl2 = nhl * nhl;
    int pt1, pt2, pt3;

    double alt, dz, dwdz, dwdz2;
    double rscale;
    double AT0, AT1, AT2;
    double a01x, a01y, a01z, a02x, a02y, a02z, a03x, a03y, a03z;
    double a12x, a12y, a12z, a23x, a23y, a23z;
    double meanwl, meanwt, meanwl2, meanwt2;
    // double o3 = 1.0 / 3.0;
    double lapx, lapy, lapz, lapr;
    double lap1, lap2;
    double sdiff;
    double funcx, funcy, funcz;
    int    jp1, jp2;

    /////////////////////////////////////////
    __shared__ double a_s[3 * (NX + 2) * (NY + 2)];
    __shared__ double wl_s[(NX + 2) * (NY + 2)];
    __shared__ double wt_s[(NX + 2) * (NY + 2)];
    /////////////////////////////////////////

    int ir = 0;
    int ir2, id;

    bool pent_ind = false; //
    int  ig;               // index in global mem

    int igh = 0; // global index in halo


    bool load_halo = compute_mem_idx(maps_d, nhl, nhl2, ig, igh, ir, ir2, pent_ind);
    id             = ig;

    if (laststep)
        for (int k = 0; k < 3; k++)
            a_s[ir * 3 + k] = divg_Mh_d[ig * nv * 3 + lev * 3 + k];
    else {
        for (int k = 0; k < 3; k++)
            a_s[ir * 3 + k] = Mh_d[ig * 3 * nv + lev * 3 + k];
        wl_s[ir] = Wh_d[ig * (nv + 1) + lev];
        wt_s[ir] = Wh_d[ig * (nv + 1) + lev + 1];
    }

    ///////////////////////////////
    //////////// Halo /////////////
    ///////////////////////////////
    if (load_halo) {
        if (igh >= 0) {
            if (laststep)
                for (int k = 0; k < 3; k++)
                    a_s[ir2 * 3 + k] = divg_Mh_d[igh * nv * 3 + lev * 3 + k];
            else {
                for (int k = 0; k < 3; k++)
                    a_s[ir2 * 3 + k] = Mh_d[igh * 3 * nv + lev * 3 + k];
                wl_s[ir2] = Wh_d[igh * (nv + 1) + lev];
                wt_s[ir2] = Wh_d[igh * (nv + 1) + lev + 1];
            }
        }
        else {
            if (laststep)
                for (int k = 0; k < 3; k++)
                    a_s[ir2 * 3 + k] = 0.0;
            else {
                for (int k = 0; k < 3; k++)
                    a_s[ir2 * 3 + k] = 0.0;
                wl_s[ir2] = 0.0;
                wt_s[ir2] = 0.0;
            }
        }
    }

    __syncthreads();
    //////////////////////////////////////////////

    if (DeepModel) {
        alt    = Altitude_d[lev];
        rscale = A / (alt + A);
    }
    else
        rscale = 1.0;

    AT0 = 0.5 * rscale;

    lapx = 0.0;
    lapy = 0.0;
    lapz = 0.0;

    if (!laststep)
        dz = 1.0 / (Altitudeh_d[lev] - Altitudeh_d[lev + 1]);

    for (int j = 0; j < 6; j++) {
        jp1 = (pent_ind) * (j == 4 ? 0 : j + 1) + (!pent_ind) * (j == 5 ? 0 : j + 1);
        jp2 = (pent_ind) * (jp1 == 4 ? 0 : jp1 + 1) + (!pent_ind) * (jp1 == 5 ? 0 : jp1 + 1);

        if (j == 0) {
            pt1 = (y + 2) * nhl + x + 1;
            pt2 = (y + 2) * nhl + x + 2;
            pt3 = (y + 1) * nhl + x + 2;
        }
        else if (j == 1) {
            pt1 = (y + 2) * nhl + x + 2;
            pt2 = (y + 1) * nhl + x + 2;
            pt3 = (y)*nhl + x + 1;
        }
        else if (j == 2) {
            pt1 = (y + 1) * nhl + x + 2;
            pt2 = (y)*nhl + x + 1;
            pt3 = pent_ind * ((y + 1) * nhl + x) + !pent_ind * ((y)*nhl + x);
        }
        else if (j == 3) {
            pt1 = (y)*nhl + x + 1;
            pt2 = pent_ind * ((y + 1) * nhl + x) + !pent_ind * ((y)*nhl + x);
            pt3 = pent_ind * ((y + 2) * nhl + x + 1) + !pent_ind * ((y + 1) * nhl + x);
        }
        else if (j == 4) {
            pt1 = pent_ind * ((y + 1) * nhl + x) + !pent_ind * ((y)*nhl + x);
            pt2 = pent_ind * ((y + 2) * nhl + x + 1) + !pent_ind * ((y + 1) * nhl + x);
            pt3 = pent_ind * ((y + 2) * nhl + x + 2) + !pent_ind * ((y + 2) * nhl + x + 1);
        }
        else {
            pt1 = (y + 1) * nhl + x;
            pt2 = (y + 2) * nhl + x + 1;
            pt3 = (y + 2) * nhl + x + 2;
        }

        if (!laststep) {
            // meanwl  = (wl_s[ir] + wl_s[pt1] + wl_s[pt2]) * o3;
            // meanwt  = (wt_s[ir] + wt_s[pt1] + wt_s[pt2]) * o3;
            // meanwl2 = (wl_s[ir] + wl_s[pt2] + wl_s[pt3]) * o3;
            // meanwt2 = (wt_s[ir] + wt_s[pt2] + wt_s[pt3]) * o3;

            // do tri-linear interpolation to get vertical mom. at corners
            meanwl = (areas_d[id * 6 * 3 + j * 3 + 0] * wl_s[ir]
                      + areas_d[id * 6 * 3 + j * 3 + 1] * wl_s[pt1]
                      + areas_d[id * 6 * 3 + j * 3 + 2] * wl_s[pt2])
                     / areasTr_d[id * 6 + j];
            meanwt = (areas_d[id * 6 * 3 + j * 3 + 0] * wt_s[ir]
                      + areas_d[id * 6 * 3 + j * 3 + 1] * wt_s[pt1]
                      + areas_d[id * 6 * 3 + j * 3 + 2] * wt_s[pt2])
                     / areasTr_d[id * 6 + j];
            meanwl2 = (areas_d[id * 6 * 3 + jp1 * 3 + 0] * wl_s[ir]
                       + areas_d[id * 6 * 3 + jp1 * 3 + 1] * wl_s[pt2]
                       + areas_d[id * 6 * 3 + jp1 * 3 + 2] * wl_s[pt3])
                      / areasTr_d[id * 6 + jp1];
            meanwt2 = (areas_d[id * 6 * 3 + jp1 * 3 + 0] * wt_s[ir]
                       + areas_d[id * 6 * 3 + jp1 * 3 + 1] * wt_s[pt2]
                       + areas_d[id * 6 * 3 + jp1 * 3 + 2] * wt_s[pt3])
                      / areasTr_d[id * 6 + jp1];

            dwdz  = (1 - (pent_ind && j == 5)) * (meanwl - meanwt) * dz;
            dwdz2 = (1 - (pent_ind && j == 5)) * (meanwl2 - meanwt2) * dz;
        }
        else {
            dwdz  = 0.0;
            dwdz2 = 0.0;
        }

        if (j == 0) {
            a01x = (a_s[ir * 3 + 0] + a_s[pt1 * 3 + 0]);
            a01y = (a_s[ir * 3 + 1] + a_s[pt1 * 3 + 1]);
            a01z = (a_s[ir * 3 + 2] + a_s[pt1 * 3 + 2]);
            a12x = (a_s[pt1 * 3 + 0] + a_s[pt2 * 3 + 0]);
            a12y = (a_s[pt1 * 3 + 1] + a_s[pt2 * 3 + 1]);
            a12z = (a_s[pt1 * 3 + 2] + a_s[pt2 * 3 + 2]);
            AT1  = 0.5 * rscale / areasTr_d[id * 6 + j];
        }


        AT2  = 0.5 * rscale / areasTr_d[id * 6 + jp1];
        a02x = (a_s[ir * 3 + 0] + a_s[pt2 * 3 + 0]);
        a02y = (a_s[ir * 3 + 1] + a_s[pt2 * 3 + 1]);
        a02z = (a_s[ir * 3 + 2] + a_s[pt2 * 3 + 2]);
        a23x = (a_s[pt2 * 3 + 0] + a_s[pt3 * 3 + 0]);
        a23y = (a_s[pt2 * 3 + 1] + a_s[pt3 * 3 + 1]);
        a23z = (a_s[pt2 * 3 + 2] + a_s[pt3 * 3 + 2]);
        a03x = (a_s[ir * 3 + 0] + a_s[pt3 * 3 + 0]);
        a03y = (a_s[ir * 3 + 1] + a_s[pt3 * 3 + 1]);
        a03z = (a_s[ir * 3 + 2] + a_s[pt3 * 3 + 2]);

        if (j == 0) {
            lap1 =
                (-a01x * nvecti_d[id * 6 * 3 + j * 3 + 0] - a01y * nvecti_d[id * 6 * 3 + j * 3 + 1]
                 - a01z * nvecti_d[id * 6 * 3 + j * 3 + 2] + a12x * nvecte_d[id * 6 * 3 + j * 3 + 0]
                 + a12y * nvecte_d[id * 6 * 3 + j * 3 + 1] + a12z * nvecte_d[id * 6 * 3 + j * 3 + 2]
                 + a02x * nvecti_d[id * 6 * 3 + jp1 * 3 + 0]
                 + a02y * nvecti_d[id * 6 * 3 + jp1 * 3 + 1]
                 + a02z * nvecti_d[id * 6 * 3 + jp1 * 3 + 2])
                * AT1;

            lap2 = (-a02x * nvecti_d[id * 6 * 3 + jp1 * 3 + 0]
                    - a02y * nvecti_d[id * 6 * 3 + jp1 * 3 + 1]
                    - a02z * nvecti_d[id * 6 * 3 + jp1 * 3 + 2]
                    + a23x * nvecte_d[id * 6 * 3 + jp1 * 3 + 0]
                    + a23y * nvecte_d[id * 6 * 3 + jp1 * 3 + 1]
                    + a23z * nvecte_d[id * 6 * 3 + jp1 * 3 + 2]
                    + a03x * nvecti_d[id * 6 * 3 + jp2 * 3 + 0]
                    + a03y * nvecti_d[id * 6 * 3 + jp2 * 3 + 1]
                    + a03z * nvecti_d[id * 6 * 3 + jp2 * 3 + 2])
                   * AT2;
        }
        else {
            lap1 = lap2;
            lap2 = (-a02x * nvecti_d[id * 6 * 3 + jp1 * 3 + 0]
                    - a02y * nvecti_d[id * 6 * 3 + jp1 * 3 + 1]
                    - a02z * nvecti_d[id * 6 * 3 + jp1 * 3 + 2]
                    + a23x * nvecte_d[id * 6 * 3 + jp1 * 3 + 0]
                    + a23y * nvecte_d[id * 6 * 3 + jp1 * 3 + 1]
                    + a23z * nvecte_d[id * 6 * 3 + jp1 * 3 + 2]
                    + a03x * nvecti_d[id * 6 * 3 + jp2 * 3 + 0]
                    + a03y * nvecti_d[id * 6 * 3 + jp2 * 3 + 1]
                    + a03z * nvecti_d[id * 6 * 3 + jp2 * 3 + 2])
                   * AT2;
        }

        lapx += (lap1 + dwdz + lap2 + dwdz2) * nvecoa_d[id * 6 * 3 + j * 3 + 0] * AT0;
        lapy += (lap1 + dwdz + lap2 + dwdz2) * nvecoa_d[id * 6 * 3 + j * 3 + 1] * AT0;
        lapz += (lap1 + dwdz + lap2 + dwdz2) * nvecoa_d[id * 6 * 3 + j * 3 + 2] * AT0;

        if (pent_ind && j == 4)
            break;
    }

    if (laststep) {
        sdiff = K_d[lev];

        funcx = func_r_d[id * 3 + 0];
        funcy = func_r_d[id * 3 + 1];
        funcz = func_r_d[id * 3 + 2];

        lapr = funcx * lapx + funcy * lapy + funcz * lapz;

        lapx = lapx - funcx * lapr;
        lapy = lapy - funcy * lapr;
        lapz = lapz - funcz * lapr;

        DivM_d[id * nv * 3 + lev * 3 + 0] = -lapx * sdiff;
        DivM_d[id * nv * 3 + lev * 3 + 1] = -lapy * sdiff;
        DivM_d[id * nv * 3 + lev * 3 + 2] = -lapz * sdiff;
    }
    else {
        divg_Mh_d[id * nv * 3 + lev * 3 + 0] = lapx;
        divg_Mh_d[id * nv * 3 + lev * 3 + 1] = lapy;
        divg_Mh_d[id * nv * 3 + lev * 3 + 2] = lapz;
    }
}

template<int NN>
__global__ void DivM_Op_Poles(double* DivM_d,
                              double* divg_Mh_d,
                              double* Mh_d,
                              double* Wh_d,
                              double* K_d,
                              double* areasTr_d,
                              double* areas_d,
                              double* nvecoa_d,
                              double* nvecti_d,
                              double* nvecte_d,
                              double* func_r_d,
                              double* Altitudeh_d,
                              double* Altitude_d,
                              double  A,
                              int*    local_d,
                              int     num,
                              bool    laststep,
                              bool    DeepModel) {


    int id = blockIdx.x * blockDim.x + threadIdx.x;
    id += num - 2; // Poles

    int nv  = gridDim.y;
    int lev = blockIdx.y;

    double alt, sdiff;
    double dwdz, dwdz2, dz;
    double meanwl, meanwl2;
    double meanwt, meanwt2;
    double rscale;
    // double o3 = 1.0 / 3.0;
    double lapx, lapy, lapz, lapr, lap1, lap2;
    double a01x, a01y, a01z, a02x, a02y, a02z, a03x, a03y, a03z;
    double a12x, a12y, a12z, a23x, a23y, a23z;
    double AT0, AT1, AT2;
    double funcx, funcy, funcz;

    /////////////////////////////////////////
    __shared__ double a_p[(NN + 1) * 3];
    __shared__ double wl_p[NN + 1];
    __shared__ double wt_p[NN + 1];

    __shared__ double areasTr_p[NN];
    __shared__ double nvecoa_p[NN * 3];
    __shared__ double nvecti_p[NN * 3];
    __shared__ double nvecte_p[NN * 3];
    __shared__ int    local_p[NN];
    __shared__ double areas_p[NN * 3];

    /////////////////////////////////////////

    for (int i = 0; i < 5; i++)
        local_p[i] = local_d[id * 6 + i];
    for (int i = 0; i < 5; i++)
        areasTr_p[i] = areasTr_d[id * 6 + i];
    for (int i = 0; i < 5; i++)
        for (int k = 0; k < 3; k++)
            nvecoa_p[i * 3 + k] = nvecoa_d[id * 6 * 3 + i * 3 + k];
    for (int i = 0; i < 5; i++)
        for (int k = 0; k < 3; k++)
            nvecti_p[i * 3 + k] = nvecti_d[id * 6 * 3 + i * 3 + k];
    for (int i = 0; i < 5; i++)
        for (int k = 0; k < 3; k++)
            nvecte_p[i * 3 + k] = nvecte_d[id * 6 * 3 + i * 3 + k];
    for (int i = 0; i < 5; i++)
        for (int k = 0; k < 3; k++)
            areas_p[i * 3 + k] = areas_d[id * 6 * 3 + i * 3 + k];

    if (laststep)
        sdiff = K_d[lev];

    alt = Altitude_d[lev];

    if (laststep) {
        for (int k = 0; k < 3; k++)
            a_p[0 + k] = divg_Mh_d[id * nv * 3 + lev * 3 + k];
        for (int i = 1; i < 6; i++)
            for (int k = 0; k < 3; k++)
                a_p[i * 3 + k] = divg_Mh_d[local_p[i - 1] * nv * 3 + lev * 3 + k];
    }
    else {
        for (int k = 0; k < 3; k++)
            a_p[0 * 3 + k] = Mh_d[id * 3 * nv + lev * 3 + k];
        for (int i = 1; i < 6; i++)
            for (int k = 0; k < 3; k++)
                a_p[i * 3 + k] = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + k];

        wl_p[0] = Wh_d[id * (nv + 1) + lev];
        for (int i = 1; i < 6; i++)
            wl_p[i] = Wh_d[local_p[i - 1] * (nv + 1) + lev];

        wt_p[0] = Wh_d[id * (nv + 1) + lev + 1];
        for (int i = 1; i < 6; i++)
            wt_p[i] = Wh_d[local_p[i - 1] * (nv + 1) + lev + 1];
    }

    if (DeepModel) {
        alt    = Altitude_d[lev];
        rscale = A / (alt + A);
    }
    else
        rscale = 1.0;

    AT0 = 0.5 * rscale;

    lapx = 0.0;
    lapy = 0.0;
    lapz = 0.0;

    if (!laststep)
        dz = 1.0 / (Altitudeh_d[lev] - Altitudeh_d[lev + 1]);

    for (int k = 0; k < 5; k++) {
        int j   = k + 1;
        int jp1 = (j == 5 ? 1 : j + 1);
        int jp2 = (jp1 == 5 ? 1 : jp1 + 1);
        int kp1 = (k == 4 ? 0 : k + 1);
        int kp2 = (kp1 == 4 ? 0 : kp1 + 1);

        if (laststep) {
            dwdz  = 0.0;
            dwdz2 = 0.0;
        }
        else {
            // interpolates vertical velocity to middle of triangle
            // currently uses mean but could be more accurate with tri-linear interp
            // meanwl  = (wl_p[0] + wl_p[j] + wl_p[jp1]) * o3;
            // meanwt  = (wt_p[0] + wt_p[j] + wt_p[jp1]) * o3;
            // meanwl2 = (wl_p[0] + wl_p[j] + wl_p[jp2]) * o3;
            // meanwt2 = (wt_p[0] + wt_p[j] + wt_p[jp2]) * o3;

            // use tri-linear interpolation
            meanwl = (areas_p[k * 3 + 0] * wl_p[0] + areas_p[k * 3 + 1] * wl_p[j]
                      + areas_p[k * 3 + 1] * wl_p[jp1])
                     / areasTr_p[k];
            meanwt = (areas_p[k * 3 + 0] * wt_p[0] + areas_p[k * 3 + 1] * wt_p[j]
                      + areas_p[k * 3 + 1] * wt_p[jp1])
                     / areasTr_p[k];
            meanwl2 = (areas_p[kp1 * 3 + 0] * wl_p[0] + areas_p[kp1 * 3 + 1] * wl_p[jp1]
                       + areas_p[kp1 * 3 + 1] * wl_p[jp2])
                      / areasTr_p[kp1];
            meanwt2 = (areas_p[kp1 * 3 + 0] * wt_p[0] + areas_p[kp1 * 3 + 1] * wt_p[jp1]
                       + areas_p[kp1 * 3 + 1] * wt_p[jp2])
                      / areasTr_p[kp1];

            dwdz  = (meanwl - meanwt) * dz;
            dwdz2 = (meanwl2 - meanwt2) * dz;
        }

        // set up terms of divergence at each corner of two triangles
        if (k == 0) {
            a01x = (a_p[0 * 3 + 0] + a_p[j * 3 + 0]);
            a01y = (a_p[0 * 3 + 1] + a_p[j * 3 + 1]);
            a01z = (a_p[0 * 3 + 2] + a_p[j * 3 + 2]);
            a12x = (a_p[j * 3 + 0] + a_p[jp1 * 3 + 0]);
            a12y = (a_p[j * 3 + 1] + a_p[jp1 * 3 + 1]);
            a12z = (a_p[j * 3 + 2] + a_p[jp1 * 3 + 2]);
            AT1  = 0.5 * rscale / (areasTr_p[k]);
        }

        a02x = (a_p[0 * 3 + 0] + a_p[jp1 * 3 + 0]);
        a02y = (a_p[0 * 3 + 1] + a_p[jp1 * 3 + 1]);
        a02z = (a_p[0 * 3 + 2] + a_p[jp1 * 3 + 2]);
        a23x = (a_p[jp1 * 3 + 0] + a_p[jp2 * 3 + 0]);
        a23y = (a_p[jp1 * 3 + 1] + a_p[jp2 * 3 + 1]);
        a23z = (a_p[jp1 * 3 + 2] + a_p[jp2 * 3 + 2]);
        a03x = (a_p[0 * 3 + 0] + a_p[jp2 * 3 + 0]);
        a03y = (a_p[0 * 3 + 1] + a_p[jp2 * 3 + 1]);
        a03z = (a_p[0 * 3 + 2] + a_p[jp2 * 3 + 2]);
        AT2  = 0.5 * rscale / (areasTr_p[kp1]);

        // calculate divergence of two adjacent triangles
        if (k == 0) {
            lap1 = (-a01x * nvecti_p[k * 3 + 0] - a01y * nvecti_p[k * 3 + 1]
                    - a01z * nvecti_p[k * 3 + 2] + a12x * nvecte_p[k * 3 + 0]
                    + a12y * nvecte_p[k * 3 + 1] + a12z * nvecte_p[k * 3 + 2]
                    + a02x * nvecti_p[kp1 * 3 + 0] + a02y * nvecti_p[kp1 * 3 + 1]
                    + a02z * nvecti_p[kp1 * 3 + 2])
                   * AT1;

            lap2 = (-a02x * nvecti_p[kp1 * 3 + 0] - a02y * nvecti_p[kp1 * 3 + 1]
                    - a02z * nvecti_p[kp1 * 3 + 2] + a23x * nvecte_p[kp1 * 3 + 0]
                    + a23y * nvecte_p[kp1 * 3 + 1] + a23z * nvecte_p[kp1 * 3 + 2]
                    + a03x * nvecti_p[kp2 * 3 + 0] + a03y * nvecti_p[kp2 * 3 + 1]
                    + a03z * nvecti_p[kp2 * 3 + 2])
                   * AT2;
        }
        else {
            lap1 = lap2;

            lap2 = (-a02x * nvecti_p[kp1 * 3 + 0] - a02y * nvecti_p[kp1 * 3 + 1]
                    - a02z * nvecti_p[kp1 * 3 + 2] + a23x * nvecte_p[kp1 * 3 + 0]
                    + a23y * nvecte_p[kp1 * 3 + 1] + a23z * nvecte_p[kp1 * 3 + 2]
                    + a03x * nvecti_p[kp2 * 3 + 0] + a03y * nvecti_p[kp2 * 3 + 1]
                    + a03z * nvecti_p[kp2 * 3 + 2])
                   * AT2;
        }
        // gradient terms
        lapx += (lap1 + dwdz + lap2 + dwdz2) * nvecoa_p[k * 3 + 0] * AT0;
        lapy += (lap1 + dwdz + lap2 + dwdz2) * nvecoa_p[k * 3 + 1] * AT0;
        lapz += (lap1 + dwdz + lap2 + dwdz2) * nvecoa_p[k * 3 + 2] * AT0;
    }

    if (laststep) {
        sdiff = K_d[lev];
        funcx = func_r_d[id * 3 + 0];
        funcy = func_r_d[id * 3 + 1];
        funcz = func_r_d[id * 3 + 2];

        lapr = funcx * lapx + funcy * lapy + funcz * lapz;

        lapx = lapx - funcx * lapr;
        lapy = lapy - funcy * lapr;
        lapz = lapz - funcz * lapr;

        DivM_d[id * nv * 3 + lev * 3 + 0] = -lapx * sdiff;
        DivM_d[id * nv * 3 + lev * 3 + 1] = -lapy * sdiff;
        DivM_d[id * nv * 3 + lev * 3 + 2] = -lapz * sdiff;
    }
    else {
        divg_Mh_d[id * nv * 3 + lev * 3 + 0] = lapx;
        divg_Mh_d[id * nv * 3 + lev * 3 + 1] = lapy;
        divg_Mh_d[id * nv * 3 + lev * 3 + 2] = lapz;
    }
}
