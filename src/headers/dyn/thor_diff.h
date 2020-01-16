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
// Description: 4th order hyper-diffusion (see equations 43 to 48 from Mendonca et al. 2016)
//
//
// Method: Laplacian: first computes the gradient in the primary grid
//        (triangles), which places the calculated vectors at the corners
//         of the control volume avoiding an extra interpolation step.
//         Then using these vectors we apply Eq. 6 to the control volume
//         to obtain the horizontal laplacian.
//
// Known limitations: None
//
//
// Known issues: None
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
__global__ void Diffusion_Op(double* diffmh_d,
                             double* diffw_d,
                             double* diffrh_d,
                             double* diffpr_d,
                             double* diff_d,
                             double* Mh_d,
                             double* Rho_d,
                             double* temperature_d,
                             double* W_d,
                             double* areasTr_d,
                             double* nvecoa_d,
                             double* nvecti_d,
                             double* nvecte_d,
                             double* func_r_d,
                             double* K_d,
                             double* Altitude_d,
                             double  A,
                             double* Rd_d,
                             int*    maps_d,
                             int     nl_region,
                             bool    laststep,
                             bool    DeepModel) {

    int x = threadIdx.x;
    int y = threadIdx.y;
    //int ib  = blockIdx.x;
    int nv  = gridDim.y;
    int lev = blockIdx.y;
    int var = blockIdx.z;

    int nhl  = nl_region + 2;
    int nhl2 = nhl * nhl;
    int pt1, pt2, pt3;

    double alt;
    double rscale;
    double sdiff, vdiff;
    // double funcx, funcy, funcz;
    double AT1, AT2;
    double o3 = 1.0 / 3.0;
    double o6 = 1.0 / 6.0;
    double lap, lap1, lap2, lap3, lap4, lap5, lap6;
    double lapx1, lapx2, lapy1, lapy2, lapz1, lapz2;
    // double dmhz, dmhr;


    /////////////////////////////////////////
    __shared__ double a_s[(NX + 2) * (NY + 2)];
    __shared__ double Rho_s[(NX + 2) * (NY + 2)];
    /////////////////////////////////////////

    int ir = 0;

    int ir2, id, jp1, jp2;


    bool pent_ind = false; //
    int  ig;               // index in global mem

    int igh = 0; // global index in halo

    // Load shared memory


    bool load_halo = compute_mem_idx(maps_d, nhl, nhl2, ig, igh, ir, ir2, pent_ind);
    id             = ig;


    int n_faces = (pent_ind ? 5 : 6);

    Rho_s[ir] = Rho_d[ig * nv + lev];

    ///////////////////////////////
    //////////// Halo /////////////
    ///////////////////////////////
    if (load_halo) {
        if (igh >= 0) {
            Rho_s[ir2] = Rho_d[igh * nv + lev];
        }
        else {


            Rho_s[ir2] = 0.0;
        }
    }
    __syncthreads();
    //////////////////////////////////////////////
    // These set the arguments that the Laplacian will operate on
    if (laststep)
        // second time thru, arg is previous result
        a_s[ir] = diff_d[id * nv * 6 + lev * 6 + var];
    else {
        // first time thru, arg is fluid property
        if (var == 0)
            a_s[ir] = Rho_d[id * nv + lev];
        else if (var == 1)
            a_s[ir] = Mh_d[id * nv * 3 + lev * 3 + 0];
        else if (var == 2)
            a_s[ir] = Mh_d[id * nv * 3 + lev * 3 + 1];
        else if (var == 3)
            a_s[ir] = Mh_d[id * nv * 3 + lev * 3 + 2];
        else if (var == 4)
            a_s[ir] = W_d[id * nv + lev];
        else if (var == 5)
            a_s[ir] = temperature_d[id * nv + lev];
    }
    // we want the velocities, not the momentum
    if (var >= 1 && var <= 4 && !laststep)
        a_s[ir] = a_s[ir] / Rho_s[ir];

    ///////////////////////////////
    //////////// Halo /////////////
    ///////////////////////////////
    if (load_halo) {
        if (igh >= 0) {
            if (laststep)
                a_s[ir2] = diff_d[igh * nv * 6 + lev * 6 + var];
            else {
                if (var == 0)
                    a_s[ir2] = Rho_d[igh * nv + lev];
                else if (var == 1)
                    a_s[ir2] = Mh_d[igh * nv * 3 + lev * 3 + 0];
                else if (var == 2)
                    a_s[ir2] = Mh_d[igh * nv * 3 + lev * 3 + 1];
                else if (var == 3)
                    a_s[ir2] = Mh_d[igh * nv * 3 + lev * 3 + 2];
                else if (var == 4)
                    a_s[ir2] = W_d[igh * nv + lev];
                else if (var == 5)
                    a_s[ir2] = temperature_d[igh * nv + lev];
            }
            if (var >= 1 && var <= 4 && !laststep)
                a_s[ir2] = a_s[ir2] / Rho_s[ir2];
        }
        else
            a_s[ir2] = 0.0;
    }


    __syncthreads();
    //////////////////////////////////////////////

    if (DeepModel) {
        alt    = Altitude_d[lev];
        rscale = A / (alt + A);
    }
    else
        rscale = 1.0;

    if (laststep) {
        if (var < 5)
            sdiff = -K_d[lev];
        else
            sdiff = -K_d[lev] * Rd_d[id * nv + lev]; // multiply by gas constant in temperature eqn
    }

    lap = 0.0;

    for (int j = 0; j < 6; j++) {

        jp1 = (j + 1) % n_faces;
        jp2 = (j + 2) % n_faces;

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

        if (laststep) {
            if (var == 0) {
                vdiff = 0.5 * rscale * sdiff;
            }
            else {
                vdiff = 0.5 * rscale
                        * (2.0 * Rho_s[ir] + Rho_s[pt1] + 2.0 * Rho_s[pt2] + Rho_s[pt3]) * o6
                        * sdiff;
            }
        }
        else
            vdiff = 0.5 * rscale;

        if (j == 0) {
            // lapi are the values 'a_s' the gradient will operate on and INCLUDE the
            // correction terms from the gradient calculation, where the
            // value at the center of the triangle is interpolated as 1/3*(a_s1+a_s2+a_s3)
            AT1  = rscale / areasTr_d[id * 6 + j];
            lap1 = o6 * (a_s[ir] + a_s[pt1]) - o3 * a_s[pt2];
            lap2 = o6 * (a_s[pt1] + a_s[pt2]) - o3 * a_s[ir];
            lap3 = o6 * (a_s[ir] + a_s[pt2]) - o3 * a_s[pt1];
        }
        AT2  = rscale / areasTr_d[id * 6 + jp1];
        lap4 = o6 * (a_s[ir] + a_s[pt2]) - o3 * a_s[pt3];
        lap5 = o6 * (a_s[pt2] + a_s[pt3]) - o3 * a_s[ir];
        lap6 = o6 * (a_s[ir] + a_s[pt3]) - o3 * a_s[pt2];

        if (j == 0) {
            // gradient calculation around triangle centered at two corners of control volume
            lapx1 =
                (-lap1 * nvecti_d[id * 6 * 3 + j * 3 + 0] + lap2 * nvecte_d[id * 6 * 3 + j * 3 + 0]
                 + lap3 * nvecti_d[id * 6 * 3 + jp1 * 3 + 0])
                * AT1;

            lapx2 = (-lap4 * nvecti_d[id * 6 * 3 + jp1 * 3 + 0]
                     + lap5 * nvecte_d[id * 6 * 3 + jp1 * 3 + 0]
                     + lap6 * nvecti_d[id * 6 * 3 + jp2 * 3 + 0])
                    * AT2;

            lapy1 =
                (-lap1 * nvecti_d[id * 6 * 3 + j * 3 + 1] + lap2 * nvecte_d[id * 6 * 3 + j * 3 + 1]
                 + lap3 * nvecti_d[id * 6 * 3 + jp1 * 3 + 1])
                * AT1;

            lapy2 = (-lap4 * nvecti_d[id * 6 * 3 + jp1 * 3 + 1]
                     + lap5 * nvecte_d[id * 6 * 3 + jp1 * 3 + 1]
                     + lap6 * nvecti_d[id * 6 * 3 + jp2 * 3 + 1])
                    * AT2;

            lapz1 =
                (-lap1 * nvecti_d[id * 6 * 3 + j * 3 + 2] + lap2 * nvecte_d[id * 6 * 3 + j * 3 + 2]
                 + lap3 * nvecti_d[id * 6 * 3 + jp1 * 3 + 2])
                * AT1;

            lapz2 = (-lap4 * nvecti_d[id * 6 * 3 + jp1 * 3 + 2]
                     + lap5 * nvecte_d[id * 6 * 3 + jp1 * 3 + 2]
                     + lap6 * nvecti_d[id * 6 * 3 + jp2 * 3 + 2])
                    * AT2;
        }
        else {
            // gradient calculation at next two corners of control volume
            lapx1 = lapx2;

            lapx2 = (-lap4 * nvecti_d[id * 6 * 3 + jp1 * 3 + 0]
                     + lap5 * nvecte_d[id * 6 * 3 + jp1 * 3 + 0]
                     + lap6 * nvecti_d[id * 6 * 3 + jp2 * 3 + 0])
                    * AT2;

            lapy1 = lapy2;

            lapy2 = (-lap4 * nvecti_d[id * 6 * 3 + jp1 * 3 + 1]
                     + lap5 * nvecte_d[id * 6 * 3 + jp1 * 3 + 1]
                     + lap6 * nvecti_d[id * 6 * 3 + jp2 * 3 + 1])
                    * AT2;

            lapz1 = lapz2;

            lapz2 = (-lap4 * nvecti_d[id * 6 * 3 + jp1 * 3 + 2]
                     + lap5 * nvecte_d[id * 6 * 3 + jp1 * 3 + 2]
                     + lap6 * nvecti_d[id * 6 * 3 + jp2 * 3 + 2])
                    * AT2;
        }

        // divergence of gradient (sum over j)
        lap += ((lapx1 + lapx2) * nvecoa_d[id * 6 * 3 + j * 3 + 0]
                + (lapy1 + lapy2) * nvecoa_d[id * 6 * 3 + j * 3 + 1]
                + (lapz1 + lapz2) * nvecoa_d[id * 6 * 3 + j * 3 + 2])
               * vdiff;

        if (pent_ind && j == 4)
            break;
    }

    if (laststep) {
        if (var == 0)
            diffrh_d[id * nv + lev] = lap;
        if (var == 1)
            diffmh_d[id * nv * 3 + lev * 3 + 0] = lap;
        if (var == 2)
            diffmh_d[id * nv * 3 + lev * 3 + 1] = lap;
        if (var == 3) { //zero out radial component
            // funcx = func_r_d[id * 3 + 0];
            // funcy = func_r_d[id * 3 + 1];
            // funcz = func_r_d[id * 3 + 2];
            // dmhz  = lap;
            // //does diffmh_d always get updated before this step?? need to verify!!
            // dmhr = funcx * diffmh_d[id * nv * 3 + lev * 3 + 0]
            //        + funcy * diffmh_d[id * nv * 3 + lev * 3 + 1] + funcz * dmhz;
            // diffmh_d[id * nv * 3 + lev * 3 + 0] += -funcx * dmhr;
            // diffmh_d[id * nv * 3 + lev * 3 + 1] += -funcy * dmhr;
            // diffmh_d[id * nv * 3 + lev * 3 + 2] = dmhz - funcz * dmhr;
            diffmh_d[id * nv * 3 + lev * 3 + 2] = lap;
        }
        if (var == 4)
            diffw_d[id * nv + lev] = lap;
        if (var == 5)
            diffpr_d[id * nv + lev] = lap;
    }
    else
        diff_d[id * nv * 6 + lev * 6 + var] = lap;
}

template<int NN>
__global__ void Diffusion_Op_Poles(double* diffmh_d,
                                   double* diffw_d,
                                   double* diffrh_d,
                                   double* diffpr_d,
                                   double* diff_d,
                                   double* Mh_d,
                                   double* Rho_d,
                                   double* temperature_d,
                                   double* W_d,
                                   double* func_r_d,
                                   double* areasTr_d,
                                   double* nvecoa_d,
                                   double* nvecti_d,
                                   double* nvecte_d,
                                   double* K_d,
                                   double* Altitude_d,
                                   double* Altitudeh_d,
                                   double  A,
                                   double* Rd_d,
                                   int*    local_d,
                                   int     num,
                                   bool    laststep,
                                   bool    DeepModel) {

    int id = blockIdx.x * blockDim.x + threadIdx.x;
    id += num - 2; // Poles

    int nv  = gridDim.y;
    int lev = blockIdx.y;
    int var = blockIdx.z;

    double alt, sdiff, vdiff;
    double rscale;
    int    jp1, jp2;
    int    kp1, kp2;
    double o3 = 1.0 / 3.0;
    double o6 = 1.0 / 6.0;
    double AT1, AT2;
    // double dmhz, dmhr;
    // double funcx, funcy, funcz;
    double lap, lap1, lap2, lap3, lap4, lap5, lap6;
    double lapx1, lapx2, lapy1, lapy2, lapz1, lapz2;
    /////////////////////////////////////////
    __shared__ double a_p[NN + 1];
    __shared__ double Rho_p[NN + 1];

    __shared__ double areasTr_p[NN];
    __shared__ double nvecoa_p[NN * 3];
    __shared__ double nvecti_p[NN * 3];
    __shared__ double nvecte_p[NN * 3];
    __shared__ int    local_p[NN];
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


    alt = Altitude_d[lev];

    Rho_p[0] = Rho_d[id * nv + lev];
    for (int i = 1; i < 6; i++)
        Rho_p[i] = Rho_d[local_p[i - 1] * nv + lev];

    if (laststep) {
        if (var < 5)
            sdiff = -K_d[lev];
        else
            sdiff = -K_d[lev] * Rd_d[id * nv + lev];
    }

    if (laststep) {
        a_p[0] = diff_d[id * nv * 6 + lev * 6 + var];
        for (int i = 1; i < 6; i++)
            a_p[i] = diff_d[local_p[i - 1] * nv * 6 + lev * 6 + var];
    }
    else {
        if (var == 0) {
            a_p[0] = Rho_p[0];
            for (int i = 1; i < 6; i++)
                a_p[i] = Rho_p[i];
        }
        if (var == 1) {
            a_p[0] = Mh_d[id * 3 * nv + lev * 3 + 0] / Rho_p[0];
            for (int i = 1; i < 6; i++)
                a_p[i] = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 0] / Rho_p[i];
        }
        if (var == 2) {
            a_p[0] = Mh_d[id * 3 * nv + lev * 3 + 1] / Rho_p[0];
            for (int i = 1; i < 6; i++)
                a_p[i] = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 1] / Rho_p[i];
        }
        if (var == 3) {
            a_p[0] = Mh_d[id * 3 * nv + lev * 3 + 2] / Rho_p[0];
            for (int i = 1; i < 6; i++)
                a_p[i] = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 2] / Rho_p[i];
        }
        if (var == 4) {
            a_p[0] = W_d[id * nv + lev] / Rho_p[0];
            for (int i = 1; i < 6; i++)
                a_p[i] = W_d[local_p[i - 1] * nv + lev] / Rho_p[i];
        }
        if (var == 5) {
            a_p[0] = temperature_d[id * nv + lev];
            for (int i = 1; i < 6; i++)
                a_p[i] = temperature_d[local_p[i - 1] * nv + lev];
        }


        if (DeepModel)
            rscale = A / (alt + A);
        else
            rscale = 1.0;

        lap = 0.0;

        for (int k = 0; k < 5; k++) {
            int j = k + 1;
            jp1   = (j + 1) % 5;
            jp2   = (j + 2) % 5;
            kp1   = (k + 1) % 5;
            kp2   = (k + 2) % 5;

            if (laststep) {
                if (var == 0) {
                    vdiff = 0.5 * sdiff * rscale;
                }
                else {
                    vdiff = 0.5 * sdiff * rscale
                            * (2.0 * Rho_p[0] + Rho_p[j] + 2.0 * Rho_p[jp1] + Rho_p[jp2]) * o6;
                }
            }
            else
                vdiff = 0.5 * rscale;

            if (k == 0) {
                lap1 = (o6 * (a_p[0] + a_p[j]) - o3 * a_p[jp1]);
                lap2 = (o6 * (a_p[j] + a_p[jp1]) - o3 * a_p[0]);
                lap3 = (o6 * (a_p[0] + a_p[jp1]) - o3 * a_p[j]);
                AT1  = rscale / areasTr_p[k];
            }
            AT2  = rscale / areasTr_p[kp1];
            lap4 = (o6 * (a_p[0] + a_p[jp1]) - o3 * a_p[jp2]);
            lap5 = (o6 * (a_p[jp1] + a_p[jp2]) - o3 * a_p[0]);
            lap6 = (o6 * (a_p[0] + a_p[jp2]) - o3 * a_p[jp1]);

            if (k == 0) {
                lapx1 = (-lap1 * nvecti_p[k * 3 + 0] + lap2 * nvecte_p[k * 3 + 0]
                         + lap3 * nvecti_p[kp1 * 3 + 0])
                        * AT1;
                lapx2 = (-lap4 * nvecti_p[kp1 * 3 + 0] + lap5 * nvecte_p[kp1 * 3 + 0]
                         + lap6 * nvecti_p[kp2 * 3 + 0])
                        * AT2;

                lapy1 = (-lap1 * nvecti_p[k * 3 + 1] + lap2 * nvecte_p[k * 3 + 1]
                         + lap3 * nvecti_p[kp1 * 3 + 1])
                        * AT1;
                lapy2 = (-lap4 * nvecti_p[kp1 * 3 + 1] + lap5 * nvecte_p[kp1 * 3 + 1]
                         + lap6 * nvecti_p[kp2 * 3 + 1])
                        * AT2;

                lapz1 = (-lap1 * nvecti_p[k * 3 + 2] + lap2 * nvecte_p[k * 3 + 2]
                         + lap3 * nvecti_p[kp1 * 3 + 2])
                        * AT1;
                lapz2 = (-lap4 * nvecti_p[kp1 * 3 + 2] + lap5 * nvecte_p[kp1 * 3 + 2]
                         + lap6 * nvecti_p[kp2 * 3 + 2])
                        * AT2;
            }
            else {
                lapx1 = lapx2;
                lapx2 = (-lap4 * nvecti_p[kp1 * 3 + 0] + lap5 * nvecte_p[kp1 * 3 + 0]
                         + lap6 * nvecti_p[kp2 * 3 + 0])
                        * AT2;

                lapy1 = lapy2;
                lapy2 = (-lap4 * nvecti_p[kp1 * 3 + 1] + lap5 * nvecte_p[kp1 * 3 + 1]
                         + lap6 * nvecti_p[kp2 * 3 + 1])
                        * AT2;

                lapz1 = lapz2;
                lapz2 = (-lap4 * nvecti_p[kp1 * 3 + 2] + lap5 * nvecte_p[kp1 * 3 + 2]
                         + lap6 * nvecti_p[kp2 * 3 + 2])
                        * AT2;
            }
            lap += ((lapx1 + lapx2) * nvecoa_p[k * 3 + 2] + (lapy1 + lapy2) * nvecoa_p[k * 3 + 1]
                    + (lapz1 + lapz2) * nvecoa_p[k * 3 + 2])
                   * vdiff;
        }
        if (laststep) {
            if (var == 0)
                diffrh_d[id * nv + lev] = lap;
            if (var == 1)
                diffmh_d[id * 3 * nv + lev * 3 + 0] = lap;
            if (var == 2)
                diffmh_d[id * 3 * nv + lev * 3 + 1] = lap;
            if (var == 3) {
                // funcx = func_r_d[id * 3 + 0];
                // funcy = func_r_d[id * 3 + 1];
                // funcz = func_r_d[id * 3 + 2];
                // dmhz  = lap;
                // dmhr  = funcx * diffmh_d[id * nv * 3 + lev * 3 + 0]
                //        + funcy * diffmh_d[id * nv * 3 + lev * 3 + 1] + funcz * dmhz;
                // diffmh_d[id * nv * 3 + lev * 3 + 0] += -funcx * dmhr;
                // diffmh_d[id * nv * 3 + lev * 3 + 1] += -funcy * dmhr;
                // diffmh_d[id * nv * 3 + lev * 3 + 2] = dmhz - funcz * dmhr;
                diffmh_d[id * 3 * nv + lev * 3 + 2] = lap;
            }
            if (var == 4)
                diffw_d[id * nv + lev] = lap;
            if (var == 5)
                diffpr_d[id * nv + lev] = lap;
        }
        else
            diff_d[id * nv * 6 + lev * 6 + var] = lap;
    }
}

__global__ void Correct_Horizontal(double* diffmh_d, double* diffmv_d, double* func_r_d, int num) {
    //this function removes any spurious vertical component in horizontal momentum diffusion
    int id  = blockIdx.x * blockDim.x + threadIdx.x;
    int nv  = gridDim.y;
    int lev = blockIdx.y;

    if (id < num) {
        double dmhr, dmvr;
        dmhr = func_r_d[id * 3 + 0] * diffmh_d[id * nv * 3 + lev * 3 + 0]
               + func_r_d[id * 3 + 1] * diffmh_d[id * nv * 3 + lev * 3 + 1]
               + func_r_d[id * 3 + 2] * diffmh_d[id * nv * 3 + lev * 3 + 2];

        diffmh_d[id * nv * 3 + lev * 3 + 0] += -func_r_d[id * 3 + 0] * dmhr;
        diffmh_d[id * nv * 3 + lev * 3 + 1] += -func_r_d[id * 3 + 1] * dmhr;
        diffmh_d[id * nv * 3 + lev * 3 + 2] += -func_r_d[id * 3 + 2] * dmhr;

        dmvr = func_r_d[id * 3 + 0] * diffmv_d[id * nv * 3 + lev * 3 + 0]
               + func_r_d[id * 3 + 1] * diffmv_d[id * nv * 3 + lev * 3 + 1]
               + func_r_d[id * 3 + 2] * diffmv_d[id * nv * 3 + lev * 3 + 2];

        diffmv_d[id * nv * 3 + lev * 3 + 0] += -func_r_d[id * 3 + 0] * dmvr;
        diffmv_d[id * nv * 3 + lev * 3 + 1] += -func_r_d[id * 3 + 1] * dmvr;
        diffmv_d[id * nv * 3 + lev * 3 + 2] += -func_r_d[id * 3 + 2] * dmvr;
    }
}
