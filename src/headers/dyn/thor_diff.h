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
                             double* areas_d,
                             double* nvecoa_d,
                             double* nvecti_d,
                             double* nvecte_d,
                             double* func_r_d,
                             double* K_d,
                             double* Altitude_d,
                             double  A,
                             double* Rd_d,
                             double* Cp_d,
                             int*    maps_d,
                             int     nl_region,
                             bool    laststep,
                             bool    DeepModel,
                             bool    DiffSponge,
                             int     order_diff_sponge,
                             double* Kdh2_d,
                             double* boundary_flux_d,
                             bool    energy_equation) {

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
    double sdiff, vdiff, astar;
    // double funcx, funcy, funcz;
    double AT1, AT2;
    // double o3 = 1.0 / 3.0;
    // double o6 = 1.0 / 6.0;
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
    if (laststep) {
        // second time thru, arg is previous result
        if (var == 0) {
            a_s[ir] = diff_d[id * nv * 6 + lev * 6 + var];
        }
        else {
            a_s[ir] = diff_d[id * nv * 6 + lev * 6 + var] * Rho_s[ir];
        }
    }
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
        else if (var == 5) {
            if (energy_equation) {
                a_s[ir] =
                    (Cp_d[id * nv + lev] - Rd_d[id * nv + lev]) * temperature_d[id * nv + lev];
            }
            else {
                a_s[ir] = Rd_d[id * nv + lev] * temperature_d[id * nv + lev];
            }
        }
    }
    // we want the velocities, not the momentum
    if (var >= 1 && var <= 4 && !laststep)
        a_s[ir] = a_s[ir] / Rho_s[ir];

    ///////////////////////////////
    //////////// Halo /////////////
    ///////////////////////////////
    if (load_halo) {
        if (igh >= 0) {
            if (laststep) {
                if (var == 0) {
                    a_s[ir2] = diff_d[igh * nv * 6 + lev * 6 + var];
                }
                else {
                    a_s[ir2] = diff_d[igh * nv * 6 + lev * 6 + var] * Rho_s[ir2];
                }
            }
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
                else if (var == 5) {
                    if (energy_equation) {
                        a_s[ir2] = (Cp_d[igh * nv + lev] - Rd_d[igh * nv + lev])
                                   * temperature_d[igh * nv + lev];
                    }
                    else {
                        a_s[ir2] = Rd_d[igh * nv + lev] * temperature_d[igh * nv + lev];
                    }
                }
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
        sdiff = -K_d[lev];
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
            // if (var == 0) {
            vdiff = 0.5 * rscale * sdiff;
            // }
            // else {
            //     vdiff = 0.5 * rscale
            //             * (2.0 * Rho_s[ir] + Rho_s[pt1] + 2.0 * Rho_s[pt2] + Rho_s[pt3]) * o6
            //             * sdiff;
            // }
        }
        else
            vdiff = 0.5 * rscale;

        if (j == 0) {
            // lapi are the values 'a_s' the gradient will operate on and INCLUDE the
            // correction terms from the gradient calculation
            AT1   = rscale / areasTr_d[id * 6 + j];
            astar = (areas_d[id * 6 * 3 + j * 3 + 0] * a_s[ir]
                     + areas_d[id * 6 * 3 + j * 3 + 1] * a_s[pt1]
                     + areas_d[id * 6 * 3 + j * 3 + 2] * a_s[pt2])
                    / areasTr_d[id * 6 + j];
            lap1 = (0.5 * (a_s[pt1] + a_s[ir]) - astar);
            lap2 = (0.5 * (a_s[pt1] + a_s[pt2]) - astar);
            lap3 = (0.5 * (a_s[ir] + a_s[pt2]) - astar);
        }
        AT2   = rscale / areasTr_d[id * 6 + jp1];
        astar = (areas_d[id * 6 * 3 + jp1 * 3 + 0] * a_s[ir]
                 + areas_d[id * 6 * 3 + jp1 * 3 + 1] * a_s[pt2]
                 + areas_d[id * 6 * 3 + jp1 * 3 + 2] * a_s[pt3])
                / areasTr_d[id * 6 + jp1];
        lap4 = (0.5 * (a_s[ir] + a_s[pt2]) - astar);
        lap5 = (0.5 * (a_s[pt2] + a_s[pt3]) - astar);
        lap6 = (0.5 * (a_s[ir] + a_s[pt3]) - astar);

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
        if (var == 0) {
            diffrh_d[id * nv + lev] = lap;
        }
        if (var == 1) {
            diffmh_d[id * nv * 3 + lev * 3 + 0] = lap;
            if (DiffSponge && order_diff_sponge == 2) {
                diffmh_d[id * nv * 3 + lev * 3 + 0] +=
                    Rho_d[id * nv + lev] * Kdh2_d[lev] * diff_d[id * nv * 6 + lev * 6 + var];
            }
        }
        if (var == 2) {
            diffmh_d[id * nv * 3 + lev * 3 + 1] = lap;
            if (DiffSponge && order_diff_sponge == 2) {
                diffmh_d[id * nv * 3 + lev * 3 + 1] +=
                    Rho_d[id * nv + lev] * Kdh2_d[lev] * diff_d[id * nv * 6 + lev * 6 + var];
            }
        }
        if (var == 3) {
            diffmh_d[id * nv * 3 + lev * 3 + 2] = lap;
            if (DiffSponge && order_diff_sponge == 2) {
                diffmh_d[id * nv * 3 + lev * 3 + 2] +=
                    Rho_d[id * nv + lev] * Kdh2_d[lev] * diff_d[id * nv * 6 + lev * 6 + var];
            }
        }
        if (var == 4) {
            diffw_d[id * nv + lev] = lap;
            if (DiffSponge && order_diff_sponge == 2) {
                diffw_d[id * nv + lev] +=
                    Rho_d[id * nv + lev] * Kdh2_d[lev] * diff_d[id * nv * 6 + lev * 6 + var];
            }
        }
        if (var == 5)
            diffpr_d[id * nv + lev] = lap;
    }
    else {
        diff_d[id * nv * 6 + lev * 6 + var] = lap;
    }
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
                                   double* areas_d,
                                   double* nvecoa_d,
                                   double* nvecti_d,
                                   double* nvecte_d,
                                   double* K_d,
                                   double* Altitude_d,
                                   double* Altitudeh_d,
                                   double  A,
                                   double* Rd_d,
                                   double* Cp_d,
                                   int*    local_d,
                                   int     num,
                                   bool    laststep,
                                   bool    DeepModel,
                                   bool    DiffSponge,
                                   int     order_diff_sponge,
                                   double* Kdh2_d,
                                   double* boundary_flux_d,
                                   bool    energy_equation) {

    int id = blockIdx.x * blockDim.x + threadIdx.x;
    id += num - 2; // Poles

    int nv  = gridDim.y;
    int lev = blockIdx.y;
    int var = blockIdx.z;

    double alt, sdiff, vdiff, astar;
    double rscale;
    int    jp1, jp2;
    int    kp1, kp2;
    // double o6 = 1.0 / 6.0;
    double AT1, AT2;
    double lap, lap1, lap2, lap3, lap4, lap5, lap6;
    double lapx1, lapx2, lapy1, lapy2, lapz1, lapz2;
    /////////////////////////////////////////
    __shared__ double a_p[NN + 1];
    __shared__ double Rho_p[NN + 1];

    __shared__ double areasTr_p[NN];
    __shared__ double areas_p[NN * 3];
    __shared__ double nvecoa_p[NN * 3];
    __shared__ double nvecti_p[NN * 3];
    __shared__ double nvecte_p[NN * 3];
    __shared__ int    local_p[NN];
    /////////////////////////////////////////

    for (int i = 0; i < 5; i++) //why do we do this loop five times in sequence??
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

    alt = Altitude_d[lev];

    Rho_p[0] = Rho_d[id * nv + lev];
    for (int i = 1; i < 6; i++)
        Rho_p[i] = Rho_d[local_p[i - 1] * nv + lev];

    if (laststep) {
        sdiff = -K_d[lev];
    }

    if (laststep) {
        if (var == 0) {
            a_p[0] = diff_d[id * nv * 6 + lev * 6 + var];
        }
        else {
            a_p[0] = diff_d[id * nv * 6 + lev * 6 + var] * Rho_p[0];
        }
        for (int i = 1; i < 6; i++) {
            if (var == 0) {
                a_p[i] = diff_d[local_p[i - 1] * nv * 6 + lev * 6 + var];
            }
            else {
                a_p[i] = diff_d[local_p[i - 1] * nv * 6 + lev * 6 + var] * Rho_p[i];
            }
        }
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
            if (energy_equation) {
                a_p[0] = (Cp_d[id * nv + lev] - Rd_d[id * nv + lev]) * temperature_d[id * nv + lev];
                for (int i = 1; i < 6; i++)
                    a_p[i] = (Cp_d[local_p[i - 1] * nv + lev] - Rd_d[local_p[i - 1] * nv + lev])
                             * temperature_d[local_p[i - 1] * nv + lev];
            }
            else {
                a_p[0] = Rd_d[id * nv + lev] * temperature_d[id * nv + lev];
                for (int i = 1; i < 6; i++)
                    a_p[i] =
                        Rd_d[local_p[i - 1] * nv + lev] * temperature_d[local_p[i - 1] * nv + lev];
            }
        }
    }


    if (DeepModel)
        rscale = A / (alt + A);
    else
        rscale = 1.0;

    lap = 0.0;

    for (int k = 0; k < 5; k++) {
        int j = k + 1;
        jp1   = (j) % 5 + 1;
        jp2   = (j + 1) % 5 + 1;
        kp1   = (k + 1) % 5;
        kp2   = (k + 2) % 5;

        if (laststep) {
            // if (var == 0) {
            vdiff = 0.5 * sdiff * rscale;
            // }
            // else {
            //     vdiff = 0.5 * sdiff * rscale
            //             * (2.0 * Rho_p[0] + Rho_p[j] + 2.0 * Rho_p[jp1] + Rho_p[jp2]) * o6;
            // }
        }
        else
            vdiff = 0.5 * rscale;

        if (k == 0) {
            AT1   = rscale / areasTr_p[k];
            astar = (areas_p[k * 3 + 0] * a_p[0] + areas_p[k * 3 + 1] * a_p[j]
                     + areas_p[k * 3 + 2] * a_p[jp1])
                    / areasTr_p[k];
            lap1 = (0.5 * (a_p[j] + a_p[0]) - astar);
            lap2 = (0.5 * (a_p[j] + a_p[jp1]) - astar);
            lap3 = (0.5 * (a_p[0] + a_p[jp1]) - astar);
        }
        AT2   = rscale / areasTr_p[kp1];
        astar = (areas_p[kp1 * 3 + 0] * a_p[0] + areas_p[kp1 * 3 + 1] * a_p[jp1]
                 + areas_p[kp1 * 3 + 2] * a_p[jp2])
                / areasTr_p[kp1];
        lap4 = (0.5 * (a_p[jp1] + a_p[0]) - astar);
        lap5 = (0.5 * (a_p[jp2] + a_p[jp1]) - astar);
        lap6 = (0.5 * (a_p[0] + a_p[jp2]) - astar);

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

        lap += ((lapx1 + lapx2) * nvecoa_p[k * 3 + 0] + (lapy1 + lapy2) * nvecoa_p[k * 3 + 1]
                + (lapz1 + lapz2) * nvecoa_p[k * 3 + 2])
               * vdiff;
    }
    if (laststep) {
        if (var == 0)
            diffrh_d[id * nv + lev] = lap;
        if (var == 1) {
            diffmh_d[id * 3 * nv + lev * 3 + 0] = lap;
            if (DiffSponge && order_diff_sponge == 2) {
                diffmh_d[id * nv * 3 + lev * 3 + 0] +=
                    Rho_d[id * nv + lev] * Kdh2_d[lev] * diff_d[id * nv * 6 + lev * 6 + var];
            }
        }
        if (var == 2) {
            diffmh_d[id * 3 * nv + lev * 3 + 1] = lap;
            if (DiffSponge && order_diff_sponge == 2) {
                diffmh_d[id * nv * 3 + lev * 3 + 1] +=
                    Rho_d[id * nv + lev] * Kdh2_d[lev] * diff_d[id * nv * 6 + lev * 6 + var];
            }
        }
        if (var == 3) {
            diffmh_d[id * 3 * nv + lev * 3 + 2] = lap;
            if (DiffSponge && order_diff_sponge == 2) {
                diffmh_d[id * nv * 3 + lev * 3 + 2] +=
                    Rho_d[id * nv + lev] * Kdh2_d[lev] * diff_d[id * nv * 6 + lev * 6 + var];
            }
        }
        if (var == 4) {
            diffw_d[id * nv + lev] = lap;
            if (DiffSponge && order_diff_sponge == 2) {
                diffw_d[id * nv + lev] +=
                    Rho_d[id * nv + lev] * Kdh2_d[lev] * diff_d[id * nv * 6 + lev * 6 + var];
            }
        }
        if (var == 5)
            diffpr_d[id * nv + lev] = lap;
    }
    else {
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

__global__ void Diffusion_Op_Vert(double* diffmv_d,
                                  double* diffwv_d,
                                  double* diffrv_d,
                                  double* diffprv_d,
                                  double* diffv_d1,
                                  double* diffv_d2,
                                  double* Mh_d,
                                  double* Rho_d,
                                  double* temperature_d,
                                  double* W_d,
                                  double* func_r_d,
                                  double* Kv_d,
                                  double* Altitude_d,
                                  double  A,
                                  double  Rd,
                                  int     num,
                                  int     step_num,
                                  bool    DeepModel) {

    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id < num) {
        int nv  = gridDim.y;
        int lev = blockIdx.y;
        int var = blockIdx.z;

        /////////////////////////////////////////
        double a_s[5];
        double Rho_s[5];
        double r_s[5];
        /////////////////////////////////////////

        double dz, lapl;
        int    ilev;
        int    pad, nterms;

        if (lev == 0) { //current lev is ilev = 0
            pad    = 0;
            nterms = 3;
        }
        else if (lev == 1 || lev == nv - 2) { //current lev is ilev = 1
            pad    = 1;
            nterms = 3;
        }
        else if (lev == nv - 1) { //current lev is ilev = 2
            pad    = 2;
            nterms = 3;
        }
        else { //current lev is ilev = 2
            pad    = 2;
            nterms = 5;
        }

        for (ilev = 0; ilev < nterms; ilev++) {
            Rho_s[ilev] = Rho_d[id * nv + lev - pad + ilev];
            r_s[ilev]   = A + Altitude_d[lev - pad + ilev];

            if (step_num == 0) {
                if (var == 0)
                    a_s[ilev] = Rho_d[id * nv + lev - pad + ilev];
                else if (var == 1)
                    a_s[ilev] = Mh_d[id * nv * 3 + (lev - pad + ilev) * 3 + 0] / Rho_s[ilev];
                else if (var == 2)
                    a_s[ilev] = Mh_d[id * nv * 3 + (lev - pad + ilev) * 3 + 1] / Rho_s[ilev];
                else if (var == 3)
                    a_s[ilev] = Mh_d[id * nv * 3 + (lev - pad + ilev) * 3 + 2] / Rho_s[ilev];
                else if (var == 4)
                    a_s[ilev] = W_d[id * nv + lev - pad + ilev] / Rho_s[ilev];
                else if (var == 5)
                    a_s[ilev] = temperature_d[id * nv + lev - pad + ilev];
            }
            else if (step_num == 1) {
                a_s[ilev] = diffv_d1[id * nv * 6 + (lev - pad + ilev) * 6 + var];
            }
            else if (step_num == 2) {
                a_s[ilev] = -Kv_d[lev] * diffv_d2[id * nv * 6 + (lev - pad + ilev) * 6 + var];
                if (var > 0)
                    a_s[ilev] *= Rho_s[ilev];
                if (var == 5)
                    a_s[ilev] *= Rd;
            }
            if (DeepModel)
                a_s[ilev] *= r_s[ilev];
        }

        // grid spacing in vertical
        dz = Altitude_d[1] - Altitude_d[0];

        // now, calculate vertical laplacian term
        if (nterms == 3) {                                        //second order accuracy
            lapl = (a_s[2] - 2.0 * a_s[1] + a_s[0]) / pow(dz, 2); //
        }
        else if (nterms == 5) { //fourth order accuracy
            lapl = (-a_s[4] + 16.0 * a_s[3] - 30.0 * a_s[2] + 16.0 * a_s[1] - a_s[0]) / pow(dz, 2)
                   / 12.0;
        }
        if (DeepModel)
            lapl *= 1.0 / r_s[pad];

        // if (id == 0 && var == 0 && lev == 5) {
        //     // printf("argghhh1!\n"); //stupid place to put a stupid break point
        //     lapl *= 1.0;
        // }

        if (step_num == 0) {
            diffv_d1[id * nv * 6 + lev * 6 + var] = lapl; //
        }
        else if (step_num == 1) {
            diffv_d2[id * nv * 6 + lev * 6 + var] = lapl; //
        }
        else if (step_num == 2) {
            if (var == 0)
                diffrv_d[id * nv + lev] = lapl;
            if (var == 1)
                diffmv_d[id * nv * 3 + lev * 3 + 0] = lapl;
            if (var == 2)
                diffmv_d[id * nv * 3 + lev * 3 + 1] = lapl;
            if (var == 3)
                diffmv_d[id * nv * 3 + lev * 3 + 2] = lapl;
            if (var == 4)
                diffwv_d[id * nv + lev] = lapl;
            if (var == 5)
                diffprv_d[id * nv + lev] = lapl;
        }
    }
}


template<int nv>
__global__ void vertical_diff_joao(double* diffmv_d,      //
                                   double* diffwv_d,      //
                                   double* diffrv_d,      //
                                   double* diffprv_d,     //
                                   double* Mh_d,          //
                                   double* Rho_d,         //
                                   double* Temperature_d, //
                                   double* W_d,           //
                                   double* Kv_d,          //
                                   double* Altitude_d,    //
                                   double* Altitudeh_d,
                                   double  A,  //
                                   double  Rd, //
                                   bool    DeepModel,
                                   int     num) { //

    //for reference: original code by joao mendonca

    int id  = blockIdx.x * blockDim.x + threadIdx.x;
    int var = blockIdx.z;

    double ai[nv + 2];
    double af[nv + 2];

    double dz[nv];
    double dzh[nv + 1];

    double rhoh[nv + 1];

    if (id < num) {
        for (int lev = 0; lev < nv + 2; lev++) {
            if (lev == 0) {
                double r_d = Rho_d[id * nv + 0];
                if (var == 0)
                    ai[lev] = r_d;
                if (var == 1)
                    ai[lev] = Mh_d[id * nv * 3 + 0] / r_d;
                if (var == 2)
                    ai[lev] = Mh_d[id * nv * 3 + 1] / r_d;
                if (var == 3)
                    ai[lev] = Mh_d[id * nv * 3 + 2] / r_d;
                if (var == 4)
                    ai[lev] = W_d[id * nv + 0] / r_d;
                if (var == 5)
                    ai[lev] = Temperature_d[id * nv + 0];
            }
            else if (lev == nv + 1) {
                double r_d = Rho_d[id * nv + nv - 1];
                if (var == 0)
                    ai[lev] = r_d;
                if (var == 1)
                    ai[lev] = Mh_d[id * nv * 3 + (nv - 1) * 3 + 0] / r_d;
                if (var == 2)
                    ai[lev] = Mh_d[id * nv * 3 + (nv - 1) * 3 + 1] / r_d;
                if (var == 3)
                    ai[lev] = Mh_d[id * nv * 3 + (nv - 1) * 3 + 2] / r_d;
                if (var == 4)
                    ai[lev] = W_d[id * nv + (nv - 1)] / r_d;
                if (var == 5)
                    ai[lev] = Temperature_d[id * nv + nv - 1];
            }
            else {
                double r_d = Rho_d[id * nv + lev - 1];
                if (var == 0)
                    ai[lev] = r_d;
                if (var == 1)
                    ai[lev] = Mh_d[id * nv * 3 + (lev - 1) * 3 + 0] / r_d;
                if (var == 2)
                    ai[lev] = Mh_d[id * nv * 3 + (lev - 1) * 3 + 1] / r_d;
                if (var == 3)
                    ai[lev] = Mh_d[id * nv * 3 + (lev - 1) * 3 + 2] / r_d;
                if (var == 4)
                    ai[lev] = W_d[id * nv + (lev - 1)] / r_d;
                if (var == 5)
                    ai[lev] = Temperature_d[id * nv + lev - 1];
            }

            if (lev < nv)
                dz[lev] = 1.0 / (Altitudeh_d[lev + 1] - Altitudeh_d[lev]);
            if (lev == 0)
                dzh[lev] = 1.0 / (2.0 * Altitude_d[0]);
            else if (lev == nv)
                dzh[lev] = 1.0 / (2 * (Altitudeh_d[nv] - Altitude_d[nv - 1]));
            else if (lev < nv)
                dzh[lev] = 1.0 / (Altitude_d[lev] - Altitude_d[lev - 1]);
        }

        for (int l = 0; l < 2; l++) {
            for (int lev = 0; lev < nv; lev++)
                af[lev] = ((ai[lev + 2] - ai[lev + 1]) * dzh[lev + 1]
                           - (ai[lev + 1] - ai[lev]) * dzh[lev])
                          * dz[lev];
            if (l == 0) {
                ai[0]      = af[0];
                ai[nv + 1] = af[nv - 1];
                for (int lev = 0; lev < nv; lev++)
                    ai[lev + 1] = af[lev];
            }
            else {
                ai[0]      = af[0];
                ai[nv + 1] = af[nv - 1];
                for (int lev = 0; lev < nv; lev++)
                    ai[lev + 1] = af[lev];
            }
        }

        for (int lev = 0; lev < nv + 1; lev++) {
            if (lev == 0)
                rhoh[0] = Rho_d[id * nv] * Kv_d[0];
            else if (lev == nv)
                rhoh[nv] = Rho_d[id * nv + nv - 1] * Kv_d[nv - 1];
            else {
                double rl   = Rho_d[id * nv + lev - 1] * Kv_d[lev - 1];
                double rp   = Rho_d[id * nv + lev] * Kv_d[lev];
                double xi   = Altitudeh_d[lev];
                double xim  = Altitude_d[lev - 1];
                double xip  = Altitude_d[lev];
                double intl = (xi - xip) / (xim - xip);
                double intp = (xi - xim) / (xip - xim);
                rhoh[lev]   = (rl * intl + rp * intp);
            }
        }

        for (int lev = 0; lev < nv; lev++) {
            if (var != 5) {
                af[lev] = ((ai[lev + 2] - ai[lev + 1]) * dzh[lev + 1] * rhoh[lev + 1]
                           - (ai[lev + 1] - ai[lev]) * dzh[lev] * rhoh[lev])
                          * dz[lev];
            }
            else {
                af[lev] = ((ai[lev + 2] - ai[lev + 1]) * dzh[lev + 1] * rhoh[lev + 1]
                           - (ai[lev + 1] - ai[lev]) * dzh[lev] * rhoh[lev])
                          * dz[lev] * Rd;
            }
        }

        for (int lev = 0; lev < nv; lev++) {
            if (var == 0)
                diffrv_d[id * nv + lev] = af[lev];
            if (var == 1)
                diffmv_d[id * nv * 3 + lev * 3 + 0] = af[lev];
            if (var == 2)
                diffmv_d[id * nv * 3 + lev * 3 + 1] = af[lev];
            if (var == 3)
                diffmv_d[id * nv * 3 + lev * 3 + 2] = af[lev];
            if (var == 4)
                diffwv_d[id * nv + lev] = af[lev];
            if (var == 5)
                diffprv_d[id * nv + lev] = af[lev];
        }
    }
}

__global__ void vertical_diff(double* diffmv_d,  //
                              double* diffwv_d,  //
                              double* diffrv_d,  //
                              double* diffprv_d, //
                              double* diff_d,
                              double* diff2_d,
                              double* Mh_d,          //
                              double* Rho_d,         //
                              double* Temperature_d, //
                              double* W_d,           //
                              double* Kv_d,          //
                              double* Altitude_d,    //
                              double* Altitudeh_d,
                              double  A,  //
                              double  Rd, //
                              bool    DeepModel,
                              int     num,
                              int     nv) { //

    // based on joao's code, with updates and corrections by russell

    int id  = blockIdx.x * blockDim.x + threadIdx.x;
    int var = blockIdx.z;

    double sign        = 1.0;
    bool   include_rho = false;
    //order of diffusion: can be changed (minimum = 4, even numbers only) but
    //the calculation of Kv needs to be adjusted in esp_initial.cu
    int order = 6;

    if (id < num) {
        double dzi, dzup, dzlow;
        double ri, rup, rlow;
        for (int lev = 0; lev < nv + 2; lev++) {
            if (lev == 0) {
                double r_d = Rho_d[id * nv + 0];
                if (var == 0)
                    diff2_d[id * (nv + 2) * 6 + lev * 6 + var] = r_d;
                if (var == 1)
                    diff2_d[id * (nv + 2) * 6 + lev * 6 + var] = Mh_d[id * nv * 3 + 0] / r_d;
                if (var == 2)
                    diff2_d[id * (nv + 2) * 6 + lev * 6 + var] = Mh_d[id * nv * 3 + 1] / r_d;
                if (var == 3)
                    diff2_d[id * (nv + 2) * 6 + lev * 6 + var] = Mh_d[id * nv * 3 + 2] / r_d;
                if (var == 4)
                    diff2_d[id * (nv + 2) * 6 + lev * 6 + var] = W_d[id * nv + 0] / r_d;
                if (var == 5)
                    diff2_d[id * (nv + 2) * 6 + lev * 6 + var] = Temperature_d[id * nv + 0];
            }
            else if (lev == nv + 1) {
                double r_d = Rho_d[id * nv + nv - 1];
                if (var == 0)
                    diff2_d[id * (nv + 2) * 6 + lev * 6 + var] = r_d;
                if (var == 1)
                    diff2_d[id * (nv + 2) * 6 + lev * 6 + var] =
                        Mh_d[id * nv * 3 + (nv - 1) * 3 + 0] / r_d;
                if (var == 2)
                    diff2_d[id * (nv + 2) * 6 + lev * 6 + var] =
                        Mh_d[id * nv * 3 + (nv - 1) * 3 + 1] / r_d;
                if (var == 3)
                    diff2_d[id * (nv + 2) * 6 + lev * 6 + var] =
                        Mh_d[id * nv * 3 + (nv - 1) * 3 + 2] / r_d;
                if (var == 4)
                    diff2_d[id * (nv + 2) * 6 + lev * 6 + var] = W_d[id * nv + (nv - 1)] / r_d;
                if (var == 5)
                    diff2_d[id * (nv + 2) * 6 + lev * 6 + var] = Temperature_d[id * nv + nv - 1];
            }
            else {
                double r_d = Rho_d[id * nv + lev - 1];
                if (var == 0)
                    diff2_d[id * (nv + 2) * 6 + lev * 6 + var] = r_d;
                if (var == 1)
                    diff2_d[id * (nv + 2) * 6 + lev * 6 + var] =
                        Mh_d[id * nv * 3 + (lev - 1) * 3 + 0] / r_d;
                if (var == 2)
                    diff2_d[id * (nv + 2) * 6 + lev * 6 + var] =
                        Mh_d[id * nv * 3 + (lev - 1) * 3 + 1] / r_d;
                if (var == 3)
                    diff2_d[id * (nv + 2) * 6 + lev * 6 + var] =
                        Mh_d[id * nv * 3 + (lev - 1) * 3 + 2] / r_d;
                if (var == 4)
                    diff2_d[id * (nv + 2) * 6 + lev * 6 + var] = W_d[id * nv + (lev - 1)] / r_d;
                if (var == 5)
                    diff2_d[id * (nv + 2) * 6 + lev * 6 + var] = Temperature_d[id * nv + lev - 1];
            }
        }

        int max_l = order / 2 - 1;
        for (int l = 0; l < max_l; l++) {
            for (int lev = 0; lev < nv; lev++) {
                dzi = 1.0 / (Altitudeh_d[lev + 1] - Altitudeh_d[lev]);
                if (DeepModel) {
                    ri   = A + Altitude_d[lev];
                    rup  = A + Altitudeh_d[lev + 1];
                    rlow = A + Altitudeh_d[lev];
                }
                else {
                    ri   = 1.0;
                    rup  = 1.0;
                    rlow = 1.0;
                }
                if (lev == 0) {
                    dzup  = 1.0 / (Altitude_d[1] - Altitude_d[0]);
                    dzlow = 1.0 / (2.0 * Altitude_d[0]);
                }
                else if (lev == nv - 1) {
                    dzup  = 1.0 / (2 * (Altitudeh_d[lev + 1] - Altitude_d[lev]));
                    dzlow = 1.0 / (Altitude_d[lev] - Altitude_d[lev - 1]);
                }
                else {
                    dzup  = 1.0 / (Altitude_d[lev + 1] - Altitude_d[lev]);
                    dzlow = 1.0 / (Altitude_d[lev] - Altitude_d[lev - 1]);
                }
                diff_d[id * nv * 6 + lev * 6 + var] =
                    (pow(rup, 2)
                         * (diff2_d[id * (nv + 2) * 6 + (lev + 2) * 6 + var]
                            - diff2_d[id * (nv + 2) * 6 + (lev + 1) * 6 + var])
                         * dzup
                     - pow(rlow, 2)
                           * (diff2_d[id * (nv + 2) * 6 + (lev + 1) * 6 + var]
                              - diff2_d[id * (nv + 2) * 6 + lev * 6 + var])
                           * dzlow)
                    * dzi / pow(ri, 2);
            }
            if (l < max_l - 1) {
                diff2_d[id * (nv + 2) * 6 + 0 * 6 + var] = diff_d[id * nv * 6 + 0 * 6 + var];
                diff2_d[id * (nv + 2) * 6 + (nv + 1) * 6 + var] =
                    diff_d[id * nv * 6 + (nv - 1) * 6 + var];
                for (int lev = 0; lev < nv; lev++)
                    diff2_d[id * (nv + 2) * 6 + (lev + 1) * 6 + var] =
                        diff_d[id * nv * 6 + lev * 6 + var];
            }
            else {
                if (var == 0) {
                    diff2_d[id * (nv + 2) * 6 + 0 * 6 + var] =
                        sign * Kv_d[0] * diff_d[id * nv * 6 + 0 * 6 + var];
                    diff2_d[id * (nv + 2) * 6 + (nv + 1) * 6 + var] =
                        sign * Kv_d[nv - 1] * diff_d[id * nv * 6 + (nv - 1) * 6 + var];
                    for (int lev = 0; lev < nv; lev++)
                        diff2_d[id * (nv + 2) * 6 + (lev + 1) * 6 + var] =
                            sign * Kv_d[lev] * diff_d[id * nv * 6 + lev * 6 + var];
                }
                else {
                    diff2_d[id * (nv + 2) * 6 + 0 * 6 + var] =
                        sign * Rho_d[id * nv + 0] * Kv_d[0] * diff_d[id * nv * 6 + 0 * 6 + var];
                    diff2_d[id * (nv + 2) * 6 + (nv + 1) * 6 + var] =
                        sign * Rho_d[id * nv + nv - 1] * Kv_d[nv - 1]
                        * diff_d[id * nv * 6 + (nv - 1) * 6 + var];
                    for (int lev = 0; lev < nv; lev++)
                        diff2_d[id * (nv + 2) * 6 + (lev + 1) * 6 + var] =
                            sign * Rho_d[id * nv + lev] * Kv_d[lev]
                            * diff_d[id * nv * 6 + lev * 6 + var];
                }
            }
        }

        for (int lev = 0; lev < nv; lev++) {
            dzi = 1.0 / (Altitudeh_d[lev + 1] - Altitudeh_d[lev]);
            if (DeepModel) {
                ri   = A + Altitude_d[lev];
                rup  = A + Altitudeh_d[lev + 1];
                rlow = A + Altitudeh_d[lev];
            }
            else {
                ri   = 1.0;
                rup  = 1.0;
                rlow = 1.0;
            }
            if (lev == 0) {
                dzup  = 1.0 / (Altitude_d[1] - Altitude_d[0]);
                dzlow = 1.0 / (2.0 * Altitude_d[0]);
            }
            else if (lev == nv - 1) {
                dzup  = 1.0 / (2 * (Altitudeh_d[lev + 1] - Altitude_d[lev]));
                dzlow = 1.0 / (Altitude_d[lev] - Altitude_d[lev - 1]);
            }
            else {
                dzup  = 1.0 / (Altitude_d[lev + 1] - Altitude_d[lev]);
                dzlow = 1.0 / (Altitude_d[lev] - Altitude_d[lev - 1]);
            }
            if (var != 5) {
                diff_d[id * nv * 6 + lev * 6 + var] =
                    (pow(rup, 2)
                         * (diff2_d[id * (nv + 2) * 6 + (lev + 2) * 6 + var]
                            - diff2_d[id * (nv + 2) * 6 + (lev + 1) * 6 + var])
                         * dzup
                     - pow(rlow, 2)
                           * (diff2_d[id * (nv + 2) * 6 + (lev + 1) * 6 + var]
                              - diff2_d[id * (nv + 2) * 6 + lev * 6 + var])
                           * dzlow)
                    * dzi / pow(ri, 2);
            }
            else {
                diff_d[id * nv * 6 + lev * 6 + var] =
                    (pow(rup, 2)
                         * (diff2_d[id * (nv + 2) * 6 + (lev + 2) * 6 + var]
                            - diff2_d[id * (nv + 2) * 6 + (lev + 1) * 6 + var])
                         * dzup
                     - pow(rlow, 2)
                           * (diff2_d[id * (nv + 2) * 6 + (lev + 1) * 6 + var]
                              - diff2_d[id * (nv + 2) * 6 + lev * 6 + var])
                           * dzlow)
                    * dzi * Rd / pow(ri, 2);
            }
        }

        for (int lev = 0; lev < nv; lev++) {
            if (var == 0) {
                if (include_rho) {
                    diffrv_d[id * nv + lev] = diff_d[id * nv * 6 + lev * 6 + var];
                }
                else {
                    diffrv_d[id * nv + lev] = 0.0;
                }
            }
            if (var == 1)
                diffmv_d[id * nv * 3 + lev * 3 + 0] = diff_d[id * nv * 6 + lev * 6 + var];
            if (var == 2)
                diffmv_d[id * nv * 3 + lev * 3 + 1] = diff_d[id * nv * 6 + lev * 6 + var];
            if (var == 3)
                diffmv_d[id * nv * 3 + lev * 3 + 2] = diff_d[id * nv * 6 + lev * 6 + var];
            if (var == 4)
                diffwv_d[id * nv + lev] = diff_d[id * nv * 6 + lev * 6 + var];
            if (var == 5)
                diffprv_d[id * nv + lev] = diff_d[id * nv * 6 + lev * 6 + var];
        }
    }
}
