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
// Description: Compute the 3D advection and coriolis terms.
//
//
// Method: Central finite-volume in 3D cartesian coordinates.
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

#include "kernel_halo_helpers.h"


// Computes the 3D advection terms and the velocities in cartesians.
// (no poles)
template<int NX, int NY>
__global__ void Compute_Advec_Cori1(
    double3* Adv_d,      // Advection term.
    double3* v_d,        // 3D velocity (cartesians).
    double3* Mh_d,       // Horizontal momentum.
    double3* div_d,      // Divergence operator.
    double*  W_d,        // Vertical momentum.
    double*  Rho_d,      // Density.
    double*  Altitude_d, // Altitude.
    double   A,          // Radius.
    double3* func_r_d,   // Unit vectors normal to the spherical surface.
    int*     maps_d,     // Global indexes.
    int      nl_region,  // Length of the side of each rhombi (from the sphere decomposition).
    bool DeepModel) {    // Switches on and off the deep atmosphere solution (see headers/define.h).

    int x = threadIdx.x;
    int y = threadIdx.y;
    //int ib  = blockIdx.x;
    int nv  = gridDim.y;
    int lev = blockIdx.y;

    int     pt1, pt2, pt3, pt4, pt5, pt6;
    double3 div0, div1, div2, div3, div4, div5, div6;
    int     nhl  = nl_region + 2;
    int     nhl2 = nhl * nhl;

    double alt;
    double rho;
    double rscale;

    double w;


    /////////////////////////////////////////
    __shared__ double3 v_s[(NX + 2) * (NY + 2)];
    __shared__ double3 M_s[(NX + 2) * (NY + 2)];
    __shared__ double3 nflx_s[NX * NY];
    /////////////////////////////////////////

    int ir = 0; // index in region
    int iri, ir2, id;


    bool pent_ind = false; //
    int  ig;               // index in global mem

    int igh = 0; // global index in halo

    // Load shared memory


    bool load_halo = compute_mem_idx(maps_d, nhl, nhl2, ig, igh, ir, ir2, pent_ind);
    id             = ig;


    M_s[ir]        = Mh_d[ig * nv + lev];
    rho            = Rho_d[ig * nv + lev];
    w              = W_d[ig * nv + lev];
    double3 func_r = func_r_d[ig];

    v_s[ir].x = M_s[ir].x / rho + (w / rho) * func_r.x;
    v_s[ir].y = M_s[ir].y / rho + (w / rho) * func_r.y;
    v_s[ir].z = M_s[ir].z / rho + (w / rho) * func_r.z;

    ///////////////////////////////
    //////////// Halo /////////////
    ///////////////////////////////

    if (load_halo) {
        if (igh >= 0) {
            M_s[ir2]       = Mh_d[igh * nv + lev];
            rho            = Rho_d[igh * nv + lev];
            w              = W_d[igh * nv + lev];
            double3 func_r = func_r_d[igh];

            v_s[ir2].x = M_s[ir2].x / rho + (w / rho) * func_r.x;
            v_s[ir2].y = M_s[ir2].y / rho + (w / rho) * func_r.y;
            v_s[ir2].z = M_s[ir2].z / rho + (w / rho) * func_r.z;
        }
        // we can probably ignore this as it's never used and set to 0 by divergence for pentagons.
        //   either set it to zero automagically or use mappings
        else {
            // set to 0 if it is a pentagon index
            v_s[ir2].x = 0.0;
            v_s[ir2].y = 0.0;
            v_s[ir2].z = 0.0;
            M_s[ir2].x = 0.0;
            M_s[ir2].y = 0.0;
            M_s[ir2].z = 0.0;
        }
    }

    __syncthreads();

    //////////////////////////////////////////////
    // Inner index
    iri = y * nl_region + x;

    // Neighbours.
    pt1 = (y + 2) * nhl + x + 1;
    pt2 = (y + 2) * nhl + x + 2;
    pt3 = (y + 1) * nhl + x + 2;
    pt4 = (y)*nhl + x + 1;
    pt5 = pent_ind * ((y + 1) * nhl + x) + (!pent_ind) * ((y)*nhl + x);
    pt6 = (y + 1) * nhl + x;

    // Sets the weight needed for the deep atmosphere solution.
    if (DeepModel) {
        alt    = Altitude_d[lev];
        rscale = A / (alt + A);
    }
    else
        rscale = 1.0;

    // Initialize fluxes.
    nflx_s[iri].x = 0.0;
    nflx_s[iri].y = 0.0;
    nflx_s[iri].z = 0.0;

    // Calculate divergence.
    {
        div0 = div_d[id * 7 + 0];
        div1 = div_d[id * 7 + 1];
        div2 = div_d[id * 7 + 2];
        div3 = div_d[id * 7 + 3];
        div4 = div_d[id * 7 + 4];
        div5 = div_d[id * 7 + 5];
        div6 = div_d[id * 7 + 6]; // For pent_ind = 1 div6 is equal to 0.

        nflx_s[iri].x += rscale
                         * (div0.x * v_s[ir].x * M_s[ir].x + div1.x * v_s[pt1].x * M_s[pt1].x
                            + div2.x * v_s[pt2].x * M_s[pt2].x + div3.x * v_s[pt3].x * M_s[pt3].x
                            + div4.x * v_s[pt4].x * M_s[pt4].x + div5.x * v_s[pt5].x * M_s[pt5].x
                            + div6.x * v_s[pt6].x * M_s[pt6].x);

        nflx_s[iri].y += rscale
                         * (div0.x * v_s[ir].y * M_s[ir].x + div1.x * v_s[pt1].y * M_s[pt1].x
                            + div2.x * v_s[pt2].y * M_s[pt2].x + div3.x * v_s[pt3].y * M_s[pt3].x
                            + div4.x * v_s[pt4].y * M_s[pt4].x + div5.x * v_s[pt5].y * M_s[pt5].x
                            + div6.x * v_s[pt6].y * M_s[pt6].x);

        nflx_s[iri].z += rscale
                         * (div0.x * v_s[ir].z * M_s[ir].x + div1.x * v_s[pt1].z * M_s[pt1].x
                            + div2.x * v_s[pt2].z * M_s[pt2].x + div3.x * v_s[pt3].z * M_s[pt3].x
                            + div4.x * v_s[pt4].z * M_s[pt4].x + div5.x * v_s[pt5].z * M_s[pt5].x
                            + div6.x * v_s[pt6].z * M_s[pt6].x);

        nflx_s[iri].x += rscale
                         * (div0.y * v_s[ir].x * M_s[ir].y + div1.y * v_s[pt1].x * M_s[pt1].y
                            + div2.y * v_s[pt2].x * M_s[pt2].y + div3.y * v_s[pt3].x * M_s[pt3].y
                            + div4.y * v_s[pt4].x * M_s[pt4].y + div5.y * v_s[pt5].x * M_s[pt5].y
                            + div6.y * v_s[pt6].x * M_s[pt6].y);

        nflx_s[iri].y += rscale
                         * (div0.y * v_s[ir].y * M_s[ir].y + div1.y * v_s[pt1].y * M_s[pt1].y
                            + div2.y * v_s[pt2].y * M_s[pt2].y + div3.y * v_s[pt3].y * M_s[pt3].y
                            + div4.y * v_s[pt4].y * M_s[pt4].y + div5.y * v_s[pt5].y * M_s[pt5].y
                            + div6.y * v_s[pt6].y * M_s[pt6].y);

        nflx_s[iri].z += rscale
                         * (div0.y * v_s[ir].z * M_s[ir].y + div1.y * v_s[pt1].z * M_s[pt1].y
                            + div2.y * v_s[pt2].z * M_s[pt2].y + div3.y * v_s[pt3].z * M_s[pt3].y
                            + div4.y * v_s[pt4].z * M_s[pt4].y + div5.y * v_s[pt5].z * M_s[pt5].y
                            + div6.y * v_s[pt6].z * M_s[pt6].y);

        nflx_s[iri].x += rscale
                         * (div0.z * v_s[ir].x * M_s[ir].z + div1.z * v_s[pt1].x * M_s[pt1].z
                            + div2.z * v_s[pt2].x * M_s[pt2].z + div3.z * v_s[pt3].x * M_s[pt3].z
                            + div4.z * v_s[pt4].x * M_s[pt4].z + div5.z * v_s[pt5].x * M_s[pt5].z
                            + div6.z * v_s[pt6].x * M_s[pt6].z);

        nflx_s[iri].y += rscale
                         * (div0.z * v_s[ir].y * M_s[ir].z + div1.z * v_s[pt1].y * M_s[pt1].z
                            + div2.z * v_s[pt2].y * M_s[pt2].z + div3.z * v_s[pt3].y * M_s[pt3].z
                            + div4.z * v_s[pt4].y * M_s[pt4].z + div5.z * v_s[pt5].y * M_s[pt5].z
                            + div6.z * v_s[pt6].y * M_s[pt6].z);
        nflx_s[iri].z += rscale
                         * (div0.z * v_s[ir].z * M_s[ir].z + div1.z * v_s[pt1].z * M_s[pt1].z
                            + div2.z * v_s[pt2].z * M_s[pt2].z + div3.z * v_s[pt3].z * M_s[pt3].z
                            + div4.z * v_s[pt4].z * M_s[pt4].z + div5.z * v_s[pt5].z * M_s[pt5].z
                            + div6.z * v_s[pt6].z * M_s[pt6].z);
    }
    // Return values (3D advection term and velocities in cartesians).
    Adv_d[id * nv + lev] = nflx_s[iri];
    v_d[id * nv + lev]   = v_s[ir];
}


// Computes the 3D advection terms and the velocities in cartesians.
// (poles)
template<int NN>
__global__ void Compute_Advec_Cori_Poles(double* Adv_d,         //
                                         double* v_d,           //
                                         double* Mh_d,          //
                                         double* div_d,         //
                                         double* W_d,           //
                                         double* Rho_d,         //
                                         double* Altitude_d,    //
                                         double  A,             //
                                         double* func_r_d,      //
                                         int*    point_local_d, //
                                         int     num,           //
                                         int     nv,            //
                                         bool    DeepModel) {      //

    int id = blockIdx.x * blockDim.x + threadIdx.x;
    id += num - 2; // Poles

    double alt, w;
    double rscale;

    double rho;

    /////////////////////////////////////////
    __shared__ double v_p[NN * 3];
    __shared__ double M_p[NN * 3];
    __shared__ double func_r_p[3 * NN];
    __shared__ double div_p[3 * 7];
    __shared__ int    local_p[NN];

    double nflxx;
    double nflxy;
    double nflxz;

    /////////////////////////////////////////
    if (id < num) {
        for (int i = 0; i < 5; i++)
            local_p[i] = point_local_d[id * 6 + i];
        func_r_p[0] = func_r_d[id * 3 + 0];
        func_r_p[1] = func_r_d[id * 3 + 1];
        func_r_p[2] = func_r_d[id * 3 + 2];
        for (int i = 1; i < 6; i++) {
            func_r_p[i * 3 + 0] = func_r_d[local_p[i - 1] * 3 + 0];
            func_r_p[i * 3 + 1] = func_r_d[local_p[i - 1] * 3 + 1];
            func_r_p[i * 3 + 2] = func_r_d[local_p[i - 1] * 3 + 2];
        }
        for (int i = 0; i < 7; i++)
            for (int k = 0; k < 3; k++)
                div_p[i * 3 + k] = div_d[id * 7 * 3 + i * 3 + k];

        for (int lev = 0; lev < nv; lev++) {
            v_p[0] = Mh_d[id * 3 * nv + lev * 3 + 0];
            v_p[1] = Mh_d[id * 3 * nv + lev * 3 + 1];
            v_p[2] = Mh_d[id * 3 * nv + lev * 3 + 2];
            M_p[0] = v_p[0];
            M_p[1] = v_p[1];
            M_p[2] = v_p[2];
            rho    = Rho_d[id * nv + lev];
            w      = W_d[id * nv + lev];
            v_p[0] = v_p[0] / rho + (w / rho) * func_r_p[0];
            v_p[1] = v_p[1] / rho + (w / rho) * func_r_p[1];
            v_p[2] = v_p[2] / rho + (w / rho) * func_r_p[2];
            for (int i = 1; i < 6; i++) {
                v_p[i * 3 + 0] = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 0];
                v_p[i * 3 + 1] = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 1];
                v_p[i * 3 + 2] = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 2];
                M_p[i * 3 + 0] = v_p[i * 3 + 0];
                M_p[i * 3 + 1] = v_p[i * 3 + 1];
                M_p[i * 3 + 2] = v_p[i * 3 + 2];
                rho            = Rho_d[local_p[i - 1] * nv + lev];
                w              = W_d[local_p[i - 1] * nv + lev];

                v_p[i * 3 + 0] = v_p[i * 3 + 0] / rho + (w / rho) * func_r_p[i * 3 + 0];
                v_p[i * 3 + 1] = v_p[i * 3 + 1] / rho + (w / rho) * func_r_p[i * 3 + 1];
                v_p[i * 3 + 2] = v_p[i * 3 + 2] / rho + (w / rho) * func_r_p[i * 3 + 2];
            }

            if (DeepModel) {
                alt    = Altitude_d[lev];
                rscale = A / (alt + A);
            }
            else
                rscale = 1.0;

            nflxx = 0.0;
            nflxy = 0.0;
            nflxz = 0.0;

            for (int k = 0; k < 3; k++) {
                nflxx += rscale
                         * (div_p[3 * 0 + k] * v_p[0 * 3 + 0] * M_p[0 * 3 + k]
                            + div_p[3 * 1 + k] * v_p[1 * 3 + 0] * M_p[1 * 3 + k]
                            + div_p[3 * 2 + k] * v_p[2 * 3 + 0] * M_p[2 * 3 + k]
                            + div_p[3 * 3 + k] * v_p[3 * 3 + 0] * M_p[3 * 3 + k]
                            + div_p[3 * 4 + k] * v_p[4 * 3 + 0] * M_p[4 * 3 + k]
                            + div_p[3 * 5 + k] * v_p[5 * 3 + 0] * M_p[5 * 3 + k]);

                nflxy += rscale
                         * (div_p[3 * 0 + k] * v_p[0 * 3 + 1] * M_p[0 * 3 + k]
                            + div_p[3 * 1 + k] * v_p[1 * 3 + 1] * M_p[1 * 3 + k]
                            + div_p[3 * 2 + k] * v_p[2 * 3 + 1] * M_p[2 * 3 + k]
                            + div_p[3 * 3 + k] * v_p[3 * 3 + 1] * M_p[3 * 3 + k]
                            + div_p[3 * 4 + k] * v_p[4 * 3 + 1] * M_p[4 * 3 + k]
                            + div_p[3 * 5 + k] * v_p[5 * 3 + 1] * M_p[5 * 3 + k]);

                nflxz += rscale
                         * (div_p[3 * 0 + k] * v_p[0 * 3 + 2] * M_p[0 * 3 + k]
                            + div_p[3 * 1 + k] * v_p[1 * 3 + 2] * M_p[1 * 3 + k]
                            + div_p[3 * 2 + k] * v_p[2 * 3 + 2] * M_p[2 * 3 + k]
                            + div_p[3 * 3 + k] * v_p[3 * 3 + 2] * M_p[3 * 3 + k]
                            + div_p[3 * 4 + k] * v_p[4 * 3 + 2] * M_p[4 * 3 + k]
                            + div_p[3 * 5 + k] * v_p[5 * 3 + 2] * M_p[5 * 3 + k]);
            }
            // Return values (3D advection term and velocities in cartesians).
            Adv_d[id * 3 * nv + lev * 3 + 0] = nflxx;
            Adv_d[id * 3 * nv + lev * 3 + 1] = nflxy;
            Adv_d[id * 3 * nv + lev * 3 + 2] = nflxz;
            v_d[id * 3 * nv + lev * 3 + 0]   = v_p[0];
            v_d[id * 3 * nv + lev * 3 + 1]   = v_p[1];
            v_d[id * 3 * nv + lev * 3 + 2]   = v_p[2];
        }
    }
}

// Projects the advection terms on the spherical surface and
// also adds the vertical contribution.
__global__ void Compute_Advec_Cori2(double* Adv_d,       //
                                    double* v_d,         //
                                    double* Wh_d,        //
                                    double* Rho_d,       //
                                    double* Altitude_d,  //
                                    double* Altitudeh_d, //
                                    double* func_r_d,    //
                                    double  Omega,       //
                                    double  A,           //
                                    int     nv,          //
                                    int     num,         //
                                    bool    DeepModel) {    //

    int id = blockIdx.x * blockDim.x + threadIdx.x;

    double altht, althl;
    double altt, alt, altl;
    double r2p, r2m, r2l;

    double wht, whl;
    double vxl, vx, vxt;
    double vyl, vy, vyt;
    double vzl, vz, vzt;

    double intt, intl;
    double dvwxl, dvwyl, dvwzl;
    double dvwxt, dvwyt, dvwzt;
    double davx, davy, davz;
    double rho;

    double xi, xim, xip, dz;
    if (id < num) {
        for (int lev = 0; lev < nv; lev++) {
            if (lev == 0) {
                dvwxl = 0.0;
                dvwyl = 0.0;
                dvwzl = 0.0;

                altht = Altitudeh_d[lev + 1];
                althl = Altitudeh_d[lev];
                alt   = Altitude_d[lev];
                altt  = Altitude_d[lev + 1];
                vx    = v_d[id * 3 * nv + lev * 3 + 0];
                vxt   = v_d[id * 3 * nv + (lev + 1) * 3 + 0];
                vy    = v_d[id * 3 * nv + lev * 3 + 1];
                vyt   = v_d[id * 3 * nv + (lev + 1) * 3 + 1];
                vz    = v_d[id * 3 * nv + lev * 3 + 2];
                vzt   = v_d[id * 3 * nv + (lev + 1) * 3 + 2];
                wht   = Wh_d[id * (nv + 1) + lev + 1];

                // Linear interpolation.
                xi  = altht;
                xim = alt;
                xip = altt;

                intt = (xi - xip) / (xim - xip);
                intl = (xi - xim) / (xip - xim);

                dvwxt = (vx * intt + vxt * intl) * wht;
                dvwyt = (vy * intt + vyt * intl) * wht;
                dvwzt = (vz * intt + vzt * intl) * wht;
            }
            else if (lev == nv - 1) {
                dvwxt = 0.0;
                dvwyt = 0.0;
                dvwzt = 0.0;

                xi  = althl;
                xim = altl;
                xip = alt;

                intt = (xi - xip) / (xim - xip);
                intl = (xi - xim) / (xip - xim);

                dvwxl = (vxl * intt + vx * intl) * whl;
                dvwyl = (vyl * intt + vy * intl) * whl;
                dvwzl = (vzl * intt + vz * intl) * whl;
            }
            else {
                xi  = althl;
                xim = altl;
                xip = alt;

                intt = (xi - xip) / (xim - xip);
                intl = (xi - xim) / (xip - xim);

                dvwxl = (vxl * intt + vx * intl) * whl;
                dvwyl = (vyl * intt + vy * intl) * whl;
                dvwzl = (vzl * intt + vz * intl) * whl;

                xi  = altht;
                xim = alt;
                xip = altt;

                intt = (xi - xip) / (xim - xip);
                intl = (xi - xim) / (xip - xim);

                dvwxt = (vx * intt + vxt * intl) * wht;
                dvwyt = (vy * intt + vyt * intl) * wht;
                dvwzt = (vz * intt + vzt * intl) * wht;
            }

            if (DeepModel) {
                r2p = pow(altht + A, 2.0);
                r2m = pow(alt + A, 2.0);
                r2l = pow(althl + A, 2.0);
            }
            else {
                r2p = 1.0;
                r2m = 1.0;
                r2l = 1.0;
            }

            dz   = 1.0 / ((althl - altht) * r2m);
            davx = (dvwxl * r2l - dvwxt * r2p) * dz;
            davy = (dvwyl * r2l - dvwyt * r2p) * dz;
            davz = (dvwzl * r2l - dvwzt * r2p) * dz;
            rho  = Rho_d[id * nv + lev];

            // Advection + Coriolis.
            double Cx, Cy, Cz;
            if (DeepModel) {
                Cx = -2.0 * Omega * vy * rho;
                Cy = 2.0 * Omega * vx * rho;
                Cz = 0.0;
            }
            else {
                Cx = 2 * Omega
                     * (vz * func_r_d[id * 3 + 2] * func_r_d[id * 3 + 1]
                        - vy * pow(func_r_d[id * 3 + 2], 2))
                     * rho;
                Cy = -2 * Omega
                     * (vz * func_r_d[id * 3 + 2] * func_r_d[id * 3 + 0]
                        - vx * pow(func_r_d[id * 3 + 2], 2))
                     * rho;
                Cz = 2 * Omega
                     * (vy * func_r_d[id * 3 + 2] * func_r_d[id * 3 + 0]
                        - vx * func_r_d[id * 3 + 2] * func_r_d[id * 3 + 1])
                     * rho;
            }

            Adv_d[id * 3 * nv + lev * 3 + 0] += davx + Cx;
            Adv_d[id * 3 * nv + lev * 3 + 1] += davy + Cy;
            Adv_d[id * 3 * nv + lev * 3 + 2] += davz + Cz;

            if (lev < nv - 1) {
                althl = altht;
                altht = Altitudeh_d[lev + 2];
                altl  = alt;
                alt   = altt;
                if (lev < nv - 2)
                    altt = Altitude_d[lev + 2];
                whl = wht;
                if (lev < nv - 2)
                    wht = Wh_d[id * (nv + 1) + lev + 2];
                vxl = vx;
                vx  = vxt;
                if (lev < nv - 2)
                    vxt = v_d[id * 3 * nv + (lev + 2) * 3 + 0];
                vyl = vy;
                vy  = vyt;
                if (lev < nv - 2)
                    vyt = v_d[id * 3 * nv + (lev + 2) * 3 + 1];
                vzl = vz;
                vz  = vzt;
                if (lev < nv - 2)
                    vzt = v_d[id * 3 * nv + (lev + 2) * 3 + 2];
            }
        }
    }
}
