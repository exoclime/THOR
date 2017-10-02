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
// Current Code Owner: Joao Mendonca, EEG. joao.mendonca@csh.unibe.ch
//
// History:
// Version Date       Comment
// ======= ====       =======
//
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////
// Computes the 3D advection terms and the velocities in cartesians.
// (no poles)
template <int NX, int NY>
__global__ void Compute_Advec_Cori1(double * Adv_d     , // Advection term.
                                    double * v_d       , // 3D velocity (cartesians).
                                    double * Mh_d      , // Horizontal momentum.
                                    double * div_d     , // Divergence operator.
                                    double * W_d       , // Vertical momentum.
                                    double * Rho_d     , // Density.
                                    double * Altitude_d, // Altitude.
                                    double   A         , // Radius.
                                    double * func_r_d  , // Unit vectors normal to the spherical surface.
                                    int * maps_d       , // Global indexes. 
                                    int nl_region      , // Length of the side of each rhombi (from the sphere decomposition).
                                    bool DeepModel     ){// Switches on and off the deep atmosphere solution (see headers/define.h).

    int x  = threadIdx.x;
    int y  = threadIdx.y;
    int ib = blockIdx.x;
    int nv = gridDim.y ;
    int lev= blockIdx.y;

    int pt1, pt2, pt3, pt4, pt5, pt6;
    double div0, div1, div2, div3, div4, div5, div6;
    int nhl = nl_region + 2;
    int nhl2 = nhl*nhl;

    double alt;
    double rho;
    double rscale;

    double w;

    int ir = (y + 1)*nhl + x + 1;   // Region index
    int iri, ir2, twot, id;

    /////////////////////////////////////////
    __shared__ double v_s[(NX + 2)*(NY + 2) * 3];
    __shared__ double M_s[(NX + 2)*(NY + 2) * 3];
    __shared__ double nflxx_s[NX*NY];
    __shared__ double nflxy_s[NX*NY];
    __shared__ double nflxz_s[NX*NY];
    /////////////////////////////////////////

    bool pent_ind = 0; // 
    int ig;

    // Load shared memory
    ig = maps_d[ib*nhl2 + ir];
    id = ig;
    if (x == 0 && y == 0) if (maps_d[ib * nhl2] == -1) pent_ind = 1;

    v_s[ir * 3 + 0] = Mh_d[ig * 3 * nv + lev * 3 + 0];
    v_s[ir * 3 + 1] = Mh_d[ig * 3 * nv + lev * 3 + 1];
    v_s[ir * 3 + 2] = Mh_d[ig * 3 * nv + lev * 3 + 2];
    M_s[ir * 3 + 0] = v_s[ir * 3 + 0];
    M_s[ir * 3 + 1] = v_s[ir * 3 + 1];
    M_s[ir * 3 + 2] = v_s[ir * 3 + 2];
    rho = Rho_d[ig*nv + lev];
    w   = W_d[ig*nv + lev];

    v_s[ir * 3 + 0] = v_s[ir * 3 + 0] / rho + (w / rho) * func_r_d[ig * 3 + 0];
    v_s[ir * 3 + 1] = v_s[ir * 3 + 1] / rho + (w / rho) * func_r_d[ig * 3 + 1];
    v_s[ir * 3 + 2] = v_s[ir * 3 + 2] / rho + (w / rho) * func_r_d[ig * 3 + 2];

    ///////////////////////////////
    //////////// Halo /////////////
    ///////////////////////////////
    if (x == 0) {
        ir2 = (y + 1) * nhl + x;
        ig = maps_d[ib * nhl2 + ir2];
        v_s[ir2 * 3 + 0] = Mh_d[ig * 3 * nv + lev * 3 + 0];
        v_s[ir2 * 3 + 1] = Mh_d[ig * 3 * nv + lev * 3 + 1];
        v_s[ir2 * 3 + 2] = Mh_d[ig * 3 * nv + lev * 3 + 2];
        M_s[ir2 * 3 + 0] = v_s[ir2 * 3 + 0];
        M_s[ir2 * 3 + 1] = v_s[ir2 * 3 + 1];
        M_s[ir2 * 3 + 2] = v_s[ir2 * 3 + 2];
        rho = Rho_d[ig*nv + lev];
        w = W_d[ig*nv + lev];
        v_s[ir2 * 3 + 0] = v_s[ir2 * 3 + 0] / rho + (w / rho) * func_r_d[ig * 3 + 0];
        v_s[ir2 * 3 + 1] = v_s[ir2 * 3 + 1] / rho + (w / rho) * func_r_d[ig * 3 + 1];
        v_s[ir2 * 3 + 2] = v_s[ir2 * 3 + 2] / rho + (w / rho) * func_r_d[ig * 3 + 2];
    }
    if (x == nhl - 3){
        ir2 = (y + 1) * nhl + x + 2;
        ig = maps_d[ib * nhl2 + ir2];
        v_s[ir2 * 3 + 0] = Mh_d[ig * 3 * nv + lev * 3 + 0];
        v_s[ir2 * 3 + 1] = Mh_d[ig * 3 * nv + lev * 3 + 1];
        v_s[ir2 * 3 + 2] = Mh_d[ig * 3 * nv + lev * 3 + 2];
        M_s[ir2 * 3 + 0] = v_s[ir2 * 3 + 0];
        M_s[ir2 * 3 + 1] = v_s[ir2 * 3 + 1];
        M_s[ir2 * 3 + 2] = v_s[ir2 * 3 + 2];
        rho = Rho_d[ig*nv + lev];
        w = W_d[ig*nv + lev];
        v_s[ir2 * 3 + 0] = v_s[ir2 * 3 + 0] / rho + (w / rho) * func_r_d[ig * 3 + 0];
        v_s[ir2 * 3 + 1] = v_s[ir2 * 3 + 1] / rho + (w / rho) * func_r_d[ig * 3 + 1];
        v_s[ir2 * 3 + 2] = v_s[ir2 * 3 + 2] / rho + (w / rho) * func_r_d[ig * 3 + 2];
    }
    if (y == 0){
        twot = 1;
        ir2 = y * nhl + (x + 1);
        if (x == 0) twot = 2;

        for (int k = 0; k < twot; k++){
            if (k == 1) ir2 = y * nhl + x;
            ig = maps_d[ib * nhl2 + ir2];
            if (ig >= 0){
                v_s[ir2 * 3 + 0] = Mh_d[ig * 3 * nv + lev * 3 + 0];
                v_s[ir2 * 3 + 1] = Mh_d[ig * 3 * nv + lev * 3 + 1];
                v_s[ir2 * 3 + 2] = Mh_d[ig * 3 * nv + lev * 3 + 2];
                M_s[ir2 * 3 + 0] = v_s[ir2 * 3 + 0];
                M_s[ir2 * 3 + 1] = v_s[ir2 * 3 + 1];
                M_s[ir2 * 3 + 2] = v_s[ir2 * 3 + 2];
                rho = Rho_d[ig * nv + lev];
                w = W_d[ig * nv + lev];
                v_s[ir2 * 3 + 0] = v_s[ir2 * 3 + 0] / rho + (w / rho) * func_r_d[ig * 3 + 0];
                v_s[ir2 * 3 + 1] = v_s[ir2 * 3 + 1] / rho + (w / rho) * func_r_d[ig * 3 + 1];
                v_s[ir2 * 3 + 2] = v_s[ir2 * 3 + 2] / rho + (w / rho) * func_r_d[ig * 3 + 2];
            }
            else{
                v_s[ir2 * 3 + 0] = 0.0;
                v_s[ir2 * 3 + 1] = 0.0;
                v_s[ir2 * 3 + 2] = 0.0;
                M_s[ir2 * 3 + 0] = 0.0;
                M_s[ir2 * 3 + 1] = 0.0;
                M_s[ir2 * 3 + 2] = 0.0;
            }
        }
    }
    if (y == nhl - 3) {
        twot = 1;
        ir2 = (y + 2) * nhl + (x + 1);
        if (x == nhl - 3) twot = 2;
        for (int k = 0; k < twot; k++){
            if (k == 1) ir2 = (y + 2) * nhl + (x + 2);
            ig = maps_d[ib * nhl2 + ir2];
            v_s[ir2 * 3 + 0] = Mh_d[ig * 3 * nv + lev * 3 + 0];
            v_s[ir2 * 3 + 1] = Mh_d[ig * 3 * nv + lev * 3 + 1];
            v_s[ir2 * 3 + 2] = Mh_d[ig * 3 * nv + lev * 3 + 2];
            M_s[ir2 * 3 + 0] = v_s[ir2 * 3 + 0];
            M_s[ir2 * 3 + 1] = v_s[ir2 * 3 + 1];
            M_s[ir2 * 3 + 2] = v_s[ir2 * 3 + 2];
            rho = Rho_d[ig * nv + lev];
            w = W_d[ig * nv + lev];
            v_s[ir2 * 3 + 0] = v_s[ir2 * 3 + 0] / rho + (w / rho) * func_r_d[ig * 3 + 0];
            v_s[ir2 * 3 + 1] = v_s[ir2 * 3 + 1] / rho + (w / rho) * func_r_d[ig * 3 + 1];
            v_s[ir2 * 3 + 2] = v_s[ir2 * 3 + 2] / rho + (w / rho) * func_r_d[ig * 3 + 2];
        }
    }
    __syncthreads();

    //////////////////////////////////////////////
    // Inner index
    iri = y * nl_region + x;
    
    // Neighbours.
    pt1 = (y + 2)*nhl + x + 1;
    pt2 = (y + 2)*nhl + x + 2;
    pt3 = (y + 1)*nhl + x + 2;
    pt4 = (y    )*nhl + x + 1;
    pt5 = pent_ind*((y + 1)*nhl + x) + (!pent_ind)*((y   )*nhl + x);
    pt6 = (y + 1)*nhl + x    ;

    // Sets the weight needed for the deep atmosphere solution.
    if (DeepModel){
        alt = Altitude_d[lev];
        rscale = A / (alt + A);
    }
    else rscale = 1.0;

    // Initialize fluxes.
    nflxx_s[iri] = 0.0;
    nflxy_s[iri] = 0.0;
    nflxz_s[iri] = 0.0;

    // Calculate divergence.
    for (int k = 0; k < 3; k++){
        div0 = div_d[id * 7 * 3 + 3 * 0 + k];
        div1 = div_d[id * 7 * 3 + 3 * 1 + k];
        div2 = div_d[id * 7 * 3 + 3 * 2 + k];
        div3 = div_d[id * 7 * 3 + 3 * 3 + k];
        div4 = div_d[id * 7 * 3 + 3 * 4 + k];
        div5 = div_d[id * 7 * 3 + 3 * 5 + k];
        div6 = div_d[id * 7 * 3 + 3 * 6 + k]; // For pent_ind = 1 div6 is equal to 0.

        nflxx_s[iri] += rscale*(div0 * v_s[ir * 3 + 0]  * M_s[ir * 3 + k] +
                                div1 * v_s[pt1 * 3 + 0] * M_s[pt1 * 3 + k] +
                                div2 * v_s[pt2 * 3 + 0] * M_s[pt2 * 3 + k] +
                                div3 * v_s[pt3 * 3 + 0] * M_s[pt3 * 3 + k] +
                                div4 * v_s[pt4 * 3 + 0] * M_s[pt4 * 3 + k] +
                                div5 * v_s[pt5 * 3 + 0] * M_s[pt5 * 3 + k] +
                                div6 * v_s[pt6 * 3 + 0] * M_s[pt6 * 3 + k]);

        nflxy_s[iri] += rscale*(div0 * v_s[ir * 3 + 1]  * M_s[ir * 3 + k] +
                                div1 * v_s[pt1 * 3 + 1] * M_s[pt1 * 3 + k] +
                                div2 * v_s[pt2 * 3 + 1] * M_s[pt2 * 3 + k] +
                                div3 * v_s[pt3 * 3 + 1] * M_s[pt3 * 3 + k] +
                                div4 * v_s[pt4 * 3 + 1] * M_s[pt4 * 3 + k] +
                                div5 * v_s[pt5 * 3 + 1] * M_s[pt5 * 3 + k] +
                                div6 * v_s[pt6 * 3 + 1] * M_s[pt6 * 3 + k]);

        nflxz_s[iri] += rscale*(div0 * v_s[ir * 3 + 2]  * M_s[ir * 3 + k]  +
                                div1 * v_s[pt1 * 3 + 2] * M_s[pt1 * 3 + k] +
                                div2 * v_s[pt2 * 3 + 2] * M_s[pt2 * 3 + k] +
                                div3 * v_s[pt3 * 3 + 2] * M_s[pt3 * 3 + k] +
                                div4 * v_s[pt4 * 3 + 2] * M_s[pt4 * 3 + k] +
                                div5 * v_s[pt5 * 3 + 2] * M_s[pt5 * 3 + k] +
                                div6 * v_s[pt6 * 3 + 2] * M_s[pt6 * 3 + k]);
    }
    // Return values (3D advection term and velocities in cartesians).
    Adv_d[id * 3 * nv + lev * 3 + 0] = nflxx_s[iri];
    Adv_d[id * 3 * nv + lev * 3 + 1] = nflxy_s[iri];
    Adv_d[id * 3 * nv + lev * 3 + 2] = nflxz_s[iri];
    v_d[id * 3 * nv + lev * 3 + 0] = v_s[ir * 3 + 0];
    v_d[id * 3 * nv + lev * 3 + 1] = v_s[ir * 3 + 1];
    v_d[id * 3 * nv + lev * 3 + 2] = v_s[ir * 3 + 2];
}

// Computes the 3D advection terms and the velocities in cartesians.
// (poles)
template <int NN>
__global__ void Compute_Advec_Cori_Poles(double * Adv_d        , //
                                         double * v_d          , //
                                         double * Mh_d         , //
                                         double * div_d        , //
                                         double * W_d          , //
                                         double * Rho_d        , //
                                         double * Altitude_d   , //
                                         double   A            , //
                                         double * func_r_d     , //
                                         int    * point_local_d, //
                                         int      num          , //
                                         int      nv           , //
                                         bool     DeepModel    ){//

    int id = blockIdx.x * blockDim.x + threadIdx.x;
    id += num - 2;     // Poles

    double alt, w;
    double rscale;

    double rho;
    
    /////////////////////////////////////////
    __shared__ double v_p[NN * 3];
    __shared__ double M_p[NN * 3];
    __shared__ double func_r_p[3 * NN];
    __shared__ double div_p[3 * 7];
    __shared__ int local_p[NN];

    double nflxx;
    double nflxy;
    double nflxz;

    /////////////////////////////////////////
    if (id < num){
        for (int i = 0; i < 5; i++)local_p[i] = point_local_d[id * 6 + i];
        func_r_p[0] = func_r_d[id * 3 + 0];
        func_r_p[1] = func_r_d[id * 3 + 1];
        func_r_p[2] = func_r_d[id * 3 + 2];
        for (int i = 1; i < 6; i++){
            func_r_p[i * 3 + 0] = func_r_d[local_p[i - 1] * 3 + 0];
            func_r_p[i * 3 + 1] = func_r_d[local_p[i - 1] * 3 + 1];
            func_r_p[i * 3 + 2] = func_r_d[local_p[i - 1] * 3 + 2];
        }
        for (int i = 0; i < 7; i++) for (int k = 0; k < 3; k++) div_p[i * 3 + k] = div_d[id * 7 * 3 + i * 3 + k];

        for (int lev = 0; lev < nv; lev++){
            v_p[0] = Mh_d[id * 3 * nv + lev * 3 + 0];
            v_p[1] = Mh_d[id * 3 * nv + lev * 3 + 1];
            v_p[2] = Mh_d[id * 3 * nv + lev * 3 + 2];
            M_p[0] = v_p[0];
            M_p[1] = v_p[1];
            M_p[2] = v_p[2];
            rho = Rho_d[id*nv + lev];
            w = W_d[id*nv + lev];
            v_p[0] = v_p[0] / rho + (w / rho) * func_r_p[0];
            v_p[1] = v_p[1] / rho + (w / rho) * func_r_p[1];
            v_p[2] = v_p[2] / rho + (w / rho) * func_r_p[2];
            for (int i = 1; i < 6; i++){
                v_p[i * 3 + 0] = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 0];
                v_p[i * 3 + 1] = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 1];
                v_p[i * 3 + 2] = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 2];
                M_p[i * 3 + 0] = v_p[i * 3 + 0];
                M_p[i * 3 + 1] = v_p[i * 3 + 1];
                M_p[i * 3 + 2] = v_p[i * 3 + 2];
                rho = Rho_d[local_p[i - 1] * nv + lev];
                w = W_d[local_p[i - 1] * nv + lev];

                v_p[i * 3 + 0] = v_p[i * 3 + 0] / rho + (w / rho) * func_r_p[i * 3 + 0];
                v_p[i * 3 + 1] = v_p[i * 3 + 1] / rho + (w / rho) * func_r_p[i * 3 + 1];
                v_p[i * 3 + 2] = v_p[i * 3 + 2] / rho + (w / rho) * func_r_p[i * 3 + 2];
            }

            if (DeepModel){
                alt = Altitude_d[lev];
                rscale = A / (alt + A);
            }
            else rscale = 1.0;

            nflxx = 0.0;
            nflxy = 0.0;
            nflxz = 0.0;

            for (int k = 0; k < 3; k++){
                nflxx   += rscale*(div_p[3 * 0 + k] * v_p[0 * 3 + 0] * M_p[0 * 3 + k] +
                                   div_p[3 * 1 + k] * v_p[1 * 3 + 0] * M_p[1 * 3 + k] +
                                   div_p[3 * 2 + k] * v_p[2 * 3 + 0] * M_p[2 * 3 + k] +
                                   div_p[3 * 3 + k] * v_p[3 * 3 + 0] * M_p[3 * 3 + k] +
                                   div_p[3 * 4 + k] * v_p[4 * 3 + 0] * M_p[4 * 3 + k] +
                                   div_p[3 * 5 + k] * v_p[5 * 3 + 0] * M_p[5 * 3 + k]);

                nflxy   += rscale*(div_p[3 * 0 + k] * v_p[0 * 3 + 1] * M_p[0 * 3 + k] +
                                   div_p[3 * 1 + k] * v_p[1 * 3 + 1] * M_p[1 * 3 + k] +
                                   div_p[3 * 2 + k] * v_p[2 * 3 + 1] * M_p[2 * 3 + k] +
                                   div_p[3 * 3 + k] * v_p[3 * 3 + 1] * M_p[3 * 3 + k] +
                                   div_p[3 * 4 + k] * v_p[4 * 3 + 1] * M_p[4 * 3 + k] +
                                   div_p[3 * 5 + k] * v_p[5 * 3 + 1] * M_p[5 * 3 + k]);

                nflxz   += rscale*(div_p[3 * 0 + k] * v_p[0 * 3 + 2] * M_p[0 * 3 + k] +
                                   div_p[3 * 1 + k] * v_p[1 * 3 + 2] * M_p[1 * 3 + k] +
                                   div_p[3 * 2 + k] * v_p[2 * 3 + 2] * M_p[2 * 3 + k] +
                                   div_p[3 * 3 + k] * v_p[3 * 3 + 2] * M_p[3 * 3 + k] +
                                   div_p[3 * 4 + k] * v_p[4 * 3 + 2] * M_p[4 * 3 + k] +
                                   div_p[3 * 5 + k] * v_p[5 * 3 + 2] * M_p[5 * 3 + k]);
            }
            // Return values (3D advection term and velocities in cartesians).
            Adv_d[id * 3 * nv + lev * 3 + 0] = nflxx;
            Adv_d[id * 3 * nv + lev * 3 + 1] = nflxy;
            Adv_d[id * 3 * nv + lev * 3 + 2] = nflxz;
            v_d[id * 3 * nv + lev * 3 + 0] = v_p[0];
            v_d[id * 3 * nv + lev * 3 + 1] = v_p[1];
            v_d[id * 3 * nv + lev * 3 + 2] = v_p[2];
        }
    }
}

// Projects the advection terms on the spherical surface and
// also adds the vertical contribution.
__global__ void Compute_Advec_Cori2(double * Adv_d      , //
                                    double * v_d        , //
                                    double * Wh_d       , //
                                    double * Rho_d      , //
                                    double * Altitude_d , //
                                    double * Altitudeh_d, //
                                    double   Omega      , //
                                    double   A          , //
                                    int      nv         , //
                                    int      num        , //
                                    bool     DeepModel  ){// 

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
    if (id < num){
        for (int lev = 0; lev < nv; lev++){
            if (lev == 0){
                dvwxl = 0.0; dvwyl = 0.0; dvwzl = 0.0;

                altht = Altitudeh_d[lev + 1];
                althl = Altitudeh_d[lev];
                alt   = Altitude_d[lev];
                altt  = Altitude_d[lev + 1];
                vx = v_d[id * 3 * nv + lev * 3 + 0];
                vxt= v_d[id * 3 * nv + (lev+1) * 3 + 0];
                vy = v_d[id * 3 * nv + lev * 3 + 1];
                vyt= v_d[id * 3 * nv + (lev + 1) * 3 + 1];
                vz = v_d[id * 3 * nv + lev * 3 + 2];
                vzt= v_d[id * 3 * nv + (lev + 1) * 3 + 2];
                wht = Wh_d[id*(nv + 1) + lev + 1];

                // Linear interpolation.
                xi = altht;
                xim = alt;
                xip = altt;

                intt = (xi - xip) / (xim - xip);
                intl = (xi - xim) / (xip - xim);
                
                dvwxt = (vx * intt + vxt * intl) * wht;
                dvwyt = (vy * intt + vyt * intl) * wht;
                dvwzt = (vz * intt + vzt * intl) * wht;
            }
            else if (lev == nv - 1){
                dvwxt = 0.0; dvwyt = 0.0; dvwzt = 0.0;

                xi = althl;
                xim = altl;
                xip = alt;

                intt = (xi - xip) / (xim - xip);
                intl = (xi - xim) / (xip - xim);

                dvwxl = (vxl * intt + vx * intl) * whl;
                dvwyl = (vyl * intt + vy * intl) * whl;
                dvwzl = (vzl * intt + vz * intl) * whl;
            }
            else{
                xi = althl;
                xim = altl;
                xip = alt;

                intt = (xi - xip) / (xim - xip);
                intl = (xi - xim) / (xip - xim);

                dvwxl = (vxl * intt + vx * intl) * whl;
                dvwyl = (vyl * intt + vy * intl) * whl;
                dvwzl = (vzl * intt + vz * intl) * whl;

                xi = altht;
                xim = alt;
                xip = altt;

                intt = (xi - xip) / (xim - xip);
                intl = (xi - xim) / (xip - xim);

                dvwxt = (vx * intt + vxt * intl) * wht;
                dvwyt = (vy * intt + vyt * intl) * wht;
                dvwzt = (vz * intt + vzt * intl) * wht;
            }

            if (DeepModel){
                r2p = pow(altht + A, 2.0);
                r2m = pow(alt + A, 2.0);
                r2l = pow(althl + A, 2.0);
            }
            else{
                r2p = 1.0;
                r2m = 1.0;
                r2l = 1.0;
            }

            dz = 1.0 / ((althl - altht)*r2m);
            davx = (dvwxl*r2l - dvwxt*r2p) * dz;
            davy = (dvwyl*r2l - dvwyt*r2p) * dz;
            davz = (dvwzl*r2l - dvwzt*r2p) * dz;
            rho  = Rho_d[id * nv + lev];

            // Advection + Coriolis.
            Adv_d[id * 3 * nv + lev * 3 + 0] += davx - 2.0 * Omega * vy * rho;
            Adv_d[id * 3 * nv + lev * 3 + 1] += davy + 2.0 * Omega * vx * rho;
            Adv_d[id * 3 * nv + lev * 3 + 2] += davz;

            if (lev < nv - 1){
                althl= altht;
                altht= Altitudeh_d[lev + 2];
                altl = alt;
                alt  = altt;
                if(lev < nv - 2) altt = Altitude_d[lev + 2];
                whl  = wht;
                if (lev < nv - 2) wht = Wh_d[id*(nv + 1) + lev + 2];
                vxl  = vx;
                vx   = vxt;
                if (lev < nv - 2) vxt = v_d[id * 3 * nv + (lev + 2) * 3 + 0];
                vyl  = vy;
                vy   = vyt;
                if (lev < nv - 2) vyt = v_d[id * 3 * nv + (lev + 2) * 3 + 1];
                vzl  = vz;
                vz   = vzt;
                if (lev < nv - 2) vzt = v_d[id * 3 * nv + (lev + 2) * 3 + 2];
            }
        }
    }
}
