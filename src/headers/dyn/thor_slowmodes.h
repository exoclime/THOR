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
// Description: Computes the slow modes (see equations 33 to 35 from Mendonca et al. 2016)
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

template<int NX, int NY>
__global__ void Compute_Slow_Modes(double *SlowMh_d,
                                   double *SlowWh_d,
                                   double *SlowRho_d,
                                   double *Slowpressure_d,
                                   double *Mh_d,
                                   double *Wh_d,
                                   double *Rho_d,
                                   double *Adv_d,
                                   double *DivM_d,
                                   double *diffmh_d,
                                   double *diffw_d,
                                   double *diffrh_d,
                                   double *diffpr_d,
                                   double *diffmv_d,
                                   double *diffwv_d,
                                   double *diffrv_d,
                                   double *diffprv_d,
                                   double *pressure_d,
                                   double *h_d,
                                   double *hh_d,
                                   double *gtil_d,
                                   double *grad_d,
                                   double *div_d,
                                   double *Altitude_d,
                                   double *Altitudeh_d,
                                   double  A,
                                   double  Gravit,
                                   double  Cp,
                                   double  Rd,
                                   double *func_r_d,
                                   int *   maps_d,
                                   int     nl_region,
                                   bool    DeepModel,
                                   bool    NonHydro) {

    int x = threadIdx.x;
    int y = threadIdx.y;
    //    int ib  = blockIdx.x;
    int nv  = gridDim.y;
    int lev = blockIdx.y;

    int    pt1, pt2, pt3, pt4, pt5, pt6;
    double div0, div1, div2, div3, div4, div5, div6;
    int    nhl  = nl_region + 2;
    int    nhl2 = nhl * nhl;

    double altht, althl;
    double alt, altl;
    double r2p, r2m, r2l;
    double rscale;

    double advr, advrl;
    double advx, advy, advz;
    double advxl, advyl, advzl;
    double wht, whl;

    double rhol, rho;

    double intt, intl;

    double xi, xim, xip;
    double dwdz, dz, dpr;
    double swr, dpdz, vgp, r;

    double hhl, hht, dwhdz;
    double pressurel;
    double Cv = Cp - Rd;
    double funcx, funcy, funcz;

    /////////////////////////////////////////
    __shared__ double v_s[(NX + 2) * (NY + 2) * 3];
    __shared__ double pressure_s[(NX + 2) * (NY + 2)];
    __shared__ double h_s[(NX + 2) * (NY + 2)];
    __shared__ double nflxr_s[NX * NY];
    __shared__ double nflxv_s[3 * NX * NY];
    __shared__ double nflxp_s[NX * NY];
    /////////////////////////////////////////

    int ir = 0; // index in region
    int iri, ir2, id;


    bool pent_ind = false; //
    int  ig;               // index in global mem

    int igh = 0; // global index in halo

    // Load shared memory
    bool load_halo = compute_mem_idx(maps_d, nhl, nhl2, ig, igh, ir, ir2, pent_ind);
    id             = ig;

    //hack
    // if (ig == 0) {
    //     printf("%d %d %e %e %e %e %e\n",
    //            ig,
    //            lev,
    //            Rho_d[ig * nv + lev],
    //            pressure_d[ig * nv + lev] / Rd / Rho_d[ig * nv + lev],
    //            diffrv_d[ig * nv + lev],
    //            diffprv_d[ig * nv + lev],
    //            diffwv_d[ig * nv + lev]);
    // }

    v_s[ir * 3 + 0] = Mh_d[ig * 3 * nv + lev * 3 + 0];
    v_s[ir * 3 + 1] = Mh_d[ig * 3 * nv + lev * 3 + 1];
    v_s[ir * 3 + 2] = Mh_d[ig * 3 * nv + lev * 3 + 2];

    h_s[ir]        = h_d[ig * nv + lev];
    pressure_s[ir] = pressure_d[ig * nv + lev];

    ///////////////////////////////
    //////////// Halo /////////////
    ///////////////////////////////
    if (load_halo) {
        if (igh >= 0) {
            v_s[ir2 * 3 + 0] = Mh_d[igh * 3 * nv + lev * 3 + 0];
            v_s[ir2 * 3 + 1] = Mh_d[igh * 3 * nv + lev * 3 + 1];
            v_s[ir2 * 3 + 2] = Mh_d[igh * 3 * nv + lev * 3 + 2];
            h_s[ir2]         = h_d[igh * nv + lev];
            pressure_s[ir2]  = pressure_d[igh * nv + lev];
        }
        else {
            v_s[ir2 * 3 + 0] = 0.0;
            v_s[ir2 * 3 + 1] = 0.0;
            v_s[ir2 * 3 + 2] = 0.0;

            h_s[ir2]        = 0.0;
            pressure_s[ir2] = 0.0;
        }
    }
    __syncthreads();
    //////////////////////////////////////////////

    iri = y * nl_region + x;

    funcx = func_r_d[id * 3 + 0];
    funcy = func_r_d[id * 3 + 1];
    funcz = func_r_d[id * 3 + 2];

    advx = Adv_d[id * 3 * nv + lev * 3 + 0];
    advy = Adv_d[id * 3 * nv + lev * 3 + 1];
    advz = Adv_d[id * 3 * nv + lev * 3 + 2];

    rho = Rho_d[id * nv + lev];

    if (lev != 0) {
        advxl     = Adv_d[id * 3 * nv + (lev - 1) * 3 + 0];
        advyl     = Adv_d[id * 3 * nv + (lev - 1) * 3 + 1];
        advzl     = Adv_d[id * 3 * nv + (lev - 1) * 3 + 2];
        advrl     = advxl * funcx + advyl * funcy + advzl * funcz;
        rhol      = Rho_d[id * nv + lev - 1];
        altl      = Altitude_d[lev - 1];
        pressurel = pressure_d[id * nv + lev - 1];
    }

    advr = advx * funcx + advy * funcy + advz * funcz;

    advx += -advr * funcx;
    advy += -advr * funcy;
    advz += -advr * funcz;

    // Neighbours
    pt1 = (y + 2) * nhl + x + 1;
    pt2 = (y + 2) * nhl + x + 2;
    pt3 = (y + 1) * nhl + x + 2;
    pt4 = (y)*nhl + x + 1;
    pt5 = (pent_ind) * ((y + 1) * nhl + x) + (!pent_ind) * ((y)*nhl + x);
    pt6 = (y + 1) * nhl + x;

    altht = Altitudeh_d[lev + 1];
    althl = Altitudeh_d[lev];
    alt   = Altitude_d[lev];

    if (DeepModel) {
        r2p    = pow(altht + A, 2.0);
        r2m    = pow(alt + A, 2.0);
        r2l    = pow(althl + A, 2.0);
        rscale = A / (alt + A);
    }
    else {
        r2p    = 1.0;
        r2m    = 1.0;
        r2l    = 1.0;
        rscale = 1.0;
    }

    nflxr_s[iri] = 0.0;
    nflxp_s[iri] = 0.0;

    for (int k = 0; k < 3; k++) {

        div0 = div_d[id * 7 * 3 + 3 * 0 + k];
        div1 = div_d[id * 7 * 3 + 3 * 1 + k];
        div2 = div_d[id * 7 * 3 + 3 * 2 + k];
        div3 = div_d[id * 7 * 3 + 3 * 3 + k];
        div4 = div_d[id * 7 * 3 + 3 * 4 + k];
        div5 = div_d[id * 7 * 3 + 3 * 5 + k];
        div6 = div_d[id * 7 * 3 + 3 * 6 + k];

        nflxv_s[iri * 3 + k] = rscale
                               * (grad_d[id * 7 * 3 + 3 * 0 + k] * pressure_s[ir]
                                  + grad_d[id * 7 * 3 + 3 * 1 + k] * pressure_s[pt1]
                                  + grad_d[id * 7 * 3 + 3 * 2 + k] * pressure_s[pt2]
                                  + grad_d[id * 7 * 3 + 3 * 3 + k] * pressure_s[pt3]
                                  + grad_d[id * 7 * 3 + 3 * 4 + k] * pressure_s[pt4]
                                  + grad_d[id * 7 * 3 + 3 * 5 + k] * pressure_s[pt5]
                                  + grad_d[id * 7 * 3 + 3 * 6 + k] * pressure_s[pt6]);

        nflxr_s[iri] +=
            rscale
            * (div0 * v_s[ir * 3 + k] + div1 * v_s[pt1 * 3 + k] + div2 * v_s[pt2 * 3 + k]
               + div3 * v_s[pt3 * 3 + k] + div4 * v_s[pt4 * 3 + k] + div5 * v_s[pt5 * 3 + k]
               + div6 * v_s[pt6 * 3 + k]);

        nflxp_s[iri] += rscale
                        * (div0 * v_s[ir * 3 + k] * h_s[ir] + div1 * v_s[pt1 * 3 + k] * h_s[pt1]
                           + div2 * v_s[pt2 * 3 + k] * h_s[pt2] + div3 * v_s[pt3 * 3 + k] * h_s[pt3]
                           + div4 * v_s[pt4 * 3 + k] * h_s[pt4] + div5 * v_s[pt5 * 3 + k] * h_s[pt5]
                           + div6 * v_s[pt6 * 3 + k] * h_s[pt6]);
    }

    if (lev == 0) {
        whl = 0.0;
        wht = Wh_d[id * (nv + 1) + lev + 1];
        hhl = 0.0;
        hht = hh_d[id * (nv + 1) + lev + 1];
    }
    else {
        whl = Wh_d[id * (nv + 1) + lev];
        wht = Wh_d[id * (nv + 1) + lev + 1];
        hhl = hh_d[id * (nv + 1) + lev];
        hht = hh_d[id * (nv + 1) + lev + 1];
    }
    dz    = altht - althl;
    dwdz  = (wht * r2p - whl * r2l) / (dz * r2m);
    dwhdz = (wht * r2p * hht - whl * r2l * hhl) / (dz * r2m);

    //Mh
    dpr =
        nflxv_s[iri * 3 + 0] * funcx + nflxv_s[iri * 3 + 1] * funcy + nflxv_s[iri * 3 + 2] * funcz;

    SlowMh_d[id * 3 * nv + lev * 3 + 0] =
        -(nflxv_s[iri * 3 + 0] - dpr * funcx) - advx + DivM_d[id * 3 * nv + lev * 3 + 0]
        + diffmh_d[id * 3 * nv + lev * 3 + 0] + diffmv_d[id * 3 * nv + lev * 3 + 0];
    SlowMh_d[id * 3 * nv + lev * 3 + 1] =
        -(nflxv_s[iri * 3 + 1] - dpr * funcy) - advy + DivM_d[id * 3 * nv + lev * 3 + 1]
        + diffmh_d[id * 3 * nv + lev * 3 + 1] + diffmv_d[id * 3 * nv + lev * 3 + 1];
    SlowMh_d[id * 3 * nv + lev * 3 + 2] =
        -(nflxv_s[iri * 3 + 2] - dpr * funcz) - advz + DivM_d[id * 3 * nv + lev * 3 + 2]
        + diffmh_d[id * 3 * nv + lev * 3 + 2] + diffmv_d[id * 3 * nv + lev * 3 + 2];

    // Wh
    if (lev == 0) {
        SlowWh_d[id * (nv + 1) + 0]  = 0.0;
        SlowWh_d[id * (nv + 1) + nv] = 0.0;
    }
    else {
        xi  = althl;
        xim = altl;
        xip = alt;

        intt = (xi - xip) / (xim - xip);
        intl = (xi - xim) / (xip - xim);

        if (NonHydro)
            swr =
                (-diffw_d[id * nv + lev - 1] - diffwv_d[id * nv + lev - 1] + advrl + rhol * Gravit)
                    * intt
                + (-diffw_d[id * nv + lev] - diffwv_d[id * nv + lev] + advr + rho * Gravit) * intl;
        else
            swr = (-diffw_d[id * nv + lev - 1] - diffwv_d[id * nv + lev - 1] + rhol * Gravit) * intt
                  + (-diffw_d[id * nv + lev] - diffwv_d[id * nv + lev] + rho * Gravit) * intl;

        dz   = alt - altl;
        dpdz = (pressure_s[ir] - pressurel) / dz;
        swr += dpdz;

        SlowWh_d[id * (nv + 1) + lev] = -swr;
    }

    // Rho
    nflxr_s[iri] += dwdz;
    SlowRho_d[id * nv + lev] = -nflxr_s[iri] + diffrh_d[id * nv + lev] + diffrv_d[id * nv + lev];

    // pressure
    r   = 1.0 / rho;
    vgp = (nflxv_s[iri * 3 + 0] - dpr * funcx) * v_s[ir * 3 + 0] * r
          + (nflxv_s[iri * 3 + 1] - dpr * funcy) * v_s[ir * 3 + 1] * r
          + (nflxv_s[iri * 3 + 2] - dpr * funcz) * v_s[ir * 3 + 2] * r;

    Slowpressure_d[id * nv + lev] =
        (Rd / Cv) * (-nflxp_s[iri] - gtil_d[id * nv + lev] - dwhdz + vgp) + diffpr_d[id * nv + lev]
        + diffprv_d[id * nv + lev];

    // if (diffprv_d[id * nv + lev] > 0) {
    //     printf("hi!\n"); //
    // }
}

template<int NN>
__global__ void Compute_Slow_Modes_Poles(double *SlowMh_d,
                                         double *SlowWh_d,
                                         double *SlowRho_d,
                                         double *Slowpressure_d,
                                         double *Mh_d,
                                         double *Wh_d,
                                         double *Rho_d,
                                         double *Adv_d,
                                         double *DivM_d,
                                         double *diffmh_d,
                                         double *diffw_d,
                                         double *diffrh_d,
                                         double *diffpr_d,
                                         double *diffmv_d,
                                         double *diffwv_d,
                                         double *diffrv_d,
                                         double *diffprv_d,
                                         double *pressure_d,
                                         double *h_d,
                                         double *hh_d,
                                         double *gtil_d,
                                         double *grad_d,
                                         double *div_d,
                                         double *Altitude_d,
                                         double *Altitudeh_d,
                                         double  A,
                                         double  Gravit,
                                         double  Cp,
                                         double  Rd,
                                         double *func_r_d,
                                         int *   point_local_d,
                                         int     nv,
                                         int     num,
                                         bool    DeepModel,
                                         bool    NonHydro) {

    int id = blockIdx.x * blockDim.x + threadIdx.x;
    id += num - 2; // Poles

    double altht, althl;
    double alt, altl;
    double r2p, r2m, r2l;
    double rscale;

    double advrt, advrl;
    double advx, advy, advz;
    double wht, whl;

    double rhol, rhot;

    double intt, intl;

    double xi, xim, xip;
    double dwdz, dz, dpr, r;
    double swr, dpdz, vgp;

    double hhl, hht, dwhdz;
    double Cv = Cp - Rd;
    double pressurel;

    /////////////////////////////////////////

    __shared__ double v_p[NN * 3];
    __shared__ double pressure_p[NN];
    __shared__ double h_p[NN];
    __shared__ double func_r_p[3];
    __shared__ double div_p[3 * 7];
    __shared__ double grad_p[3 * 7];
    __shared__ double nflxv_p[3];
    __shared__ int    local_p[6];

    double nflxr_p;
    double nflxp_p;

    /////////////////////////////////////////
    if (id < num) {
        for (int i = 0; i < 5; i++) local_p[i] = point_local_d[id * 6 + i];
        func_r_p[0] = func_r_d[id * 3 + 0];
        func_r_p[1] = func_r_d[id * 3 + 1];
        func_r_p[2] = func_r_d[id * 3 + 2];
        for (int i = 0; i < 7; i++)
            for (int k = 0; k < 3; k++) div_p[i * 3 + k] = div_d[id * 7 * 3 + i * 3 + k];
        for (int i = 0; i < 7; i++)
            for (int k = 0; k < 3; k++) grad_p[i * 3 + k] = grad_d[id * 7 * 3 + i * 3 + k];

        for (int lev = 0; lev < nv; lev++) {

            v_p[0]        = Mh_d[id * 3 * nv + lev * 3 + 0];
            v_p[1]        = Mh_d[id * 3 * nv + lev * 3 + 1];
            v_p[2]        = Mh_d[id * 3 * nv + lev * 3 + 2];
            pressure_p[0] = pressure_d[id * nv + lev];
            h_p[0]        = h_d[id * nv + lev];
            for (int i = 1; i < 6; i++) {
                v_p[i * 3 + 0] = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 0];
                v_p[i * 3 + 1] = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 1];
                v_p[i * 3 + 2] = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 2];

                pressure_p[i] = pressure_d[local_p[i - 1] * nv + lev];
                h_p[i]        = h_d[local_p[i - 1] * nv + lev];
            }

            advx = Adv_d[id * 3 * nv + lev * 3 + 0];
            advy = Adv_d[id * 3 * nv + lev * 3 + 1];
            advz = Adv_d[id * 3 * nv + lev * 3 + 2];

            if (lev > 0) advrl = advrt;
            advrt = advx * func_r_p[0] + advy * func_r_p[1] + advz * func_r_p[2];

            advx += -advrt * func_r_p[0];
            advy += -advrt * func_r_p[1];
            advz += -advrt * func_r_p[2];

            if (lev == 0) {
                altht = Altitudeh_d[lev + 1];
                althl = Altitudeh_d[lev];
            }
            if (lev > 0) altl = alt;
            alt = Altitude_d[lev];

            if (lev > 0) rhol = rhot;
            rhot = Rho_d[id * nv + lev];

            if (DeepModel) {
                r2p    = pow(altht + A, 2.0);
                r2m    = pow(alt + A, 2.0);
                r2l    = pow(althl + A, 2.0);
                rscale = A / (alt + A);
            }
            else {
                r2p    = 1.0;
                r2m    = 1.0;
                r2l    = 1.0;
                rscale = 1.0;
            }

            nflxr_p = 0.0;
            nflxp_p = 0.0;

            for (int k = 0; k < 3; k++) {
                nflxv_p[k] =
                    rscale
                    * (grad_p[3 * 0 + k] * pressure_p[0] + grad_p[3 * 1 + k] * pressure_p[1]
                       + grad_p[3 * 2 + k] * pressure_p[2] + grad_p[3 * 3 + k] * pressure_p[3]
                       + grad_p[3 * 4 + k] * pressure_p[4] + grad_p[3 * 5 + k] * pressure_p[5]);

                nflxr_p +=
                    rscale
                    * (div_p[3 * 0 + k] * v_p[0 * 3 + k] + div_p[3 * 1 + k] * v_p[1 * 3 + k]
                       + div_p[3 * 2 + k] * v_p[2 * 3 + k] + div_p[3 * 3 + k] * v_p[3 * 3 + k]
                       + div_p[3 * 4 + k] * v_p[4 * 3 + k] + div_p[3 * 5 + k] * v_p[5 * 3 + k]);

                nflxp_p += rscale
                           * (div_p[3 * 0 + k] * v_p[0 * 3 + k] * h_p[0]
                              + div_p[3 * 1 + k] * v_p[1 * 3 + k] * h_p[1]
                              + div_p[3 * 2 + k] * v_p[2 * 3 + k] * h_p[2]
                              + div_p[3 * 3 + k] * v_p[3 * 3 + k] * h_p[3]
                              + div_p[3 * 4 + k] * v_p[4 * 3 + k] * h_p[4]
                              + div_p[3 * 5 + k] * v_p[5 * 3 + k] * h_p[5]);
            }

            if (lev == 0) {
                whl = 0.0;
                wht = Wh_d[id * (nv + 1) + lev + 1];
                hhl = 0.0;
                hht = hh_d[id * (nv + 1) + lev + 1];
            }

            dz    = altht - althl;
            dwdz  = (wht * r2p - whl * r2l) / (dz * r2m);
            dwhdz = (wht * r2p * hht - whl * r2l * hhl) / (dz * r2m);

            // Mh
            dpr = nflxv_p[0] * func_r_p[0] + nflxv_p[1] * func_r_p[1] + nflxv_p[2] * func_r_p[2];

            SlowMh_d[id * 3 * nv + lev * 3 + 0] =
                -(nflxv_p[0] - dpr * func_r_p[0]) - advx + DivM_d[id * 3 * nv + lev * 3 + 0]
                + diffmh_d[id * 3 * nv + lev * 3 + 0] + diffmv_d[id * 3 * nv + lev * 3 + 0];
            SlowMh_d[id * 3 * nv + lev * 3 + 1] =
                -(nflxv_p[1] - dpr * func_r_p[1]) - advy + DivM_d[id * 3 * nv + lev * 3 + 1]
                + diffmh_d[id * 3 * nv + lev * 3 + 1] + diffmv_d[id * 3 * nv + lev * 3 + 1];
            SlowMh_d[id * 3 * nv + lev * 3 + 2] =
                -(nflxv_p[2] - dpr * func_r_p[2]) - advz + DivM_d[id * 3 * nv + lev * 3 + 2]
                + diffmh_d[id * 3 * nv + lev * 3 + 2] + diffmv_d[id * 3 * nv + lev * 3 + 2];
            // Wh
            if (lev == 0) {
                SlowWh_d[id * (nv + 1)]      = 0.0;
                SlowWh_d[id * (nv + 1) + nv] = 0.0;
            }
            else {
                xi  = althl;
                xim = altl;
                xip = alt;

                intt = (xi - xip) / (xim - xip);
                intl = (xi - xim) / (xip - xim);
                if (NonHydro)
                    swr = (-diffw_d[id * nv + lev - 1] - diffwv_d[id * nv + lev - 1] + advrl
                           + rhol * Gravit)
                              * intt
                          + (-diffw_d[id * nv + lev] - diffwv_d[id * nv + lev] + advrt
                             + rhot * Gravit)
                                * intl;
                else
                    swr =
                        (-diffw_d[id * nv + lev - 1] - diffwv_d[id * nv + lev - 1] + rhol * Gravit)
                            * intt
                        + (-diffw_d[id * nv + lev] - diffwv_d[id * nv + lev] + rhot * Gravit)
                              * intl;

                dz   = alt - altl;
                dpdz = (pressure_p[0] - pressurel) / dz;
                swr  = swr + dpdz;

                SlowWh_d[id * (nv + 1) + lev] = -swr;
            }

            // Rho
            nflxr_p += dwdz;
            SlowRho_d[id * nv + lev] = -nflxr_p + diffrh_d[id * nv + lev] + diffrv_d[id * nv + lev];

            // pressure
            r   = 1.0 / rhot;
            vgp = (nflxv_p[0] - dpr * func_r_p[0]) * v_p[0] * r
                  + (nflxv_p[1] - dpr * func_r_p[1]) * v_p[1] * r
                  + (nflxv_p[2] - dpr * func_r_p[2]) * v_p[2] * r;

            Slowpressure_d[id * nv + lev] =
                (Rd / Cv) * (-nflxp_p - gtil_d[id * nv + lev] - dwhdz + vgp)
                + diffpr_d[id * nv + lev] + diffprv_d[id * nv + lev];

            if (lev < nv - 1) {
                pressurel = pressure_p[0];

                pressure_p[0] = pressure_d[id * nv + lev + 1];
                for (int i = 1; i < 6; i++)
                    pressure_p[i] = pressure_d[local_p[i - 1] * nv + lev + 1];
            }

            if (lev != nv - 1) {
                althl = altht;
                altht = Altitudeh_d[lev + 2];
                whl   = wht;
                wht   = Wh_d[id * (nv + 1) + lev + 2];
                hhl   = hht;
                hht   = hh_d[id * (nv + 1) + lev + 2];
            }
        }
        SlowWh_d[id * (nv + 1) + nv] = 0.0;
    }
}
