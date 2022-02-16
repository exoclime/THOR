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
// Description: Computes the fast modes (see equations 33 to 35 from Mendonca et al. 2016)
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
#include "kernel_halo_helpers.h"

#include "debug.h"
#include "diagnostics.h"

#define REALLY_SMALL 1e-300

template<int NX, int NY>
__global__ void Momentum_Eq(double *M_d,
                            double *pressure_d,
                            double *SlowMh_d,
                            double *grad_d,
                            double *Altitude_d,
                            double *DivM_d,
                            double  A,
                            double *func_r_d,
                            double  dt,
                            int *   maps_d,
                            int     nl_region,
                            bool    DeepModel) {

    int x = threadIdx.x;
    int y = threadIdx.y;
    //int ib  = blockIdx.x;
    int nv  = gridDim.y;
    int lev = blockIdx.y;

    int pt1, pt2, pt3, pt4, pt5, pt6;
    int nhl  = nl_region + 2;
    int nhl2 = nhl * nhl;

    int iri;

    __shared__ double nflxv_s[3 * NX * NY];
    __shared__ double pressure_s[(NX + 2) * (NY + 2)];


    int ir = 0; // index in region

    int    ir2, id;
    double vr;
    double alt, rscale;
    double Mx, My, Mz;
    double funcx, funcy, funcz;


    bool pent_ind = false; //
    int  ig;               // index in global mem

    int igh = 0; // global index in halo

    // Load shared memory
    bool load_halo = compute_mem_idx(maps_d, nhl, nhl2, ig, igh, ir, ir2, pent_ind);
    id             = ig;
    pressure_s[ir] = pressure_d[ig * nv + lev];

    ///////////////////////////////
    //////////// Halo /////////////
    ///////////////////////////////
    if (load_halo) {
        if (igh >= 0) {
            pressure_s[ir2] = pressure_d[igh * nv + lev];
        }
        else
            pressure_s[ir2] = 0.0;
    }
    __syncthreads();
    //////////////////////////////////////////////

    iri = y * nl_region + x;

    funcx = func_r_d[id * 3 + 0];
    funcy = func_r_d[id * 3 + 1];
    funcz = func_r_d[id * 3 + 2];

    pt1 = (y + 2) * nhl + x + 1;
    pt2 = (y + 2) * nhl + x + 2;
    pt3 = (y + 1) * nhl + x + 2;
    pt4 = (y)*nhl + x + 1;
    pt5 = (pent_ind) * ((y + 1) * nhl + x) + (!pent_ind) * ((y)*nhl + x);
    pt6 = (y + 1) * nhl + x;

    if (DeepModel) {
        alt    = Altitude_d[lev];
        rscale = A / (alt + A);
    }
    else
        rscale = 1.0;

    for (int k = 0; k < 3; k++) {
        nflxv_s[iri * 3 + k] = rscale
                               * (grad_d[id * 7 * 3 + 3 * 0 + k] * pressure_s[ir]
                                  + grad_d[id * 7 * 3 + 3 * 1 + k] * pressure_s[pt1]
                                  + grad_d[id * 7 * 3 + 3 * 2 + k] * pressure_s[pt2]
                                  + grad_d[id * 7 * 3 + 3 * 3 + k] * pressure_s[pt3]
                                  + grad_d[id * 7 * 3 + 3 * 4 + k] * pressure_s[pt4]
                                  + grad_d[id * 7 * 3 + 3 * 5 + k] * pressure_s[pt5]
                                  + grad_d[id * 7 * 3 + 3 * 6 + k] * pressure_s[pt6]);
    }

    Mx = (-nflxv_s[iri * 3 + 0] + SlowMh_d[id * 3 * nv + lev * 3 + 0]
          + DivM_d[id * 3 * nv + lev * 3 + 0])
         * dt;
    My = (-nflxv_s[iri * 3 + 1] + SlowMh_d[id * 3 * nv + lev * 3 + 1]
          + DivM_d[id * 3 * nv + lev * 3 + 1])
         * dt;
    Mz = (-nflxv_s[iri * 3 + 2] + SlowMh_d[id * 3 * nv + lev * 3 + 2]
          + DivM_d[id * 3 * nv + lev * 3 + 2])
         * dt;

    vr = Mx * funcx + My * funcy + Mz * funcz;

    Mx += -vr * funcx;
    My += -vr * funcy;
    Mz += -vr * funcz;

    // Updates momenta
    M_d[id * nv * 3 + lev * 3 + 0] += Mx;
    M_d[id * nv * 3 + lev * 3 + 1] += My;
    M_d[id * nv * 3 + lev * 3 + 2] += Mz;
}

template<int NN>
__global__ void Momentum_Eq_Poles(double *M_d,
                                  double *pressure_d,
                                  double *SlowMh_d,
                                  double *grad_d,
                                  double *Altitude_d,
                                  double *DivM_d,
                                  double  A,
                                  double *func_r_d,
                                  double  dt,
                                  int *   point_local_d,
                                  int     nv,
                                  int     num,
                                  bool    DeepModel) {

    int id = blockIdx.x * blockDim.x + threadIdx.x;
    id += num - 2; // Poles

    __shared__ double grad_p[3 * 7];
    __shared__ double func_r_p[3];
    __shared__ double nflxv_p[3];
    __shared__ double pressure_p[NN];
    __shared__ int    local_p[NN];

    double vr;
    double alt, rscale;
    double Mx, My, Mz;

    if (id < num) {
        for (int i = 0; i < 5; i++)
            local_p[i] = point_local_d[id * 6 + i];
        func_r_p[0] = func_r_d[id * 3 + 0];
        func_r_p[1] = func_r_d[id * 3 + 1];
        func_r_p[2] = func_r_d[id * 3 + 2];
        for (int i = 0; i < 7; i++)
            for (int k = 0; k < 3; k++)
                grad_p[i * 3 + k] = grad_d[id * 7 * 3 + i * 3 + k];

        for (int lev = 0; lev < nv; lev++) {

            pressure_p[0] = pressure_d[id * nv + lev];
            for (int i = 1; i < 6; i++)
                pressure_p[i] = pressure_d[local_p[i - 1] * nv + lev];

            alt = Altitude_d[lev];

            if (DeepModel)
                rscale = A / (alt + A);
            else
                rscale = 1.0;

            for (int k = 0; k < 3; k++) {
                // clang-format off
                nflxv_p[k] = rscale * (grad_p[3 * 0 + k] * pressure_p[0]
                                     + grad_p[3 * 1 + k] * pressure_p[1]
                                     + grad_p[3 * 2 + k] * pressure_p[2]
                                     + grad_p[3 * 3 + k] * pressure_p[3]
                                     + grad_p[3 * 4 + k] * pressure_p[4]
                                     + grad_p[3 * 5 + k] * pressure_p[5]);
                // clang-format on
            }

            Mx = (-nflxv_p[0] + SlowMh_d[id * 3 * nv + lev * 3 + 0]
                  + DivM_d[id * 3 * nv + lev * 3 + 0])
                 * dt;
            My = (-nflxv_p[1] + SlowMh_d[id * 3 * nv + lev * 3 + 1]
                  + DivM_d[id * 3 * nv + lev * 3 + 1])
                 * dt;
            Mz = (-nflxv_p[2] + SlowMh_d[id * 3 * nv + lev * 3 + 2]
                  + DivM_d[id * 3 * nv + lev * 3 + 2])
                 * dt;

            vr = Mx * func_r_p[0] + My * func_r_p[1] + Mz * func_r_p[2];

            Mx += -vr * func_r_p[0];
            My += -vr * func_r_p[1];
            Mz += -vr * func_r_p[2];

            // Updates momenta
            M_d[id * nv * 3 + lev * 3 + 0] += Mx;
            M_d[id * nv * 3 + lev * 3 + 1] += My;
            M_d[id * nv * 3 + lev * 3 + 2] += Mz;
        }
    }
}

template<int NX, int NY>
__global__ void Density_Pressure_Eqs(double *      pressure_d,
                                     double *      pressurek_d,
                                     double *      Rho_d,
                                     double *      Rhok_d,
                                     double *      Mh_d,
                                     double *      Mhk_d,
                                     double *      Wh_d,
                                     double *      Whk_d,
                                     double *      pt_d,
                                     double *      pth_d,
                                     double *      epotential_d,
                                     double *      epotentialh_d,
                                     double *      ekinetic_d,
                                     double *      ekinetich_d,
                                     double *      Etotal_tau_d,
                                     double *      h_d,
                                     double *      hh_d,
                                     double *      SlowRho_d,
                                     double *      profx_Qheat_d,
                                     double *      diffpr_d,
                                     double *      diffprv_d,
                                     double *      div_d,
                                     double *      Altitude_d,
                                     double *      Altitudeh_d,
                                     double *      Cp_d,
                                     double *      Rd_d,
                                     double        A,
                                     double        P_Ref,
                                     double        Gravit,
                                     double        dt,
                                     int *         maps_d,
                                     int           nl_region,
                                     bool          DeepModel,
                                     bool          energy_equation,
                                     unsigned int *diagnostics_flag,
                                     diag_data *   diagnostics_data) {
    // This function is one of the tricky part in the algorithm.
    // It uses the output of the vertical solver using a thomas algorithm
    int x = threadIdx.x;
    int y = threadIdx.y;
    //    int ib  = blockIdx.x;
    int nv  = gridDim.y;
    int lev = blockIdx.y;

    int               pt1, pt2, pt3, pt4, pt5, pt6;
    double            div0, div1, div2, div3, div4, div5, div6;
    int               nhl  = nl_region + 2;
    int               nhl2 = nhl * nhl;
    __shared__ double nflxr_s[NX * NY];
    __shared__ double nflxpt_s[NX * NY];
    __shared__ double v_s[3 * (NX + 2) * (NY + 2)];
    __shared__ double v1_s[3 * (NX + 2) * (NY + 2)];
    __shared__ double pt_s[(NX + 2) * (NY + 2)];

    double alt, r2p, r2m, r2l, rscale;
    double wht, whl;
    double wht2, whl2;
    double pht, phl;
    double dz, dwdz;
    double dwptdz;
    double aux, r, p;
    double altht, althl;
    double Cv;
    double dr2dz;


    int ir = 0; // index in region
    int iri, ir2, id;


    bool pent_ind = false; //
    int  ig;               // index in global mem

    int igh = 0; // global index in halo

    // Load shared memory


    bool load_halo = compute_mem_idx(maps_d, nhl, nhl2, ig, igh, ir, ir2, pent_ind);
    id             = ig;

    Cv = Cp_d[id * nv + lev] - Rd_d[id * nv + lev];

    // Load shared memory

    // horizontal momentum deviation
    v_s[ir * 3 + 0] = Mh_d[ig * 3 * nv + lev * 3 + 0];
    v_s[ir * 3 + 1] = Mh_d[ig * 3 * nv + lev * 3 + 1];
    v_s[ir * 3 + 2] = Mh_d[ig * 3 * nv + lev * 3 + 2];

    // total horizontal momentum
    v1_s[ir * 3 + 0] = v_s[ir * 3 + 0] + Mhk_d[ig * 3 * nv + lev * 3 + 0];
    v1_s[ir * 3 + 1] = v_s[ir * 3 + 1] + Mhk_d[ig * 3 * nv + lev * 3 + 1];
    v1_s[ir * 3 + 2] = v_s[ir * 3 + 2] + Mhk_d[ig * 3 * nv + lev * 3 + 2];

    if (energy_equation) {
        //use potential temp array to store sum of enthalpy, potential, kinetic
        pt_s[ir] = (h_d[ig * nv + lev] + epotential_d[ig * nv + lev] + ekinetic_d[ig * nv + lev]);
    }
    else {
        pt_s[ir] = pt_d[ig * nv + lev];
    }
    ///////////////////////////////
    //////////// Halo /////////////
    ///////////////////////////////
    if (load_halo) {
        if (igh >= 0) {
            v_s[ir2 * 3 + 0]  = Mh_d[igh * 3 * nv + lev * 3 + 0];
            v_s[ir2 * 3 + 1]  = Mh_d[igh * 3 * nv + lev * 3 + 1];
            v_s[ir2 * 3 + 2]  = Mh_d[igh * 3 * nv + lev * 3 + 2];
            v1_s[ir2 * 3 + 0] = v_s[ir2 * 3 + 0] + Mhk_d[igh * 3 * nv + lev * 3 + 0];
            v1_s[ir2 * 3 + 1] = v_s[ir2 * 3 + 1] + Mhk_d[igh * 3 * nv + lev * 3 + 1];
            v1_s[ir2 * 3 + 2] = v_s[ir2 * 3 + 2] + Mhk_d[igh * 3 * nv + lev * 3 + 2];
            if (energy_equation) {
                pt_s[ir2] = (h_d[igh * nv + lev] + epotential_d[igh * nv + lev]
                             + ekinetic_d[igh * nv + lev]);
            }
            else {
                pt_s[ir2] = pt_d[igh * nv + lev];
            }
        }
        else {
            v_s[ir2 * 3 + 0]  = 0.0;
            v_s[ir2 * 3 + 1]  = 0.0;
            v_s[ir2 * 3 + 2]  = 0.0;
            v1_s[ir2 * 3 + 0] = 0.0;
            v1_s[ir2 * 3 + 1] = 0.0;
            v1_s[ir2 * 3 + 2] = 0.0;
            pt_s[ir2]         = 0.0;
        }
    }

    __syncthreads();
    //////////////////////////////////////////////

    iri = (y)*nl_region + x;

    pt1 = (y + 2) * nhl + x + 1;
    pt2 = (y + 2) * nhl + x + 2;
    pt3 = (y + 1) * nhl + x + 2;
    pt4 = (y)*nhl + x + 1;
    pt5 = (pent_ind) * ((y + 1) * nhl + x) + (!pent_ind) * ((y)*nhl + x);
    pt6 = (y + 1) * nhl + x;

    altht = Altitudeh_d[lev + 1];
    althl = Altitudeh_d[lev];

    dz = altht - althl;
    if (DeepModel) {
        alt    = Altitude_d[lev];
        r2p    = pow(altht + A, 2.0);
        r2m    = pow(alt + A, 2.0);
        r2l    = pow(althl + A, 2.0);
        rscale = A / (alt + A);
        dr2dz  = (pow(altht + A, 3.0) - pow(althl + A, 3.0)) / 3.0;
    }
    else {
        r2p    = 1.0;
        r2m    = 1.0;
        r2l    = 1.0;
        rscale = 1.0;
        dr2dz  = r2m * dz;
    }

    nflxr_s[iri]  = 0.0;
    nflxpt_s[iri] = 0.0;

    for (int k = 0; k < 3; k++) {

        div0 = div_d[id * 7 * 3 + 3 * 0 + k];
        div1 = div_d[id * 7 * 3 + 3 * 1 + k];
        div2 = div_d[id * 7 * 3 + 3 * 2 + k];
        div3 = div_d[id * 7 * 3 + 3 * 3 + k];
        div4 = div_d[id * 7 * 3 + 3 * 4 + k];
        div5 = div_d[id * 7 * 3 + 3 * 5 + k];
        div6 = div_d[id * 7 * 3 + 3 * 6 + k];

        // clang-format off
        nflxr_s[iri] += rscale * (div0 * v_s[ir * 3 + k]
                                + div1 * v_s[pt1 * 3 + k]
                                + div2 * v_s[pt2 * 3 + k]
                                + div3 * v_s[pt3 * 3 + k]
                                + div4 * v_s[pt4 * 3 + k]
                                + div5 * v_s[pt5 * 3 + k]
                                + div6 * v_s[pt6 * 3 + k]);

        nflxpt_s[iri] += rscale * (div0 * v1_s[ir * 3 + k] * pt_s[ir]
                                 + div1 * v1_s[pt1 * 3 + k] * pt_s[pt1]
                                 + div2 * v1_s[pt2 * 3 + k] * pt_s[pt2]
                                 + div3 * v1_s[pt3 * 3 + k] * pt_s[pt3]
                                 + div4 * v1_s[pt4 * 3 + k] * pt_s[pt4]
                                 + div5 * v1_s[pt5 * 3 + k] * pt_s[pt5]
                                 + div6 * v1_s[pt6 * 3 + k] * pt_s[pt6]);
        // clang-format on
    }

    if (lev == 0) {
        whl  = 0.0;
        wht  = Wh_d[id * (nv + 1) + lev + 1];
        whl2 = 0.0;
        wht2 = Wh_d[id * (nv + 1) + lev + 1] + Whk_d[id * (nv + 1) + lev + 1];
        phl  = 0.0;
        if (energy_equation) {
            pht = (hh_d[id * (nv + 1) + lev + 1] + epotentialh_d[id * (nv + 1) + lev + 1]
                   + ekinetich_d[id * (nv + 1) + lev + 1]);
        }
        else {
            pht = pth_d[id * (nv + 1) + lev + 1];
        }
    }
    else {
        whl  = Wh_d[id * (nv + 1) + lev];
        wht  = Wh_d[id * (nv + 1) + lev + 1];
        whl2 = Wh_d[id * (nv + 1) + lev] + Whk_d[id * (nv + 1) + lev];
        wht2 = Wh_d[id * (nv + 1) + lev + 1] + Whk_d[id * (nv + 1) + lev + 1];
        if (energy_equation) {
            phl = (hh_d[id * (nv + 1) + lev] + epotentialh_d[id * (nv + 1) + lev]
                   + ekinetich_d[id * (nv + 1) + lev]);
            pht = (hh_d[id * (nv + 1) + lev + 1] + epotentialh_d[id * (nv + 1) + lev + 1]
                   + ekinetich_d[id * (nv + 1) + lev + 1]);
        }
        else {
            phl = pth_d[id * (nv + 1) + lev];
            pht = pth_d[id * (nv + 1) + lev + 1];
        }
    }

    // uses thomas algorithm computed vertical velocities in dyn/thor_vertical_int.h :: vertical_eq
    // wht/whl and wht/whl2 -> can have big error

    dwdz   = (wht * r2p - whl * r2l) / (dr2dz);
    dwptdz = (wht2 * pht * r2p - whl2 * phl * r2l) / (dz * r2m);

    // the errors can sometimes make aux become negative and cause issues later
    aux = -(nflxpt_s[iri] + dwptdz) * dt;               //advection terms in thermo eqn
    r   = Rhok_d[id * nv + lev] + Rho_d[id * nv + lev]; //density at time tau


    // Updates density
    nflxr_s[iri] += dwdz; //hack to test mass conservation
    Rho_d[id * nv + lev] +=
        (SlowRho_d[id * nv + lev] - nflxr_s[iri]) * dt; //density deviation at time tau+dtau

    // back to thermo equation
    if (energy_equation) {
        double Epot_kin, w;
        aux += Etotal_tau_d[id * nv + lev]
               + (diffpr_d[id * nv + lev] + diffprv_d[id * nv + lev]) * dt
               + profx_Qheat_d[id * nv + lev] * dt; // new total energy

        Etotal_tau_d[id * nv + lev] = aux; //store Etotal for next small step

        // vert mom at tau + dtau, center of cell
        w = (whl2 * (Altitudeh_d[lev + 1] - Altitude_d[lev])
             + wht2 * (Altitude_d[lev] - Altitudeh_d[lev]))
            / (Altitudeh_d[lev + 1] - Altitudeh_d[lev]);

        //new potential + kinetic energy using newest momenta and density
        Epot_kin = (Rho_d[id * nv + lev] + Rhok_d[id * nv + lev]) * Gravit * Altitude_d[lev]
                   + 0.5
                         * (pow(v1_s[ir * 3 + 0], 2) + pow(v1_s[ir * 3 + 1], 2)
                            + pow(v1_s[ir * 3 + 2], 2) + pow(w, 2))
                         / (Rho_d[id * nv + lev] + Rhok_d[id * nv + lev]);
        //pressure perturbation
        pressure_d[id * nv + lev] =
            Rd_d[id * nv + lev] / Cv * (aux - Epot_kin) - pressurek_d[id * nv + lev];
    }
    else {
        double pt; //, pt_p;
        pt = (P_Ref / (Rd_d[id * nv + lev] * r))
             * pow((pressure_d[id * nv + lev] + pressurek_d[id * nv + lev]) / P_Ref,
                   Cv / Cp_d[id * nv + lev]);

        // negative aux can get wrongly negative here and cause NaN in fractional power computation.
        // r: current density
        aux += pt * r;

#if defined(DIAG_CHECK_DENSITY_PRESSURE_EQ_AUX) || defined(DIAG_CHECK_DENSITY_PRESSURE_EQ_P_NAN)
        // clear diagnostics memory
        diagnostics_data[id * nv + lev].flag = 0;
        diagnostics_data[id * nv + lev].data = make_double4(0.0, 0.0, 0.0, 0.0);
#endif // defined(DIAG_CHECK_DENSITY_PRESSURE_EQ_AUX) || defined(DIAG_CHECK_DENSITY_PRESSURE_EQ_P_NAN)


#ifdef DIAG_CHECK_DENSITY_PRESSURE_EQ_AUX
        if (aux < 0.0) {
            // set global diagnostics flag atomically.
            atomicOr(diagnostics_flag, NEGATIVE_VALUE);
            // add the element diagnostics flag to diag array
            diagnostics_data[id * nv + lev].flag |= NEGATIVE_VALUE;
            // store aux value in diagnostics data array
            diagnostics_data[id * nv + lev].data.x = aux;
        }
#endif // DIAG_CHECK_DENSITY_PRESSURE_EQ_AUX

        // fractional power that returns a NaN
        p = P_Ref * pow(Rd_d[id * nv + lev] * aux / P_Ref, Cp_d[id * nv + lev] / Cv);
#ifdef DIAG_CHECK_DENSITY_PRESSURE_EQ_P_NAN
        if (isnan(p)) {
            atomicOr(diagnostics_flag, NAN_VALUE);
            diagnostics_data[id * nv + lev].flag |= NAN_VALUE;
        }
#endif // DIAG_CHECK_DENSITY_PRESSURE_EQ_P_NAN
        pressure_d[id * nv + lev] = p - pressurek_d[id * nv + lev]
                                    + (diffpr_d[id * nv + lev] + diffprv_d[id * nv + lev]) * dt
                                    + Rd_d[id * nv + lev] / Cv * profx_Qheat_d[id * nv + lev] * dt;
    }
}


template<int NN>
__global__ void Density_Pressure_Eqs_Poles(double *      pressure_d,
                                           double *      pressurek_d,
                                           double *      Rho_d,
                                           double *      Rhok_d,
                                           double *      Mh_d,
                                           double *      Mhk_d,
                                           double *      Wh_d,
                                           double *      Whk_d,
                                           double *      pt_d,
                                           double *      pth_d,
                                           double *      epotential_d,
                                           double *      epotentialh_d,
                                           double *      ekinetic_d,
                                           double *      ekinetich_d,
                                           double *      Etotal_tau_d,
                                           double *      h_d,
                                           double *      hh_d,
                                           double *      SlowRho_d,
                                           double *      profx_Qheat_d,
                                           double *      diffpr_d,
                                           double *      diffprv_d,
                                           double *      div_d,
                                           double *      Altitude_d,
                                           double *      Altitudeh_d,
                                           double *      Cp_d,
                                           double *      Rd_d,
                                           double        A,
                                           double        P_Ref,
                                           double        Gravit,
                                           double        dt,
                                           int *         point_local_d,
                                           int           num,
                                           int           nv,
                                           bool          DeepModel,
                                           bool          energy_equation,
                                           unsigned int *diagnostics_flag,
                                           diag_data *   diagnostics_data) {

    int id = blockIdx.x * blockDim.x + threadIdx.x;
    id += num - 2; // Poles

    __shared__ double div_p[3 * 7];
    __shared__ double v_p[3 * NN];
    __shared__ double v1_p[3 * NN];
    __shared__ double pt_p[NN];
    __shared__ int    local_p[NN];

    double nflxr_p;
    double nflxpt_p;

    double alt, r2p, r2m, r2l, rscale;
    double wht, whl;
    double wht2, whl2;
    double pht, phl;
    double dz, dwdz;
    double dwptdz;
    double aux, r, p;
    double altht, althl;
    double Cv;

    if (id < num) {
        for (int i = 0; i < 5; i++)
            local_p[i] = point_local_d[id * 6 + i];
        for (int i = 0; i < 7; i++)
            for (int k = 0; k < 3; k++)
                div_p[i * 3 + k] = div_d[id * 7 * 3 + i * 3 + k];

        for (int lev = 0; lev < nv; lev++) {
            Cv = Cp_d[id * nv + lev] - Rd_d[id * nv + lev];
            // horizontal momentum deviation
            v_p[0] = Mh_d[id * 3 * nv + lev * 3 + 0];
            v_p[1] = Mh_d[id * 3 * nv + lev * 3 + 1];
            v_p[2] = Mh_d[id * 3 * nv + lev * 3 + 2];
            // total horizontal momentum
            v1_p[0] = Mh_d[id * 3 * nv + lev * 3 + 0] + Mhk_d[id * 3 * nv + lev * 3 + 0];
            v1_p[1] = Mh_d[id * 3 * nv + lev * 3 + 1] + Mhk_d[id * 3 * nv + lev * 3 + 1];
            v1_p[2] = Mh_d[id * 3 * nv + lev * 3 + 2] + Mhk_d[id * 3 * nv + lev * 3 + 2];
            if (energy_equation) {
                pt_p[0] =
                    (h_d[id * nv + lev] + epotential_d[id * nv + lev] + ekinetic_d[id * nv + lev]);
            }
            else {
                pt_p[0] = pt_d[id * nv + lev];
            }
            for (int i = 1; i < 6; i++) {
                v_p[i * 3 + 0]  = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 0];
                v_p[i * 3 + 1]  = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 1];
                v_p[i * 3 + 2]  = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 2];
                v1_p[i * 3 + 0] = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 0]
                                  + Mhk_d[local_p[i - 1] * 3 * nv + lev * 3 + 0];
                v1_p[i * 3 + 1] = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 1]
                                  + Mhk_d[local_p[i - 1] * 3 * nv + lev * 3 + 1];
                v1_p[i * 3 + 2] = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 2]
                                  + Mhk_d[local_p[i - 1] * 3 * nv + lev * 3 + 2];
                if (energy_equation) {
                    pt_p[i] =
                        (h_d[local_p[i - 1] * nv + lev] + epotential_d[local_p[i - 1] * nv + lev]
                         + ekinetic_d[local_p[i - 1] * nv + lev]);
                }
                else {
                    pt_p[i] = pt_d[local_p[i - 1] * nv + lev];
                }
            }

            if (lev == 0) {
                altht = Altitudeh_d[lev + 1];
                althl = Altitudeh_d[lev];
            }

            alt = Altitude_d[lev];

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

            nflxr_p  = 0.0;
            nflxpt_p = 0.0;

            for (int k = 0; k < 3; k++) {
                // clang-format off
                nflxr_p += rscale * (div_p[3 * 0 + k] * v_p[0 * 3 + k]
                                  + div_p[3 * 1 + k] * v_p[1 * 3 + k]
                                  + div_p[3 * 2 + k] * v_p[2 * 3 + k]
                                  + div_p[3 * 3 + k] * v_p[3 * 3 + k]
                                  + div_p[3 * 4 + k] * v_p[4 * 3 + k]
                                  + div_p[3 * 5 + k] * v_p[5 * 3 + k]);

                nflxpt_p += rscale
                              * (div_p[3 * 0 + k] * v1_p[0 * 3 + k] * pt_p[0]
                               + div_p[3 * 1 + k] * v1_p[1 * 3 + k] * pt_p[1]
                               + div_p[3 * 2 + k] * v1_p[2 * 3 + k] * pt_p[2]
                               + div_p[3 * 3 + k] * v1_p[3 * 3 + k] * pt_p[3]
                               + div_p[3 * 4 + k] * v1_p[4 * 3 + k] * pt_p[4]
                               + div_p[3 * 5 + k] * v1_p[5 * 3 + k] * pt_p[5]);
                // clang-format on
            }

            if (lev == 0) {
                whl  = 0.0;
                wht  = Wh_d[id * (nv + 1) + lev + 1];
                whl2 = 0.0;
                wht2 = Wh_d[id * (nv + 1) + lev + 1] + Whk_d[id * (nv + 1) + lev + 1];
                phl  = 0.0;
                if (energy_equation) {
                    pht = (hh_d[id * (nv + 1) + lev + 1] + epotentialh_d[id * (nv + 1) + lev + 1]
                           + ekinetich_d[id * (nv + 1) + lev + 1]);
                }
                else {
                    pht = pth_d[id * (nv + 1) + lev + 1];
                }
            }

            dz     = altht - althl;
            dwdz   = (wht * r2p - whl * r2l) / (dz * r2m);
            dwptdz = (wht2 * pht * r2p - whl2 * phl * r2l) / (dz * r2m);

            aux = -(nflxpt_p + dwptdz) * dt;                    //advection terms in thermo eqn
            r   = Rhok_d[id * nv + lev] + Rho_d[id * nv + lev]; //density at time tau

            // Updates density
            nflxr_p += dwdz;
            Rho_d[id * nv + lev] +=
                (SlowRho_d[id * nv + lev] - nflxr_p) * dt; //density at time tau+dtau

            //back to thermo equation
            if (energy_equation) {
                double Epot_kin, w;
                aux += Etotal_tau_d[id * nv + lev]
                       + (diffpr_d[id * nv + lev] + diffprv_d[id * nv + lev]) * dt
                       + profx_Qheat_d[id * nv + lev] * dt; // new total energy

                Etotal_tau_d[id * nv + lev] = aux; //store Etotal for next small step

                // vert mom at tau + dtau, center of cell
                w = (whl2 * (Altitudeh_d[lev + 1] - Altitude_d[lev])
                     + wht2 * (Altitude_d[lev] - Altitudeh_d[lev]))
                    / (Altitudeh_d[lev + 1] - Altitudeh_d[lev]);

                //new potential + kinetic energy using newest momenta and density
                Epot_kin = (Rho_d[id * nv + lev] + Rhok_d[id * nv + lev]) * Gravit * Altitude_d[lev]
                           + 0.5
                                 * (pow(v1_p[0 * 3 + 0], 2) + pow(v1_p[0 * 3 + 1], 2)
                                    + pow(v1_p[0 * 3 + 2], 2) + pow(w, 2))
                                 / (Rho_d[id * nv + lev] + Rhok_d[id * nv + lev]);
                //pressure perturbation
                pressure_d[id * nv + lev] =
                    Rd_d[id * nv + lev] / Cv * (aux - Epot_kin) - pressurek_d[id * nv + lev];
            }
            else {
                double pt, ptmp;
                ptmp = pressure_d[id * nv + lev];
                pt   = (P_Ref / (Rd_d[id * nv + lev] * r))
                     * pow((ptmp + pressurek_d[id * nv + lev]) / P_Ref, Cv / Cp_d[id * nv + lev]);

                aux += pt * r;
#if defined(DIAG_CHECK_DENSITY_PRESSURE_EQ_AUX) || defined(DIAG_CHECK_DENSITY_PRESSURE_EQ_P_NAN)
                // clear diagnostics memory
                diagnostics_data[id * nv + lev].flag = 0;
                diagnostics_data[id * nv + lev].data = make_double4(0.0, 0.0, 0.0, 0.0);
#endif

#ifdef DIAG_CHECK_DENSITY_PRESSURE_EQ_AUX

                if (aux < 0.0) {
                    // set flags. Use binary OR operator to enable the flag
                    // without overwriting other flags that could be present
                    atomicOr(diagnostics_flag, NEGATIVE_VALUE);
                    diagnostics_data[id * nv + lev].flag |= NEGATIVE_VALUE;
                    diagnostics_data[id * nv + lev].data.x = aux;
                }
#endif // DIAG_CHECK_DENSITY_PRESSURE_EQ_AUX

                p = P_Ref * pow(Rd_d[id * nv + lev] * aux / P_Ref, Cp_d[id * nv + lev] / Cv); //
#ifdef DIAG_CHECK_DENSITY_PRESSURE_EQ_P_NAN
                if (isnan(p)) {
                    // set our second flag
                    atomicOr(diagnostics_flag, NAN_VALUE);
                    diagnostics_data[id * nv + lev].flag |= NAN_VALUE;
                }
#endif // DIAG_CHECK_DENSITY_PRESSURE_EQ_P_NAN


                // Updates pressure
                pressure_d[id * nv + lev] =
                    p - pressurek_d[id * nv + lev]
                    + (diffpr_d[id * nv + lev] + diffprv_d[id * nv + lev]) * dt
                    + Rd_d[id * nv + lev] / Cv * profx_Qheat_d[id * nv + lev] * dt;
            }


            if (lev != nv - 1) {
                althl = altht;
                altht = Altitudeh_d[lev + 2];
                whl   = wht;
                wht   = Wh_d[id * (nv + 1) + lev + 2];
                whl2  = wht2;
                wht2  = Wh_d[id * (nv + 1) + lev + 2] + Whk_d[id * (nv + 1) + lev + 2];
                phl   = pht;
                if (energy_equation) {
                    pht = (hh_d[id * (nv + 1) + lev + 2] + epotentialh_d[id * (nv + 1) + lev + 2]
                           + ekinetich_d[id * (nv + 1) + lev + 2]);
                }
                else {
                    pht = pth_d[id * (nv + 1) + lev + 2];
                }
            }
        }
    }
}
