//
//
//
//
// Description:
//
//
// Method:
//   Calls to
//
// Known limitations:
//
//
// Known issues:
//
//
// Current Code Owner:
//
// History:
// Version Date       Comment
// ======= ====       =======
//
// 1.0
//
////////////////////////////////////////////////////////////////////////

#include "debug_helpers.h"

__global__ void annual_insol(double *insol_ann_d, double *insol_d, int nstep, int num) {

    int id = blockIdx.x * blockDim.x + threadIdx.x;

    if (id < num) {
        insol_ann_d[id] = insol_ann_d[id] * (nstep - 1) / nstep + insol_d[id] / nstep;
    }
}

__device__ void radcsw(double *phtemp,
                       double  coszrs,
                       double  r_orb,
                       double *dtemp,
                       double *tau_d,
                       double *fsw_up_d,
                       double *fsw_dn_d,
                       double *Altitude_d,
                       double *Altitudeh_d,
                       double  incflx,
                       double  alb,
                       double  kappa_sw,
                       double  ps0,
                       double  gravit,
                       double  A,
                       int     id,
                       int     nv,
                       double *insol_d,
                       bool    DeepModel) {

    //  Calculate upward, downward, and net flux.
    //  Downward Directed Radiation

    // double gocp;
    //double tau      = (tausw / ps0) * (phtemp[id * (nv + 1) + nv]);
    double tau      = (kappa_sw / gravit) * (phtemp[id * (nv + 1) + nv]);
    insol_d[id]     = incflx * pow(r_orb, -2) * coszrs;
    double flux_top = insol_d[id] * (1.0 - alb);
    double rup, rlow;

    // Extra layer to avoid over heating at the top.
    fsw_dn_d[id * (nv + 1) + nv] = flux_top * exp(-(1.0 / coszrs) * tau);

    // Normal integration
    for (int lev = nv; lev >= 1; lev--)
        fsw_dn_d[id * (nv + 1) + lev - 1] =
            fsw_dn_d[id * (nv + 1) + lev]
            * exp(-(1.0 / coszrs) * tau_d[id * nv * 2 + (lev - 1) * 2]);
    for (int lev = 0; lev <= nv; lev++)
        fsw_up_d[id * (nv + 1) + lev] = 0.0;

    // Update temperature rates.
    for (int lev = 0; lev < nv; lev++) {
        if (DeepModel) { //this seems to cause strange problems at TOA, set both factors to 1 for now
            // rup =
            //     (Altitudeh_d[lev + 1] + A) / (Altitude_d[lev] + A); //vertical scaling in divergence
            // rlow = (Altitudeh_d[lev] + A) / (Altitude_d[lev] + A);
            rup  = 1.0;
            rlow = 1.0;
        }
        else {
            rup  = 1.0;
            rlow = 1.0;
        }
        dtemp[id * nv + lev] =
            -(pow(rlow, 2) * (fsw_up_d[id * (nv + 1) + lev] - fsw_dn_d[id * (nv + 1) + lev])
              - pow(rup, 2)
                    * (fsw_up_d[id * (nv + 1) + lev + 1] - fsw_dn_d[id * (nv + 1) + lev + 1]))
            / ((Altitudeh_d[lev] - Altitudeh_d[lev + 1]));
        // gocp = gravit / Cp;
        // dtemp[id * nv + lev] =
        //     gocp
        //     * ((fsw_up_d[id * (nv + 1) + lev] - fsw_dn_d[id * (nv + 1) + lev])
        //        - (fsw_up_d[id * (nv + 1) + lev + 1] - fsw_dn_d[id * (nv + 1) + lev + 1]))
        //     / (phtemp[id * (nv + 1) + lev] - phtemp[id * (nv + 1) + lev + 1]);

        // printf("%d %e\n", lev, (dtemp2 - dtemp[id * nv + lev]) / dtemp2);
        if (isnan(dtemp[id * nv + lev])) {
            printf("stop here");
        }
    }
}

__device__ int factorial_num(int number) {

    int temp;
    int n1;
    temp = 1;
    n1   = number;
    while (n1 != 0) {
        temp = temp * n1;
        n1   = n1 - 1;
    }
    return temp;
}

__device__ double source_func_lin(double bb, double bl, double bt, double tau, double diff_ang) {

    double e1 = 0.0;
    double e2 = 0.0;
    double e  = 0.0;

    if (tau >= 1e-10) {
        e1 = bb - bl
             + (bl - (diff_ang / (tau / 2.0)) * (bb - bl)) * (1.0 - exp(-(tau / 2.0) / diff_ang));
        e2 = bl - bt
             + (bt - (diff_ang / (tau / 2.0)) * (bl - bt)) * (1.0 - exp(-(tau / 2.0) / diff_ang));
        e = e1 + e2 * exp(-(tau / 2.0) / diff_ang);
    }
    else {
        for (int i = 0; i < 5; i++) {
            int fac = factorial_num(i + 2);
            e1      = e1
                 + (pow(-1.0, i + 2.0)) * ((bl + (i + 1) * bl) / fac)
                       * pow((tau / 2.0) / diff_ang, i + 1.0);
            e2 = e2
                 + (pow(-1.0, i + 2.0)) * ((bt + (i + 1) * bb) / fac)
                       * pow((tau / 2.0) / diff_ang, i + 1.0);
        }
        e = e1 + e2 * exp(-(tau / 2.0) / diff_ang);
    }

    return e;
}

__device__ void radclw(double *phtemp,
                       double *ttemp,
                       double *thtemp,
                       double *dtemp,
                       double *tau_d,
                       double *flw_up_d,
                       double *flw_dn_d,
                       double *Altitude_d,
                       double *Altitudeh_d,
                       double  diff_ang,
                       double  tint,
                       double  gravit,
                       bool    surface,
                       double  Tsurface,
                       double  A,
                       int     id,
                       int     nv,
                       bool    DeepModel) {

    // double gocp = gravit / Cp;
    double tb, tl, tt;
    double bb, bl, bt;
    double rup, rlow;

    double bc = 5.677036E-8; //Stefan–Boltzmann constant W⋅m−2⋅K−4

    //
    //  Calculate upward, downward, and net flux.
    //  Downward Directed Radiation
    //
    flw_dn_d[id * (nv + 1) + nv] = 0.0; // Upper boundary
    for (int lev = nv - 1; lev >= 0; lev--) {
        double ed = 0.0;
        if (tau_d[id * nv * 2 + 2 * lev + 1] < 0.0)
            tau_d[id * nv * 2 + 2 * lev + 1] = 0.0;

        tb = thtemp[id * (nv + 1) + lev];
        tl = ttemp[id * nv + lev];
        tt = thtemp[id * (nv + 1) + lev + 1];

        bb = bc * tb * tb * tb * tb;
        bl = bc * tl * tl * tl * tl;
        bt = bc * tt * tt * tt * tt;

        ed = source_func_lin(bb, bl, bt, tau_d[id * nv * 2 + 2 * lev + 1], diff_ang);

        flw_dn_d[id * (nv + 1) + lev] =
            ed
            + flw_dn_d[id * (nv + 1) + lev + 1]
                  * exp(-(1. / diff_ang) * tau_d[id * nv * 2 + 2 * lev + 1]);
    }
    //
    //  Upward Directed Radiation
    //
    flw_up_d[id * (nv + 1) + 0] = bc * tint * tint * tint * tint; // Lower boundary;
    if (surface == true) {
        flw_up_d[id * (nv + 1) + 0] += bc * Tsurface * Tsurface * Tsurface * Tsurface;
    }
    else {
        // reflecting boundary
        flw_up_d[id * (nv + 1) + 0] += flw_dn_d[id * (nv + 1) + 0];
    }
    for (int lev = 1; lev <= nv; lev++) {

        double eu = 0.0;

        if (tau_d[id * nv * 2 + 2 * (lev - 1) + 1] < 0.0)
            tau_d[id * nv * 2 + 2 * (lev - 1) + 1] = 0.0;

        tb = thtemp[id * (nv + 1) + lev - 1];
        tl = ttemp[id * nv + lev - 1];
        tt = thtemp[id * (nv + 1) + lev];

        bb = bc * tb * tb * tb * tb;
        bl = bc * tl * tl * tl * tl;
        bt = bc * tt * tt * tt * tt;

        eu = source_func_lin(bt, bl, bb, tau_d[id * nv * 2 + 2 * (lev - 1) + 1], diff_ang);

        flw_up_d[id * (nv + 1) + lev] =
            eu
            + flw_up_d[id * (nv + 1) + lev - 1]
                  * exp(-(1. / diff_ang) * tau_d[id * nv * 2 + 2 * (lev - 1) + 1]);
    }

    for (int lev = 0; lev < nv; lev++) {
        if (DeepModel) { //this seems to cause strange problems at TOA, set both factors to 1 for now
            // rup =
            //     (Altitudeh_d[lev + 1] + A) / (Altitude_d[lev] + A); //vertical scaling in divergence
            // rlow = (Altitudeh_d[lev] + A) / (Altitude_d[lev] + A);
            rup  = 1.0;
            rlow = 1.0;
        }
        else {
            rup  = 1.0;
            rlow = 1.0;
        }
        dtemp[id * nv + lev] =
            dtemp[id * nv + lev]
            - (pow(rlow, 2) * (flw_up_d[id * (nv + 1) + lev] - flw_dn_d[id * (nv + 1) + lev])
               - pow(rup, 2)
                     * (flw_up_d[id * (nv + 1) + lev + 1] - flw_dn_d[id * (nv + 1) + lev + 1]))
                  / ((Altitudeh_d[lev] - Altitudeh_d[lev + 1]));
        // dtemp[id * nv + lev] =
        //     dtemp[id * nv + lev]
        //     + gocp
        //           * ((flw_up_d[id * (nv + 1) + lev] - flw_dn_d[id * (nv + 1) + lev])
        //              - (flw_up_d[id * (nv + 1) + lev + 1] - flw_dn_d[id * (nv + 1) + lev + 1]))
        //           / (phtemp[id * (nv + 1) + lev] - phtemp[id * (nv + 1) + lev + 1]);
        if (isnan(dtemp[id * nv + lev])) {
            printf("stop here");
        }
    }
}


__device__ void computetau(double *tau_d,
                           double *pressure_d,
                           double *Rho_d,
                           double *Altitudeh_d,
                           double  kappa_sw,
                           double  kappa_lw,
                           double  n_sw,
                           double  n_lw,
                           double  f_lw,
                           double  ps0,
                           int     id,
                           int     nv) {

    // need to update this to use scaling for non-uniform mixing
    for (int lev = 0; lev < nv; lev++) {
        //***shortwave optical depth across layer***
        // tau_d[id * 2 * nv + lev * 2] =
        //     (tausw / pow(ps0, n_sw))
        //     * (pow(phtemp[id * (nv + 1) + lev], n_sw) - pow(phtemp[id * (nv + 1) + lev + 1], n_sw));
        // the above relation breaks if pressure is not monotonic
        tau_d[id * 2 * nv + lev * 2] =
            n_sw * kappa_sw * pow(pressure_d[id * nv + lev] / ps0, n_sw - 1) * Rho_d[id * nv + lev]
            * (Altitudeh_d[lev + 1] - Altitudeh_d[lev]);

        //***longwave optical depth across layer***
        // tau_d[id * 2 * nv + lev * 2 + 1] =
        //     (taulw * f_lw / ps0) * (phtemp[id * (nv + 1) + lev] - phtemp[id * (nv + 1) + lev + 1])
        //     + (taulw * (1 - f_lw) / pow(ps0, n_lw))
        //           * (pow(phtemp[id * (nv + 1) + lev], n_lw)
        //              - pow(phtemp[id * (nv + 1) + lev + 1], n_lw));
        // above also breaks if pressure is not monotonic
        tau_d[id * 2 * nv + lev * 2 + 1] =
            kappa_lw
            * (1 + n_lw * (1 - f_lw) / f_lw * pow(pressure_d[id * nv + lev] / ps0, n_lw - 1))
            * Rho_d[id * nv + lev] * (Altitudeh_d[lev + 1] - Altitudeh_d[lev]);
    }
}

__global__ void rtm_dual_band(double *pressure_d,
                              double *Rho_d,
                              double *temperature_d,
                              double *flw_up_d,
                              double *flw_dn_d,
                              double *fsw_up_d,
                              double *fsw_dn_d,
                              double *tau_d,
                              double  gravit,
                              double *Cp_d,
                              double *lonlat_d,
                              double *Altitude_d,
                              double *Altitudeh_d,
                              double *phtemp,
                              double *dtemp,
                              double *ttemp,
                              double *thtemp,
                              double  timestep,
                              double  tstar,
                              double  planet_star_dist,
                              double  radius_star,
                              double  diff_ang,
                              double  tint,
                              double  alb,
                              double  kappa_sw,
                              double  kappa_lw,
                              bool    latf_lw_mod,
                              double  kappa_lw_pole,
                              double  n_sw,
                              double  n_lw,
                              double  f_lw,
                              double  incflx,
                              double  ps0,
                              int     num,
                              int     nv,
                              int     nvi,
                              double  A,
                              double  r_orb,
                              double *zenith_angles,
                              double *insol_d,
                              bool    surface,
                              double  Csurf,
                              double *Tsurface_d,
                              double *dTsurf_dt_d,
                              double *surf_flux_d,
                              double *areasT_d,
                              double *ASR_d,
                              double *OLR_d,
                              double *profx_Qheat_d,
                              double *DG_Qheat_d, // internal qheat for debugging
                              double *Rd_d,
                              double  Qheat_scaling,
                              bool    gcm_off,
                              bool    rt1Dmode,
                              bool    DeepModel) {


    //
    //  Description:
    //
    //
    //
    //  Input: .
    //
    //  Output:
    //

    int id = blockIdx.x * blockDim.x + threadIdx.x;

    double coszrs;
    double ps, psm;
    double pp, ptop;

    double xi, xip, xim, a, b;

    if (id < num) {

        for (int lev = 0; lev < nv; lev++) {
            dtemp[id * nv + lev] = 0.0;
        }
        // Calculate pressures and temperatures at interfaces
        for (int lev = 0; lev <= nv; lev++) {
            if (lev == 0) {
                psm = pressure_d[id * nv + 1]
                      - Rho_d[id * nv + 0] * gravit * (-Altitude_d[0] - Altitude_d[1]);
                ps = 0.5 * (pressure_d[id * nv + 0] + psm);

                phtemp[id * nvi + 0] = ps;
                ttemp[id * nv + 0]   = temperature_d[id * nv + 0];
                thtemp[id * nvi + 0] = ttemp[id * nv + 0];
            }
            else if (lev == nv) {
                // pp = pressure_d[id*nv + nv-2] - Rho_d[id*nv + nv-1] * gravit * (2*Altitudeh_d[nv]-Altitude_d[nv-1]-Altitude_d[nv-2]);
                pp = pressure_d[id * nv + nv - 2]
                     + (pressure_d[id * nv + nv - 1] - pressure_d[id * nv + nv - 2])
                           / (Altitude_d[nv - 1] - Altitude_d[nv - 2])
                           * (2 * Altitudeh_d[nv] - Altitude_d[nv - 1] - Altitude_d[nv - 2]);
                if (pp < 0)
                    pp = 0; //prevents pressure at the top from becoming negative
                ptop = 0.5 * (pressure_d[id * nv + nv - 1] + pp);

                phtemp[id * nvi + nv] = ptop;
                thtemp[id * nvi + nv] = temperature_d[id * nv + nv - 1];
            }
            else {
                ttemp[id * nv + lev] = temperature_d[id * nv + lev];
                xi                   = Altitudeh_d[lev];
                xim                  = Altitude_d[lev - 1];
                xip                  = Altitude_d[lev];
                a                    = (xi - xip) / (xim - xip);
                b                    = (xi - xim) / (xip - xim);

                phtemp[id * nvi + lev] =
                    pressure_d[id * nv + lev - 1] * a + pressure_d[id * nv + lev] * b;
                thtemp[id * nvi + lev] =
                    temperature_d[id * nv + lev - 1] * a + ttemp[id * nv + lev] * b;
            }
        }

        // zenith angle
        if (rt1Dmode) {
            coszrs = 0.5;
        }
        else {
            coszrs = zenith_angles[id];
        }

        //hack for daily averaged insolation
        // coszrs = 0.25 * (1 + 1.4 * 0.25 * (1 - 3 * pow(sin(lonlat_d[id * 2 + 1]), 2)));

        // Compute opacities
        double kappa_lw_lat;
        if (latf_lw_mod) {
            //latitude dependence of opacity, for e.g., earth
            kappa_lw_lat =
                kappa_lw + (kappa_lw_pole - kappa_lw) * pow(sin(lonlat_d[id * 2 + 1]), 2);
        }
        else {
            kappa_lw_lat = kappa_lw;
        }
        //computetau(tau_d, phtemp, coszrs, tausw, taulw_lat, n_sw, n_lw, f_lw, ps0, id, nv);
        computetau(tau_d,
                   pressure_d,
                   Rho_d,
                   Altitudeh_d,
                   kappa_sw,
                   kappa_lw_lat,
                   n_sw,
                   n_lw,
                   f_lw,
                   ps0,
                   id,
                   nv);

        for (int lev = 0; lev <= nv; lev++) {
            fsw_up_d[id * nvi + lev] = 0.0;
        }
        for (int lev = 0; lev <= nv; lev++) {
            fsw_dn_d[id * nvi + lev] = 0.0;
        }

        if (coszrs > 0.0) {
            radcsw(phtemp,
                   coszrs,
                   r_orb,
                   dtemp,
                   tau_d,
                   fsw_up_d,
                   fsw_dn_d,
                   Altitude_d,
                   Altitudeh_d,
                   incflx,
                   alb,
                   kappa_sw,
                   ps0,
                   gravit,
                   A,
                   id,
                   nv,
                   insol_d,
                   DeepModel);
        }
        else {
            insol_d[id] = 0;
        }

        if (surface == true) {
            surf_flux_d[id] = fsw_dn_d[id * nvi + 0] - fsw_up_d[id * nvi + 0];
        }

        //calculate ASR for this point
        double rscale;
        if (DeepModel) {
            rscale = (A + Altitudeh_d[nv]) / A;
        }
        else {
            rscale = 1.0;
        }
        ASR_d[id] = fsw_dn_d[id * nvi + nv] * areasT_d[id] * pow(rscale, 2);

        for (int lev = 0; lev <= nv; lev++) {
            flw_up_d[id * nvi + lev] = 0.0;
        }
        for (int lev = 0; lev <= nv; lev++) {
            flw_dn_d[id * nvi + lev] = 0.0;
        }

        radclw(phtemp,
               ttemp,
               thtemp,
               dtemp,
               tau_d,
               flw_up_d,
               flw_dn_d,
               Altitude_d,
               Altitudeh_d,
               diff_ang,
               tint,
               gravit,
               surface,
               Tsurface_d[id],
               A,
               id,
               nv,
               DeepModel);

        if (surface == true) {
            surf_flux_d[id] += flw_dn_d[id * nvi + 0] - flw_up_d[id * nvi + 0];
            Tsurface_d[id] += surf_flux_d[id] * timestep / Csurf
                              + dTsurf_dt_d[id] * timestep; // put dTsurf here temporarily
            if (Tsurface_d[id] < 0)
                Tsurface_d[id] = 0;
        }

        //calculate OLR for this point
        OLR_d[id] = flw_up_d[id * nvi + nv] * areasT_d[id] * pow(rscale, 2);

        for (int lev = 0; lev < nv; lev++) {
            // if (gcm_off) {
            //     temperature_d[id * nv + lev] = ttemp[id * nv + lev]
            //                                    + 1.0 / (Cp_d[id * nv + lev] - Rd_d[id * nv + lev])
            //                                          * dtemp[id * nv + lev] / Rho_d[id * nv + lev]
            //                                          * timestep;
            //     if (temperature_d[id * nv + lev] < 0)
            //         temperature_d[id * nv + lev] = 0.0;
            //     profx_Qheat_d[id * nv + lev] += Qheat_scaling * dtemp[id * nv + lev];
            // }
            // else {
            if (pressure_d[id * nv + lev]
                    + Rd_d[id * nv + lev] / (Cp_d[id * nv + lev] - Rd_d[id * nv + lev])
                          * dtemp[id * nv + lev] * timestep
                < 0) {
                //trying to prevent too much cooling resulting in negative pressure in dyn core
                dtemp[id * nv + lev] = -pressure_d[id * nv + lev] / timestep;
            }
            DG_Qheat_d[id * nv + lev] = dtemp[id * nv + lev];
            profx_Qheat_d[id * nv + lev] += Qheat_scaling * dtemp[id * nv + lev];
            if (isnan(profx_Qheat_d[id * nv + lev])) {
                printf("stop here");
            }
            // }
        }
    }
}



//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


// Calculates the IR band Rosseland mean opacity (local T) according to the
// Freedman et al. (2014) fit and coefficents

__device__ void kernel_k_Ross_Freedman(double Tin, double Pin, double met, double &k_IR) {
    // dependcies
    //// powl from math
    //// log10l from math        
    //// atan from math
    //// onedivpi -> namespace constants::onedivpi

    // Input:
    // T - Local gas temperature [K]
    // P - Local gas pressure [pa]
    // met - Local metallicity [M/H] (log10l from solar, solar [M/H] = 0.0)

    // Call by reference (Input&Output):
    // k_IR - IR band Rosseland mean opacity [m2 kg-1]

    const double pi = atan((double)(1)) * 4;
    const double onedivpi = 1.0 / pi;

    // Coefficent parameters for Freedman et al. (2014) table fit
    double c1 = 10.602;
    double c2 = 2.882;
    double c3 = 6.09e-15;
    double c4 = 2.954;
    double c5 = -2.526;
    double c6 = 0.843;
    double c7 = -5.490;
    double c8_l = -14.051, c8_h = 82.241;
    double c9_l = 3.055, c9_h = -55.456;
    double c10_l = 0.024, c10_h = 8.754;
    double c11_l = 1.877, c11_h = 0.7048;
    double c12_l = -0.445, c12_h = -0.0414;
    double c13_l = 0.8321, c13_h = 0.8321;

    // work variables
    double k_lowP;
    double k_hiP;
    double Tl10;
    double Pl10;

    // start operations



    Tl10 = log10((double)(Tin));
    Pl10 = log10(Pin * 10.0); // Convert to dyne cm-2 and log 10

    // Low pressure expression
    k_lowP = c1 * atan((double)(Tl10 - c2)) -
        (c3 / (Pl10 + c4)) * exp((double)(pow((double)(Tl10 - c5), 2.0))) + c6 * met + c7;

    // De log10l
    k_lowP = pow((double)(10.0), k_lowP);

    // Temperature split for coefficents = 800 K
    if (Tin <= 800.0)
    {
        k_hiP = c8_l + c9_l * Tl10 + c10_l * pow((double)(Tl10), 2.0) +
            Pl10 * (c11_l + c12_l * Tl10) +
            c13_l * met * (0.5 + onedivpi * atan((double)((Tl10 - ((double)2.5)) / (double)0.2)));
    }
    else
    {
        k_hiP = c8_h + c9_h * Tl10 +
            c10_h * pow((double)(Tl10), 2.0) + Pl10 * (c11_h + c12_h * Tl10) +
            c13_h * met * (0.5 + onedivpi * atan((double)((Tl10 - ((double)2.5)) / (double)0.2)));
    }

    // De log10l
    k_hiP = pow((double)(10.0), k_hiP);

    // Total Rosseland mean opacity - converted to m2 kg-1
    k_IR = (k_lowP + k_hiP) / ((double)10.0);

    // Avoid divergence in fit for large values
    if (k_IR > 1.0e10)
    {
        k_IR = 1.0e10;
    }
}



///////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////




__device__  void linear_log_interp(int id,
    int i,
    int nlay,
    int nlay1,
    double *Altitude_d,
    double *Altitudeh_d, 
    double *Tl,
    double *Te__df_e) {
    // dependcies
    //// pow from math
    //// log10 from math

    // work variables
    double eAlt;
    double lTl1;
    double lTl2;
    double lAlt1;
    double lAlt2;
    double norm;

    // start operations
    eAlt = log10((double)(Altitudeh_d[i +1])) ;
    lAlt1 = log10((double)(Altitude_d[ i ])) ;
    lAlt2 = log10((double)(Altitude_d[ i + 1])) ;
    lTl1 = log10((double)(Tl[id * nlay + i + 1])) ;
    lTl2 = log10((double)(Tl[id * nlay + i])) ;

    norm = (1.0) / (lAlt2 - lAlt1);

    Te__df_e[id * nlay1 + i + 1] = pow((double)(10.0), ( (lTl1 * (lAlt2 - eAlt) + lTl2 * (eAlt - lAlt1)) * norm) );
}

///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/*
__device__  void linear_interp(int id,
    int i,
    int nlay,
    int nlay1,
    double *Altitude_d,
    double *Altitudeh_d, 
    double *Tl,
    double *Te__df_e) {
    // dependcies
    //// pow from math
    //// log10 from math

    // work variables
    
    double ContributionFactorFromBelow;
    double ContributionFactorFromAbove;
    


    // start operations
    ContributionFactorFromBelow = (Altitudeh_d[lev] - Altitude_d[lev]) / (Altitude_d[lev - 1] - Altitude_d[lev]);
    ContributionFactorFromAbove =  (Altitudeh_d[lev] - Altitude_d[lev - 1]) / (Altitude_d[lev] - Altitude_d[lev - 1]);

    Te__df_e[id * nlay1 + i] = Tl[id * nlay + i-1] * ContributionFactorFromBelow    +   Tl[id * nlay + i] * ContributionFactorFromAbove;

}

*/

///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

__device__ void tau_struct(int id,
    int nlev,
    double *Rho_d,
    double *Altitudeh_d,
    double *kRoss,
    int nchan,
    int channel,
    double *tau_struc_e) {

    // work variables
    double tau_sum;
    double tau_lay;
    double delPdelAlt;
    int level;

    // running sum of optical depth
    tau_sum = 0.0;

    // start operations
    //  Upper most tau_struc is given by some low pressure value (here 1e-9 bar = 1e-4 pa)
    //dP = (p_half(1) - 1e-4)
    //tau_lay = (kRoss(1) * dP) / grav
    //tau_sum = tau_sum + tau_lay
    tau_struc_e[id*(nlev+1) + nlev+1] = tau_sum;

    // Integrate from top to bottom    

    for (level = nlev-1; level > -1; level--) 
    {
        // Pressure difference between layer edges
        delPdelAlt = (Altitudeh_d[level + 1] - Altitudeh_d[level + 0]);

        if (delPdelAlt < 0.0000001)
        {
            delPdelAlt = 0.0000001;
        }

        // Optical depth of layer assuming hydrostatic equilibirum  //////////////////// kappa * rho * delta height  (old: kappa *delP/gravity)
        tau_lay = kRoss[id*nlev*nchan + channel * nlev + level] * delPdelAlt * Rho_d[id*nlev  + level];

        // Add to running sum
        tau_sum = tau_sum + tau_lay;

        // Optical depth structure is running sum
        tau_struc_e[id*(nlev+1) + level] = tau_sum;
    }

}

///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

__device__  void sw_grey_down(int id,
    int nlay1,
    double solar,
    double *solar_tau,
    double *sw_down__df_e,
    double *mu) {
    // dependencies
    //// expl -> math

    // start operations
    for (int i = nlay1-1; i >-1; i--)
    {
        sw_down__df_e[id * nlay1 + i] = solar * mu[id] * exp(-solar_tau[id * nlay1 + i] / mu[id]);
    }

}

///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

__device__  void lw_grey_updown_linear(int id,
    int nlay,
    int nlay1,
    double *be__df_e,
    double *tau_IRe__df_e,
    double *lw_up__df_e,
    double *lw_down__df_e,
    double *dtau__dff_l,
    double *del__dff_l,
    double *edel__dff_l,
    double *e0i__dff_l,
    double *e1i__dff_l,
    double *Am__dff_l,
    double *Bm__dff_l,
    double *lw_up_g__dff_e,
    double *lw_down_g__dff_e) {
    // dependencies
    //// expll -> math
    //// atan -> math

    const double pi = atan((double)(1)) * 4;
    const double  twopi = 2.0 * pi;

    // Work variables and arrays
    int k, g;



    //Gauss quadrature variables
    const int gauss_ng = 2;
    double uarr[gauss_ng];
    double w[gauss_ng];

    uarr[0] = 0.21132487;
    uarr[1] = 0.78867513;
    w[0] = 0.5;
    w[1] = 0.5;

    for (k = nlay-1; k >-1; k--)
    {
        dtau__dff_l[id*nlay + k ] = (tau_IRe__df_e[id*nlay1 + k] - tau_IRe__df_e[id*nlay1 + k + 1]);
    }

    // Zero the flux arrays
    for (k = 0; k < nlay1; k++)
    {
        lw_down__df_e[id*nlay1+k] = 0.0;
        lw_up__df_e[id*nlay1 + k] = 0.0;
    }

    // Start loops to integrate in mu space
    for (g = 0; g < gauss_ng; g++)
    {
        // Prepare loop
        for (k = 0; k < nlay; k++)
        {
            // Olson & Kunasz (1987) parameters
            del__dff_l[id*nlay + k] = dtau__dff_l[id * nlay + k] / uarr[g];
            edel__dff_l[id * nlay + k] = exp((double)(-del__dff_l[id * nlay + k]));
            e0i__dff_l[id * nlay + k] = 1.0 - edel__dff_l[id * nlay + k];
            e1i__dff_l[id * nlay + k] = del__dff_l[id * nlay + k] - e0i__dff_l[id * nlay + k];

            Am__dff_l[id * nlay + k] = e0i__dff_l[id * nlay + k] - e1i__dff_l[id * nlay + k] / del__dff_l[id * nlay + k]; // Am[k] = Gp[k], just indexed differently
            Bm__dff_l[id * nlay + k] = e1i__dff_l[id * nlay + k] / del__dff_l[id * nlay + k]; // Bm[k] = Bp[k], just indexed differently

        }

        // Peform downward loop first
        // Top boundary condition
        lw_down_g__dff_e[id * nlay1 +  nlay+1] = 0.0;
        for (k = nlay-1; k > -1; k--)
        {
            lw_down_g__dff_e[id * nlay1 +  k] = lw_down_g__dff_e[id * nlay1 +  k + 1] * edel__dff_l[id * nlay + k] + 
            Am__dff_l[id * nlay + k] * be__df_e[id * nlay1 + k+1] + Bm__dff_l[id * nlay + k] * be__df_e[id * nlay1 + k]; // TS intensity
        }


        // Peform upward loop
        // Lower boundary condition
        
        lw_up_g__dff_e[id * nlay1 + 0  ] = be__df_e[id * nlay1 + 0];
        for (k = 0; k < nlay; k++)
        {
            lw_up_g__dff_e[id * nlay1 + k+1] = lw_up_g__dff_e[id * nlay1 + k] * edel__dff_l[id * nlay + k] +
                Bm__dff_l[id * nlay + k] * be__df_e[id * nlay1 + k + 1] + Am__dff_l[id * nlay + k] * be__df_e[id * nlay1 + k]; // TS intensity
        }

        // Sum up flux arrays with Gauss weights and points
        for (k = nlay1-1; k > -1; k--)
        {
            lw_down__df_e[id * nlay1 + k] = lw_down__df_e[id * nlay1 + k ] + lw_down_g__dff_e[id * nlay1 +  k] * w[g] * uarr[g];
            lw_up__df_e[id * nlay1 + k ] = lw_up__df_e[id * nlay1 + k ] + lw_up_g__dff_e[id * nlay1 + k] * w[g] * uarr[g];
        }
    }

    for (k = 0; k < nlay1; k++)
    {
        lw_down__df_e[id * nlay1 + k] = twopi * lw_down__df_e[id * nlay1 + k];
        lw_up__df_e[id * nlay1 + k] = twopi * lw_up__df_e[id * nlay1 + k];
    }

}


///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////



    __device__ void Kitzmann_TS_noscatt(int id,
        const int nlay,
        const int nlay1,
        double *Altitude_d,
        double *Altitudeh_d,
        double *Rho_d,
        double *Tl,
        double *pl,
        double *pe,
        double *k_V_3_nv_d,
        double *k_IR_2_nv_d,
        double *Beta_V_3_d,
        double *Beta_2_d,
        double *net_F_nvi_d,
        double *mu_s,
        double Finc,
        double Fint,
        double grav,
        double AB_d,
        //Kitzman working variables
        double *tau_Ve__df_e,
        double *tau_IRe__df_e,
        double *Te__df_e,
        double *be__df_e, 
        double *sw_down__df_e,
        double *sw_down_b__df_e,
        double *sw_up__df_e,
        double *lw_down__df_e,
        double *lw_down_b__df_e,
        double *lw_up__df_e,
        double *lw_up_b__df_e,
        double *lw_net__df_e,
        double *sw_net__df_e,
        // lw_grey_updown_linear working variables
        double *dtau__dff_l,
        double *del__dff_l, 
        double *edel__dff_l,
        double *e0i__dff_l,
        double *e1i__dff_l,
        double *Am__dff_l,
        double *Bm__dff_l,
        double *lw_up_g__dff_e,
        double *lw_down_g__dff_e) 
    {
        // dependcies
        //// powll -> include math
        //// log10f -> include math
        //// nlay -> layers
        //// nlay1 -> layers +1
        //// linear_log_interp -> function
        //// tau_struct -> function
        //// sw_grey_down -> function
        //// lw_grey_updown_linear -> function

        const double pi = atan(1.0) * 4;
        const double twopi = 2.0 * pi;
        const double StBC = 5.670374419e-8;

        // work variables
        double Finc_B;
       
        // start operation
        double out;
                
        /*
        // Find temperature at layer edges through linear interpolation and extrapolation
        for (int i = nlay - 2; i > -1; i--)
        {

            
            linear_log_interp(id,
                i,
                nlay,
                nlay1,
                Altitude_d,
                Altitudeh_d,
                Tl,
                Te__df_e);

            //printf("---------\n");
            if (pe[id*nlay1 + i + 1] == 0.0)
            {                
                    printf("pe[id*nlay1 + i + 1] is zer0\n");
            }
             if (pe[id*nlay1 + i + 1] < 0.0)
            {                
                    printf("pe[id*nlay1 + i + 1] is negative\n");
            }

            if (pe[id*nlay1 + i + 1] > pe[id * nlay1 + i ])
            {                
                    printf("pe[id*nlay1 + i + 1] > pe[id * nlay1 + i ] \n");
            }
            if (pl[id*nlay + i + 1] > pl[id * nlay + i ])
            {                
                    printf("pl[id*nlay + i + 1] > pl[id * nlay + i ] \n");
            }

            if (pl[id * nlay + i + 1] == 0.0)
            {                
                    printf("pl[id * nlay + i + 1] is zer0\n");
            }
            if (pl[id * nlay + i + 1] < 0.0)
            {                
                    printf("pl[id * nlay + i + 1] is negative\n");
            }
            if (pl[id * nlay + i ] == 0.0)
            {                
                    printf("pl[id * nlay + i ] is zer0\n");
            }
            if (pl[id * nlay + i ] < 0.0)
            {                
                    printf("pl[id * nlay + i ] is negative\n");
            }
            if (Tl[id * nlay + i+1] == 0.0)
            {                
                    printf("Tl[id * nlay + i+1] is zer0\n");
            }
            if (Tl[id * nlay + i+1] < 0.0)
            {                
                    printf("Tl[id * nlay + i+1] is negative\n");
            }
            if (Tl[id * nlay + i] == 0.0)
            {                
                    printf("Tl[id * nlay + i] is zer0\n");
            }
            if (Tl[id * nlay + i] < 0.0)
            {                
                    printf("Tl[id * nlay + i] is negative\n");
            }

            if (Te__df_e[id * nlay1 + i] == 0.0 && i != 0)
            {                
                   
                    printf("Te__df_e[id * nlay1 + i] is zer0 at level %d\n", i);
            }
            if (Te__df_e[id * nlay1 + i] < 0.0)
            {                
                    printf("Te__df_e[id * nlay1 + i] is negative\n");
            }

            //printf("level:%d   \n", i);
            
            //printf("---------\n");
        }
        
        Te__df_e[id * nlay1 + nlay1] = Tl[id * nlay + nlay] + (Altitudeh_d[ nlay1-1] - Altitudeh_d[ nlay1 ]) / 
            ( Altitude_d[ nlay1 ] - Altitudeh_d[ nlay-1] ) * (Tl[id * nlay + nlay] - Te__df_e[id * nlay1 + nlay]);

        Te__df_e[id * nlay1 + 0] = Tl[id * nlay + 0] + (Altitudeh_d[ 1] - Altitudeh_d[ 0]) /
            (Altitude_d[ 1] - Altitude_d[ 0]) *
            (Tl[id * nlay + 0] - Te__df_e[id * nlay1 + 1]);

        */

        // Shortwave fluxes
        for (int i = 0; i < nlay1; i++)
        {
            sw_down__df_e[id * nlay1 + i] = 0.0;
            sw_up__df_e[id * nlay1 + i] = 0.0;
        }
        for (int channel = 0; channel < 3; channel++)
        {
           

            // Incident flux in band
            if (mu_s[id]>0.0)
            {
                // Find the opacity structure
                tau_struct(id,
                    nlay,
                    Rho_d,
                    Altitudeh_d,
                    k_V_3_nv_d,
                    3,
                    channel,
                    tau_Ve__df_e);
                
                //printf("tau_struct finished\n");

                for (int i = 0; i < nlay1; i++)
                {
                    if (tau_Ve__df_e[id * nlay1 + i] < 0.0)
                    {                
                        printf("tau_Ve__df_e[id * nlay1 + i] is negative at level %d and in channel %d\n", i, channel);
                    }
                    if (isnan(tau_Ve__df_e[id * nlay1 + i]))
                    {                
                        printf("tau_Ve__df_e[id * nlay1 + i] is negativeat level %d and in channel %d\n", i, channel);
                    }
                }
                for (int i = 0; i < nlay; i++)
                {
                    if (Rho_d[id * nlay + i] < 0.0)
                    {                
                        printf("Rho_d[id * nlay + i] is negative at level %d and in channel %d\n", i, channel);
                    }
                    if (isnan(Rho_d[id * nlay + i]))
                    {                
                        printf("Rho_d[id * nlay + i] is negativeat level %d and in channel %d\n", i, channel);
                    }
                }
                

                Finc_B = Finc * Beta_V_3_d[id * 3 + channel];

                
                sw_grey_down(id,
                    nlay1,
                    Finc_B,
                    tau_Ve__df_e,
                    sw_down_b__df_e,
                    mu_s);

                //printf("sw_grey_down finished\n");            

            } else
            {
                printf("sw darkside\n");

                for (int i = nlay1-1; i >-1; i--)
                {
                    sw_down_b__df_e[id * nlay1 + i] = 0.0;
                    tau_Ve__df_e[id * nlay1 + i] = 1.0;
                }
                
            }
            

            // Sum all bands
            for (int i = 0; i < nlay1; i++)
            {
                sw_down__df_e[id * nlay1 + i] = sw_down__df_e[id * nlay1 + i] + sw_down_b__df_e[id * nlay1 + i];
            }
        }

        // Long wave two-stream fluxes
        for (int i = 0; i < nlay1; i++)
        {
            lw_down__df_e[id * nlay1 + i] = 0.0;
            lw_up__df_e[id * nlay1 + i] = 0.0;
            
            lw_up_b__df_e[id * nlay1 + i] = 0.0;
            lw_down_b__df_e[id * nlay1 + i] = 0.0;
        }
        for (int channel = 0; channel < 2; channel++)
        {
            // Find the opacity structure
            tau_struct(id,
            nlay,
            Rho_d,
            Altitudeh_d,
            k_IR_2_nv_d,
            2,
            channel,
            tau_IRe__df_e);

            for (int i = 0; i < nlay1; i++)
                {
                    if (tau_Ve__df_e[id * nlay1 + i] < 0.0)
                    {                
                        printf("tau_Ve__df_e[id * nlay1 + i] is negative at level %d and in channel %d\n", i, channel);
                    }
                    if (isnan(tau_Ve__df_e[id * nlay1 + i]))
                    {                
                        printf("tau_Ve__df_e[id * nlay1 + i] is negativeat level %d and in channel %d\n", i, channel);
                    }
                }

            //printf("tau_struct finished\n");

            // Blackbody fluxes (note divide by pi for correct units)
            for (int i = 0; i < nlay1; i++)
            {
                be__df_e[id * nlay1 + i] = StBC * pow((Te__df_e[id * nlay1 + i]), 4.0) / pi * Beta_2_d[id * 2 + channel];

            }

            // Calculate lw flux
            lw_grey_updown_linear(id,
                nlay,
                nlay1,
                be__df_e,
                tau_IRe__df_e,
                lw_up_b__df_e,
                lw_down_b__df_e,
                dtau__dff_l,
                del__dff_l,
                edel__dff_l,
                e0i__dff_l,
                e1i__dff_l,
                Am__dff_l,
                Bm__dff_l,
                lw_up_g__dff_e,
                lw_down_g__dff_e);

            //printf("lw_grey_updown_linear finished\n");
            


            // Sum all bands
            for (int i = 0; i < nlay1; i++)
            {
                lw_up__df_e[id * nlay1 + i] = lw_up__df_e[id * nlay1 + i] + lw_up_b__df_e[id * nlay1 + i];
                lw_down__df_e[id * nlay1 + i] = lw_down__df_e[id * nlay1 + i] + lw_down_b__df_e[id * nlay1 + i];
            }

        }

        // Net fluxes
        for (int i = 0; i < nlay1; i++)
        {
            lw_net__df_e[id * nlay1 + i] = lw_up__df_e[id * nlay1 + i] - lw_down__df_e[id * nlay1 + i];
            sw_net__df_e[id * nlay1 + i] = sw_up__df_e[id * nlay1 + i] - sw_down__df_e[id * nlay1 + i];
            net_F_nvi_d[id * nlay1 + i] = lw_net__df_e[id * nlay1 + i] + sw_net__df_e[id * nlay1 + i];
        }

        net_F_nvi_d[id * nlay1 + 0] = Fint;

        //printf("Kitzmann finished\n");



    }


///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


__global__ void rtm_picket_fence(double *pressure_d,
                              double *temperature_d,
                              double *Rho_d,
                              double  gravit,
                              double *Cp_d,
                              double *lonlat_d,
                              double *Altitude_d,
                              double *Altitudeh_d,
                              double r_rob,
                              double *phtemp,
                              double *dtemp,
                              double  timestep,
                              double  tint,
                              double  alb,
                              double  kappa_lw,
                              bool    latf_lw_mod,
                              double  kappa_lw_pole,              
                              double  incflx,
                              int     num,
                              int     nv,
                              int     nvi,
                              double  A,
                              double *zenith_angles,
                              double *insol_d,
                              bool    surface,
                              double *Tsurface_d, // ??? relevant
                              double *surf_flux_d,  // ??? relevant
                              double *areasT_d,
                              double *ASR_d, // ??? relevant
                              double *OLR_d, // ??? relevant
                              double *profx_Qheat_d,
                              double *DG_Qheat_d, // internal qheat for debugging
                              double *Rd_d,
                              double  Qheat_scaling,
                              double  met,
                              double *k_IR_2_nv_d,
                              double *k_V_3_nv_d,
                              double *gam_V_3_d,
                              double *gam_1_d,
                              double *gam_2_d,
                              double *Beta_V_3_d,
                              double *Beta_2_d,
                              double *net_F_nvi_d,    
                              double *AB_d,
                              //Kitzman working variables
                              double *tau_Ve__df_e,
                              double *tau_IRe__df_e,
                              double *Te__df_e,
                              double *be__df_e, 
                              double *sw_down__df_e,
                              double *sw_down_b__df_e,
                              double *sw_up__df_e,
                              double *lw_down__df_e,
                              double *lw_down_b__df_e,
                              double *lw_up__df_e,
                              double *lw_up_b__df_e,
                              double *lw_net__df_e,
                              double *sw_net__df_e,
                              // lw_grey_updown_linear working variables
                              double *dtau__dff_l,
                              double *del__dff_l,
                              double *edel__dff_l,
                              double *e0i__dff_l,
                              double *e1i__dff_l,
                              double *Am__dff_l,
                              double *Bm__dff_l,
                              double *lw_up_g__dff_e, 
                              double *lw_down_g__dff_e,
                              //general model parameters
                              bool    rt1Dmode,
                              bool    DeepModel) {


    //
    //  Description:
    //
    //
    //
    //  Input: .
    //
    //  Output:
    //

    int id = blockIdx.x * blockDim.x + threadIdx.x;

    double coszrs;
    double ps, psm;
    double pp, ptop;

    
    double ContributionFactorFromBelow;
    double ContributionFactorFromAbove;

    if (id < num) {

        for (int channel = 0; channel < 3; channel++) {
            if (gam_V_3_d[id*3 + channel]==0.0){
                printf("gam_V_3_d contains 0 in the channel\n",  channel);
                
                for (int chan = 0; chan < 3; chan++) {
                    //printf("gam_V_3_d contains 0 in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d channel:%d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, chan, &gam_V_3_d[id*3 + chan]);
                }
               // __threadfence();         // ensure store issued before trap
                //asm("trap;");            // kill kernel with error           
            }

            if (isnan(gam_V_3_d[id*3 + channel])){
                printf("gam_V_3_d contains NaN in the channel\n",  channel);
                for (int chan = 0; chan < 3; chan++) {
                    //printf("gam_V_3_d contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d channel:%d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, chan, &gam_V_3_d[id*3 + chan]);
                }
                __threadfence();         // ensure store issued before trap
                asm("trap;");            // kill kernel with error           
            }
        }
        
           

        for (int lev = 0; lev < nv; lev++) {
            /*
            if (isnan(pressure_d[id*nv + lev]) ){
                for (int level = 0; level < nv; level++) {
                    printf("pressure_d contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d lev:%d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, level, &pressure_d[id*nv + level]);
                }
                __threadfence();         // ensure store issued before trap
                asm("trap;");            // kill kernel with error           
            }
            
            if (isnan(temperature_d[id*nv + lev]) && lev !=0){
                printf("temperature_d contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d lev:%d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, lev, &temperature_d[id*nv + lev]);
                __threadfence();         // ensure store issued before trap
                asm("trap;");            // kill kernel with error           
            }
            
            if (temperature_d[id*nv + lev]==0 && lev !=0){
                printf("temperature_d contains 0 in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d lev:%d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, lev, &temperature_d[id*nv + lev]);
                __threadfence();         // ensure store issued before trap
                asm("trap;");            // kill kernel with error           
            }
            
            if (pressure_d[id*nv + lev]==0 && lev !=0){
                for (int lev = 0; lev < nv; lev++) {
                    printf("pressure_d contains 0 in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d lev:%d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, lev, &pressure_d[id*nv + lev]);
                }
                __threadfence();         // ensure store issued before trap
                asm("trap;");            // kill kernel with error           
            }
            if (pressure_d[id*nv + lev] < 0){
                for (int lev = 0; lev < nv; lev++) {
                    printf("pressure_d contains negative pressure in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d lev:%d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, lev, &pressure_d[id*nv + lev]);
                }
                __threadfence();         // ensure store issued before trap
                asm("trap;");            // kill kernel with error           
            }

            
             
            
            */
            dtemp[id * nv + lev] = 0.0;


        }
        
        // Calculate pressures and temperatures at interfaces
        /*
        for (int lev = 0; lev <= nv; lev++) {
            if (lev == 0) {
                psm = pressure_d[id * nv + 1]
                      - Rho_d[id * nv + 0] * gravit * (-Altitude_d[0] - Altitude_d[1]);
                ps = 0.5 * (pressure_d[id * nv + 0] + psm);

                phtemp[id * nvi + 0] = ps;
                
            }
            else if (lev == nv) {
                pp = pressure_d[id * nv + nv - 2]
                     + (pressure_d[id * nv + nv - 1] - pressure_d[id * nv + nv - 2])
                           / (Altitude_d[nv - 1] - Altitude_d[nv - 2])
                           * (2 * Altitudeh_d[nv] - Altitude_d[nv - 1] - Altitude_d[nv - 2]);
                if (pp < 0)
                    pp = 0; //prevents pressure at the top from becoming negative
                ptop = 0.5 * (pressure_d[id * nv + nv - 1] + pp);

                phtemp[id * nvi + nv] = ptop;
               
            }
            else {
                
                xi                   = Altitudeh_d[lev];
                xim                  = Altitude_d[lev - 1];
                xip                  = Altitude_d[lev];
                a                    = (xi - xip) / (xim - xip);
                b                    = (xi - xim) / (xip - xim);

                phtemp[id * nvi + lev] = pressure_d[id * nv + lev - 1] * a + pressure_d[id * nv + lev] * b;
                
            }
        }
        */

       double ContributionFactorFromBelow;
       double ContributionFactorFromAbove;

        for (int lev = 0; lev <= nv; lev++) {
            if (lev == 0) {
                psm = pressure_d[id * nv + 1]
                      - Rho_d[id * nv + 0] * gravit * (-Altitude_d[0] - Altitude_d[1]);
                ps = 0.5 * (pressure_d[id * nv + 0] + psm);

                phtemp[id * nvi + 0] = ps;

                Te__df_e[id * nvi + 0] = tint;
                
            }
            else if (lev == nv) {
                pp = pressure_d[id * nv + nv - 2]
                     + (pressure_d[id * nv + nv - 1] - pressure_d[id * nv + nv - 2])
                           / (Altitude_d[nv - 1] - Altitude_d[nv - 2])
                           * (2 * Altitudeh_d[nv] - Altitude_d[nv - 1] - Altitude_d[nv - 2]);
                if (pp < 0)
                    pp = 0; //prevents pressure at the top from becoming negative
                ptop = 0.5 * (pressure_d[id * nv + nv - 1] + pp);

                phtemp[id * nvi + nv] = ptop;

                pp = temperature_d[id * nv + nv - 2]
                     + (temperature_d[id * nv + nv - 1] - temperature_d[id * nv + nv - 2])
                           / (Altitude_d[nv - 1] - Altitude_d[nv - 2])
                           * (2 * Altitudeh_d[nv] - Altitude_d[nv - 1] - Altitude_d[nv - 2]);
                if (pp < 0)
                    pp = 0; //prevents temperature at the top from becoming negative
                ptop = 0.5 * (temperature_d[id * nv + nv - 1] + pp);

                Te__df_e[id * nvi + nv] = ptop;
               
            }
            else if (lev == nvi) {
                pp = pressure_d[id * nv + nv - 2]
                     + (pressure_d[id * nv + nv - 1] - pressure_d[id * nv + nv - 2])
                           / (Altitude_d[nv - 1] - Altitude_d[nv - 2])
                           * (2 * Altitudeh_d[nvi] - Altitude_d[nv - 1] - Altitude_d[nv - 2]);
                if (pp < 0)
                    pp = 0; //prevents pressure at the top from becoming negative
                ptop = 0.5 * (pressure_d[id * nv + nv]  + pp);

                phtemp[id * nvi + nvi] = ptop;

                pp = temperature_d[id * nv + nv - 2]
                     + (temperature_d[id * nv + nv - 1] - temperature_d[id * nv + nv - 2])
                           / (Altitude_d[nv - 1] - Altitude_d[nv - 2])
                           * (2 * Altitudeh_d[nvi] - Altitude_d[nv - 1] - Altitude_d[nv - 2]);
                if (pp < 0)
                    pp = 0; //prevents temperature at the top from becoming negative
                ptop = 0.5 * (temperature_d[id * nv + nv ] + pp);

                Te__df_e[id * nvi + nvi] = ptop;
               
            }
            else {
                
                
                ContributionFactorFromBelow = (Altitudeh_d[lev] - Altitude_d[lev]) / (Altitude_d[lev - 1] - Altitude_d[lev]);
                ContributionFactorFromAbove =  (Altitudeh_d[lev] - Altitude_d[lev - 1]) / (Altitude_d[lev] - Altitude_d[lev - 1]);

                phtemp[id * nvi + lev] = pressure_d[id * nv + lev - 1] * ContributionFactorFromBelow + pressure_d[id * nv + lev] * ContributionFactorFromAbove;

                Te__df_e[id * nvi + lev] = temperature_d[id * nv + lev - 1] * ContributionFactorFromBelow + temperature_d[id * nv + lev] * ContributionFactorFromAbove;
                
            }
        }

        for (int lev = 0; lev < nvi; lev++)
        {
            if (Te__df_e[id * nvi + lev] == 0.0)
            {
                printf("Te__df_e has a zero at the level:%u\n", lev);
                __threadfence();         // ensure store issued before trap
                asm("trap;");            // kill kernel with error
            }
            if (isnan(Te__df_e[id * nvi + lev]))
            {
                printf("Te__df_e has a NaN at the level:%u\n", lev);
                __threadfence();         // ensure store issued before trap
                asm("trap;");            // kill kernel with error
            }
            if (phtemp[id * nvi + lev] == 0.0)
            {
                printf("phtemp has a zero at the level:%u\n", lev);
                __threadfence();         // ensure store issued before trap
                asm("trap;");            // kill kernel with error
            }
            if (isnan(phtemp[id * nvi + lev]))
            {
                printf("phtemp has a NaN at the level:%u\n", lev);
                __threadfence();         // ensure store issued before trap
                asm("trap;");            // kill kernel with error
            }
            
        }
        for (int lev = 0; lev < nv; lev++)
        {
            if (temperature_d[id * nv + lev] == 0.0)
            {
                printf("temperature_d has a zero at the level:%u\n", lev);
                __threadfence();         // ensure store issued before trap
                asm("trap;");            // kill kernel with error
            }
            if (isnan(temperature_d[id * nv + lev]))
            {
                printf("temperature_d has a NaN at the level:%u\n", lev);
                __threadfence();         // ensure store issued before trap
                asm("trap;");            // kill kernel with error
            }
            if (pressure_d[id * nv + lev] == 0.0)
            {
                printf("pressure_d has a zero at the level:%u\n", lev);
                __threadfence();         // ensure store issued before trap
                asm("trap;");            // kill kernel with error
            }
            if (isnan(pressure_d[id * nv + lev]))
            {
                printf("pressure_d has a NaN at the level:%u\n", lev);
                __threadfence();         // ensure store issued before trap
                asm("trap;");            // kill kernel with error
            }
            
        }
        

        // zenith angle
        if (rt1Dmode) {
            coszrs = 0.5;
        }
        else {
            coszrs = zenith_angles[id];
        }

        // kappa calculation loop here if using non-constant kappa
        for (int level = 0; level < nv; level++)
        {

            kernel_k_Ross_Freedman(temperature_d[id * nv + level],
                pressure_d[id * nv + level],
                met,
                k_IR_2_nv_d[id * nv * 2 + 0 * nv + level]);

            // Find the visual Rosseland mean opacity from gam_V_3_d


            for (int channel = 0; channel < 3; channel++)
            {
                k_V_3_nv_d[id * nv * 3 + channel * nv + level] = k_IR_2_nv_d[id * nv * 2 + 0 * nv + level] * gam_V_3_d[id * 3 + channel];
            }


            // Find the IR Rosseland mean opacity in each IR picket fence band
            // Note: 2nd band done first here to avoid overwrite

            if (latf_lw_mod) {

                //latitude dependence of opacity, for e.g., earth           

                k_IR_2_nv_d[id * nv * 2 + 1 * nv + level] = k_IR_2_nv_d[id * nv * 2 + 0 * nv + level] * gam_2_d[id];
                k_IR_2_nv_d[id * nv * 2 + 0 * nv + level] = k_IR_2_nv_d[id * nv * 2 + 0 * nv + level] * gam_1_d[id];

                k_IR_2_nv_d[id * nv * 2 + 1 * nv + level] = k_IR_2_nv_d[id * nv * 2 + 1 * nv + level] + 
                    k_IR_2_nv_d[id * nv * 2 + 1 * nv + level]*(kappa_lw_pole / kappa_lw) * pow(sin(lonlat_d[id * 2 + 1]), 2);
                k_IR_2_nv_d[id * nv * 2 + 0 * nv + level] = k_IR_2_nv_d[id * nv * 2 + 0 * nv + level] +
                    k_IR_2_nv_d[id * nv * 2 + 0 * nv + level]*(kappa_lw_pole / kappa_lw) * pow(sin(lonlat_d[id * 2 + 1]), 2);
            }
            else {
                k_IR_2_nv_d[id * nv * 2 + 1 * nv + level] = k_IR_2_nv_d[id * nv * 2 + 0 * nv + level] * gam_2_d[id];
                k_IR_2_nv_d[id * nv * 2 + 0 * nv + level] = k_IR_2_nv_d[id * nv * 2 + 0 * nv + level] * gam_1_d[id];
            }

            
        }

        // !! Radiation - Comment in what scheme you want to use - Heng model won't work!
        
        if (zenith_angles[id] > 0.0) {
            insol_d[id] = incflx * pow(r_rob,-2);
            double flux_top = insol_d[id] * (1-alb);

            Kitzmann_TS_noscatt(id,
                nv,
                nvi,
                Altitude_d,
                Altitudeh_d,
                Rho_d,
                temperature_d,
                pressure_d,
                phtemp,
                k_V_3_nv_d,
                k_IR_2_nv_d,
                Beta_V_3_d,
                Beta_2_d,
                net_F_nvi_d,
                zenith_angles,
                flux_top,
                tint,
                gravit,
                AB_d[id],
                //Kitzman working variables
                tau_Ve__df_e,
                tau_IRe__df_e,
                Te__df_e,
                be__df_e, 
                sw_down__df_e,
                sw_down_b__df_e,
                sw_up__df_e,
                lw_down__df_e,
                lw_down_b__df_e,
                lw_up__df_e,
                lw_up_b__df_e,
                lw_net__df_e,
                sw_net__df_e,
                // lw_grey_updown_linear working variables
                dtau__dff_l,
                del__dff_l, 
                edel__dff_l,
                e0i__dff_l,
                e1i__dff_l,
                Am__dff_l,
                Bm__dff_l,
                lw_up_g__dff_e,
                lw_down_g__dff_e);

            // ?????? Is insol_d needed ?
            insol_d[id] =  insol_d[id] * coszrs* (1-alb); 
        }
        else {
            insol_d[id] = 0;

            Kitzmann_TS_noscatt(id,
                nv,
                nvi,
                Altitude_d,
                Altitudeh_d,
                Rho_d,
                temperature_d,
                pressure_d,
                phtemp,
                k_V_3_nv_d,
                k_IR_2_nv_d,
                Beta_V_3_d,
                Beta_2_d,
                net_F_nvi_d,
                zenith_angles,
                insol_d[id],
                tint,
                gravit,
                AB_d[id],
                //Kitzman working variables
                tau_Ve__df_e,
                tau_IRe__df_e,
                Te__df_e,
                be__df_e, 
                sw_down__df_e,
                sw_down_b__df_e,
                sw_up__df_e,
                lw_down__df_e,
                lw_down_b__df_e,
                lw_up__df_e,
                lw_up_b__df_e,
                lw_net__df_e,
                sw_net__df_e,
                // lw_grey_updown_linear working variables
                dtau__dff_l,
                del__dff_l, 
                edel__dff_l,
                e0i__dff_l,
                e1i__dff_l,
                Am__dff_l,
                Bm__dff_l,
                lw_up_g__dff_e,
                lw_down_g__dff_e);
        }

        printf("Kitzmann really finished\n");

        /*
        
        if (isnan(zenith_angles[id] ) ) {
            
           
            printf("zenith_angles contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d --value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, &zenith_angles[id]);
                        //temperature_d[id * nv + level] = id * nv + level;
                        //__threadfence();         // ensure store issued before trap
                        //asm("trap;");            // kill kernel with error
            
        } 

        for (int channel = 0; channel < 2 ; channel++)
        {
            if (isnan(Beta_2_d[id*2 + channel] ) ) {
                for (int chan = 0; chan < 2; chan++)
                {
                    printf("Beta_2_d contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d channel:%d --value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, chan, &Beta_2_d[id*2 + chan]);
                        //temperature_d[id * nv + level] = id * nv + level;
                        __threadfence();         // ensure store issued before trap
                        asm("trap;");            // kill kernel with error
                }
            }

            if (Beta_2_d[id*2 + channel] ==0 ) {
                for (int chan = 0; chan < 2; chan++)
                {
                    printf("Beta_2_d contains zeroes in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d channel:%d --value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, chan, &Beta_2_d[id*2 + chan]);
                        //temperature_d[id * nv + level] = id * nv + level;
                        __threadfence();         // ensure store issued before trap
                        asm("trap;");            // kill kernel with error
                }
            }
            
        }
        printf("Beta_2_d print out\n");      
        
        for (int channel = 0; channel < 3; channel++)
            {
            
                if (gam_V_3_d[id * 3 + channel]  != gam_V_3_d[id * 3 + channel] ) {
                        printf("gam_V_3_d contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d  --value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, &gam_V_3_d[id*3 + channel]);
                        //temperature_d[id * nv + level] = id * nv + level;
                        //__threadfence();         // ensure store issued before trap
                        //asm("trap;");            // kill kernel with error
                }
            }
        
        printf("gam_V_3_d print out\n");      
            
        for (int level = 0; level < nv; level++)
        {
            if (isnan(temperature_d[id * nv + level]) ) {
                for (int lev = 0; lev < nv; lev++)
                {
                    printf("temperature_d contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d -- level: %d --value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, lev, &temperature_d[id*nv + lev]);
                   
                    __threadfence();         // ensure store issued before trap
                    asm("trap;");            // kill kernel with error
                }
            }
            if (isnan(temperature_d[id * nv + level]) ) {
                for (int lev = 0; lev < nv; lev++)
                {
                    printf("temperature_d contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d  level: %d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, lev, &temperature_d[id*nv + lev]);
                    
                    __threadfence();         // ensure store issued before trap
                    asm("trap;");            // kill kernel with error
                }
            }
            if (isnan(pressure_d[id * nv + level]) ) {
                for (int lev = 0; lev < nv; lev++)
                {
                    printf("pressure_d contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d -- level: %d --value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, lev, &pressure_d[id*nv + lev]);
                   
                    __threadfence();         // ensure store issued before trap
                    asm("trap;");            // kill kernel with error
                }
            }
            if (isnan(pressure_d[id * nv + level]) ) {
                for (int lev = 0; lev < nv; lev++)
                {
                    printf("pressure_d contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d  level: %d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, lev, &pressure_d[id*nv + lev]);
                    
                    __threadfence();         // ensure store issued before trap
                    asm("trap;");            // kill kernel with error
                }
            }

        }
        */

        //printf("temperature_d & pressure_d print out\n");  
        /*
        for (int level = 0; level < nvi; level++)
        {
            
            if (isnan(tau_IRe__df_e[id * nvi + level]) ) {
                for (int lev = 0; lev < nvi; lev++)
                {
                    printf("tau_IRe__df_e contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d -- level: %d --value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, lev, &tau_IRe__df_e[id*nvi + lev]);
                   
                    __threadfence();         // ensure store issued before trap
                    asm("trap;");            // kill kernel with error
                }
            }
            if (isnan(tau_Ve__df_e[id * nvi + level]) ) {
                for (int lev = 0; lev < nvi; lev++)
                {
                    printf("tau_Ve__df_e contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d  level: %d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, lev, &tau_Ve__df_e[id*nvi + lev]);
                    
                    __threadfence();         // ensure store issued before trap
                    asm("trap;");            // kill kernel with error
                }
            }
            

           
            if (phtemp[id * nvi + level]==0 ) {
                for (int lev = 0; lev < nvi; lev++)
                {
                    printf("phtemp contains zeroes in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d  level: %d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, lev, &phtemp[id*nvi + lev]);
                    
                    __threadfence();         // ensure store issued before trap
                    asm("trap;");            // kill kernel with error
                }
            }

             if (isnan(phtemp[id * nvi + level]) ) {
                for (int lev = 0; lev < nvi; lev++)
                {
                    printf("phtemp contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d  level: %d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, lev, &phtemp[id*nvi + lev]);
                    
                    __threadfence();         // ensure store issued before trap
                    asm("trap;");            // kill kernel with error
                }
            }

            
            

            
            if (isnan(Te__df_e[id * nvi + level]) ) {
                for (int lev = 0; lev < nvi; lev++)
                {
                    printf("Te__df_e contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d  level: %d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, lev, &Te__df_e[id*nvi + lev]);
                    
                    //__threadfence();         // ensure store issued before trap
                    //asm("trap;");            // kill kernel with error
                }
            }
            if (Te__df_e[id * nvi + level]==0 ) {
                for (int lev = 0; lev < nvi; lev++)
                {
                    printf("Te__df_e contains zeroes in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d  level: %d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, lev, &Te__df_e[id*nvi + lev]);
                    
                    //__threadfence();         // ensure store issued before trap
                    //asm("trap;");            // kill kernel with error
                }
            }
            

            if (isnan(be__df_e[id * nvi + level]) ) {
                for (int lev = 0; lev < nvi; lev++)
                {
                    printf("be__df_e contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d  level: %d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, lev, &be__df_e[id*nvi + lev]);
                    
                    //__threadfence();         // ensure store issued before trap
                    //asm("trap;");            // kill kernel with error
                }
            }  

             if (be__df_e[id * nvi + level]==0 ) {
                for (int lev = 0; lev < nvi; lev++)
                {
                    printf("be__df_e contains zeroes in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d  level: %d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, lev, &be__df_e[id*nvi + lev]);
                    
                    //__threadfence();         // ensure store issued before trap
                    //asm("trap;");            // kill kernel with error
                }
            }

            
        
        }
         */  
        
        for (int level = 0; level < nv; level++)
        {
            if ( isnan(temperature_d[id * nv + level] ) ) {
                    printf("temperature_d contains NaNs at level %d\n", level);
                    
                    //printf("temperature_d contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d --value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, &temperature_d[id*nv + level]);
                    //temperature_d[id * nv + level] = id * nv + level;
                    //__threadfence();         // ensure store issued before trap
                    //asm("trap;");            // kill kernel with error
            }
            
            for (int channel = 0; channel < 3; channel++)
            {
            
                if ( isnan(k_V_3_nv_d[id * nv * 3 + channel * nv + level] ) ) {

                    printf("k_V_3_nv_d contains NaNs at level: %d \n",level);
                    for (int lev = 0; lev < nv; lev++)
                    {
                        
                        
                        printf("k_V_3_nv_d contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d  level: %d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id,lev, &k_V_3_nv_d[id*nv*3 +channel*nv + lev]);
                        //temperature_d[id * nv + level] = id * nv + level;
                        //__threadfence();         // ensure store issued before trap
                        //asm("trap;");            // kill kernel with error
                    }
                }
            }

            for (int channel = 0; channel < 2; channel++)
            {
                if ( isnan(k_IR_2_nv_d[id * nv * 2 + channel * nv + level] ) ) {
                    printf("k_IR_2_nv_d contains NaNs at level:%d \n", level);
                    for (int lev = 0; lev < nv; lev++)
                    {
                        
                        //printf("k_IR_2_nv_d contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d level:%d --value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, lev, &k_IR_2_nv_d[id*nv*2 +channel*nv + lev]);

                    }
                    

                }
                
            }              
        }

        //printf("k_V_3_nv_d & k_IR_2_nv_d print out\n");
        
        
        
        
        
        
        for (int level = 0; level < nv; level++)
        {
            dtemp[id * nv + level] = 1* //(gravit / Cp_d) *
                (net_F_nvi_d[id * nvi + level+1] - net_F_nvi_d[id * nvi + level]) / 
                (Altitudeh_d[level+1] - Altitudeh_d[level]);

            if (dtemp[id * nv + level]<0.0)
            {
                //printf("dtemp contains a negative value at level:%d \n",  level);                    
            }

            if (dtemp[id * nv + level]==0.0)
            {
                //printf("dtemp contains a 0 at level:%d \n",  level);                    
            }

            if (isnan(dtemp[id * nv + level]))
            {
                printf("dtemp contains a 0 at level:%d \n",  level);                    
            }

        }

        printf("dtemp is computed\n");
        
        
        

        

        if (surface == true) {
            Tsurface_d[id] = timestep * dtemp[id * nv + 0];
            surf_flux_d[id] =  dtemp[id * nv + 0];
            
            if (Tsurface_d[id] < 0)
                Tsurface_d[id] = 0;
        }

        //printf("Tsurface_d is computed\n");

        //calculate ASR for this point
        double rscale;
        if (DeepModel) {
            rscale = (A + Altitudeh_d[nv]) / A;
        }
        else {
            rscale = 1.0;
        }

        //printf("rscale is computed\n");

        
        
        
        
        
        
        for (int level = 0; level < nvi; level++)
        {
                                    
            if (isnan(lw_up__df_e[id * nvi + level]) ) {
                for (int lev = 0; lev < nvi; lev++)
                {
                    printf("lw_up__df_e has NaNs at the level:%u\n", lev);
                    //printf("lw_up__df_e contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d timestep:%d level:%d \n", blockIdx.x, blockDim.x, threadIdx.x, id, timestep, lev);
                    //printf("lw_up__df_e has the value:%u\n",  &lw_up__df_e[id*nvi + lev]);
                }
                //sw_net__df_e[id * nv + level] = id * nv + level;
                //__threadfence();         // ensure store issued before trap
                //asm("trap;");            // kill kernel with error
            }
            if (lw_up__df_e[id * nvi + level]<0.0 ) {
                printf("lw_up__df_e has a negative value at the level:%u\n", level);
                //printf("lw_up__df_e has a negative value in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d timestep:%d level:%d \n", blockIdx.x, blockDim.x, threadIdx.x, id, timestep, level);
                    
                for (int lev = 0; lev < nvi; lev++)
                {
                    //printf("lw_up__df_e has 0 at the level:%u\n", lev);
                    //printf("lw_up__df_e contains 0 in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d timestep:%d level:%d \n", blockIdx.x, blockDim.x, threadIdx.x, id, timestep, lev);
                    //printf("lw_up__df_e has the value:%u\n",  &lw_up__df_e[id*nvi + lev]);
                }
                //sw_net__df_e[id * nv + level] = id * nv + level;
                //__threadfence();         // ensure store issued before trap
                //asm("trap;");            // kill kernel with error
            }

            if (lw_up__df_e[id * nvi + level] == 0.0 ) {
                printf("lw_up__df_e has 0 at the level:%u\n", level);
                //printf("lw_up__df_e contains 0 in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d timestep:%d level:%d \n", blockIdx.x, blockDim.x, threadIdx.x, id, timestep, level);
                    
                for (int lev = 0; lev < nvi; lev++)
                {
                    //printf("lw_up__df_e has 0 at the level:%u\n", lev);
                    //printf("lw_up__df_e contains 0 in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d timestep:%d level:%d \n", blockIdx.x, blockDim.x, threadIdx.x, id, timestep, lev);
                    //printf("lw_up__df_e has the value:%u\n",  &lw_up__df_e[id*nvi + lev]);
                }
                //sw_net__df_e[id * nv + level] = id * nv + level;
                //__threadfence();         // ensure store issued before trap
                //asm("trap;");            // kill kernel with error
            }

            if (lw_down__df_e[id * nvi + level] == 0.0  ) {
                printf("lw_down__df_e has 0 at the level:%u\n", level);
                //printf("lw_down__df_e contains 0 in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d timestep:%d level:%d \n", blockIdx.x, blockDim.x, threadIdx.x, id, timestep, level);
                   
                for (int lev = 0; lev < nvi; lev++)
                {
                    //printf("lw_down__df_e has 0 at the level:%u\n", lev);
                    //printf("lw_down__df_e contains 0 in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d timestep:%d level:%d \n", blockIdx.x, blockDim.x, threadIdx.x, id, timestep, lev);
                    //printf("lw_down__df_e has the value:%u\n",  &lw_down__df_e[id*nvi + level]);
                }
                //sw_net__df_e[id * nv + level] = id * nv + level;
                //__threadfence();         // ensure store issued before trap
                //asm("trap;");            // kill kernel with error
            }

            if (lw_down__df_e[id * nvi + level] < 0.0  ) {
                //printf("lw_down__df_e has a negative value at the level:%u\n", level);
                //printf("lw_down__df_e has a negative value in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d timestep:%d level:%d \n", blockIdx.x, blockDim.x, threadIdx.x, id, timestep, level);
                   
                for (int lev = 0; lev < nvi; lev++)
                {
                    //printf("lw_down__df_e has 0 at the level:%u\n", lev);
                    //printf("lw_down__df_e contains 0 in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d timestep:%d level:%d \n", blockIdx.x, blockDim.x, threadIdx.x, id, timestep, lev);
                    //printf("lw_down__df_e has the value:%u\n",  &lw_down__df_e[id*nvi + level]);
                }
                //sw_net__df_e[id * nv + level] = id * nv + level;
                //__threadfence();         // ensure store issued before trap
                //asm("trap;");            // kill kernel with error
            }
            if (isnan(lw_down__df_e[id * nvi + level] ) ) {
                printf("lw_down__df_e has a NaN at the level:%u\n", level);
                for (int lev = 0; lev < nvi; lev++)
                {
                    //printf("lw_down__df_e has NaNs at the level:%u\n", lev);
                    //printf("lw_down__df_e contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d timestep:%d level:%d \n", blockIdx.x, blockDim.x, threadIdx.x, id, timestep, lev);
                    //printf("lw_down__df_e has the value:%u\n",  &lw_down__df_e[id*nvi + level]);
                }
                //sw_net__df_e[id * nv + level] = id * nv + level;
                //__threadfence();         // ensure store issued before trap
                //asm("trap;");            // kill kernel with error
            }
            
            if (isnan(sw_up__df_e[id * nvi + level] )  ) {
                printf("sw_up__df_e has a NaN at the level:%u\n", level);
                for (int lev = 0; lev < nvi; lev++)
                {
                    //printf("sw_up__df_e has NaNs at the level:%u\n", lev);
                    //printf("sw_up__df_e contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d timestep:%d level: %d \n", blockIdx.x, blockDim.x, threadIdx.x, id, timestep, lev);
                    //printf("sw_up__df_e has the value:%u\n", &sw_up__df_e[id*nvi+lev]);
                }
               
               // __threadfence();         // ensure store issued before trap
                //asm("trap;");            // kill kernel with error
            }
            if (isnan(sw_down__df_e[id * nvi + level] )  ) {
                printf("sw_down__df_e has a NaN at the level:%u\n", level);
                  
                for (int lev = 0; lev < nvi; lev++)
                {
                    //printf("sw_down__df_e has NaNs at the level:%u\n", lev);
                    //printf("sw_down__df_e contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d timestep:%d level:%d \n", blockIdx.x, blockDim.x, threadIdx.x, id, timestep, lev);
                    //printf("sw_down__df_e has the value:%u\n", &sw_down__df_e[id*nvi + lev]);
                }
                //sw_net__df_e[id * nv + level] = id * nv + level;
                //__threadfence();         // ensure store issued before trap
                //asm("trap;");            // kill kernel with error
            }

            if (sw_down__df_e[id * nvi + level] < 0.0  ) {

                //printf("sw_down__df_e has a negative value at the level:%u\n", level);
                  
                for (int lev = 0; lev < nvi; lev++)
                {
                    //printf("sw_down__df_e has NaNs at the level:%u\n", lev);
                    //printf("sw_down__df_e has a negative value in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d timestep:%d level:%d \n", blockIdx.x, blockDim.x, threadIdx.x, id, timestep, lev);
                    //printf("sw_down__df_e has the value:%u\n", &sw_down__df_e[id*nvi + lev]);
                }
                //sw_net__df_e[id * nv + level] = id * nv + level;
                //__threadfence();         // ensure store issued before trap
                //asm("trap;");            // kill kernel with error
            }

            if (sw_down__df_e[id * nvi + level] == 0.0  ) {

                //printf("sw_down__df_e has a zero at the level:%u\n", level);                  
                for (int lev = 0; lev < nvi; lev++)
                {
                    //printf("sw_down__df_e has NaNs at the level:%u\n", lev);
                    //printf("sw_down__df_e contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d timestep:%d level:%d \n", blockIdx.x, blockDim.x, threadIdx.x, id, timestep, lev);
                    //printf("sw_down__df_e has the value:%u\n", &sw_down__df_e[id*nvi + lev]);
                }
                //sw_net__df_e[id * nv + level] = id * nv + level;
                //__threadfence();         // ensure store issued before trap
                //asm("trap;");            // kill kernel with error
            }
            
            
            if (isnan(sw_net__df_e[id * nvi + level])  ) {
                printf("sw_net__df_e contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d level:%d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, level, &sw_net__df_e[id*nvi + level]);
                //sw_net__df_e[id * nv + level] = id * nv + level;
                //__threadfence();         // ensure store issued before trap
                //asm("trap;");            // kill kernel with error
            }


            //calculate OLR for this point
            
             
            if (isnan(lw_net__df_e[id * nvi + level])) {
                printf("lw_net__df_e contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d  level:%d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, level, &lw_net__df_e[id*nvi + level]);
                  
                //__threadfence();          // ensure store issued before trap
                asm("trap;");            // kill kernel with error            
            }
            

            
            if ( lw_net__df_e[id * nvi + level] ==0.0  ) {
                    printf("lw_net__df_e contains 0 in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d  level:%d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, level, &lw_net__df_e[id*nvi + level]);
                    //lw_net__df_e[id] = id;
                    //__threadfence();         // ensure store issued before trap
                    //asm("trap;");            // kill kernel with error
            }            

        }

        if (zenith_angles[id]>0.0)
        {
            printf("zenith_angles is larger than 0\n");
        } else
        {
            printf("zenith_angles is zero or smaller\n");
        }
        if (zenith_angles[id]==0.0)
        {
            printf("zenith_angles is 0\n");
        }

        if (zenith_angles[id]>1)
        {
             printf("zenith_angles is larger than 1\n");

        }
        
        

        printf("before ASR_d\n");
        
        ASR_d[id] = sw_down__df_e[id * nvi + nv]; // * areasT_d[id] * pow(rscale, 2);
        printf("ASR_d is computed\n");
            
            
            if (isnan(ASR_d[id] )) {
                printf("ASR_d contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, &ASR_d[id]);
                //ASR_d[id] = id;
                __threadfence();         // ensure store issued before trap
                asm("trap;");            // kill kernel with error
            }

            if (ASR_d[id] == 0.0 ) {
                //printf("my message:ASR_d contains 0 in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, &ASR_d[id]);
                //OLR_d[id] = id;
                //__threadfence();         // ensure store issued before trap
                //asm("trap;");            // kill kernel with error
            }
             if (ASR_d[id] < 0.0 ) {
                printf("my message:ASR_d contains a negative value in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, &ASR_d[id]);
                //OLR_d[id] = id;
                //__threadfence();         // ensure store issued before trap
                //asm("trap;");            // kill kernel with error
            }

            

            
        
        OLR_d[id] = lw_up__df_e[id * nvi + nv] ; //*areasT_d[id] * pow(rscale, 2);
        printf("OLR_d is computed\n");

        if (isnan(OLR_d[id] )) {
                printf("OLR_d contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, &OLR_d[id]);
                //ASR_d[id] = id;
                __threadfence();         // ensure store issued before trap
                asm("trap;");            // kill kernel with error
            }

        if (OLR_d[id] < 0.0 ) {
                printf("my message:ASR_d contains a negative value in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, &OLR_d[id]);
                //OLR_d[id] = id;
                //__threadfence();         // ensure store issued before trap
                //asm("trap;");            // kill kernel with error
            }

        if (isnan(areasT_d[id] ) ) {
                printf("areasT_d[id] contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, &areasT_d[id]);
                //OLR_d[id] = id;
                //__threadfence();         // ensure store issued before trap
                //asm("trap;");            // kill kernel with error
            }

        if (areasT_d[id] <1.0 ) {
                printf("areasT_d[id] small than 1 in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, &areasT_d[id]);
                //OLR_d[id] = id;
                //__threadfence();         // ensure store issued before trap
                //asm("trap;");            // kill kernel with error
            }
        if (areasT_d[id] <0.01 ) {
                printf("areasT_d[id] small than 0.01 in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, &areasT_d[id]);
                //OLR_d[id] = id;
                //__threadfence();         // ensure store issued before trap
                //asm("trap;");            // kill kernel with error
            }
        
        if (areasT_d[id] <0.0001 ) {
                printf("areasT_d[id] small than 0.0001 in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, &areasT_d[id]);
                //OLR_d[id] = id;
                //__threadfence();         // ensure store issued before trap
                //asm("trap;");            // kill kernel with error
            }
        if (areasT_d[id] <0.000001 ) {
                printf("areasT_d[id] small than 0.000001 in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, &areasT_d[id]);
                //OLR_d[id] = id;
                //__threadfence();         // ensure store issued before trap
                //asm("trap;");            // kill kernel with error
            }

            if (areasT_d[id] == 0.0 ) {
                printf("areasT_d[id] contains 0 in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, &areasT_d[id]);
                //OLR_d[id] = id;
                //__threadfence();         // ensure store issued before trap
                //asm("trap;");            // kill kernel with error
            }

            if (rscale == 0.0 ) {
                printf("rscale contains 0 in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d value:%d\n", blockIdx.x, blockDim.x, threadIdx.x, id, rscale);
                //OLR_d[id] = id;
                //__threadfence();         // ensure store issued before trap
                //asm("trap;");            // kill kernel with error
            }
             if (isnan(rscale) ) {
                printf("rscale contains NaN in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d value:%d\n", blockIdx.x, blockDim.x, threadIdx.x, id, rscale);
                //OLR_d[id] = id;
                //__threadfence();         // ensure store issued before trap
                //asm("trap;");            // kill kernel with error
            }
            
            
            if (isnan(OLR_d[id] ) ) {
                printf("my message:OLR_d contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, &OLR_d[id]);
                //OLR_d[id] = id;
                //__threadfence();         // ensure store issued before trap
                //asm("trap;");            // kill kernel with error
            }
            if (OLR_d[id] == 0 ) {
                printf("my message:OLR_d contains 0 in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, &OLR_d[id]);
                //OLR_d[id] = id;
                //__threadfence();         // ensure store issued before trap
                //asm("trap;");            // kill kernel with error
            }
            
            
        
        
        
        for (int lev = 0; lev < nv; lev++) {
        
            /*
        
            if ( isnan(dtemp[id * nv + lev]) ) {
                for (int lev = 0; lev < nv; lev++) {
                
                    printf("dtemp contains NaNs in blockIdx.x:%d * blockDim.x:%d + threadIdx.x:%d = globalThreadId:%d timestep:%d level:%d value:%u\n", blockIdx.x, blockDim.x, threadIdx.x, id, timestep, lev, &dtemp[id * nv + lev]);
                }
                //dtemp[id] = id;
                //__threadfence();         // ensure store issued before trap
                //asm("trap;");            // kill kernel with error
            }
            
            
            
            if (isnan(Rd_d[id] )  ) {
                //__threadfence();         // ensure store issued before trap
                //asm("trap;");            // kill kernel with error
            }
            if (isnan(Cp_d[id] ) ) {
                //__threadfence();         // ensure store issued before trap
                //asm("trap;");            // kill kernel with error
            }
            
            */
            
            
            
            
            
            if (pressure_d[id * nv + lev]
                    + Rd_d[id * nv + lev] / (Cp_d[id * nv + lev] - Rd_d[id * nv + lev])
                          * dtemp[id * nv + lev] * timestep
                < 0) {
                //trying to prevent too much cooling resulting in negative pressure in dyn core
                dtemp[id * nv + lev] = -pressure_d[id * nv + lev] / timestep;
            }
            DG_Qheat_d[id * nv + lev] = dtemp[id * nv + lev];
            profx_Qheat_d[id * nv + lev] += Qheat_scaling * dtemp[id * nv + lev];
            if (isnan(profx_Qheat_d[id * nv + lev])) {
                printf("profx_Qheat_d has NaNs - stop here");
            }
            // }
        }

        printf("Column complete");


    }
}
