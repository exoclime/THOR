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
        if (DeepModel) {
            rup =
                (Altitudeh_d[lev + 1] + A) / (Altitude_d[lev] + A); //vertical scaling in divergence
            rlow = (Altitudeh_d[lev] + A) / (Altitude_d[lev] + A);
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
        if (DeepModel) {
            rup =
                (Altitudeh_d[lev + 1] + A) / (Altitude_d[lev] + A); //vertical scaling in divergence
            rlow = (Altitudeh_d[lev] + A) / (Altitude_d[lev] + A);
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
                              bool    latf_lw,
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

        // Compute opacities
        double kappa_lw_lat;
        if (latf_lw) {
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
            if (gcm_off) {
                temperature_d[id * nv + lev] = ttemp[id * nv + lev]
                                               + 1.0 / (Cp_d[id * nv + lev] - Rd_d[id * nv + lev])
                                                     * dtemp[id * nv + lev] / Rho_d[id * nv + lev]
                                                     * timestep;
                if (temperature_d[id * nv + lev] < 0)
                    temperature_d[id * nv + lev] = 0;
            }
            else {
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
            }
        }
    }
}
