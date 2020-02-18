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
// Description: Solves the vertical momentum (HE-VI scheme).
//
//
// Method:
//   Uses an implicit scheme.
//   Tridiagonal matrix solved numerically using a Thomas algorithm.
//
// Known limitations: The method to solve the tridiagonal matrix can be further optimized.
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
__global__ void Vertical_Eq(double *Whs_d,
                            double *Ws_d,
                            double *pressures_d,
                            double *h_d,
                            double *hh_d,
                            double *Rhos_d,
                            double *gtil_d,
                            double *gtilh_d,
                            double *Sp_d,
                            double *Sd_d,
                            double *Srh_d,
                            double *Cp_d,
                            double *Rd_d,
                            double  deltat,
                            double  Gravit,
                            double *Altitude_d,
                            double *Altitudeh_d,
                            double  A,
                            bool    NonHydro,
                            int     num,
                            int     nv,
                            int     nvi,
                            bool    DeepModel) {

    //
    //  Integration-> y'' = c3 y' + c2 y + c1
    //
    int id = blockIdx.x * blockDim.x + threadIdx.x;


    // using arrays created here either goes into local memory, which is slow and off device,
    // or uses a lot of registers and limits the number of blocks that can be run concurrently.
    // Replace it with shared memory which is on the device
    //double *cc, *dd;
    //cc = new double[nvi];
    //dd = new double[nvi];

    // define shared memory as external, because size is not known here
    extern __shared__ double mem_shared[];


    double *cc = (double *)mem_shared;
    double *dd = (double *)&mem_shared[blockDim.x * nvi];

    // double Cv = Cp - Rd;
    double C0;
    double xi, xim, xip;
    double intt, intl, inttm, intlm;
    double dSpdz, dPdz;
    double rhohs;
    double aa, bb;
    double t;
    double althl, alth, altht;
    double alt, altl;
    double Sp, Spl, p, pl, r, rl;
    double Sdl, Sd, Sdh;
    double hp, h, hm;
    double gp, g, gm;
    double whl, wht;
    double GCoRl, GCoRu, GCoR;
    double CRddl, CRddu, CRdd;
    double or2;
    double tor3;

    if (id < num) {
        for (int lev = 1; lev < nv; lev++) {
            if (lev == 1) { // Fetch data for first layer
                altht = Altitudeh_d[lev + 1];
                althl = Altitudeh_d[lev - 1];
                alth  = Altitudeh_d[lev];
                alt   = Altitude_d[lev];
                altl  = Altitude_d[lev - 1];
                hp    = hh_d[id * nvi + lev + 1];
                h     = hh_d[id * nvi + lev];
                hm    = hh_d[id * nvi + lev - 1];
                gp    = gtilh_d[id * nvi + lev + 1];
                g     = gtilh_d[id * nvi + lev];
                gm    = gtilh_d[id * nvi + lev - 1];
                Sp    = Sp_d[id * nv + lev];
                Spl   = Sp_d[id * nv + lev - 1];
                p     = pressures_d[id * nv + lev];
                pl    = pressures_d[id * nv + lev - 1];
                rl    = Rhos_d[id * nv + lev - 1];
                r     = Rhos_d[id * nv + lev];
                Sdl   = Sd_d[id * nv + lev - 1];
                Sd    = Sd_d[id * nv + lev];
                CRddl = (Cp_d[id * nv + lev - 1] - Rd_d[id * nv + lev - 1])
                        / (Rd_d[id * nv + lev - 1] * deltat * deltat);
                CRddu = (Cp_d[id * nv + lev] - Rd_d[id * nv + lev])
                        / (Rd_d[id * nv + lev] * deltat * deltat);
                GCoRl = Gravit * (Cp_d[id * nv + lev - 1] - Rd_d[id * nv + lev - 1])
                        / Rd_d[id * nv + lev - 1];
                GCoRu = Gravit * (Cp_d[id * nv + lev] - Rd_d[id * nv + lev]) / Rd_d[id * nv + lev];
            }
            if (DeepModel) {
                double dzp  = 1.0 / (altht - alth);
                double dzh  = 1.0 / (alt - altl);
                double dzm  = 1.0 / (alth - althl);
                double dzph = dzp * dzh;
                double dzmh = dzm * dzh;

                hp   = hp * pow(altht + A, 2);
                h    = h * pow(alth + A, 2);
                hm   = hm * pow(althl + A, 2);
                or2  = 1 / pow(alth + A, 2);
                tor3 = 2 / pow(alth + A, 3);

                xi  = alt;
                xim = alth;
                xip = altht;

                intt = -(xi - xip) * dzp * dzh;
                intl = (xi - xim) * dzp * dzh;

                xi  = altl;
                xim = althl;
                xip = alth;

                inttm = -(xi - xip) * dzm * dzh;
                intlm = (xi - xim) * dzm * dzh;

                // get g*Cv/Rd and Cv/Rd/dt^2 at the current interface
                CRdd = (CRddl * (alt - alth) + CRddu * (alth - altl)) / (alt - altl);
                GCoR = (GCoRl * (alt - alth) + GCoRu * (alth - altl)) / (alt - altl);

                cc[threadIdx.x * nvi + lev] = -dzph * or2 * hp - intt * (gp + GCoR - tor3 * hp);

                if (NonHydro)
                    bb = CRdd + (dzph + dzmh) * or2 * h + (intl - inttm) * (g + GCoR - tor3 * h);
                else
                    bb = (dzph + dzmh) * or2 * h + (intl - inttm) * (g + GCoR - tor3 * h);

                aa = -dzmh * or2 * hm + intlm * (gm + GCoR - tor3 * hm);

                dSpdz = (Sp - Spl) * dzh;
                dPdz  = (p - pl) * dzh;

                xi  = alth;
                xim = altl;
                xip = alt;

                intt = (xi - xip) / (xim - xip);
                intl = (xi - xim) / (xip - xim);

                rhohs = rl * intt + r * intl;
                Sdh   = Sdl * intt + Sd * intl;

                if (!NonHydro) {
                    C0 =
                        (pow(deltat, 2.0) * dSpdz + pow(deltat, 2.0) * Gravit * Sdh
                         + Gravit * deltat * rhohs + deltat * dPdz - deltat * Srh_d[id * nvi + lev])
                        * (CRdd);
                }
                else {
                    C0 = (-Whs_d[id * nvi + lev] + pow(deltat, 2.0) * dSpdz
                          + pow(deltat, 2.0) * Gravit * Sdh + Gravit * deltat * rhohs
                          + deltat * dPdz - deltat * Srh_d[id * nvi + lev])
                         * (CRdd);
                }

                if (lev < nv - 1) { // Fetch the data for next layer
                    althl = alth;
                    alth  = altht;
                    altht = Altitudeh_d[lev + 2];
                    altl  = alt;
                    alt   = Altitude_d[lev + 1];
                    hm    = hh_d[id * nvi + lev];
                    h     = hh_d[id * nvi + lev + 1];
                    hp    = hh_d[id * nvi + lev + 2];
                    gm    = g;
                    g     = gp;
                    gp    = gtilh_d[id * nvi + lev + 2];
                    Spl   = Sp;
                    Sp    = Sp_d[id * nv + lev + 1];
                    pl    = p;
                    p     = pressures_d[id * nv + lev + 1];
                    rl    = r;
                    r     = Rhos_d[id * nv + lev + 1];
                    Sdl   = Sd;
                    Sd    = Sd_d[id * nv + lev + 1];
                    CRddl = CRddu;
                    CRddu = (Cp_d[id * nv + lev + 1] - Rd_d[id * nv + lev + 1])
                            / (Rd_d[id * nv + lev + 1] * deltat * deltat);
                    GCoRl = GCoRu;
                    GCoRu = Gravit * (Cp_d[id * nv + lev + 1] - Rd_d[id * nv + lev + 1])
                            / Rd_d[id * nv + lev + 1];
                }
            }
            else { // DeepModel == false
                double dzp  = 1.0 / (altht - alth);
                double dzh  = 1.0 / (alt - altl);
                double dzm  = 1.0 / (alth - althl);
                double dzph = dzp * dzh;
                double dzmh = dzm * dzh;

                xi  = alt;
                xim = alth;
                xip = altht;

                intt = -(xi - xip) * dzp * dzh;
                intl = (xi - xim) * dzp * dzh;

                xi  = altl;
                xim = althl;
                xip = alth;

                inttm = -(xi - xip) * dzm * dzh;
                intlm = (xi - xim) * dzm * dzh;

                // get g*Cv/Rd and Cv/Rd/dt^2 at the current interface
                CRdd = (CRddl * (alt - alth) + CRddu * (alth - altl)) / (alt - altl);
                GCoR = (GCoRl * (alt - alth) + GCoRu * (alth - altl)) / (alt - altl);

                cc[threadIdx.x * nvi + lev] = -dzph * hp - intt * (gp + GCoR);

                if (NonHydro)
                    bb = CRdd + (dzph + dzmh) * h + (intl - inttm) * (g + GCoR);
                else
                    bb = (dzph + dzmh) * h + (intl - inttm) * (g + GCoR);

                aa = -dzmh * hm + intlm * (gm + GCoR);

                dSpdz = (Sp - Spl) * dzh;
                dPdz  = (p - pl) * dzh;

                xi  = alth;
                xim = altl;
                xip = alt;

                intt = (xi - xip) / (xim - xip);
                intl = (xi - xim) / (xip - xim);

                rhohs = rl * intt + r * intl;
                Sdh   = Sdl * intt + Sd * intl;

                if (!NonHydro) {
                    C0 =
                        (pow(deltat, 2.0) * dSpdz + pow(deltat, 2.0) * Gravit * Sdh
                         + Gravit * deltat * rhohs + deltat * dPdz - deltat * Srh_d[id * nvi + lev])
                        * (CRdd);
                }
                else {
                    C0 = (-Whs_d[id * nvi + lev] + pow(deltat, 2.0) * dSpdz
                          + pow(deltat, 2.0) * Gravit * Sdh + Gravit * deltat * rhohs
                          + deltat * dPdz - deltat * Srh_d[id * nvi + lev])
                         * (CRdd);
                }
                if (lev < nv - 1) { // Fetch data for the next layer
                    althl = alth;
                    alth  = altht;
                    altht = Altitudeh_d[lev + 2];
                    altl  = alt;
                    alt   = Altitude_d[lev + 1];
                    hm    = h;
                    h     = hp;
                    hp    = hh_d[id * nvi + lev + 2];
                    gm    = g;
                    g     = gp;
                    gp    = gtilh_d[id * nvi + lev + 2];
                    Spl   = Sp;
                    Sp    = Sp_d[id * nv + lev + 1];
                    pl    = p;
                    p     = pressures_d[id * nv + lev + 1];
                    rl    = r;
                    r     = Rhos_d[id * nv + lev + 1];
                    Sdl   = Sd;
                    Sd    = Sd_d[id * nv + lev + 1];
                    CRddl = CRddu;
                    CRddu = (Cp_d[id * nv + lev + 1] - Rd_d[id * nv + lev + 1])
                            / (Rd_d[id * nv + lev + 1] * deltat * deltat);
                    GCoRl = GCoRu;
                    GCoRu = Gravit * (Cp_d[id * nv + lev + 1] - Rd_d[id * nv + lev + 1])
                            / Rd_d[id * nv + lev + 1];
                }
            } // End of if (DeepModel) physics computation


            dd[threadIdx.x * nvi + lev] = -C0;
            if (lev == 1) {
                cc[threadIdx.x * nvi + 1] = cc[threadIdx.x * nvi + 1] / bb;
                dd[threadIdx.x * nvi + 1] = dd[threadIdx.x * nvi + 1] / bb;
            }
            else {
                t = 1.0 / (bb - cc[threadIdx.x * nvi + lev - 1] * aa);
                cc[threadIdx.x * nvi + lev] *= t;
                dd[threadIdx.x * nvi + lev] =
                    (dd[threadIdx.x * nvi + lev] - dd[threadIdx.x * nvi + lev - 1] * aa) * t;
            }
        } //end of loop over levels
        Whs_d[id * nvi + nv]     = 0.0;
        Whs_d[id * nvi]          = 0.0;
        Whs_d[id * nvi + nv - 1] = dd[threadIdx.x * nvi + nv - 1];
        // Whs_d[id * nvi + nv - 1] = 0.0;
        // Updates vertical momentum
        for (int lev = nvi - 2; lev > 0; lev--)
            Whs_d[id * nvi + lev] = (-cc[threadIdx.x * nvi + lev] * Whs_d[id * nvi + lev + 1]
                                     + dd[threadIdx.x * nvi + lev]);
        // Whs_d[id * nvi + lev] = 0.0;

        for (int lev = 0; lev < nv; lev++) {

            if (lev == 0) {
                althl = 0;
                altht = Altitudeh_d[1];
                whl   = 0.0;
                wht   = Whs_d[id * (nv + 1) + 1];
            }

            xi  = Altitude_d[lev];
            xim = althl;
            xip = altht;

            intt = (xi - xip) / (xim - xip);
            intl = (xi - xim) / (xip - xim);

            Ws_d[id * nv + lev] = whl * intt + wht * intl;
            if (lev < nv - 1) {
                althl = altht;
                altht = Altitudeh_d[lev + 2];
                whl   = wht;
                wht   = Whs_d[id * (nv + 1) + lev + 2];
            }
        }
    }
}


template<int NX, int NY>
__global__ void Prepare_Implicit_Vertical(double *Mh_d,
                                          double *h_d,
                                          double *div_d,
                                          double *Slowpressure_d,
                                          double *SlowRho_d,
                                          double *Sp_d,
                                          double *Sd_d,
                                          double *Altitude_d,
                                          double *Cp_d,
                                          double *Rd_d,
                                          double  A,
                                          int *   maps_d,
                                          int     nl_region,
                                          bool    DeepModel) {

    int x = threadIdx.x;
    int y = threadIdx.y;
    //int ib  = blockIdx.x;
    int nv  = gridDim.y;
    int lev = blockIdx.y;

    int    pt1, pt2, pt3, pt4, pt5, pt6;
    double div0, div1, div2, div3, div4, div5, div6;
    int    nhl  = nl_region + 2;
    int    nhl2 = nhl * nhl;


    __shared__ double nflxr_s[NX * NY];
    __shared__ double nflxp_s[NX * NY];
    __shared__ double v_s[3 * (NX + 2) * (NY + 2)];
    __shared__ double h_s[(NX + 2) * (NY + 2)];

    double alt, rscale;

    int ir = 0; // index in region
    int iri, ir2, id;


    bool pent_ind = false; //
    int  ig;               // index in global mem

    int igh = 0; // global index in halo

    // Load shared memory


    bool load_halo = compute_mem_idx(maps_d, nhl, nhl2, ig, igh, ir, ir2, pent_ind);
    id             = ig;
    double RoC     = Rd_d[id * nv + lev] / (Cp_d[id * nv + lev] - Rd_d[id * nv + lev]);


    v_s[ir * 3 + 0] = Mh_d[ig * 3 * nv + lev * 3 + 0];
    v_s[ir * 3 + 1] = Mh_d[ig * 3 * nv + lev * 3 + 1];
    v_s[ir * 3 + 2] = Mh_d[ig * 3 * nv + lev * 3 + 2];

    h_s[ir] = h_d[ig * nv + lev];

    ///////////////////////////////
    //////////// Halo /////////////
    ///////////////////////////////
    if (load_halo) {
        if (igh >= 0) {
            v_s[ir2 * 3 + 0] = Mh_d[igh * 3 * nv + lev * 3 + 0];
            v_s[ir2 * 3 + 1] = Mh_d[igh * 3 * nv + lev * 3 + 1];
            v_s[ir2 * 3 + 2] = Mh_d[igh * 3 * nv + lev * 3 + 2];
            h_s[ir2]         = h_d[igh * nv + lev];
        }
        else {
            v_s[ir2 * 3 + 0] = 0.0;
            v_s[ir2 * 3 + 1] = 0.0;
            v_s[ir2 * 3 + 2] = 0.0;
            h_s[ir2]         = 0.0;
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

    if (DeepModel) {
        alt    = Altitude_d[lev];
        rscale = A / (alt + A);
    }
    else
        rscale = 1.0;

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

        nflxr_s[iri] -=
            rscale
            * (div0 * v_s[ir * 3 + k] + div1 * v_s[pt1 * 3 + k] + div2 * v_s[pt2 * 3 + k]
               + div3 * v_s[pt3 * 3 + k] + div4 * v_s[pt4 * 3 + k] + div5 * v_s[pt5 * 3 + k]
               + div6 * v_s[pt6 * 3 + k]);

        nflxp_s[iri] -= (RoC)*rscale
                        * (div0 * v_s[ir * 3 + k] * h_s[ir] + div1 * v_s[pt1 * 3 + k] * h_s[pt1]
                           + div2 * v_s[pt2 * 3 + k] * h_s[pt2] + div3 * v_s[pt3 * 3 + k] * h_s[pt3]
                           + div4 * v_s[pt4 * 3 + k] * h_s[pt4] + div5 * v_s[pt5 * 3 + k] * h_s[pt5]
                           + div6 * v_s[pt6 * 3 + k] * h_s[pt6]);
    }

    Sp_d[id * nv + lev] = nflxp_s[iri] + Slowpressure_d[id * nv + lev];
    Sd_d[id * nv + lev] = nflxr_s[iri] + SlowRho_d[id * nv + lev];
}

template<int NN>
__global__ void Prepare_Implicit_Vertical_Poles(double *Mh_d,
                                                double *h_d,
                                                double *div_d,
                                                double *Slowpressure_d,
                                                double *SlowRho_d,
                                                double *Sp_d,
                                                double *Sd_d,
                                                double *Altitude_d,
                                                double *Cp_d,
                                                double *Rd_d,
                                                double  A,
                                                int *   point_local_d,
                                                int     num,
                                                int     nv,
                                                bool    DeepModel) {

    int id = blockIdx.x * blockDim.x + threadIdx.x;
    id += num - 2; // Poles

    __shared__ double div_p[3 * 7];
    __shared__ double v_p[3 * NN];
    __shared__ double h_p[NN];
    __shared__ int    local_p[NN];

    double nflxr_p;
    double nflxp_p;

    double alt, rscale;
    double Cv;
    if (id < num) {
        for (int i = 0; i < 5; i++)
            local_p[i] = point_local_d[id * 6 + i];
        for (int i = 0; i < 7; i++)
            for (int k = 0; k < 3; k++)
                div_p[i * 3 + k] = div_d[id * 7 * 3 + i * 3 + k];

        for (int lev = 0; lev < nv; lev++) {
            Cv     = Cp_d[id * nv + lev] - Rd_d[id * nv + lev];
            v_p[0] = Mh_d[id * 3 * nv + lev * 3 + 0];
            v_p[1] = Mh_d[id * 3 * nv + lev * 3 + 1];
            v_p[2] = Mh_d[id * 3 * nv + lev * 3 + 2];

            h_p[0] = h_d[id * nv + lev];
            for (int i = 1; i < 6; i++) {
                v_p[i * 3 + 0] = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 0];
                v_p[i * 3 + 1] = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 1];
                v_p[i * 3 + 2] = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 2];

                h_p[i] = h_d[local_p[i - 1] * nv + lev];
            }

            if (DeepModel) {
                alt    = Altitude_d[lev];
                rscale = A / (alt + A);
            }
            else
                rscale = 1.0;

            nflxr_p = 0.0;
            nflxp_p = 0.0;

            for (int k = 0; k < 3; k++) {
                nflxr_p -=
                    rscale
                    * (div_p[3 * 0 + k] * v_p[0 * 3 + k] + div_p[3 * 1 + k] * v_p[1 * 3 + k]
                       + div_p[3 * 2 + k] * v_p[2 * 3 + k] + div_p[3 * 3 + k] * v_p[3 * 3 + k]
                       + div_p[3 * 4 + k] * v_p[4 * 3 + k] + div_p[3 * 5 + k] * v_p[5 * 3 + k]);

                nflxp_p -= (Rd_d[id * nv + lev] / Cv) * rscale
                           * (div_p[3 * 0 + k] * v_p[0 * 3 + k] * h_p[0]
                              + div_p[3 * 1 + k] * v_p[1 * 3 + k] * h_p[1]
                              + div_p[3 * 2 + k] * v_p[2 * 3 + k] * h_p[2]
                              + div_p[3 * 3 + k] * v_p[3 * 3 + k] * h_p[3]
                              + div_p[3 * 4 + k] * v_p[4 * 3 + k] * h_p[4]
                              + div_p[3 * 5 + k] * v_p[5 * 3 + k] * h_p[5]);
            }
            Sp_d[id * nv + lev] = nflxp_p + Slowpressure_d[id * nv + lev];
            Sd_d[id * nv + lev] = nflxr_p + SlowRho_d[id * nv + lev];
        }
    }
}
