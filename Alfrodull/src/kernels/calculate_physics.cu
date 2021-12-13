// ==============================================================================
// This file is part of Alfrodull.
//
//     Alfrodull is free software : you can redistribute it and / or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
//
//     Alfrodull is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
//     GNU General Public License for more details.
//
//     You find a copy of the GNU General Public License in the main
//     Alfrodull directory under <license.txt>.If not, see
//     <http://www.gnu.org/licenses/>.
// ==============================================================================
//
// Kernels computing the physical quantities for each layers, used in the
// flux propagation matrix.
//
//
// Method: Helios Two Stream algorithm
//
//
// Known limitations: - Runs in a single GPU.
//
// Known issues: None
//
//
// Code contributors: Urs Schroffenegger, Matej Malik
//
// History:
// Version Date       Comment
// ======= ====       =======
// 1.0     2020-07-15 First version
//
//
////////////////////////////////////////////////////////////////////////

#include "calculate_physics.h"
#include "debug.h"
#include "physics_constants.h"
#include <stdio.h>

// fitting function for the E parameter according to "Heng, Malik & Kitzmann 2018"
__device__ double E_parameter(double w0, double g0, double i2s_transition) {
    double E;

    if (w0 > i2s_transition && g0 >= 0) {

        E = max(1.0,
                1.225 - 0.1582 * g0 - 0.1777 * w0 - 0.07465 * pow(1.0 * g0, 2.0) + 0.2351 * w0 * g0
                    - 0.05582 * pow(w0, 2.0));
    }
    else {
        E = 1.0;
    }
    return E;
}

//  calculates the transmission function
__device__ double trans_func(double epsi, double delta_tau, double w0, double g0, double E) {

    return exp(-1.0 / epsi * sqrt(E * (1.0 - w0 * g0) * (E - w0)) * delta_tau);
}

// calculates the G+ function
__device__ double
G_plus_func(double w0, double g0, double epsilon, double epsilon2, double mu_star, double E) {

    double num = w0 * (E * (1.0 - w0 * g0) + g0 * epsilon / epsilon2);

    double denom = E * pow(mu_star / epsilon, 2.0) * (E - w0) * (1.0 - w0 * g0) - 1.0;

    if (fabs(denom) < 1e-5) {
        denom = mu_star / epsilon;
        // printf("small divisor hit in G+\n");
    }

    double second_term = mu_star / epsilon + 1.0 / (E * (1.0 - w0 * g0));

    double third_term = epsilon * w0 * g0 / (epsilon2 * E * (1.0 - w0 * g0));

    double bracket = num / denom * second_term + third_term;

    double result = 0.5 * bracket;

    return result;
}

// calculates the G- function
__device__ double
G_minus_func(double w0, double g0, double epsilon, double epsilon2, double mu_star, double E) {
    double num = w0 * (E * (1.0 - w0 * g0) + g0 * epsilon / epsilon2);

    double denom = E * pow(mu_star / epsilon, 2.0) * (E - w0) * (1.0 - w0 * g0) - 1.0;

    if (fabs(denom) < 1e-5) {
        denom = mu_star / epsilon;
        // printf("small divisor hit in G-\n");
    }

    double second_term = mu_star / epsilon - 1.0 / (E * (1.0 - w0 * g0));

    double third_term = epsilon * w0 * g0 / (epsilon2 * E * (1.0 - w0 * g0));

    double bracket = num / denom * second_term - third_term;

    double result = 0.5 * bracket;

    return result;
}

__device__ double G_pm_denom(double w0, double g0, double epsi, double mu_star, double E) {
    double denom = E * pow(mu_star / epsi, 2.0) * (E - w0) * (1.0 - w0 * g0) - 1.0;

    return denom;
}


// limiting the values of the G_plus and G_minus coefficients to 1e8.
// This value is somewhat ad hoc from visual analysis. To justify, results are quite insensitive to this value.
// this is the originnal HELIOS limiter on G_pm, replaced by other limiters
__device__ double G_limiter(double G, bool debug) {

    if (abs(G) < 1e8) {
        return G;
    }
    else {
        if (debug) {
            printf("WARNING: G_functions are being artificially limited!!! \n");
        }
        return 1e8 * G / abs(G);
    }
}


// calculates the single scattering albedo w0
__device__ double single_scat_alb(double gas_scat_cross,
                                  double gas_abs,
                                  double meanmolmass,
                                  double fcloud,
                                  double cloud_scat_cross,
                                  double cloud_abs_cross,
                                  double w_0_limit) {

    return min((gas_scat_cross + fcloud * cloud_scat_cross)
                   / (gas_scat_cross + gas_abs * meanmolmass
                      + fcloud * (cloud_scat_cross + cloud_abs_cross)),
               w_0_limit);
}

// calculates the asymmetry parameter g0
__device__ double
g0_calc(double gas_scat_cross, double fcloud, double cloud_scat_cross, double g0_cloud) {
    return (fcloud * cloud_scat_cross * g0_cloud) / (gas_scat_cross + fcloud * cloud_scat_cross);
}

// calculates the two-stream coupling coefficient Zeta_minus with the scattering coefficient E
__device__ double zeta_minus(double w0, double g0, double E) {

    return 0.5 * (1.0 - sqrt((E - w0) / (E * (1.0 - w0 * g0))));
}


// calculates the two-stream coupling coefficient Zeta_plus with the scattering coefficient E
__device__ double zeta_plus(double w0, double g0, double E) {

    return 0.5 * (1.0 + sqrt((E - w0) / (E * (1.0 - w0 * g0))));
}


// calculation of transmission, w0, zeta-functions, and capital letters for the layer centers in the isothermal case
// kernel runs per wavelength bin, per wavelength sampling (?) and per layer
__global__ void trans_iso(double*       trans_wg,             // out
                          double*       delta_tau_wg,         // out
                          double*       M_term,               // out
                          double*       N_term,               // out
                          double*       P_term,               // out
                          double*       G_plus,               // out
                          double*       G_minus,              // out
                          double*       delta_colmass,        // in
                          double*       opac_wg_lay,          // in
                          double*       cloud_abs_cross_lay,  // in
                          double*       meanmolmass_lay,      // in
                          double*       scat_cross_lay,       // in
                          double*       cloud_scat_cross_lay, // in
                          double*       w_0,                  // out
                          double*       g0_wg,                // out
                          double*       g_0_cloud_lay,        // in
                          double        g_0_gas,
                          double        epsi,
                          double        epsilon2,
                          double*       zenith_angle_cols,
                          double*       mu_star_cols,
                          double        w_0_limit,
                          bool          scat,
                          int           nbin,
                          int           ny,
                          int           nlayer,
                          int           num_cols,
                          double        fcloud,
                          bool          clouds,
                          bool          scat_corr,
                          double        mu_star_wiggle_increment,
                          bool          G_pm_limiter,
                          double        G_pm_denom_limit_for_mu_star_wiggler,
                          bool          G_pm_limit_on_full_G_pm,
                          bool*         hit_G_pm_limit_global,
                          unsigned int* columns_wiggle,
                          unsigned int* columns_wiggle_request,
                          unsigned int  wiggle_request_iterator,
                          bool          debug,
                          double        i2s_transition) {
    // indices
    // wavelength bin
    int x = threadIdx.x + blockIdx.x * blockDim.x;

    // layer
    int i = threadIdx.z + blockIdx.z * blockDim.z;

    // compound column and sampling point index
    int cy = threadIdx.y + blockIdx.y * blockDim.y;

    int y = cy % ny;
    // column
    int c = cy / ny;


    if (x < nbin && y < ny && i < nlayer && c < num_cols) {

        if (columns_wiggle[c] != 0) {
            double ray_cross;
            double cloud_scat_cross;
            double g0       = g_0_gas;
            double g0_cloud = 0.0;
            double E        = 1.0;


            if (scat) {
                ray_cross = scat_cross_lay[x + nbin * i + c * nlayer * nbin];
                if (clouds)
                    cloud_scat_cross = cloud_scat_cross_lay[x];
                else
                    cloud_scat_cross = 0.0;
                // DBG: cloud_scat_cross = 0.0;
            }
            else {
                ray_cross        = 0.0;
                cloud_scat_cross = 0.0;
            }

            if (clouds) {
                g0_cloud = g_0_cloud_lay[x];
                g0       = g0_calc(ray_cross, fcloud, cloud_scat_cross, g0_cloud);
            }


            double w0 =
                single_scat_alb(ray_cross,
                                opac_wg_lay[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin],
                                meanmolmass_lay[i + c * nlayer],
                                fcloud,
                                cloud_scat_cross,
                                cloud_abs_cross_lay[x],
                                w_0_limit);


            w_0[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin]   = w0;
            g0_wg[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin] = g0;

            if (scat_corr) {
                E = E_parameter(w0, g0, i2s_transition);
            }

            delta_tau_wg[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin] =
                delta_colmass[i + c * nlayer]
                * (opac_wg_lay[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin]
                   + (ray_cross + fcloud * (cloud_abs_cross_lay[x] + cloud_scat_cross))
                         / meanmolmass_lay[i + c * nlayer]);

            double del_tau = delta_tau_wg[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin];
            trans_wg[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin] =
                trans_func(epsi, del_tau, w0, g0, E);
            double trans = trans_wg[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin];

            double zeta_min = zeta_minus(w0, g0, E);
            double zeta_pl  = zeta_plus(w0, g0, E);

            M_term[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin] =
                (zeta_min * zeta_min) * (trans * trans) - (zeta_pl * zeta_pl);
            N_term[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin] =
                zeta_pl * zeta_min * (1.0 - (trans * trans));
            P_term[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin] =
                ((zeta_min * zeta_min) - (zeta_pl * zeta_pl)) * trans;

            // DBG:
            // if (!isfinite(M_term[y + ny * x + ny * nbin * i]))
            //     printf(
            //         "abnormal M_term: %g, zeta_min: %g, trans: %g, zeta_pl: %g, "
            //         "epsi: %g, w0: %g, delta_tau: %g g0: %g, "
            //         "delta_colamss: %g, opac_wg_lay: %g, cloud_abs_cross_lay: %g, ray_cross: %g, "
            //         "cloud_scat_cross: %g, meanmolmass_lay: %g\n",
            //         M_term[y + ny * x + ny * nbin * i],
            //         zeta_min,
            //         trans,
            //         zeta_pl,
            //         epsi,
            //         w0,
            //         del_tau,
            //         g0,
            //         delta_colmass[i],
            //         opac_wg_lay[y + ny * x + ny * nbin * i],
            //         cloud_abs_cross_lay[i],
            //         ray_cross,
            //         cloud_scat_cross,
            //         meanmolmass_lay[i]);
            // if (!isfinite(N_term[y + ny * x + ny * nbin * i]))
            //     printf("abnormal N_term: %g, zeta_min: %g, trans: %g, zeta_pl: %g "
            //            "epsi: %g, w0: %g, delta_tau: %g, g0: %g\n",
            //            N_term[y + ny * x + ny * nbin * i],
            //            zeta_min,
            //            trans,
            //            zeta_pl,
            //            epsi,
            //            w0,
            //            del_tau,
            //            g0);

            double mu_star_orig = -zenith_angle_cols[c];

            // hit criteria, wiggle mu_star
            double zenith_angle_loc = acos(mu_star_orig);

            double mu_star_wiggle_factor = wiggle_request_iterator;

            double mu_star_used = cos(
                zenith_angle_loc + mu_star_wiggle_factor * mu_star_wiggle_increment / 180.0 * M_PI);

            double g_p = G_plus_func(w0, g0, epsi, epsilon2, mu_star_used, E);
            double g_m = G_minus_func(w0, g0, epsi, epsilon2, mu_star_used, E);
            /*    if (G_pm_limiter) {
                G_plus[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin] = G_limiter(g_p, debug);
                G_minus[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin] =
                    G_limiter(g_m, debug);
            }
            else {*/
            G_plus[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin]  = g_p;
            G_minus[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin] = g_m;
            /*}

            bool hit_limit = false;
            if (G_pm_limit_on_full_G_pm)
                hit_limit = (fabs(g_p) > G_pm_denom_limit_for_mu_star_wiggler)
                            || (fabs(g_m) > G_pm_denom_limit_for_mu_star_wiggler);
            else
                hit_limit = fabs(G_pm_denom(w0, g0, epsi, mu_star_used, E))
                            < G_pm_denom_limit_for_mu_star_wiggler;
            // Check G_pm criteria
            if (hit_limit) {
                // mark global iterator
                hit_G_pm_limit_global[0] = true;
                // we request that this column is recomputed at next iteration
                columns_wiggle_request[c] = 1;
                if (debug) {
                    printf(
                        "Hit G_pm denom limit, wiggle mu_star (%g) angle (%g) by %g degree to "
                        "(%g) "
                        "angle (%g) "
                        "(c: %d, l: %d, b: %d, w: %d) g_p, g_m (%g, %g)\n",
                        mu_star_orig,
                        zenith_angle_loc / M_PI * 180.0,
                        (mu_star_wiggle_factor + 1.0)
                            * mu_star_wiggle_increment, // + 1.0 because it gets applied at next iteration
                        mu_star_used,
                        zenith_angle_loc / M_PI * 180.0
                            + (mu_star_wiggle_factor + 1.0) * mu_star_wiggle_increment,
                        c,
                        i,
                        x,
                        y,
                        g_p,
                        g_m);
                }
            } */
            // debug printout when exploring NaN errors
            // if (!isfinite(g_p) || !isfinite(g_m)
            //     || !isfinite(M_term[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin])
            //     || !isfinite(N_term[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin])
            //     || !isfinite(P_term[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin])

            // ) {
            //     printf("%d, G_p: %g, G_m: %g, M: %g, N: %g, P: %g, trans: %g, z-: %g, z+: %g\n",
            //            x,
            //            g_p,
            //            g_m,
            //            M_term[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin],
            //            N_term[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin],
            //            P_term[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin],
            //            trans,
            //            zeta_min,
            //            zeta_pl);

            //     printf("%d: w0: %g, g0: %g dtau: %g gsc: %g csc: %g gabs: %g, cabs: %g cg0: %g "
            //            "mmm: %g "
            //            "dcm: "
            //            "%g mu*: %g, E: %g\n",
            //            x,
            //            w0,
            //            g0,
            //            del_tau,
            //            ray_cross,
            //            cloud_scat_cross,
            //            opac_wg_lay[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin],
            //            cloud_abs_cross_lay[x],
            //            g0_cloud,
            //            meanmolmass_lay[i + c * nlayer],
            //            delta_colmass[i + c * nlayer],
            //            mu_star_used,
            //            E);
            // }
            mu_star_cols[c] = mu_star_used;
        }
    }
}

// calculation of transmission, w0, zeta-functions, and capital letters for the non-isothermal case
__global__ void trans_noniso(double*       trans_wg_upper,
                             double*       trans_wg_lower,
                             double*       delta_tau_wg_upper,
                             double*       delta_tau_wg_lower,
                             double*       M_upper,
                             double*       M_lower,
                             double*       N_upper,
                             double*       N_lower,
                             double*       P_upper,
                             double*       P_lower,
                             double*       G_plus_upper,
                             double*       G_plus_lower,
                             double*       G_minus_upper,
                             double*       G_minus_lower,
                             double*       delta_col_upper,
                             double*       delta_col_lower,
                             double*       opac_wg_lay,
                             double*       opac_wg_int,
                             double*       cloud_abs_cross_lay,
                             double*       cloud_abs_cross_int,
                             double*       meanmolmass_lay,
                             double*       meanmolmass_int,
                             double*       scat_cross_lay,
                             double*       scat_cross_int,
                             double*       cloud_scat_cross_lay,
                             double*       cloud_scat_cross_int,
                             double*       w_0_upper,
                             double*       w_0_lower,
                             double*       g0_wg_upper,
                             double*       g0_wg_lower,
                             double*       g_0_cloud_lay,
                             double*       g_0_cloud_int,
                             double        g_0_gas,
                             double        epsi,
                             double        epsilon2,
                             double*       zenith_angle_cols,
                             double*       mu_star_cols,
                             double        w_0_limit,
                             bool          scat,
                             int           nbin,
                             int           ny,
                             int           nlayer,
                             int           num_cols,
                             double        fcloud,
                             bool          clouds,
                             bool          scat_corr,
                             double        mu_star_wiggle_increment,
                             bool          G_pm_limiter,
                             double        G_pm_denom_limit_for_mu_star_wiggler,
                             bool          G_pm_limit_on_full_G_pm,
                             bool*         hit_G_pm_limit_global,
                             unsigned int* columns_wiggle,
                             unsigned int* columns_wiggle_request,
                             unsigned int  wiggle_request_iterator,
                             bool          debug,
                             double        i2s_transition,
                             int           column_idx) {

    // wavelength bin
    int x = threadIdx.x + blockIdx.x * blockDim.x;
    // layer
    int i = threadIdx.z + blockIdx.z * blockDim.z;

    // compound column and sampling point index
    int cy = threadIdx.y + blockIdx.y * blockDim.y;

    int y = cy % ny;
    // column
    int c = cy / ny;

    if (x < nbin && y < ny && i < nlayer && c < num_cols) {
        // compute at least once
        if (columns_wiggle[c] != 0) {
            int ninterface = nlayer + 1;

            double ray_cross_up;
            double ray_cross_low;
            double cloud_scat_cross_up;
            double cloud_scat_cross_low;
            double g0_up        = g_0_gas;
            double g0_low       = g_0_gas;
            double g0_cloud_up  = 0.0;
            double g0_cloud_low = 0.0;
            double E_up         = 1.0;
            double E_low        = 1.0;

            if (scat) {
                ray_cross_up = (scat_cross_lay[x + nbin * i + c * nlayer * nbin]
                                + scat_cross_int[x + nbin * (i + 1) + c * ninterface * nbin])
                               / 2.0;
                ray_cross_low = (scat_cross_int[x + nbin * i + c * ninterface * nbin]
                                 + scat_cross_lay[x + nbin * i + c * nlayer * nbin])
                                / 2.0;
                if (clouds) {
                    cloud_scat_cross_up = (cloud_scat_cross_lay[x] + cloud_scat_cross_int[x]) / 2.0;
                    cloud_scat_cross_low =
                        (cloud_scat_cross_int[x] + cloud_scat_cross_lay[x]) / 2.0;
                }
                else {
                    cloud_scat_cross_up  = 0.0;
                    cloud_scat_cross_low = 0.0;
                }
            }
            else {
                ray_cross_up         = 0;
                ray_cross_low        = 0;
                cloud_scat_cross_up  = 0;
                cloud_scat_cross_low = 0;
            }

            if (clouds) {
                // For use with per altitude bin, need to add i*nbin indices.
                g0_cloud_up  = (g_0_cloud_lay[x] + g_0_cloud_int[x]) / 2.0;
                g0_cloud_low = (g_0_cloud_int[x] + g_0_cloud_lay[x]) / 2.0;

                g0_up  = g0_calc(ray_cross_up, fcloud, cloud_scat_cross_up, g0_cloud_up);
                g0_low = g0_calc(ray_cross_low, fcloud, cloud_scat_cross_low, g0_cloud_low);
            }


            g0_wg_upper[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin] = g0_up;
            g0_wg_lower[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin] = g0_low;

            double opac_up =
                (opac_wg_lay[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin]
                 + opac_wg_int[y + ny * x + ny * nbin * (i + 1) + c * ninterface * ny * nbin])
                / 2.0;
            double opac_low = (opac_wg_int[y + ny * x + ny * nbin * i + c * ninterface * ny * nbin]
                               + opac_wg_lay[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin])
                              / 2.0;
            double cloud_abs_cross_up  = 0.0;
            double cloud_abs_cross_low = 0.0;

            if (clouds) {
                cloud_abs_cross_up  = (cloud_abs_cross_lay[x] + cloud_abs_cross_int[x]) / 2.0;
                cloud_abs_cross_low = (cloud_abs_cross_int[x] + cloud_abs_cross_lay[x]) / 2.0;
            }

            double meanmolmass_up =
                (meanmolmass_lay[i + c * nlayer] + meanmolmass_int[i + 1 + c * ninterface]) / 2.0;
            double meanmolmass_low =
                (meanmolmass_int[i + c * ninterface] + meanmolmass_lay[i + c * nlayer]) / 2.0;

            double w_0_up = single_scat_alb(ray_cross_up,
                                            opac_up,
                                            meanmolmass_up,
                                            fcloud,
                                            cloud_scat_cross_up,
                                            cloud_abs_cross_up,
                                            w_0_limit);

            w_0_upper[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin] = w_0_up;

            double w_0_low = single_scat_alb(ray_cross_low,
                                             opac_low,
                                             meanmolmass_low,
                                             fcloud,
                                             cloud_scat_cross_low,
                                             cloud_abs_cross_low,
                                             w_0_limit);

            w_0_lower[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin] = w_0_low;


            if (scat_corr) {
                E_up  = E_parameter(w_0_up, g0_up, i2s_transition);
                E_low = E_parameter(w_0_low, g0_low, i2s_transition);
            }

            delta_tau_wg_upper[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin] =
                delta_col_upper[i + c * nlayer]
                * (opac_up
                   + (ray_cross_up + fcloud * (cloud_abs_cross_up + cloud_scat_cross_up))
                         / meanmolmass_up);
            double del_tau_up =
                delta_tau_wg_upper[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin];
            delta_tau_wg_lower[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin] =
                delta_col_lower[i + c * nlayer]
                * (opac_low
                   + (ray_cross_low + fcloud * (cloud_abs_cross_low + cloud_scat_cross_low))
                         / meanmolmass_low);
            double del_tau_low =
                delta_tau_wg_lower[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin];

            trans_wg_upper[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin] =
                trans_func(epsi, del_tau_up, w_0_up, g0_up, E_up);
            double trans_up = trans_wg_upper[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin];
            trans_wg_lower[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin] =
                trans_func(epsi, del_tau_low, w_0_low, g0_low, E_low);
            double trans_low = trans_wg_lower[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin];

            double zeta_min_up  = zeta_minus(w_0_up, g0_up, E_up);
            double zeta_min_low = zeta_minus(w_0_low, g0_low, E_low);
            double zeta_pl_up   = zeta_plus(w_0_up, g0_up, E_up);
            double zeta_pl_low  = zeta_plus(w_0_low, g0_low, E_low);

            M_upper[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin] =
                (zeta_min_up * zeta_min_up) * (trans_up * trans_up) - (zeta_pl_up * zeta_pl_up);
            M_lower[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin] =
                (zeta_min_low * zeta_min_low) * (trans_low * trans_low)
                - (zeta_pl_low * zeta_pl_low);
            N_upper[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin] =
                zeta_pl_up * zeta_min_up * (1.0 - (trans_up * trans_up));
            N_lower[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin] =
                zeta_pl_low * zeta_min_low * (1.0 - (trans_low * trans_low));
            P_upper[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin] =
                ((zeta_min_up * zeta_min_up) - (zeta_pl_up * zeta_pl_up)) * trans_up;
            P_lower[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin] =
                ((zeta_min_low * zeta_min_low) - (zeta_pl_low * zeta_pl_low)) * trans_low;

            double mu_star_orig = -zenith_angle_cols[c];

            // hit criteria, wiggle mu_star
            double zenith_angle_loc = acos(mu_star_orig);

            double mu_star_wiggle_factor = wiggle_request_iterator;

            double mu_star_used = cos(
                zenith_angle_loc + mu_star_wiggle_factor * mu_star_wiggle_increment / 180.0 * M_PI);


            double g_p_u = G_plus_func(w_0_up, g0_up, epsi, epsilon2, mu_star_used, E_up);
            double g_p_l = G_plus_func(w_0_low, g0_low, epsi, epsilon2, mu_star_used, E_low);

            double g_m_u = G_minus_func(w_0_up, g0_up, epsi, epsilon2, mu_star_used, E_up);
            double g_m_l = G_minus_func(w_0_low, g0_low, epsi, epsilon2, mu_star_used, E_low);
            /*if (G_pm_limiter) {
                G_plus_upper[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin] =
                    G_limiter(g_p_u, debug);
                G_plus_lower[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin] =
                    G_limiter(g_p_l, debug);
                G_minus_upper[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin] =
                    G_limiter(g_m_u, debug);
                G_minus_lower[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin] =
                    G_limiter(g_m_l, debug);
            }
            else {
            */
            G_plus_upper[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin]  = g_p_u;
            G_plus_lower[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin]  = g_p_l;
            G_minus_upper[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin] = g_m_u;
            G_minus_lower[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin] = g_m_l;
            /*
            }

            bool hit_limit = false;
            if (G_pm_limit_on_full_G_pm)
                hit_limit = (fabs(g_p_u) > G_pm_denom_limit_for_mu_star_wiggler)
                            || (fabs(g_p_l) > G_pm_denom_limit_for_mu_star_wiggler)
                            || (fabs(g_m_u) > G_pm_denom_limit_for_mu_star_wiggler)
                            || (fabs(g_m_l) > G_pm_denom_limit_for_mu_star_wiggler);
            else
                hit_limit = (fabs(G_pm_denom(w_0_up, g0_up, epsi, mu_star_used, E_up))
                             < G_pm_denom_limit_for_mu_star_wiggler)
                            || (fabs(G_pm_denom(w_0_low, g0_low, epsi, mu_star_used, E_low))
                                < G_pm_denom_limit_for_mu_star_wiggler);
            // Check G_pm criteria
            if (hit_limit) {
                hit_G_pm_limit_global[0] = true;
                // we request that this column is recomputed at next iteration
                columns_wiggle_request[c] = 1;

                if (debug) {
                    printf(
                        "Hit G_pm denom limit, column %d, wiggle mu_star (%g) angle (%g) by %g "
                        "degree to "
                        "(%g) "
                        "angle (%g) "
                        "(c: %d, l: %d, b: %d, w: %d) g_p_u, g_m_u (%g, %g) g_p_l, g_m_l (%g, "
                        "%g)\n",
                        column_idx + c,
                        mu_star_orig,
                        zenith_angle_loc / M_PI * 180.0,
                        (mu_star_wiggle_factor + 1.0)
                            * mu_star_wiggle_increment, // + 1.0 because it gets applied at next iteration
                        mu_star_used,
                        zenith_angle_loc / M_PI * 180.0
                            + (mu_star_wiggle_factor + 1.0) * mu_star_wiggle_increment,
                        c,
                        i,
                        x,
                        y,
                        g_p_u,
                        g_m_u,
                        g_p_l,
                        g_m_l);
                }
            } */

            // printf("%d: w0: %g, g0: %g dtau: %g gsc: %g csc: %g gabs: %g, cabs: %g cg0: %g "
            //        "mmm: %g "
            //        "dcm: "
            //        "%g mu*: %g, E: %g\n",
            //        x,
            //        w0,
            //        g0,
            //        del_tau,
            //        ray_cross,
            //        cloud_scat_cross,
            //        opac_wg_lay[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin],
            //        cloud_abs_cross_lay[x],
            //        g0_cloud,
            //        meanmolmass_lay[i + c * nlayer],
            //        delta_colmass[i + c * nlayer],
            //        mu_star_used,
            //        E);
            if (!isfinite(g_p_u) || !isfinite(g_m_u) || !isfinite(g_p_l) || !isfinite(g_m_l)
                || !isfinite(M_upper[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin])
                || !isfinite(M_lower[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin])
                || !isfinite(N_upper[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin])
                || !isfinite(N_lower[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin])
                || !isfinite(P_upper[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin])
                || !isfinite(P_lower[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin])

            )
                printf("(c: %d, l: %d, b: %d, w: %d) g_pm_u (%g, %g) g_pm_l (%g, "
                       "%g) E: (%g, %g), w0: (%g, %g), g0: (%g, %g) mu*: %g "
                       "M: (%g, %g)"
                       "N: (%g, %g)"
                       "P: (%g, %g)"
                       "tr: (%g, %g)"
                       "dt: (%g, %g)"
                       "\n",
                       c,
                       i,
                       x,
                       y,
                       g_p_u,
                       g_m_u,
                       g_p_l,
                       g_m_l,
                       E_up,
                       E_low,
                       w_0_up,
                       w_0_low,
                       g0_up,
                       g0_low,
                       mu_star_used,
                       M_upper[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin],
                       M_lower[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin],
                       N_upper[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin],
                       N_lower[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin],
                       P_upper[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin],
                       P_lower[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin],
                       trans_wg_upper[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin],
                       trans_wg_lower[y + ny * x + ny * nbin * i + c * nlayer * ny * nbin],
                       del_tau_up,
                       del_tau_low);

            mu_star_cols[c] = mu_star_used;
        }
    }
}
