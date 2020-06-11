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
// ESP -  Exoclimes Simulation Platform. (version 1.0)
//
//
//
// Method: Boundary layer (surface friction) physics module
//
//
// Known limitations: - Runs in a single GPU.
//
// Known issues: None
//
//
// If you use this code please cite the following reference:
//
//       [1] Mendonca, J.M., Grimm, S.L., Grosheintz, L., & Heng, K., ApJ, 829, 115, 2016
//
// Current Code Owner: Joao Mendonca, EEG. joao.mendonca@csh.unibe.ch
//                     Russell Deitrick (russell.deitrick@csh.unibe.ch)
//                     Urs Schroffenegger (urs.schroffenegger@csh.unibe.ch)
//
// History:
// Version Date       Comment
// ======= ====       =======
//
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////
#include "boundary_layer.h"

boundary_layer::boundary_layer() : bl_type_str("RayleighHS") {
}

boundary_layer::~boundary_layer() {
}

void boundary_layer::print_config() {
    log::printf("  Boundary layer module\n");

    // basic properties
    log::printf("    bl_type                    = %s \n", bl_type_str.c_str());
    log::printf("    surf_drag                  = %e 1/s\n", surf_drag_config);
    log::printf("    bl_sigma                   = %f \n", bl_sigma_config);

    log::printf("\n");
}

bool boundary_layer::initialise_memory(const ESP &              esp,
                                       device_RK_array_manager &phy_modules_core_arrays) {

    // cudaMalloc((void **)&dvdz_tmp, 3 * esp.nvi * esp.point_num * sizeof(double));
    cudaMalloc((void **)&d2vdz2_tmp, 3 * esp.nv * esp.point_num * sizeof(double));

    cudaMalloc((void **)&atmp, esp.nv * esp.point_num * sizeof(double));
    cudaMalloc((void **)&btmp, esp.nv * esp.point_num * sizeof(double));
    cudaMalloc((void **)&ctmp, esp.nv * esp.point_num * sizeof(double));
    cudaMalloc((void **)&cpr_tmp, esp.nv * esp.point_num * sizeof(double));
    cudaMalloc((void **)&dtmp, 3 * esp.nv * esp.point_num * sizeof(double));
    cudaMalloc((void **)&dpr_tmp, 3 * esp.nv * esp.point_num * sizeof(double));
    cudaMalloc((void **)&RiB_d, esp.nvi * esp.point_num * sizeof(double));
    cudaMalloc((void **)&KM_d, esp.nvi * esp.point_num * sizeof(double));
    cudaMalloc((void **)&KH_d, esp.nvi * esp.point_num * sizeof(double));
    cudaMalloc((void **)&CD_d, esp.point_num * sizeof(double));
    cudaMalloc((void **)&vh_lowest_d, esp.point_num * sizeof(double));
    cudaMalloc((void **)&zeta_d, esp.nvi * esp.point_num * sizeof(double));

    cudaMalloc((void **)&bl_top_lev_d, esp.point_num * sizeof(int));
    cudaMalloc((void **)&bl_top_height_d, esp.point_num * sizeof(double));

    RiB_h           = (double *)malloc(esp.nvi * esp.point_num * sizeof(double));
    KM_h            = (double *)malloc(esp.nvi * esp.point_num * sizeof(double));
    KH_h            = (double *)malloc(esp.nvi * esp.point_num * sizeof(double));
    bl_top_lev_h    = (int *)malloc(esp.point_num * sizeof(int));
    bl_top_height_h = (double *)malloc(esp.point_num * sizeof(int));
    CD_h            = (double *)malloc(esp.point_num * sizeof(double));

    cudaMemset(atmp, 0, sizeof(double) * esp.point_num * esp.nv);
    cudaMemset(btmp, 0, sizeof(double) * esp.point_num * esp.nv);
    cudaMemset(ctmp, 0, sizeof(double) * esp.point_num * esp.nv);
    cudaMemset(cpr_tmp, 0, sizeof(double) * esp.point_num * esp.nv);
    cudaMemset(dtmp, 0, sizeof(double) * 3 * esp.point_num * esp.nv);
    cudaMemset(dpr_tmp, 0, sizeof(double) * 3 * esp.point_num * esp.nv);
    cudaMemset(RiB_d, 0, sizeof(double) * esp.point_num * esp.nvi);
    cudaMemset(KM_d, 0, sizeof(double) * esp.point_num * esp.nvi);
    cudaMemset(KH_d, 0, sizeof(double) * esp.point_num * esp.nvi);
    cudaMemset(bl_top_lev_d, 0, sizeof(int) * esp.point_num);
    cudaMemset(bl_top_height_h, 0, sizeof(double) * esp.point_num);
    cudaMemset(CD_d, 0, sizeof(double) * esp.point_num);
    cudaMemset(vh_lowest_d, 0, sizeof(double) * esp.point_num);
    cudaMemset(zeta_d, 0, sizeof(double) * esp.point_num * esp.nvi);

    return true;
}


bool boundary_layer::free_memory() {

    // cudaFree(dvdz_tmp);
    cudaFree(d2vdz2_tmp);
    cudaFree(atmp);
    cudaFree(btmp);
    cudaFree(ctmp);
    cudaFree(cpr_tmp);
    cudaFree(dtmp);
    cudaFree(dpr_tmp);
    cudaFree(RiB_d);
    cudaFree(KM_d);
    cudaFree(KH_d);
    cudaFree(CD_d);
    cudaFree(vh_lowest_d);
    cudaFree(bl_top_lev_d);
    cudaFree(zeta_d);
    cudaFree(bl_top_height_d);

    return true;
}

bool boundary_layer::initial_conditions(const ESP &esp, const SimulationSetup &sim, storage *s) {
    bool config_OK = true;

    bl_type = RAYLEIGHHS;
    if (bl_type_str == "RayleighHS") {
        bl_type = RAYLEIGHHS;
        config_OK &= true;
    }
    else if (bl_type_str == "MoninObukhov" || bl_type_str == "MO") {
        bl_type = MONINOBUKHOV;
        config_OK &= true;
    }
    else if (bl_type_str == "EkmanSpiral" || bl_type_str == "Ekman") {
        bl_type = EKMANSPIRAL;
        config_OK &= true;
    }
    else {
        log::printf("bl_type config item not recognised: [%s]\n", bl_type_str.c_str());
        config_OK &= false;
    }

    if (!config_OK) {
        log::printf("Error in configuration file\n");
        exit(-1);
    }

    BLSetup(esp,
            sim,
            bl_type,
            surf_drag_config,
            bl_sigma_config,
            Ri_crit_config,
            z_rough_config,
            z_therm_config,
            f_surf_layer_config);

    return true;
}

bool boundary_layer::phy_loop(ESP &                  esp,
                              const SimulationSetup &sim,
                              int                    nstep, // Step number
                              double                 time_step) {           // Time-step [s]

    //  Number of threads per block.
    const int NTH = 256;

    //  Specify the block sizes.
    dim3 NB((esp.point_num / NTH) + 1, esp.nv, 1);
    dim3 NBi((esp.point_num / NTH) + 1, esp.nvi, 1);
    dim3 NBLEV((esp.point_num / NTH) + 1, 1, 1);

    if (bl_type == RAYLEIGHHS) {
        rayleighHS<<<NB, NTH>>>(esp.Mh_d,
                                esp.pressure_d,
                                esp.Rho_d,
                                esp.Altitude_d,
                                surf_drag,
                                bl_sigma,
                                sim.Gravit,
                                time_step,
                                esp.point_num);
    }
    else if (bl_type == MONINOBUKHOV) {
        CalcRiB<<<NBLEV, NTH>>>(esp.pressure_d,
                                esp.Rho_d,
                                esp.Mh_d,
                                esp.Tsurface_d,
                                esp.Altitude_d,
                                esp.Altitudeh_d,
                                sim.Rd,
                                sim.Cp,
                                sim.P_Ref,
                                sim.Gravit,
                                Ri_crit,
                                z_rough,
                                z_therm,
                                RiB_d,
                                bl_top_lev_d,
                                bl_top_height_d,
                                CD_d,
                                zeta_d,
                                vh_lowest_d,
                                esp.point_num,
                                esp.nv);

        CalcKM<<<NBi, NTH>>>(RiB_d,
                             zeta_d,
                             CD_d,
                             bl_top_height_d,
                             bl_top_lev_d,
                             vh_lowest_d,
                             esp.Altitude_d,
                             esp.Altitudeh_d,
                             Ri_crit,
                             z_rough,
                             z_therm,
                             f_surf_layer,
                             KM_d,
                             esp.point_num);

        // TO DO
        // calc BL height from RiB
        // adjust Tsurface for sensible heat flux
        // how to adjust pressure? adjust pt first, then compute pressure? or is there a shortcut?
        // update pressure (implicitly) here, or add to qheat?

        MomentumDiff_Impl<<<NBLEV, NTH>>>(esp.Mh_d,
                                          esp.pressure_d,
                                          esp.Rho_d,
                                          esp.Altitude_d,
                                          esp.Altitudeh_d,
                                          atmp,
                                          btmp,
                                          ctmp,
                                          cpr_tmp,
                                          dtmp,
                                          dpr_tmp,
                                          KM_d,
                                          zbl,
                                          time_step,
                                          esp.point_num,
                                          esp.nv,
                                          bl_top_lev_d);
    }
    else if (bl_type == EKMANSPIRAL) {
        // cudaMemset(dvdz_tmp, 0, sizeof(double) * 3 * esp.point_num * esp.nvi);
        // cudaMemset(d2vdz2_tmp, 0, sizeof(double) * 3 * esp.point_num * esp.nv);

        // ConstKMEkman<<<NBLEV, NTH>>>(esp.Mh_d,
        //                              esp.pressure_d,
        //                              esp.Rho_d,
        //                              esp.Altitude_d,
        //                              esp.Altitudeh_d,
        //                              d2vdz2_tmp,
        //                              KMconst,
        //                              zbl,
        //                              time_step,
        //                              esp.point_num,
        //                              esp.nv);

        // CalcRiB<<<NBLEV, NTH>>>(esp.pressure_d,
        //                         esp.Rho_d,
        //                         esp.Mh_d,
        //                         esp.Tsurface_d,
        //                         esp.Altitude_d,
        //                         esp.Altitudeh_d,
        //                         sim.Rd,
        //                         sim.Cp,
        //                         sim.P_Ref,
        //                         sim.Gravit,
        //                         Ri_crit,
        //                         RiB_d,
        //                         bl_top_lev_d,
        //                         bl_top_height_d,
        //                         esp.point_num,
        //                         esp.nv);

        // TO DO
        // need KM array, KH array, general thomas solver for KM, KH
        // calc BL height from RiB
        // adjust Tsurface for sensible heat flux
        // how to adjust pressure? adjust pt first, then compute pressure? or is there a shortcut?
        // update pressure (implicitly) here, or add to qheat?

        MomentumDiff_Impl<<<NBLEV, NTH>>>(esp.Mh_d,
                                          esp.pressure_d,
                                          esp.Rho_d,
                                          esp.Altitude_d,
                                          esp.Altitudeh_d,
                                          atmp,
                                          btmp,
                                          ctmp,
                                          cpr_tmp,
                                          dtmp,
                                          dpr_tmp,
                                          KM_d,
                                          zbl,
                                          time_step,
                                          esp.point_num,
                                          esp.nv,
                                          bl_top_lev_d);
    }


    return true;
}

bool boundary_layer::configure(config_file &config_reader) {

    config_reader.append_config_var("bl_type", bl_type_str, string(bl_type_default)); //

    // coefficient of drag strength
    config_reader.append_config_var("surf_drag", surf_drag_config, surf_drag_config);

    // percent of surface pressure where bl starts
    config_reader.append_config_var("bl_sigma", bl_sigma_config, bl_sigma_config);

    // Critical Richardson number
    config_reader.append_config_var("Ri_crit", Ri_crit_config, Ri_crit_config);

    // roughness and thermal length scaling
    config_reader.append_config_var("z_rough", z_rough_config, z_rough_config);
    config_reader.append_config_var("z_therm", z_therm_config, z_therm_config);

    // surface layer fraction
    config_reader.append_config_var("f_surf_layer", f_surf_layer_config, f_surf_layer_config);

    return true;
}

bool boundary_layer::store(const ESP &esp, storage &s) {
    if (bl_type == MONINOBUKHOV) {
        cudaMemcpy(RiB_h, RiB_d, esp.point_num * esp.nvi * sizeof(double), cudaMemcpyDeviceToHost);
        s.append_table(RiB_h, esp.point_num * esp.nvi, "/RiB", " ", "bulk Richardson number");

        cudaMemcpy(CD_h, CD_d, esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
        s.append_table(CD_h, esp.point_num, "/CD", " ", "surface drag coefficient");

        cudaMemcpy(KM_h, KM_d, esp.point_num * esp.nvi * sizeof(double), cudaMemcpyDeviceToHost);
        s.append_table(KM_h,
                       esp.point_num * esp.nvi,
                       "/KM",
                       "m^2/s",
                       "turbulent diffusion coefficient for momentum");
    }

    return true;
}

bool boundary_layer::store_init(storage &s) {
    s.append_value(surf_drag, "/surf_drag", "1/s", "surface drag coefficient");
    s.append_value(bl_sigma, "/bl_sigma", " ", "boundary layer sigma coordinate");

    s.append_value(Ri_crit, "/Ri_crit", " ", "critical Richardson number");

    s.append_value(z_rough, "/z_rough", "m", "boundary layer roughness length");

    s.append_value(z_therm, "/z_therm", "m", "boundary layer thermal length scale");

    s.append_value(f_surf_layer, "/f_surf_layer", " ", "surface layer fraction");

    return true;
}

void boundary_layer::BLSetup(const ESP &            esp,
                             const SimulationSetup &sim,
                             int                    bl_type,
                             double                 surf_drag_,
                             double                 bl_sigma_,
                             double                 Ri_crit_,
                             double                 z_rough_,
                             double                 z_therm_,
                             double                 f_surf_layer_) {
    if (bl_type == RAYLEIGHHS) {
        surf_drag = surf_drag_;
        bl_sigma  = bl_sigma_;
    }
    else if (bl_type == EKMANSPIRAL) {
        double KMconst = 12.5;

        zbl = bl_sigma_ * sim.Top_altitude;

        for (int id = 0; id < esp.point_num; id++) {
            int lev = 0;
            while (esp.Altitude_h[lev] < zbl) {
                bl_top_lev_h[id] = lev;
                lev++;
            }

            for (lev = 0; lev < esp.nvi; lev++) {
                KM_h[id * esp.nvi + lev] = KMconst;
            }
        }
        cudaMemcpy(KM_d, KM_h, esp.nvi * esp.point_num * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(bl_top_lev_d, bl_top_lev_h, esp.point_num * sizeof(int), cudaMemcpyHostToDevice);
    }
    else if (bl_type == MONINOBUKHOV) {
        //cudaMemset(KM_d, 0, esp.nvi * esp.point_num * sizeof(double));
        //cudaMemset(bl_top_lev_d, 0, esp.point_num * sizeof(int));

        //temporary (these will go in MO scheme)
        Ri_crit      = Ri_crit_;
        z_rough      = z_rough_;
        z_therm      = z_therm_;
        f_surf_layer = f_surf_layer_;
    }
    // printf("%f, %f, %d\n", zbl, esp.Altitude_h[bl_top_lev], bl_top_lev);
}


__global__ void rayleighHS(double *Mh_d,
                           double *pressure_d,
                           double *Rho_d,
                           double *Altitude_d,
                           double  surf_drag,
                           double  bl_sigma,
                           double  Gravit,
                           double  time_step,
                           int     num) {

    int id  = blockIdx.x * blockDim.x + threadIdx.x;
    int nv  = gridDim.y;
    int lev = blockIdx.y;

    if (id < num) {
        double sigma;
        double sigmab = bl_sigma;
        double kf     = surf_drag;
        double kv_hs;
        double ps, pre;
        double psm1;

        //      Calculates surface pressure
        psm1 = pressure_d[id * nv + 1]
               - Rho_d[id * nv + 0] * Gravit * (-Altitude_d[0] - Altitude_d[1]);
        ps = 0.5 * (pressure_d[id * nv + 0] + psm1);

        pre   = pressure_d[id * nv + lev];
        sigma = (pre / ps);

        //      Momentum dissipation constant.
        kv_hs = kf * max(0.0, (sigma - sigmab) / (1.0 - sigmab));

        //      Update momenta
        for (int k = 0; k < 3; k++)
            Mh_d[id * 3 * nv + lev * 3 + k] =
                Mh_d[id * 3 * nv + lev * 3 + k] / (1.0 + kv_hs * time_step);

        // Wh_d[id * (nv + 1) + lev + k] = Wh_d[id * (nv + 1) + lev + k] / (1.0 + kv_hs * time_step);
    }
}

__global__ void ConstKMEkman(double *Mh_d,
                             double *pressure_d,
                             double *Rho_d,
                             double *Altitude_d,
                             double *Altitudeh_d,
                             double *d2vdz2_tmp,
                             double  KMconst,
                             double  zbl,
                             double  time_step,
                             int     num,
                             int     nv) {

    int id = blockIdx.x * blockDim.x + threadIdx.x;
    int lev;

    if (id < num) {
        for (int k = 0; k < 3; k++) {
            // dvdz_tmp[id * 3 * (nv + 1) + 0 * 3 + k]  = 0; //boundary condition
            // dvdz_tmp[id * 3 * (nv + 1) + nv * 3 + k] = 0;
            // for (lev = 1; lev < nv; lev++) {
            //     //first derivative at interfaces (half-layers)
            //     dvdz_tmp[id * 3 * (nv + 1) + lev * 3 + k] =
            //         (Mh_d[id * 3 * nv + lev * 3 + k] / Rho_d[id * nv + lev]
            //          - Mh_d[id * 3 * nv + (lev - 1) * 3 + k] / Rho_d[id * nv + lev - 1])
            //         / (Altitude_d[lev] - Altitude_d[lev - 1]);
            // }
            // for (lev = 0; lev < nv; lev++) {
            //     d2vdz2_tmp[id * 3 * nv + lev * 3 + k] = (dvdz_tmp[id * 3 * nv + (lev + 1) * 3 + k]
            //                                              - dvdz_tmp[id * 3 * nv + lev * 3 + k])
            //                                             / (Altitudeh_d[lev + 1] - Altitudeh_d[lev]);
            //     Mh_d[id * 3 * nv + lev * 3 + k] += -Rho_d[id * nv + lev] * KMconst
            //                                        * d2vdz2_tmp[id * 3 * nv + lev * 3 + k]
            //                                        * time_step;
            // }
            for (lev = 0; lev < nv; lev++) {
                if (Altitude_d[lev] < zbl) {
                    if (lev == 0) { //lowest layer, v at lowest boundary = 0, dz0 = Altitude0
                        d2vdz2_tmp[id * 3 * nv + lev * 3 + k] =
                            ((Mh_d[id * 3 * nv + (lev + 1) * 3 + k]
                              - Mh_d[id * 3 * nv + (lev)*3 + k])
                                 / (Altitude_d[lev + 1] - Altitude_d[lev])
                             - (Mh_d[id * 3 * nv + (lev)*3 + k]) / (Altitude_d[lev]))
                            / (Altitudeh_d[lev + 1] - Altitudeh_d[lev]);
                    }
                    // else if (lev == nv - 1) { //top layer,
                    //     ((Mh_d[id * 3 * nv + (lev + 1) * 3 + k] / Rho_d[id * nv + lev + 1]
                    //       - Mh_d[id * 3 * nv + (lev)*3 + k] / Rho_d[id * nv + lev])
                    //          / (Altitude_d[lev + 1] - Altitude_d[lev])
                    //      - (Mh_d[id * 3 * nv + (lev)*3 + k] / Rho_d[id * nv + lev]
                    //         - Mh_d[id * 3 * nv + (lev - 1) * 3 + k] / Rho_d[id * nv + lev - 1])
                    //            / (Altitude_d[lev] - Altitude_d[lev - 1]))
                    //         / (Altitudeh_d[lev + 1] - Altitudeh_d[lev]);
                    // }
                    else { //might need to add a term to layer above to conserve momentum
                        d2vdz2_tmp[id * 3 * nv + lev * 3 + k] =
                            ((Mh_d[id * 3 * nv + (lev + 1) * 3 + k]
                              - Mh_d[id * 3 * nv + (lev)*3 + k])
                                 / (Altitude_d[lev + 1] - Altitude_d[lev])
                             - (Mh_d[id * 3 * nv + (lev)*3 + k]
                                - Mh_d[id * 3 * nv + (lev - 1) * 3 + k])
                                   / (Altitude_d[lev] - Altitude_d[lev - 1]))
                            / (Altitudeh_d[lev + 1] - Altitudeh_d[lev]);
                    }
                    Mh_d[id * 3 * nv + lev * 3 + k] +=
                        KMconst * d2vdz2_tmp[id * 3 * nv + lev * 3 + k] * time_step;
                }
            }
        }
    }
}

__global__ void MomentumDiff_Impl(double *Mh_d,
                                  double *pressure_d,
                                  double *Rho_d,
                                  double *Altitude_d,
                                  double *Altitudeh_d,
                                  double *atmp,
                                  double *btmp,
                                  double *ctmp,
                                  double *cpr_tmp,
                                  double *dtmp,
                                  double *dpr_tmp,
                                  double *KM_d,
                                  double  zbl,
                                  double  time_step,
                                  int     num,
                                  int     nv,
                                  int *   bl_top_lev_d) {

    //should create check on stability of thomas algorithm

    int id = blockIdx.x * blockDim.x + threadIdx.x;
    int lev;

    if (id < num) {
        for (lev = 0; lev < bl_top_lev_d[id] + 1; lev++) {
            //forward sweep
            if (lev == 0) { //lowest layer, v at lowest boundary = 0, dz0 = Altitude0
                atmp[id * nv + lev] = 0;
                btmp[id * nv + lev] =
                    -(1.0 / (Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                          * (KM_d[id * (nv + 1) + lev + 1] / (Altitude_d[lev + 1] - Altitude_d[lev])
                             + KM_d[id * (nv + 1) + lev] / Altitude_d[lev])
                      + 1.0 / time_step);
                ctmp[id * nv + lev] = KM_d[id * (nv + 1) + lev + 1]
                                      / ((Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                                         * (Altitude_d[lev + 1] - Altitude_d[lev]));
                cpr_tmp[id * nv + lev] = ctmp[id * nv + lev] / btmp[id * nv + lev];
                for (int k = 0; k < 3; k++) {
                    dtmp[id * nv * 3 + lev * 3 + k] = -Mh_d[id * 3 * nv + lev * 3 + k] / time_step;
                    dpr_tmp[id * nv * 3 + lev * 3 + k] =
                        dtmp[id * nv * 3 + lev * 3 + k] / btmp[id * nv + lev];
                }
            }
            else if (lev == bl_top_lev_d[id]) {
                atmp[id * nv + lev] = KM_d[id * (nv + 1) + lev]
                                      / ((Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                                         * (Altitude_d[lev] - Altitude_d[lev - 1]));
                btmp[id * nv + lev] =
                    -(1.0 / (Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                          * (KM_d[id * (nv + 1) + lev + 1] / (Altitude_d[lev + 1] - Altitude_d[lev])
                             + KM_d[id * (nv + 1) + lev] / (Altitude_d[lev] - Altitude_d[lev - 1]))
                      + 1.0 / time_step);
                ctmp[id * nv + lev]    = 0;
                cpr_tmp[id * nv + lev] = 0; //not used, i think
                for (int k = 0; k < 3; k++) {
                    dtmp[id * nv * 3 + lev * 3 + k] =
                        -Mh_d[id * 3 * nv + lev * 3 + k] / time_step
                        - KM_d[id * (nv + 1) + lev + 1] / (Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                              * Mh_d[id * 3 * nv + (lev + 1) * 3 + k]
                              / (Altitude_d[lev + 1] - Altitude_d[lev]);
                    dpr_tmp[id * nv * 3 + lev * 3 + k] =
                        (dtmp[id * nv * 3 + lev * 3 + k]
                         - atmp[id * nv + lev] * dpr_tmp[id * nv * 3 + (lev - 1) * 3 + k])
                        / (btmp[id * nv + lev] - atmp[id * nv + lev] * cpr_tmp[id * nv + lev - 1]);
                }
            }
            else {
                atmp[id * nv + lev] = KM_d[id * (nv + 1) + lev]
                                      / ((Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                                         * (Altitude_d[lev] - Altitude_d[lev - 1]));
                btmp[id * nv + lev] =
                    -(1.0 / (Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                          * (KM_d[id * (nv + 1) + lev + 1] / (Altitude_d[lev + 1] - Altitude_d[lev])
                             + KM_d[id * (nv + 1) + lev] / (Altitude_d[lev] - Altitude_d[lev - 1]))
                      + 1.0 / time_step);
                ctmp[id * nv + lev] = KM_d[id * (nv + 1) + lev + 1]
                                      / ((Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                                         * (Altitude_d[lev + 1] - Altitude_d[lev]));
                cpr_tmp[id * nv + lev] =
                    ctmp[id * nv + lev]
                    / (btmp[id * nv + lev] - atmp[id * nv + lev] * cpr_tmp[id * nv + lev - 1]);
                for (int k = 0; k < 3; k++) {
                    dtmp[id * nv * 3 + lev * 3 + k] = -Mh_d[id * 3 * nv + lev * 3 + k] / time_step;
                    dpr_tmp[id * nv * 3 + lev * 3 + k] =
                        (dtmp[id * nv * 3 + lev * 3 + k]
                         - atmp[id * nv + lev] * dpr_tmp[id * nv * 3 + (lev - 1) * 3 + k])
                        / (btmp[id * nv + lev] - atmp[id * nv + lev] * cpr_tmp[id * nv + lev - 1]);
                }
            }
            if (fabs(btmp[id * nv + lev])
                < (fabs(atmp[id * nv + lev]) + fabs(ctmp[id * nv + lev]))) {
                printf("Warning! Thomas algorithm in boundary layer unstable\n");
            }
        }
        // if (id == 1000) {
        //     printf("stop");
        // }

        for (lev = bl_top_lev_d[id]; lev >= 0; lev--) {
            //backward sweep
            for (int k = 0; k < 3; k++) {
                if (lev == bl_top_lev_d[id]) {
                    Mh_d[id * nv * 3 + lev * 3 + k] = dpr_tmp[id * nv * 3 + lev * 3 + k];
                }
                else {
                    Mh_d[id * nv * 3 + lev * 3 + k] =
                        (dpr_tmp[id * nv * 3 + lev * 3 + k]
                         - cpr_tmp[id * nv + lev] * Mh_d[id * nv * 3 + (lev + 1) * 3 + k]);
                }
            }
        }
        // if (id == 0) {
        //     printf("%f\n", Mh_d[id * nv * 3 + 0]);
        // }
    }
}

__global__ void CalcRiB(double *pressure_d,
                        double *Rho_d,
                        double *Mh_d,
                        double *Tsurface_d,
                        double *Altitude_d,
                        double *Altitudeh_d,
                        double  Rd,
                        double  Cp,
                        double  P_Ref,
                        double  Gravit,
                        double  Ri_crit,
                        double  z_rough,
                        double  z_therm,
                        double *RiB_d,
                        int *   bl_top_lev_d,
                        double *bl_top_height_d,
                        double *CD_d,
                        double *zeta_d,
                        double *vh_lowest_d,
                        int     num,
                        int     nv) {

    // Calculate bulk Richardson number for each level
    // The first value is defined at the midpoint between the lowest layer and the surface
    // The rest are at the interfaces between layers

    int id = blockIdx.x * blockDim.x + threadIdx.x;
    int lev;

    if (id < num) {
        double kappa = Rd / Cp;
        double p_surf, pt_surf, extrap_surf;
        double pt_layer, pt_lowest, pt_layer_below, pt_interface;
        double vh_layer, vh_layer_below, vh_interface;
        double zeta;
        for (lev = 0; lev <= nv; lev++) {
            //first find surface pressure, and calculate pt at the interfaces

            if (lev == 0) {
                //lowest level, RiB defined at midpoint between lowest and surface
                // calculate pot temp of surface
                extrap_surf = -Altitude_d[lev + 1] / (Altitude_d[lev] - Altitude_d[lev + 1]);
                p_surf =
                    pressure_d[id * nv + lev + 1]
                    + extrap_surf * (pressure_d[id * nv + lev] - pressure_d[id * nv + lev + 1]);
                pt_surf = Tsurface_d[id] * pow(p_surf / P_Ref, -kappa);

                // calculate pt and horizontal velocity of layer
                pt_layer = pow(P_Ref, kappa) * pow(pressure_d[id * nv + lev], 1.0 - kappa)
                           / (Rho_d[id * nv + lev] * Rd);
                pt_lowest = pt_layer; //will need this later
                vh_layer  = sqrt((pow(Mh_d[id * nv * 3 + lev * 3 + 0], 2)
                                 + pow(Mh_d[id * nv * 3 + lev * 3 + 1], 2)
                                 + pow(Mh_d[id * nv * 3 + lev * 3 + 2], 2)))
                           / Rho_d[id * nv + lev];

                vh_lowest_d[id] = vh_layer; //velocity of lowest layer

                if (pow(vh_layer, 2) == 0) { //zero velocity, RiB = large +number
                    RiB_d[id * nv + lev] = LARGERiB;
                }
                else { // bulk Richardson number, wrt to surface
                    RiB_d[id * (nv + 1) + lev] = Gravit * Altitude_d[lev] * (pt_layer - pt_surf)
                                                 / (pt_surf * pow(vh_layer, 2));
                }

                // convert to zeta (m-o stability param)
                zeta_d[id * (nv + 1) + lev] = RiB_2_zeta(RiB_d[id * (nv + 1) + lev],
                                                         Ri_crit,
                                                         Altitude_d[lev] / z_rough,
                                                         Altitude_d[lev] / z_therm);
                //surface drag coefficient
                if (RiB_d[id * (nv + 1) + lev] < 0) { //unstable (temporarily model as neutral)
                    CD_d[id] = pow(KVONKARMAN, 2) * pow(log(Altitude_d[lev] / z_rough), -2);
                }
                else if (RiB_d[id * (nv + 1) + lev] == 0) { //neutral
                    CD_d[id] = pow(KVONKARMAN, 2) * pow(log(Altitude_d[lev] / z_rough), -2);
                }
                else if (RiB_d[id * (nv + 1) + lev] <= Ri_crit) { //stable
                    CD_d[id] = pow(KVONKARMAN, 2)
                               * pow(log(Altitude_d[lev] / z_rough) + zeta / Ri_crit, -2);
                }
                else { //super-critical
                    CD_d[id] = 0;
                }
            }
            else if (lev == nv) {
                //what should I do at the top level??
                RiB_d[id * (nv + 1) + lev] = LARGERiB; //top level can't be incorporated into BL?
                // convert to zeta (m-o stability param)
                zeta_d[id * (nv + 1) + lev] = RiB_2_zeta(RiB_d[id * (nv + 1) + lev],
                                                         Ri_crit,
                                                         Altitude_d[lev] / z_rough,
                                                         Altitude_d[lev] / z_therm);
            }
            else {
                //potential temperatures for this layer, layer below, and interface b/w
                pt_layer_below = pt_layer;
                pt_layer       = pow(P_Ref, kappa) * pow(pressure_d[id * nv + lev], 1.0 - kappa)
                           / (Rho_d[id * nv + lev] * Rd);
                pt_interface = pt_layer_below
                               + (pt_layer - pt_layer_below)
                                     * (Altitudeh_d[lev] - Altitude_d[lev - 1])
                                     / (Altitude_d[lev] - Altitude_d[lev - 1]);

                //vh for the layers and interface
                vh_layer_below = vh_layer;
                vh_layer       = sqrt((pow(Mh_d[id * nv * 3 + lev * 3 + 0], 2)
                                 + pow(Mh_d[id * nv * 3 + lev * 3 + 1], 2)
                                 + pow(Mh_d[id * nv * 3 + lev * 3 + 2], 2)))
                           / Rho_d[id * nv + lev];
                vh_interface = vh_layer_below
                               + (vh_layer - vh_layer_below)
                                     * (Altitudeh_d[lev] - Altitude_d[lev - 1])
                                     / (Altitude_d[lev] - Altitude_d[lev - 1]);

                if (pow(vh_interface, 2) == 0) { //zero velocity, set RiB to a big +number
                    RiB_d[id * (nv + 1) + lev] = LARGERiB;
                }
                else { //bulk Ri number, wrt to lowest layer
                    RiB_d[id * (nv + 1) + lev] = Gravit * Altitudeh_d[lev]
                                                 * (pt_interface - pt_lowest)
                                                 / (pt_lowest * pow(vh_interface, 2));
                }
                // convert to zeta (m-o stability param)
                zeta_d[id * (nv + 1) + lev] = RiB_2_zeta(RiB_d[id * (nv + 1) + lev],
                                                         Ri_crit,
                                                         Altitude_d[lev] / z_rough,
                                                         Altitude_d[lev] / z_therm);
            }

            if (RiB_d[id * (nv + 1) + lev] <= Ri_crit) {
                // this layer is in the BL, if layers below are also
                if (lev == 0) {
                    bl_top_lev_d[id] = lev;
                }
                else if (bl_top_lev_d[id] == lev - 1) {
                    //make sure layer below was also in BL so it is all connected
                    bl_top_lev_d[id] = lev;
                }
            }
        }
        //take the bl height to be the height of the interface above the top bl layer
        bl_top_height_d[id] = Altitudeh_d[bl_top_lev_d[id] + 1];
        // printf("bl top = %f\n", bl_top_height_d[id]);
    }
}

__device__ double RiB_2_zeta(double RiB, double Ri_crit, double z_z0, double z_zT) {
    if (RiB <= 0) {
        return 0.0; // for now, treating unstable part as neutral
    }
    else {
        double beta = 1.0 / Ri_crit;
        return (-(2 * RiB * beta * log(z_z0) - log(z_zT))
                - sqrt(4 * RiB * beta * pow(log(z_z0), 2) - 4 * RiB * beta * log(z_z0) * log(z_zT)
                       + pow(log(z_zT), 2)))
               / (2 * (RiB * pow(beta, 2) - beta));
    }
}

__global__ void CalcKM(double *RiB_d,
                       double *zeta_d,
                       double *CD_d,
                       double *bl_top_height_d,
                       int *   bl_top_lev_d,
                       double *vh_lowest_d,
                       double *Altitude_d,
                       double *Altitudeh_d,
                       double  Ri_crit,
                       double  z_rough,
                       double  z_therm,
                       double  f_surf_layer,
                       double *KM_d,
                       int     num) {
    int id  = blockIdx.x * blockDim.x + threadIdx.x;
    int nvi = gridDim.y;
    int lev = blockIdx.y;

    if (id < num) {
        if (lev == 0) { // coefficient at surface
            KM_d[id * (nvi) + lev] = CD_d[id] * vh_lowest_d[id] * Altitude_d[lev];
        }
        else {
            // double zeta = RiB_2_zeta(RiB_d[id * (nvi) + lev],
            //                          Ri_crit,
            //                          Altitude_d[lev] / z_rough,
            //                          Altitude_d[lev] / z_therm);
            if (Altitude_d[lev] <= f_surf_layer * bl_top_height_d[id]) { //inner layer
                //condition on surface layer stability as in Frierson
                if (RiB_d[id * (nvi) + 0] < 0) { //unstable
                    KM_d[id * (nvi) + lev] =
                        KVONKARMAN * Altitudeh_d[lev] * sqrt(CD_d[id]) * vh_lowest_d[id];
                }
                else if (RiB_d[id * (nvi) + 0] == 0) { //neutral
                    KM_d[id * (nvi) + lev] =
                        KVONKARMAN * Altitudeh_d[lev] * sqrt(CD_d[id]) * vh_lowest_d[id];
                }
                else { //stable
                    KM_d[id * (nvi) + lev] = KVONKARMAN * Altitudeh_d[lev] * sqrt(CD_d[id])
                                             * vh_lowest_d[id]
                                             / (1 + zeta_d[id * (nvi) + lev] / Ri_crit);
                }
            }
            else if (Altitude_d[lev] <= bl_top_height_d[id]) { //outer layer
                double z_scale = Altitudeh_d[lev] / (f_surf_layer * bl_top_height_d[id])
                                 * pow(1
                                           - (Altitudeh_d[lev] - f_surf_layer * bl_top_height_d[id])
                                                 / ((1 - f_surf_layer) * bl_top_height_d[id]),
                                       2);
                if (RiB_d[id * (nvi) + 0] < 0) { //unstable
                    KM_d[id * (nvi) + lev] = KVONKARMAN * f_surf_layer * bl_top_height_d[id]
                                             * sqrt(CD_d[id]) * vh_lowest_d[id] * z_scale;
                }
                else if (RiB_d[id * (nvi) + 0] == 0) { //neutral
                    KM_d[id * (nvi) + lev] = KVONKARMAN * f_surf_layer * bl_top_height_d[id]
                                             * sqrt(CD_d[id]) * vh_lowest_d[id] * z_scale;
                }
                else { //stable  (zeta at f*h)
                    KM_d[id * (nvi) + lev] =
                        KVONKARMAN * f_surf_layer * bl_top_height_d[id] * sqrt(CD_d[id])
                        * vh_lowest_d[id] * z_scale
                        / (1 + zeta_d[id * (nvi) + bl_top_lev_d[id] + 1] / Ri_crit);
                }
            }
            else {
                KM_d[id * (nvi) + lev] = 0.0;
            }
        }
    }
}
