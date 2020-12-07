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
    log::printf("    surf_drag (RayleighHS)     = %e 1/s\n", surf_drag_config);
    log::printf("    bl_sigma (RayleighHS)      = %f \n", bl_sigma_config);
    log::printf("    z_rough (LocalMixL)               = %f m\n", z_rough_config);
    log::printf("    asl_transition_height (LocalMixL) = %f m\n", asl_transition_height_config);
    log::printf("    abl_asym_len (LocalMixL)          = %f m\n", abl_asym_len_config);
    log::printf("    free_asym_len (LocalMixL)         = %f m\n", abl_asym_len_config);

    log::printf("\n");
}

bool boundary_layer::initialise_memory(const ESP &              esp,
                                       device_RK_array_manager &phy_modules_core_arrays) {

    cudaMalloc((void **)&cpr_tmp, esp.nvi * esp.point_num * sizeof(double));
    cudaMalloc((void **)&dpr_tmp, 3 * esp.nvi * esp.point_num * sizeof(double));
    cudaMalloc((void **)&KM_d, esp.nvi * esp.point_num * sizeof(double));
    cudaMalloc((void **)&KH_d, esp.nvi * esp.point_num * sizeof(double));
    cudaMalloc((void **)&CM_d, esp.point_num * sizeof(double));
    cudaMalloc((void **)&CH_d, esp.point_num * sizeof(double));
    cudaMalloc((void **)&vh_lowest_d, esp.point_num * sizeof(double));
    cudaMalloc((void **)&pt_surf_d, esp.point_num * sizeof(double));
    cudaMalloc((void **)&F_sens_d, esp.point_num * sizeof(double));

    cudaMalloc((void **)&RiGrad_d, esp.nvi * esp.point_num * sizeof(double));

    cudaMalloc((void **)&bl_top_lev_d, esp.point_num * sizeof(int));
    cudaMalloc((void **)&bl_top_height_d, esp.point_num * sizeof(double));

    cudaMalloc((void **)&Rho_int_d, esp.nvi * esp.point_num * sizeof(double));
    cudaMalloc((void **)&p_int_d, esp.nvi * esp.point_num * sizeof(double));

    KM_h            = (double *)malloc(esp.nvi * esp.point_num * sizeof(double));
    KH_h            = (double *)malloc(esp.nvi * esp.point_num * sizeof(double));
    bl_top_lev_h    = (int *)malloc(esp.point_num * sizeof(int));
    bl_top_height_h = (double *)malloc(esp.point_num * sizeof(double));
    CM_h            = (double *)malloc(esp.point_num * sizeof(double));
    CH_h            = (double *)malloc(esp.point_num * sizeof(double));
    F_sens_h        = (double *)malloc(esp.point_num * sizeof(double));

    RiGrad_h = (double *)malloc(esp.nvi * esp.point_num * sizeof(double));

    cudaMemset(cpr_tmp, 0, sizeof(double) * esp.point_num * esp.nvi);
    cudaMemset(dpr_tmp, 0, sizeof(double) * 3 * esp.point_num * esp.nvi);
    cudaMemset(F_sens_d, 0, sizeof(double) * esp.point_num);
    cudaMemset(KM_d, 0, sizeof(double) * esp.point_num * esp.nvi);
    cudaMemset(KH_d, 0, sizeof(double) * esp.point_num * esp.nvi);
    cudaMemset(bl_top_lev_d, 0, sizeof(int) * esp.point_num);
    cudaMemset(bl_top_height_d, 0, sizeof(double) * esp.point_num);
    cudaMemset(CM_d, 0, sizeof(double) * esp.point_num);
    cudaMemset(CH_d, 0, sizeof(double) * esp.point_num);
    cudaMemset(vh_lowest_d, 0, sizeof(double) * esp.point_num);
    cudaMemset(pt_surf_d, 0, sizeof(double) * esp.point_num);

    cudaMemset(Rho_int_d, 0, sizeof(double) * esp.point_num * esp.nvi);
    cudaMemset(p_int_d, 0, sizeof(double) * esp.point_num * esp.nvi);

    cudaMemset(RiGrad_d, 0, sizeof(double) * esp.nvi * esp.point_num);

    return true;
}


bool boundary_layer::free_memory() {
    free(KM_h);
    free(KH_h);
    free(bl_top_lev_h);
    free(bl_top_height_h);

    free(CM_h);
    free(CH_h);
    free(RiGrad_h);
    free(F_sens_h);

    cudaFree(cpr_tmp);
    cudaFree(dpr_tmp);
    cudaFree(KM_d);
    cudaFree(KH_d);
    cudaFree(CM_d);
    cudaFree(CH_d);
    cudaFree(vh_lowest_d);
    cudaFree(pt_surf_d);
    cudaFree(bl_top_lev_d);
    cudaFree(F_sens_d);
    cudaFree(bl_top_height_d);
    cudaFree(Rho_int_d);
    cudaFree(p_int_d);
    cudaFree(RiGrad_d);

    return true;
}

bool boundary_layer::initial_conditions(const ESP &esp, const SimulationSetup &sim, storage *s) {
    bool config_OK = true;

    bl_type = RAYLEIGHHS;
    if (bl_type_str == "RayleighHS") {
        bl_type = RAYLEIGHHS;
        config_OK &= true;
    }
    else if (bl_type_str == "LocalMixL" || bl_type_str == "Local") {
        bl_type = LOCALMIXL;
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
            asl_transition_height_config,
            abl_asym_len_config,
            free_asym_len_config);

    return true;
}

bool boundary_layer::phy_loop(ESP &                  esp,
                              const SimulationSetup &sim,
                              kernel_diagnostics &   diag,
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
    else if (bl_type == LOCALMIXL) {
        CalcGradRi<<<NBi, NTH>>>(esp.pressure_d,
                                 esp.Rho_d,
                                 esp.Mh_d,
                                 esp.Tsurface_d,
                                 esp.pt_d,
                                 esp.Altitude_d,
                                 esp.Altitudeh_d,
                                 sim.Rd,
                                 sim.Cp,
                                 sim.P_Ref,
                                 sim.Gravit,
                                 Ri_crit,
                                 z_rough,
                                 RiGrad_d,
                                 pt_surf_d,
                                 p_int_d,
                                 CM_d,
                                 CH_d,
                                 Rho_int_d,
                                 KM_d,
                                 KH_d,
                                 F_sens_d,
                                 asl_transition_height,
                                 abl_asym_len,
                                 free_asym_len,
                                 esp.point_num,
                                 esp.nv);

        cudaDeviceSynchronize();

        cudaMemset(cpr_tmp, 0, sizeof(double) * esp.point_num * esp.nvi);
        cudaMemset(dpr_tmp, 0, sizeof(double) * 3 * esp.point_num * esp.nvi);

#if defined(DIAG_CHECK_BL_THOMAS_DIAG_DOM)
        //checking diagonal dominance in thomas algorithm
        diag.reset_flag();
#endif

        Momentum_Diff_Impl<<<NBLEV, NTH>>>(esp.Mh_d,
                                           esp.pressure_d,
                                           esp.Rho_d,
                                           esp.Altitude_d,
                                           esp.Altitudeh_d,
                                           cpr_tmp,
                                           dpr_tmp,
                                           KM_d,
                                           Rho_int_d,
                                           time_step,
                                           sim.A,
                                           esp.point_num,
                                           esp.nv,
                                           bl_top_lev_d, //doesn't mean anything in this case
                                           sim.DeepModel,
                                           *diag.diagnostics_global_flag,
                                           *diag.diagnostics);

        cudaDeviceSynchronize();

#if defined(DIAG_CHECK_BL_THOMAS_DIAG_DOM)
        if (diag.check_flag()) {
            diag.dump_data("bl_thomas_momentum", nstep, 0, 0, esp, esp.point_num, esp.nv);
            unsigned int flag = diag.get_flag();
            if ((flag & BL_THOMAS_NOT_DD) == BL_THOMAS_NOT_DD)
                log::printf("Non diagonally dominant matrix for thomas algorithm in boundary layer "
                            "momentum\n");
        }
#endif

        cudaMemset(cpr_tmp, 0, sizeof(double) * esp.point_num * esp.nvi);
        cudaMemset(dpr_tmp, 0, sizeof(double) * 3 * esp.point_num * esp.nvi);

#if defined(DIAG_CHECK_BL_THOMAS_DIAG_DOM)
        //checking diagonal dominance in thomas algorithm
        diag.reset_flag();
#endif

        Heat_Diff_Impl_EnergyEq<<<NBLEV, NTH>>>(esp.pt_d,
                                                esp.pressure_d,
                                                esp.temperature_d,
                                                esp.Rho_d,
                                                esp.Altitude_d,
                                                esp.Altitudeh_d,
                                                esp.Tsurface_d,
                                                cpr_tmp,
                                                dpr_tmp,
                                                KH_d,
                                                Rho_int_d,
                                                pt_surf_d,
                                                p_int_d,
                                                F_sens_d,
                                                time_step,
                                                sim.Rd,
                                                sim.Cp,
                                                sim.P_Ref,
                                                esp.Csurf,
                                                sim.A,
                                                esp.point_num,
                                                esp.nv,
                                                bl_top_lev_d,
                                                sim.DeepModel,
                                                *diag.diagnostics_global_flag,
                                                *diag.diagnostics);

        cudaDeviceSynchronize();

#if defined(DIAG_CHECK_BL_THOMAS_DIAG_DOM)
        if (diag.check_flag()) {
            diag.dump_data("bl_thomas_heateqn", nstep, 0, 0, esp, esp.point_num, esp.nv);
            unsigned int flag = diag.get_flag();
            if ((flag & BL_THOMAS_NOT_DD) == BL_THOMAS_NOT_DD)
                log::printf("Non diagonally dominant matrix for thomas algorithm in boundary layer "
                            "heat eqn\n");
        }
#endif
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

    // roughness length scaling
    config_reader.append_config_var("z_rough", z_rough_config, z_rough_config);

    return true;
}

bool boundary_layer::store(const ESP &esp, storage &s) {
    if (bl_type == LOCALMIXL) {
        cudaMemcpy(
            RiGrad_h, RiGrad_d, esp.nvi * esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
        s.append_table(
            RiGrad_h, esp.nvi * esp.point_num, "/RiGrad", " ", "gradient Richardson number");

        cudaMemcpy(CM_h, CM_d, esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
        s.append_table(CM_h, esp.point_num, "/CM", " ", "surface drag coefficient");

        cudaMemcpy(CH_h, CH_d, esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
        s.append_table(CH_h, esp.point_num, "/CH", " ", "surface heat-transfer coefficient");

        cudaMemcpy(KM_h, KM_d, esp.point_num * esp.nvi * sizeof(double), cudaMemcpyDeviceToHost);
        s.append_table(KM_h,
                       esp.point_num * esp.nvi,
                       "/KM",
                       "m^2/s",
                       "turbulent diffusion coefficient for momentum");

        cudaMemcpy(KH_h, KH_d, esp.point_num * esp.nvi * sizeof(double), cudaMemcpyDeviceToHost);
        s.append_table(KH_h,
                       esp.point_num * esp.nvi,
                       "/KH",
                       "m^2/s",
                       "turbulent diffusion coefficient for heat");

        cudaMemcpy(F_sens_h, F_sens_d, esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
        s.append_table(
            F_sens_h, esp.point_num, "/F_sens", "W/m^2", "Surface to atmos sensible heat flux");
    }

    return true;
}

bool boundary_layer::store_init(storage &s) {
    s.append_value(surf_drag, "/surf_drag", "1/s", "surface drag coefficient");
    s.append_value(bl_sigma, "/bl_sigma", " ", "boundary layer sigma coordinate");

    s.append_value(Ri_crit, "/Ri_crit", " ", "critical Richardson number");

    s.append_value(z_rough, "/z_rough", "m", "boundary layer roughness length");

    //      bl type  option
    s.append_value(int(bl_type), "/bl_type", "-", "Type of boundary layer physics");

    return true;
}

void boundary_layer::BLSetup(const ESP &            esp,
                             const SimulationSetup &sim,
                             int                    bl_type,
                             double                 surf_drag_,
                             double                 bl_sigma_,
                             double                 Ri_crit_,
                             double                 z_rough_,
                             double                 asl_transition_height_,
                             double                 abl_asym_len_,
                             double                 free_asym_len_) {
    if (bl_type == RAYLEIGHHS) {
        surf_drag = surf_drag_;
        bl_sigma  = bl_sigma_;
    }
    else if (bl_type == LOCALMIXL) {
        Ri_crit               = Ri_crit_;
        z_rough               = z_rough_;
        asl_transition_height = asl_transition_height_;
        abl_asym_len          = abl_asym_len_;
        free_asym_len         = free_asym_len_;
        for (int id = 0; id < esp.point_num; id++) {
            bl_top_lev_h[id] = esp.nv - 1; //local formulation includes entire column
        }
        cudaMemcpy(bl_top_lev_d, bl_top_lev_h, esp.point_num * sizeof(int), cudaMemcpyHostToDevice);
    }
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

__global__ void
Momentum_Diff_Impl(double *      Mh_d,
                   double *      pressure_d,
                   double *      Rho_d,
                   double *      Altitude_d,
                   double *      Altitudeh_d,
                   double *      cpr_tmp,
                   double *      dpr_tmp,
                   double *      KM_d,
                   double *      Rho_int_d,
                   double        time_step,
                   double        A,
                   int           num,
                   int           nv,
                   int *         bl_top_lev_d,
                   bool          DeepModel,
                   unsigned int *diagnostics_flag,
                   diag_data *   diagnostics_data) { //need DeepModel and planet radius here

    //should create check on stability of thomas algorithm

    int id = blockIdx.x * blockDim.x + threadIdx.x;
    int lev;

    if (id < num) {
        double rup, rlow; //radii for spherical curvature calculation
        double atmp, btmp, ctmp, dtmp;
        for (lev = 0; lev < bl_top_lev_d[id] + 1; lev++) {
            //forward sweep
            //r = Altitude_d[lev] + A;
            if (DeepModel) {
                rup  = (Altitudeh_d[lev + 1] + A) / (Altitude_d[lev] + A);
                rlow = (Altitudeh_d[lev] + A) / (Altitude_d[lev] + A);
            }
            else {
                rup  = 1.0;
                rlow = 1.0;
            }
            if (lev == 0) { //lowest layer, v at lowest boundary = 0, dz0 = Altitude0
                atmp = 0;
                btmp = -(1.0 / (Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                             * (pow(rup, 2) * Rho_int_d[id * (nv + 1) + lev + 1]
                                    * KM_d[id * (nv + 1) + lev + 1]
                                    / (Altitude_d[lev + 1] - Altitude_d[lev])
                                + pow(rlow, 2) * Rho_int_d[id * (nv + 1) + lev]
                                      * KM_d[id * (nv + 1) + lev] / Altitude_d[lev])
                         + Rho_d[id * nv + lev] / time_step);
                ctmp = pow(rup, 2) * Rho_int_d[id * (nv + 1) + lev + 1]
                       * KM_d[id * (nv + 1) + lev + 1]
                       / ((Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                          * (Altitude_d[lev + 1] - Altitude_d[lev]));
                cpr_tmp[id * (nv + 1) + lev] = ctmp / btmp;
                for (int k = 0; k < 3; k++) {
                    dtmp = -Mh_d[id * 3 * nv + lev * 3 + k] / time_step;
                    dpr_tmp[id * (nv + 1) * 3 + lev * 3 + k] = dtmp / btmp;
                }
            }
            else if (lev == bl_top_lev_d[id]) {
                atmp = pow(rlow, 2) * Rho_int_d[id * (nv + 1) + lev] * KM_d[id * (nv + 1) + lev]
                       / ((Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                          * (Altitude_d[lev] - Altitude_d[lev - 1]));
                btmp =
                    -(1.0 / (Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                          * ((pow(rlow, 2) * Rho_int_d[id * (nv + 1) + lev])
                             * KM_d[id * (nv + 1) + lev] / (Altitude_d[lev] - Altitude_d[lev - 1]))
                      + Rho_d[id * nv + lev] / time_step);
                ctmp                         = 0;
                cpr_tmp[id * (nv + 1) + lev] = 0; //not used, i think
                for (int k = 0; k < 3; k++) {
                    dtmp = -Mh_d[id * 3 * nv + lev * 3 + k] / time_step;
                    dpr_tmp[id * (nv + 1) * 3 + lev * 3 + k] =
                        (dtmp - atmp * dpr_tmp[id * (nv + 1) * 3 + (lev - 1) * 3 + k])
                        / (btmp - atmp * cpr_tmp[id * (nv + 1) + lev - 1]);
                }
            }
            else {
                atmp = pow(rlow, 2) * Rho_int_d[id * (nv + 1) + lev] * KM_d[id * (nv + 1) + lev]
                       / ((Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                          * (Altitude_d[lev] - Altitude_d[lev - 1]));
                btmp = -(1.0 / (Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                             * (pow(rup, 2) * Rho_int_d[id * (nv + 1) + lev + 1]
                                    * KM_d[id * (nv + 1) + lev + 1]
                                    / (Altitude_d[lev + 1] - Altitude_d[lev])
                                + pow(rlow, 2) * Rho_int_d[id * (nv + 1) + lev]
                                      * KM_d[id * (nv + 1) + lev]
                                      / (Altitude_d[lev] - Altitude_d[lev - 1]))
                         + Rho_d[id * nv + lev] / time_step);
                ctmp = pow(rup, 2) * Rho_int_d[id * (nv + 1) + lev + 1]
                       * KM_d[id * (nv + 1) + lev + 1]
                       / ((Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                          * (Altitude_d[lev + 1] - Altitude_d[lev]));
                cpr_tmp[id * (nv + 1) + lev] =
                    ctmp / (btmp - atmp * cpr_tmp[id * (nv + 1) + lev - 1]);
                for (int k = 0; k < 3; k++) {
                    dtmp = -Mh_d[id * 3 * nv + lev * 3 + k] / time_step;
                    dpr_tmp[id * (nv + 1) * 3 + lev * 3 + k] =
                        (dtmp - atmp * dpr_tmp[id * (nv + 1) * 3 + (lev - 1) * 3 + k])
                        / (btmp - atmp * cpr_tmp[id * (nv + 1) + lev - 1]);
                }
            }

#ifdef DIAG_CHECK_BL_THOMAS_DIAG_DOM
            //reset values
            diagnostics_data[id * nv + lev].flag = 0;
            diagnostics_data[id * nv + lev].data = make_double4(0.0, 0.0, 0.0, 0.0);
            //check if matrix diagonally dominant
            if (lev < bl_top_lev_d[id] + 1) {
                if (!(fabs(btmp) > THOMAS_DIAG_DOM_FACTOR * (fabs(atmp) + fabs(ctmp)))) {
                    atomicOr(diagnostics_flag, BL_THOMAS_NOT_DD);
                    diagnostics_data[id * nv + lev].flag   = BL_THOMAS_NOT_DD;
                    diagnostics_data[id * nv + lev].data.x = atmp;
                    diagnostics_data[id * nv + lev].data.y = btmp;
                    diagnostics_data[id * nv + lev].data.z = ctmp;
                    // printf("Warning! Thomas algorithm in boundary layer mom. equation unstable\n");
                }
            }
#endif
        }

        for (lev = bl_top_lev_d[id]; lev >= 0; lev--) {
            //backward sweep
            for (int k = 0; k < 3; k++) {
                if (lev == bl_top_lev_d[id]) {
                    Mh_d[id * nv * 3 + lev * 3 + k] =
                        dpr_tmp[id * (nv + 1) * 3 + lev * 3 + k] * Rho_d[id * nv + lev];
                }
                else {
                    Mh_d[id * nv * 3 + lev * 3 + k] =
                        (dpr_tmp[id * (nv + 1) * 3 + lev * 3 + k]
                         - cpr_tmp[id * (nv + 1) + lev] * Mh_d[id * nv * 3 + (lev + 1) * 3 + k]
                               / Rho_d[id * nv + lev + 1])
                        * Rho_d[id * nv + lev];
                }
            }
        }
    }
}

__global__ void Heat_Diff_Impl_EnergyEq(double *      pt_d,
                                        double *      pressure_d,
                                        double *      temperature_d,
                                        double *      Rho_d,
                                        double *      Altitude_d,
                                        double *      Altitudeh_d,
                                        double *      Tsurface_d,
                                        double *      cpr_tmp,
                                        double *      dpr_tmp,
                                        double *      KH_d,
                                        double *      Rho_int_d,
                                        double *      pt_surf_d,
                                        double *      p_int_d,
                                        double *      F_sens_d,
                                        double        time_step,
                                        double        Rd,
                                        double        Cp,
                                        double        P_Ref,
                                        double        Csurf,
                                        double        A,
                                        int           num,
                                        int           nv,
                                        int *         bl_top_lev_d,
                                        bool          DeepModel,
                                        unsigned int *diagnostics_flag,
                                        diag_data *   diagnostics_data) {

    //should create check on stability of thomas algorithm

    int id = blockIdx.x * blockDim.x + threadIdx.x;
    int lev;

    if (id < num) {
        double Cv    = Cp - Rd;
        double kappa = Rd / Cp;
        double rup, rlow;
        double atmp, btmp, ctmp, dtmp;
        for (lev = -1; lev < bl_top_lev_d[id] + 1; lev++) {
            //forward sweep
            if (DeepModel) {
                if (lev == -1) {
                    rup  = 1.0;
                    rlow = 1.0;
                }
                else {
                    rup  = (Altitudeh_d[lev + 1] + A) / (Altitude_d[lev] + A);
                    rlow = (Altitudeh_d[lev] + A) / (Altitude_d[lev] + A);
                }
            }
            else {
                rup  = 1.0;
                rlow = 1.0;
            }
            if (lev == -1) { //treat surface in matrix solution (need to offset abcd arrays by +1)
                atmp = 0;
                btmp = -(Cp / Csurf * Rho_int_d[id * (nv + 1)] * KH_d[id * (nv + 1) + lev + 1]
                             / Altitude_d[lev + 1]
                         + 1.0 / time_step);
                ctmp = Cp / Csurf * Rho_int_d[id * (nv + 1)] * KH_d[id * (nv + 1) + lev + 1]
                       / Altitude_d[lev + 1];
                cpr_tmp[id * (nv + 1) + lev + 1] = ctmp / btmp;
                dtmp                             = -pt_surf_d[id] / time_step;
                dpr_tmp[id * (nv + 1) + lev + 1] = dtmp / btmp;
            }
            else if (lev == 0) { //lowest layer, v at lowest boundary = 0, dz0 = Altitude0
                atmp = Cp / Cv * pow(rlow, 2) * Rho_int_d[id * (nv + 1) + lev]
                       * KH_d[id * (nv + 1) + lev] / (Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                       / Altitude_d[lev] * pow(p_int_d[id * (nv + 1) + lev] / P_Ref, kappa);
                ctmp = Cp / Cv * pow(rup, 2) * Rho_int_d[id * (nv + 1) + lev + 1]
                       * KH_d[id * (nv + 1) + lev + 1]
                       / ((Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                          * (Altitude_d[lev + 1] - Altitude_d[lev]))
                       * pow(p_int_d[id * (nv + 1) + lev + 1] / P_Ref, kappa);
                btmp = -(atmp + ctmp
                         + Rho_d[id * nv + lev] / time_step
                               * pow(pressure_d[id * nv + lev] / P_Ref, kappa));
                cpr_tmp[id * (nv + 1) + lev + 1] =
                    ctmp / (btmp - atmp * cpr_tmp[id * (nv + 1) + lev]);
                dtmp = -Rho_d[id * nv + lev] * pt_d[id * nv + lev] / time_step
                       * pow(pressure_d[id * nv + lev] / P_Ref, kappa);
                dpr_tmp[id * (nv + 1) + lev + 1] = (dtmp - atmp * dpr_tmp[id * (nv + 1) + lev])
                                                   / (btmp - atmp * cpr_tmp[id * (nv + 1) + lev]);
            }
            else if (lev == bl_top_lev_d[id]) {
                atmp = Cp / Cv * pow(rlow, 2) * Rho_int_d[id * (nv + 1) + lev]
                       * KH_d[id * (nv + 1) + lev]
                       / ((Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                          * (Altitude_d[lev] - Altitude_d[lev - 1]))
                       * pow(p_int_d[id * (nv + 1) + lev] / P_Ref, kappa);
                btmp                             = -(atmp
                         + Rho_d[id * nv + lev] / time_step
                               * pow(pressure_d[id * nv + lev] / P_Ref, kappa));
                ctmp                             = 0;
                cpr_tmp[id * (nv + 1) + lev + 1] = 0; //not used, i think
                dtmp = -Rho_d[id * nv + lev] * pt_d[id * nv + lev] / time_step
                       * pow(pressure_d[id * nv + lev] / P_Ref, kappa);
                dpr_tmp[id * (nv + 1) + lev + 1] = (dtmp - atmp * dpr_tmp[id * (nv + 1) + (lev)])
                                                   / (btmp - atmp * cpr_tmp[id * (nv + 1) + lev]);
            }
            else {
                atmp = Cp / Cv * pow(rlow, 2) * Rho_int_d[id * (nv + 1) + lev]
                       * KH_d[id * (nv + 1) + lev]
                       / ((Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                          * (Altitude_d[lev] - Altitude_d[lev - 1]))
                       * pow(p_int_d[id * (nv + 1) + lev] / P_Ref, kappa);
                ctmp = Cp / Cv * pow(rup, 2) * Rho_int_d[id * (nv + 1) + lev + 1]
                       * KH_d[id * (nv + 1) + lev + 1]
                       / ((Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                          * (Altitude_d[lev + 1] - Altitude_d[lev]))
                       * pow(p_int_d[id * (nv + 1) + lev + 1] / P_Ref, kappa);
                btmp = -(atmp + ctmp
                         + Rho_d[id * nv + lev] / time_step
                               * pow(pressure_d[id * nv + lev] / P_Ref, kappa));
                cpr_tmp[id * (nv + 1) + lev + 1] =
                    ctmp / (btmp - atmp * cpr_tmp[id * (nv + 1) + lev]);
                dtmp = -Rho_d[id * nv + lev] * pt_d[id * nv + lev] / time_step
                       * pow(pressure_d[id * nv + lev] / P_Ref, kappa);
                dpr_tmp[id * (nv + 1) + lev + 1] = (dtmp - atmp * dpr_tmp[id * (nv + 1) + (lev)])
                                                   / (btmp - atmp * cpr_tmp[id * (nv + 1) + lev]);
            }
#ifdef DIAG_CHECK_BL_THOMAS_DIAG_DOM
            //reset values
            diagnostics_data[id * nv + lev + 1].flag = 0;
            diagnostics_data[id * nv + lev + 1].data = make_double4(0.0, 0.0, 0.0, 0.0);
            //check if matrix diagonally dominant
            if (lev < bl_top_lev_d[id] + 1) {
                if (!(fabs(btmp) > THOMAS_DIAG_DOM_FACTOR * (fabs(atmp) + fabs(ctmp)))) {
                    atomicOr(diagnostics_flag, BL_THOMAS_NOT_DD);
                    diagnostics_data[id * nv + lev + 1].flag   = BL_THOMAS_NOT_DD;
                    diagnostics_data[id * nv + lev + 1].data.x = atmp;
                    diagnostics_data[id * nv + lev + 1].data.y = btmp;
                    diagnostics_data[id * nv + lev + 1].data.z = ctmp;
                }
            }
#endif
        }
        for (lev = bl_top_lev_d[id]; lev >= 0; lev--) {
            //backward sweep
            if (lev == bl_top_lev_d[id]) {
                pt_d[id * nv + lev] = dpr_tmp[id * (nv + 1) + lev + 1];
            }
            else {
                pt_d[id * nv + lev] =
                    (dpr_tmp[id * (nv + 1) + lev + 1]
                     - cpr_tmp[id * (nv + 1) + lev + 1] * pt_d[id * nv + (lev + 1)]);
            }
            temperature_d[id * nv + lev] =
                pt_d[id * nv + lev] * pow(pressure_d[id * nv + lev] / P_Ref, kappa);
            //now we can calculate new pressure
            pressure_d[id * nv + lev] = temperature_d[id * nv + lev] * (Rd * Rho_d[id * nv + lev]);
        }

        pt_surf_d[id] =
            (dpr_tmp[id * (nv + 1) + 0] - cpr_tmp[id * (nv + 1) + 0] * pt_d[id * nv + 0]);
        Tsurface_d[id] = pt_surf_d[id] * pow(p_int_d[id * (nv + 1) + 0] / P_Ref, kappa);
        F_sens_d[id] *= (pt_surf_d[id] - pt_d[id * nv + 0]); //finish calculating sensible heat flux
    }
}

//stability/scaling functions for local scheme
__device__ double stability_fM_u(double CN, double Ri, double z, double z_rough) {
    return 1.0 - 10.0 * Ri / (1 + 75.0 * CN * sqrt((z + z_rough) / z_rough * fabs(Ri)));
}

__device__ double stability_fH_u(double CN, double Ri, double z, double z_rough) {
    return 1.0 - 15.0 * Ri / (1 + 75.0 * CN * sqrt((z + z_rough) / z_rough * fabs(Ri)));
}

__device__ double stability_f_s(double Ri) {
    return 1.0 / (1.0 + 10.0 * Ri * (1 + 8.0 * Ri));
}


__global__ void CalcGradRi(double *pressure_d,
                           double *Rho_d,
                           double *Mh_d,
                           double *Tsurface_d,
                           double *pt_d,
                           double *Altitude_d,
                           double *Altitudeh_d,
                           double  Rd,
                           double  Cp,
                           double  P_Ref,
                           double  Gravit,
                           double  Ri_crit,
                           double  z_rough,
                           double *RiGrad_d,
                           double *pt_surf_d,
                           double *p_int_d,
                           double *CM_d,
                           double *CH_d,
                           double *Rho_int_d,
                           double *KM_d,
                           double *KH_d,
                           double *F_sens_d,
                           double  asl_transition_height,
                           double  abl_asym_len,
                           double  free_asym_len,
                           int     num,
                           int     nv) {

    //gradient richardson number at interfaces
    //also calculate KH, KM, CM, CH directly here
    //need to add sensible heat flux

    int id  = blockIdx.x * blockDim.x + threadIdx.x;
    int nvi = gridDim.y;
    int lev = blockIdx.y;

    if (id < num) {
        double kappa = Rd / Cp;
        double pt_below, pt_int;
        double shear2, mix_length_m, mix_length_h, asym_len_scale_m, asym_len_scale_h, e_mix_m,
            e_mix_h;

        if (lev == 0) {
            //surface drag coeffs, etc
            double extrap_surf, CN;
            extrap_surf = -Altitude_d[lev + 1] / (Altitude_d[lev] - Altitude_d[lev + 1]);
            p_int_d[id * nvi + lev] =
                pressure_d[id * nv + lev + 1]
                + extrap_surf * (pressure_d[id * nv + lev] - pressure_d[id * nv + lev + 1]);
            pt_surf_d[id] = Tsurface_d[id] * pow(p_int_d[id * nvi + lev] / P_Ref, -kappa);
            Rho_int_d[id * nvi + lev] =
                Rho_d[id * nv + lev + 1]
                + extrap_surf * (Rho_d[id * nv + lev] - Rho_d[id * nv + lev + 1]);
            pt_d[id * nv + lev] = pow(P_Ref, kappa) * pow(pressure_d[id * nv + lev], 1.0 - kappa)
                                  / (Rd * Rho_d[id * nv + lev]);
            shear2 =
                pow((Mh_d[id * nv * 3 + lev * 3 + 0]) / Rho_d[id * nv + lev] / (Altitude_d[lev]), 2)
                + pow((Mh_d[id * nv * 3 + lev * 3 + 1]) / Rho_d[id * nv + lev] / (Altitude_d[lev]),
                      2)
                + pow((Mh_d[id * nv * 3 + lev * 3 + 2]) / Rho_d[id * nv + lev] / (Altitude_d[lev]),
                      2);
            //gradient Ri number
            if (shear2 == 0) {
                RiGrad_d[id * nvi + lev] = LARGERiB;
            }
            else {
                RiGrad_d[id * nvi + lev] = Gravit / pt_d[id * nv + lev]
                                           * (pt_d[id * nv + lev] - pt_surf_d[id])
                                           / (Altitude_d[lev]) / shear2;
            }
            CN = pow(KVONKARMAN / log((Altitude_d[lev] + z_rough) / z_rough), 2);
            if (RiGrad_d[id * nvi + lev] < 0) {
                CM_d[id] =
                    CN * stability_fM_u(CN, RiGrad_d[id * nvi + lev], Altitude_d[lev], z_rough);
                CH_d[id] =
                    CN * stability_fH_u(CN, RiGrad_d[id * nvi + lev], Altitude_d[lev], z_rough);
            }
            else {
                CM_d[id] = CN * stability_f_s(RiGrad_d[id * nvi + lev]);
                CH_d[id] = CN * stability_f_s(RiGrad_d[id * nvi + lev]);
            }
            KM_d[id * nvi + lev] = CM_d[id] * sqrt(shear2) * pow(Altitude_d[lev], 2);
            KH_d[id * nvi + lev] = CH_d[id] * sqrt(shear2) * pow(Altitude_d[lev], 2);
            //calculate part of sensible heat flux (need to *pt difference when updated)
            F_sens_d[id] = Cp * pow(p_int_d[id * nvi + lev] / P_Ref, kappa)
                           * Rho_int_d[id * nvi + lev] * KH_d[id * nvi + lev] / Altitude_d[lev];
        }
        else if (lev == nv) {
            //top level of model -> K = 0
            RiGrad_d[id * nvi + lev]  = LARGERiB;
            KH_d[id * nvi + lev]      = 0.0;
            KM_d[id * nvi + lev]      = 0.0;
            Rho_int_d[id * nvi + lev] = Rho_d[id * nv + lev - 1];
            p_int_d[id * nvi + lev]   = pressure_d[id * nv + lev - 1];
        }
        else {
            //levels in between
            pt_d[id * nv + lev] = pow(P_Ref, kappa) * pow(pressure_d[id * nv + lev], 1.0 - kappa)
                                  / (Rd * Rho_d[id * nv + lev]);
            pt_below = pow(P_Ref, kappa) * pow(pressure_d[id * nv + lev - 1], 1.0 - kappa)
                       / (Rd * Rho_d[id * nv + lev - 1]);
            pt_int = pt_below
                     + (pt_d[id * nv + lev] - pt_below) * (Altitudeh_d[lev] - Altitude_d[lev - 1])
                           / (Altitude_d[lev] - Altitude_d[lev - 1]);
            Rho_int_d[id * nvi + lev] = Rho_d[id * nv + lev - 1]
                                        + (Rho_d[id * nv + lev] - Rho_d[id * nv + lev - 1])
                                              * (Altitudeh_d[lev] - Altitude_d[lev - 1])
                                              / (Altitude_d[lev] - Altitude_d[lev - 1]);
            p_int_d[id * nvi + lev] = pressure_d[id * nv + lev - 1]
                                      + (pressure_d[id * nv + lev] - pressure_d[id * nv + lev - 1])
                                            * (Altitudeh_d[lev] - Altitude_d[lev - 1])
                                            / (Altitude_d[lev] - Altitude_d[lev - 1]);
            shear2 = pow((Mh_d[id * nv * 3 + lev * 3 + 0] / Rho_d[id * nv + lev]
                          - Mh_d[id * nv * 3 + (lev - 1) * 3 + 0] / Rho_d[id * nv + lev - 1])
                             / (Altitude_d[lev] - Altitude_d[lev - 1]),
                         2)
                     + pow((Mh_d[id * nv * 3 + lev * 3 + 1] / Rho_d[id * nv + lev]
                            - Mh_d[id * nv * 3 + (lev - 1) * 3 + 1] / Rho_d[id * nv + lev - 1])
                               / (Altitude_d[lev] - Altitude_d[lev - 1]),
                           2)
                     + pow((Mh_d[id * nv * 3 + lev * 3 + 2] / Rho_d[id * nv + lev]
                            - Mh_d[id * nv * 3 + (lev - 1) * 3 + 2] / Rho_d[id * nv + lev - 1])
                               / (Altitude_d[lev] - Altitude_d[lev - 1]),
                           2);
            //gradient Ri number
            if (shear2 == 0) {
                RiGrad_d[id * nvi + lev] = LARGERiB;
            }
            else {
                RiGrad_d[id * nvi + lev] = Gravit / pt_int * (pt_d[id * nv + lev] - pt_below)
                                           / (Altitude_d[lev] - Altitude_d[lev - 1]) / shear2;
            }

            //asymptotic length scale (constant in BL, decays to lower constant in free atmosphere)
            if ((Altitudeh_d[lev] <= asl_transition_height) || (asl_transition_height < 0)) {
                asym_len_scale_m = abl_asym_len;
                asym_len_scale_h = 3 * abl_asym_len;
            }
            else {
                asym_len_scale_m = free_asym_len
                                   + (abl_asym_len - free_asym_len)
                                         * exp(1 - Altitudeh_d[lev] / asl_transition_height);
                asym_len_scale_h = free_asym_len
                                   + (3 * abl_asym_len - free_asym_len)
                                         * exp(1 - Altitudeh_d[lev] / asl_transition_height);
            }
            mix_length_m = pow(1.0 / KVONKARMAN / Altitudeh_d[lev] + 1.0 / asym_len_scale_m, -1.0);
            mix_length_h = pow(1.0 / KVONKARMAN / Altitudeh_d[lev] + 1.0 / asym_len_scale_h, -1.0);
            if (RiGrad_d[id * nvi + lev] < 0.0) {
                e_mix_m = pow(mix_length_m, 2) * shear2 * (1 - 18.0 * RiGrad_d[id * nvi + lev]);
                e_mix_h = pow(mix_length_h, 2) * shear2 * (1 - 18.0 * RiGrad_d[id * nvi + lev]);
            }
            else {
                e_mix_m = pow(mix_length_m, 2) * shear2
                          / pow((1
                                 + 10.0 * RiGrad_d[id * nvi + lev]
                                       * (1 + 8.0 * RiGrad_d[id * nvi + lev])),
                                2);
                e_mix_h = pow(mix_length_m, 2) * shear2
                          / pow((1
                                 + 10.0 * RiGrad_d[id * nvi + lev]
                                       * (1 + 8.0 * RiGrad_d[id * nvi + lev])),
                                2);
            }
            if (e_mix_m < E_MIN_MIX)
                e_mix_m = E_MIN_MIX;
            if (e_mix_h < E_MIN_MIX)
                e_mix_h = E_MIN_MIX;

            KH_d[id * nvi + lev] = mix_length_h * sqrt(e_mix_h);
            KM_d[id * nvi + lev] = mix_length_m * sqrt(e_mix_m);

            if ((isnan(KH_d[id * nvi + lev])) || (isnan(KM_d[id * nvi + lev]))) {
                printf("%f %f\n", KH_d[id * nvi + lev], KM_d[id * nvi + lev]);
            }
        }
    }
}
