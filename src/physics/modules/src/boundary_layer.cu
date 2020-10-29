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

    cudaMalloc((void **)&atmp, esp.nvi * esp.point_num * sizeof(double));
    cudaMalloc((void **)&btmp, esp.nvi * esp.point_num * sizeof(double));
    cudaMalloc((void **)&ctmp, esp.nvi * esp.point_num * sizeof(double));
    cudaMalloc((void **)&cpr_tmp, esp.nvi * esp.point_num * sizeof(double));
    cudaMalloc((void **)&dtmp, 3 * esp.nvi * esp.point_num * sizeof(double));
    cudaMalloc((void **)&dpr_tmp, 3 * esp.nvi * esp.point_num * sizeof(double));
    cudaMalloc((void **)&RiB_d, esp.point_num * sizeof(double));
    cudaMalloc((void **)&KM_d, esp.nvi * esp.point_num * sizeof(double));
    cudaMalloc((void **)&KH_d, esp.nvi * esp.point_num * sizeof(double));
    cudaMalloc((void **)&CD_d, esp.point_num * sizeof(double));
    cudaMalloc((void **)&CH_d, esp.point_num * sizeof(double));
    cudaMalloc((void **)&vh_lowest_d, esp.point_num * sizeof(double));
    cudaMalloc((void **)&pt_surf_d, esp.point_num * sizeof(double));
    // cudaMalloc((void **)&p_surf_d, esp.point_num * sizeof(double));
    //cudaMalloc((void **)&Rho_surf_d, esp.point_num * sizeof(double));
    // cudaMalloc((void **)&zeta_d, esp.nvi * esp.point_num * sizeof(double));
    cudaMalloc((void **)&L_MO_d, esp.point_num * sizeof(double));
    cudaMalloc((void **)&F_sens_d, esp.point_num * sizeof(double));
    cudaMalloc((void **)&counter_grad_d, esp.nvi * esp.point_num * sizeof(double));

    cudaMalloc((void **)&RiGrad_d, esp.nvi * esp.point_num * sizeof(double));


    cudaMalloc((void **)&bl_top_lev_d, esp.point_num * sizeof(int));
    // cudaMalloc((void **)&sl_top_lev_d, esp.point_num * sizeof(int));
    cudaMalloc((void **)&bl_top_height_d, esp.point_num * sizeof(double));

    cudaMalloc((void **)&Rho_int_d, esp.nvi * esp.point_num * sizeof(double));
    cudaMalloc((void **)&p_int_d, esp.nvi * esp.point_num * sizeof(double));

    RiB_h        = (double *)malloc(esp.point_num * sizeof(double));
    KM_h         = (double *)malloc(esp.nvi * esp.point_num * sizeof(double));
    KH_h         = (double *)malloc(esp.nvi * esp.point_num * sizeof(double));
    bl_top_lev_h = (int *)malloc(esp.point_num * sizeof(int));
    // sl_top_lev_h    = (int *)malloc(esp.point_num * sizeof(int));
    bl_top_height_h = (double *)malloc(esp.point_num * sizeof(double));
    CD_h            = (double *)malloc(esp.point_num * sizeof(double));
    CH_h            = (double *)malloc(esp.point_num * sizeof(double));

    RiGrad_h = (double *)malloc(esp.nvi * esp.point_num * sizeof(double));

    cudaMemset(atmp, 0, sizeof(double) * esp.point_num * esp.nvi);
    cudaMemset(btmp, 0, sizeof(double) * esp.point_num * esp.nvi);
    cudaMemset(ctmp, 0, sizeof(double) * esp.point_num * esp.nvi);
    cudaMemset(cpr_tmp, 0, sizeof(double) * esp.point_num * esp.nvi);
    cudaMemset(dtmp, 0, sizeof(double) * 3 * esp.point_num * esp.nvi);
    cudaMemset(dpr_tmp, 0, sizeof(double) * 3 * esp.point_num * esp.nvi);
    cudaMemset(RiB_d, 0, sizeof(double) * esp.point_num);
    cudaMemset(L_MO_d, 0, sizeof(double) * esp.point_num);
    cudaMemset(F_sens_d, 0, sizeof(double) * esp.point_num);
    cudaMemset(counter_grad_d, 0, sizeof(double) * esp.nvi * esp.point_num);
    cudaMemset(KM_d, 0, sizeof(double) * esp.point_num * esp.nvi);
    cudaMemset(KH_d, 0, sizeof(double) * esp.point_num * esp.nvi);
    cudaMemset(bl_top_lev_d, 0, sizeof(int) * esp.point_num);
    // cudaMemset(sl_top_lev_d, 0, sizeof(int) * esp.point_num);
    cudaMemset(bl_top_height_d, 0, sizeof(double) * esp.point_num);
    cudaMemset(CD_d, 0, sizeof(double) * esp.point_num);
    cudaMemset(CH_d, 0, sizeof(double) * esp.point_num);
    cudaMemset(vh_lowest_d, 0, sizeof(double) * esp.point_num);
    cudaMemset(pt_surf_d, 0, sizeof(double) * esp.point_num);
    // cudaMemset(p_surf_d, 0, sizeof(double) * esp.point_num);
    // cudaMemset(Rho_surf_d, 0, sizeof(double) * esp.point_num);
    // cudaMemset(zeta_d, 0, sizeof(double) * esp.point_num * esp.nvi);
    cudaMemset(Rho_int_d, 0, sizeof(double) * esp.point_num * esp.nvi);
    cudaMemset(p_int_d, 0, sizeof(double) * esp.point_num * esp.nvi);

    cudaMemset(RiGrad_d, 0, sizeof(double) * esp.nvi * esp.point_num);

    return true;
}


bool boundary_layer::free_memory() {
    free(RiB_h);
    free(KM_h);
    free(KH_h);
    free(bl_top_lev_h);
    free(bl_top_height_h);
    // free(sl_top_lev_h);
    free(CD_h);
    free(CH_h);
    free(RiGrad_h);

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
    cudaFree(CH_d);
    cudaFree(vh_lowest_d);
    cudaFree(pt_surf_d);
    // cudaFree(p_surf_d);
    // cudaFree(Rho_surf_d);
    cudaFree(bl_top_lev_d);
    // cudaFree(zeta_d);
    cudaFree(L_MO_d);
    cudaFree(F_sens_d);
    cudaFree(bl_top_height_d);
    // cudaFree(sl_top_lev_d);
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
    else if (bl_type_str == "MoninObukhov" || bl_type_str == "MO") {
        bl_type = MONINOBUKHOV;
        config_OK &= true;
    }
    else if (bl_type_str == "EkmanSpiral" || bl_type_str == "Ekman") {
        bl_type = EKMANSPIRAL;
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
            z_therm_config,
            f_surf_layer_config);

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
    else if (bl_type == MONINOBUKHOV) {
        // CalcRiB<<<NBLEV, NTH>>>(esp.pressure_d,
        //                         esp.Rho_d,
        //                         esp.Mh_d,
        //                         esp.Tsurface_d,
        //                         esp.pt_d,
        //                         esp.Altitude_d,
        //                         esp.Altitudeh_d,
        //                         sim.Rd,
        //                         sim.Cp,
        //                         sim.P_Ref,
        //                         sim.Gravit,
        //                         Ri_crit,
        //                         z_rough,
        //                         z_therm,
        //                         RiB_d,
        //                         bl_top_lev_d,
        //                         bl_top_height_d,
        //                         sl_top_lev_d,
        //                         f_surf_layer,
        //                         pt_surf_d,
        //                         p_surf_d,
        //                         Rho_surf_d,
        //                         CD_d,
        //                         CH_d,
        //                         zeta_d,
        //                         vh_lowest_d,
        //                         Rho_int_d,
        //                         esp.point_num,
        //                         esp.nv);

        Calc_MOlength_Cdrag_BLdepth<<<NBLEV, NTH>>>(esp.pressure_d,
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
                                                    z_therm,
                                                    RiB_d,
                                                    bl_top_lev_d,
                                                    bl_top_height_d,
                                                    f_surf_layer,
                                                    pt_surf_d,
                                                    p_int_d,
                                                    CD_d,
                                                    CH_d,
                                                    L_MO_d,
                                                    F_sens_d,
                                                    vh_lowest_d,
                                                    Rho_int_d,
                                                    esp.point_num,
                                                    esp.nv);

        cudaDeviceSynchronize();

        CalcKM_KH<<<NBi, NTH>>>(RiB_d,
                                L_MO_d,
                                CD_d,
                                CH_d,
                                bl_top_height_d,
                                bl_top_lev_d,
                                F_sens_d,
                                counter_grad_d,
                                vh_lowest_d,
                                esp.Altitude_d,
                                esp.Altitudeh_d,
                                Ri_crit,
                                z_rough,
                                z_therm,
                                f_surf_layer,
                                pt_surf_d,
                                KM_d,
                                KH_d,
                                sim.Gravit,
                                esp.point_num);

        // TO DO
        // calc BL height from RiB
        // adjust Tsurface for sensible heat flux
        // how to adjust pressure? adjust pt first, then compute pressure? or is there a shortcut?
        // update pressure (implicitly) here, or add to qheat?
        cudaDeviceSynchronize();

        cudaMemset(atmp, 0, sizeof(double) * esp.point_num * esp.nvi);
        cudaMemset(btmp, 0, sizeof(double) * esp.point_num * esp.nvi);
        cudaMemset(ctmp, 0, sizeof(double) * esp.point_num * esp.nvi);
        cudaMemset(cpr_tmp, 0, sizeof(double) * esp.point_num * esp.nvi);
        cudaMemset(dtmp, 0, sizeof(double) * 3 * esp.point_num * esp.nvi);
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
                                           atmp,
                                           btmp,
                                           ctmp,
                                           cpr_tmp,
                                           dtmp,
                                           dpr_tmp,
                                           KM_d,
                                           Rho_int_d,
                                           time_step,
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
            diag.dump_data("bl_thomas_momentum", nstep, 0, 0, esp, esp.point_num, esp.nv);
            unsigned int flag = diag.get_flag();
            if ((flag & BL_THOMAS_NOT_DD) == BL_THOMAS_NOT_DD)
                log::printf("Non diagonally dominant matrix for thomas algorithm in boundary layer "
                            "momentum\n");
        }
#endif

        cudaMemset(atmp, 0, sizeof(double) * esp.point_num * esp.nvi);
        cudaMemset(btmp, 0, sizeof(double) * esp.point_num * esp.nvi);
        cudaMemset(ctmp, 0, sizeof(double) * esp.point_num * esp.nvi);
        cudaMemset(cpr_tmp, 0, sizeof(double) * esp.point_num * esp.nvi);
        cudaMemset(dtmp, 0, sizeof(double) * 3 * esp.point_num * esp.nvi);
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
                                                atmp,
                                                btmp,
                                                ctmp,
                                                cpr_tmp,
                                                dtmp,
                                                dpr_tmp,
                                                KH_d,
                                                Rho_int_d,
                                                pt_surf_d,
                                                p_int_d,
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
        // Heat_Diff_Expl<<<NBLEV, NTH>>>(esp.pt_d,
        //                                esp.Rho_d,
        //                                esp.Altitude_d,
        //                                esp.Altitudeh_d,
        //                                esp.Tsurface_d,
        //                                esp.dTsurf_dt_d,
        //                                KH_d,
        //                                Rho_int_d,
        //                                pt_surf_d,
        //                                Rho_surf_d,
        //                                p_surf_d,
        //                                time_step,
        //                                sim.Rd,
        //                                sim.Cp,
        //                                sim.P_Ref,
        //                                esp.Csurf,
        //                                esp.profx_Qheat_d,
        //                                esp.point_num,
        //                                esp.nv,
        //                                bl_top_lev_d);
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

#if defined(DIAG_CHECK_BL_THOMAS_DIAG_DOM)
        //checking diagonal dominance in thomas algorithm
        diag.reset_flag();
#endif

        Momentum_Diff_Impl<<<NBLEV, NTH>>>(esp.Mh_d,
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
                                           Rho_int_d, ///broken!
                                           time_step,
                                           sim.A,
                                           esp.point_num,
                                           esp.nv,
                                           bl_top_lev_d,
                                           sim.DeepModel,
                                           *diag.diagnostics_global_flag,
                                           *diag.diagnostics);
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
                                 z_therm,
                                 RiGrad_d,
                                 pt_surf_d,
                                 p_int_d,
                                 CD_d,
                                 CH_d,
                                 Rho_int_d,
                                 KM_d,
                                 KH_d,
                                 esp.point_num,
                                 esp.nv);

        cudaDeviceSynchronize();

        cudaMemset(atmp, 0, sizeof(double) * esp.point_num * esp.nvi);
        cudaMemset(btmp, 0, sizeof(double) * esp.point_num * esp.nvi);
        cudaMemset(ctmp, 0, sizeof(double) * esp.point_num * esp.nvi);
        cudaMemset(cpr_tmp, 0, sizeof(double) * esp.point_num * esp.nvi);
        cudaMemset(dtmp, 0, sizeof(double) * 3 * esp.point_num * esp.nvi);
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
                                           atmp,
                                           btmp,
                                           ctmp,
                                           cpr_tmp,
                                           dtmp,
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

        cudaMemset(atmp, 0, sizeof(double) * esp.point_num * esp.nvi);
        cudaMemset(btmp, 0, sizeof(double) * esp.point_num * esp.nvi);
        cudaMemset(ctmp, 0, sizeof(double) * esp.point_num * esp.nvi);
        cudaMemset(cpr_tmp, 0, sizeof(double) * esp.point_num * esp.nvi);
        cudaMemset(dtmp, 0, sizeof(double) * 3 * esp.point_num * esp.nvi);
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
                                                atmp,
                                                btmp,
                                                ctmp,
                                                cpr_tmp,
                                                dtmp,
                                                dpr_tmp,
                                                KH_d,
                                                Rho_int_d,
                                                pt_surf_d,
                                                p_int_d,
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

    // roughness and thermal length scaling
    config_reader.append_config_var("z_rough", z_rough_config, z_rough_config);
    config_reader.append_config_var("z_therm", z_therm_config, z_therm_config);

    // surface layer fraction
    config_reader.append_config_var("f_surf_layer", f_surf_layer_config, f_surf_layer_config);

    return true;
}

bool boundary_layer::store(const ESP &esp, storage &s) {
    if (bl_type == MONINOBUKHOV) {
        cudaMemcpy(RiB_h, RiB_d, esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
        s.append_table(RiB_h, esp.point_num, "/RiB", " ", "bulk Richardson number");

        cudaMemcpy(CD_h, CD_d, esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
        s.append_table(CD_h, esp.point_num, "/CD", " ", "surface drag coefficient");

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

        cudaMemcpy(bl_top_height_h,
                   bl_top_height_d,
                   esp.point_num * sizeof(double),
                   cudaMemcpyDeviceToHost);
        s.append_table(
            bl_top_height_h, esp.point_num, "/bl_top_height", "m", "Height of boundary layer");

        cudaMemcpy(bl_top_lev_h, bl_top_lev_d, esp.point_num * sizeof(int), cudaMemcpyDeviceToHost);
        s.append_table(bl_top_lev_h,
                       esp.point_num,
                       "/bl_top_lev",
                       " ",
                       "Index of highest level in boundary layer");
    }
    else if (bl_type == LOCALMIXL) {
        cudaMemcpy(
            RiGrad_h, RiGrad_d, esp.nvi * esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
        s.append_table(
            RiGrad_h, esp.nvi * esp.point_num, "/RiGrad", " ", "gradient Richardson number");

        cudaMemcpy(CD_h, CD_d, esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
        s.append_table(CD_h, esp.point_num, "/CD", " ", "surface drag coefficient");

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
    else if (bl_type == LOCALMIXL) {
        Ri_crit = Ri_crit_;
        z_rough = z_rough_;
        z_therm = z_therm_;
        for (int id = 0; id < esp.point_num; id++) {
            bl_top_lev_h[id] = esp.nv - 1; //local formulation includes entire column
        }
        cudaMemcpy(bl_top_lev_d, bl_top_lev_h, esp.point_num * sizeof(int), cudaMemcpyHostToDevice);
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

__global__ void
Momentum_Diff_Impl(double *      Mh_d,
                   double *      pressure_d,
                   double *      Rho_d,
                   double *      Altitude_d,
                   double *      Altitudeh_d,
                   double *      atmp,
                   double *      btmp,
                   double *      ctmp,
                   double *      cpr_tmp,
                   double *      dtmp,
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
                atmp[id * (nv + 1) + lev] = 0;
                btmp[id * (nv + 1) + lev] =
                    -(1.0 / (Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                          * (pow(rup, 2) * Rho_int_d[id * (nv + 1) + lev + 1]
                                 * KM_d[id * (nv + 1) + lev + 1]
                                 / (Altitude_d[lev + 1] - Altitude_d[lev])
                             + pow(rlow, 2) * Rho_int_d[id * (nv + 1) + lev]
                                   * KM_d[id * (nv + 1) + lev] / Altitude_d[lev])
                      + Rho_d[id * nv + lev] / time_step);
                ctmp[id * (nv + 1) + lev] = pow(rup, 2) * Rho_int_d[id * (nv + 1) + lev + 1]
                                            * KM_d[id * (nv + 1) + lev + 1]
                                            / ((Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                                               * (Altitude_d[lev + 1] - Altitude_d[lev]));
                cpr_tmp[id * (nv + 1) + lev] =
                    ctmp[id * (nv + 1) + lev] / btmp[id * (nv + 1) + lev];
                for (int k = 0; k < 3; k++) {
                    dtmp[id * (nv + 1) * 3 + lev * 3 + k] =
                        -Mh_d[id * 3 * nv + lev * 3 + k] / time_step;
                    dpr_tmp[id * (nv + 1) * 3 + lev * 3 + k] =
                        dtmp[id * (nv + 1) * 3 + lev * 3 + k] / btmp[id * (nv + 1) + lev];
                }
            }
            else if (lev == bl_top_lev_d[id]) {
                atmp[id * (nv + 1) + lev] = pow(rlow, 2) * Rho_int_d[id * (nv + 1) + lev]
                                            * KM_d[id * (nv + 1) + lev]
                                            / ((Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                                               * (Altitude_d[lev] - Altitude_d[lev - 1]));
                btmp[id * (nv + 1) + lev] =
                    -(1.0 / (Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                          * (pow(rup, 2) * Rho_int_d[id * (nv + 1) + lev + 1]
                                 * KM_d[id * (nv + 1) + lev + 1]
                                 / (Altitude_d[lev + 1] - Altitude_d[lev])
                             + (pow(rlow, 2) * Rho_int_d[id * (nv + 1) + lev])
                                   * KM_d[id * (nv + 1) + lev]
                                   / (Altitude_d[lev] - Altitude_d[lev - 1]))
                      + Rho_d[id * nv + lev] / time_step);
                ctmp[id * (nv + 1) + lev]    = 0;
                cpr_tmp[id * (nv + 1) + lev] = 0; //not used, i think
                for (int k = 0; k < 3; k++) {
                    dtmp[id * (nv + 1) * 3 + lev * 3 + k] =
                        -Mh_d[id * 3 * nv + lev * 3 + k] / time_step
                        - pow(rup, 2) * Rho_int_d[id * (nv + 1) + lev + 1]
                              * KM_d[id * (nv + 1) + lev + 1]
                              / (Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                              * Mh_d[id * 3 * nv + (lev + 1) * 3 + k] / Rho_d[id * nv + lev + 1]
                              / (Altitude_d[lev + 1] - Altitude_d[lev]);
                    dpr_tmp[id * (nv + 1) * 3 + lev * 3 + k] =
                        (dtmp[id * (nv + 1) * 3 + lev * 3 + k]
                         - atmp[id * (nv + 1) + lev]
                               * dpr_tmp[id * (nv + 1) * 3 + (lev - 1) * 3 + k])
                        / (btmp[id * (nv + 1) + lev]
                           - atmp[id * (nv + 1) + lev] * cpr_tmp[id * (nv + 1) + lev - 1]);
                }
            }
            else {
                atmp[id * (nv + 1) + lev] = pow(rlow, 2) * Rho_int_d[id * (nv + 1) + lev]
                                            * KM_d[id * (nv + 1) + lev]
                                            / ((Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                                               * (Altitude_d[lev] - Altitude_d[lev - 1]));
                btmp[id * (nv + 1) + lev] =
                    -(1.0 / (Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                          * (pow(rup, 2) * Rho_int_d[id * (nv + 1) + lev + 1]
                                 * KM_d[id * (nv + 1) + lev + 1]
                                 / (Altitude_d[lev + 1] - Altitude_d[lev])
                             + pow(rlow, 2) * Rho_int_d[id * (nv + 1) + lev]
                                   * KM_d[id * (nv + 1) + lev]
                                   / (Altitude_d[lev] - Altitude_d[lev - 1]))
                      + Rho_d[id * nv + lev] / time_step);
                ctmp[id * (nv + 1) + lev] = pow(rup, 2) * Rho_int_d[id * (nv + 1) + lev + 1]
                                            * KM_d[id * (nv + 1) + lev + 1]
                                            / ((Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                                               * (Altitude_d[lev + 1] - Altitude_d[lev]));
                cpr_tmp[id * (nv + 1) + lev] =
                    ctmp[id * (nv + 1) + lev]
                    / (btmp[id * (nv + 1) + lev]
                       - atmp[id * (nv + 1) + lev] * cpr_tmp[id * (nv + 1) + lev - 1]);
                for (int k = 0; k < 3; k++) {
                    dtmp[id * (nv + 1) * 3 + lev * 3 + k] =
                        -Mh_d[id * 3 * nv + lev * 3 + k] / time_step;
                    dpr_tmp[id * (nv + 1) * 3 + lev * 3 + k] =
                        (dtmp[id * (nv + 1) * 3 + lev * 3 + k]
                         - atmp[id * (nv + 1) + lev]
                               * dpr_tmp[id * (nv + 1) * 3 + (lev - 1) * 3 + k])
                        / (btmp[id * (nv + 1) + lev]
                           - atmp[id * (nv + 1) + lev] * cpr_tmp[id * (nv + 1) + lev - 1]);
                }
            }

#ifdef DIAG_CHECK_BL_THOMAS_DIAG_DOM
            //reset values
            diagnostics_data[id * nv + lev].flag = 0;
            diagnostics_data[id * nv + lev].data = make_double4(0.0, 0.0, 0.0, 0.0);
            //check if matrix diagonally dominant
            if (lev < bl_top_lev_d[id] + 1) {
                if (!(fabs(btmp[id * (nv + 1) + lev])
                      > THOMAS_DIAG_DOM_FACTOR
                            * (fabs(atmp[id * (nv + 1) + lev])
                               + fabs(ctmp[id * (nv + 1) + lev])))) {
                    atomicOr(diagnostics_flag, BL_THOMAS_NOT_DD);
                    diagnostics_data[id * nv + lev].flag   = BL_THOMAS_NOT_DD;
                    diagnostics_data[id * nv + lev].data.x = atmp[id * (nv + 1) + lev];
                    diagnostics_data[id * nv + lev].data.y = btmp[id * (nv + 1) + lev];
                    diagnostics_data[id * nv + lev].data.z = ctmp[id * (nv + 1) + lev];
                    // printf("Warning! Thomas algorithm in boundary layer mom. equation unstable\n");
                }
            }
#endif
        }
        // if (id == 1000) {
        //     printf("stop");
        // }

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
        // if (id == 0) {
        //     printf("%f\n", Mh_d[id * nv * 3 + 0]);
        // }
    }
}

__global__ void Heat_Diff_Impl(double *      pt_d,
                               double *      pressure_d,
                               double *      Rho_d,
                               double *      Altitude_d,
                               double *      Altitudeh_d,
                               double *      Tsurface_d,
                               double *      atmp,
                               double *      btmp,
                               double *      ctmp,
                               double *      cpr_tmp,
                               double *      dtmp,
                               double *      dpr_tmp,
                               double *      KH_d,
                               double *      Rho_int_d,
                               double *      pt_surf_d,
                               double *      p_surf_d,
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
                atmp[id * nv + lev + 1] = 0;
                btmp[id * nv + lev + 1] =
                    -(Rho_int_d[id * (nv + 1)] * Cp * KH_d[id * (nv + 1) + lev + 1]
                          / Altitude_d[lev + 1]
                      + Csurf / time_step);
                ctmp[id * nv + lev + 1] = Rho_int_d[id * (nv + 1)] * Cp
                                          * KH_d[id * (nv + 1) + lev + 1] / Altitude_d[lev + 1];
                cpr_tmp[id * nv + lev + 1] = ctmp[id * nv + lev + 1] / btmp[id * nv + lev + 1];
                dtmp[id * nv + lev + 1]    = -Csurf * pt_surf_d[id] / time_step;
                dpr_tmp[id * nv + lev + 1] = dtmp[id * nv + lev + 1] / btmp[id * nv + lev + 1];
            }
            else if (lev == 0) { //lowest layer, v at lowest boundary = 0, dz0 = Altitude0
                atmp[id * nv + lev + 1] =
                    pow(rlow, 2) * Rho_int_d[id * (nv + 1) + lev] * KH_d[id * (nv + 1) + lev]
                    / (Altitudeh_d[lev + 1] - Altitudeh_d[lev]) / Altitude_d[lev];
                btmp[id * nv + lev + 1] =
                    -(1.0 / (Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                          * (pow(rup, 2) * Rho_int_d[id * (nv + 1) + lev + 1]
                                 * KH_d[id * (nv + 1) + lev + 1]
                                 / (Altitude_d[lev + 1] - Altitude_d[lev])
                             + pow(rlow, 2) * Rho_int_d[id * (nv + 1) + lev]
                                   * KH_d[id * (nv + 1) + lev] / Altitude_d[lev])
                      + Rho_d[id * nv + lev] / time_step);
                ctmp[id * nv + lev + 1] = pow(rup, 2) * Rho_int_d[id * (nv + 1) + lev + 1]
                                          * KH_d[id * (nv + 1) + lev + 1]
                                          / ((Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                                             * (Altitude_d[lev + 1] - Altitude_d[lev]));
                cpr_tmp[id * nv + lev + 1] =
                    ctmp[id * nv + lev + 1]
                    / (btmp[id * nv + lev + 1] - atmp[id * nv + lev + 1] * cpr_tmp[id * nv + lev]);
                dtmp[id * nv + lev + 1] = -Rho_d[id * nv + lev] * pt_d[id * nv + lev] / time_step;
                dpr_tmp[id * nv + lev + 1] =
                    (dtmp[id * nv + lev + 1] - atmp[id * nv + lev + 1] * dpr_tmp[id * nv + lev])
                    / (btmp[id * nv + lev + 1] - atmp[id * nv + lev + 1] * cpr_tmp[id * nv + lev]);
            }
            else if (lev == bl_top_lev_d[id]) {
                atmp[id * nv + lev + 1] = pow(rlow, 2) * Rho_int_d[id * (nv + 1) + lev]
                                          * KH_d[id * (nv + 1) + lev]
                                          / ((Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                                             * (Altitude_d[lev] - Altitude_d[lev - 1]));
                btmp[id * nv + lev + 1]    = -(1.0 / (Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                                                * (pow(rup, 2) * Rho_int_d[id * (nv + 1) + lev + 1]
                                                       * KH_d[id * (nv + 1) + lev + 1]
                                                       / (Altitude_d[lev + 1] - Altitude_d[lev])
                                                   + pow(rlow, 2) * Rho_int_d[id * (nv + 1) + lev]
                                                         * KH_d[id * (nv + 1) + lev]
                                                         / (Altitude_d[lev] - Altitude_d[lev - 1]))
                                            + Rho_d[id * nv + lev] / time_step);
                ctmp[id * nv + lev + 1]    = 0;
                cpr_tmp[id * nv + lev + 1] = 0; //not used, i think
                dtmp[id * nv + lev + 1]    = -Rho_d[id * nv + lev] * pt_d[id * nv + lev] / time_step
                                          - pow(rup, 2) * Rho_int_d[id * (nv + 1) + lev + 1]
                                                * KH_d[id * (nv + 1) + lev + 1]
                                                / (Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                                                * pt_d[id * nv + (lev + 1)]
                                                / (Altitude_d[lev + 1] - Altitude_d[lev]);
                dpr_tmp[id * nv + lev + 1] =
                    (dtmp[id * nv + lev + 1] - atmp[id * nv + lev + 1] * dpr_tmp[id * nv + (lev)])
                    / (btmp[id * nv + lev + 1] - atmp[id * nv + lev + 1] * cpr_tmp[id * nv + lev]);
            }
            else {
                atmp[id * nv + lev + 1] = pow(rlow, 2) * Rho_int_d[id * (nv + 1) + lev]
                                          * KH_d[id * (nv + 1) + lev]
                                          / ((Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                                             * (Altitude_d[lev] - Altitude_d[lev - 1]));
                btmp[id * nv + lev + 1] = -(1.0 / (Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                                                * (pow(rup, 2) * Rho_int_d[id * (nv + 1) + lev + 1]
                                                       * KH_d[id * (nv + 1) + lev + 1]
                                                       / (Altitude_d[lev + 1] - Altitude_d[lev])
                                                   + pow(rlow, 2) * Rho_int_d[id * (nv + 1) + lev]
                                                         * KH_d[id * (nv + 1) + lev]
                                                         / (Altitude_d[lev] - Altitude_d[lev - 1]))
                                            + Rho_d[id * nv + lev] / time_step);
                ctmp[id * nv + lev + 1] = pow(rup, 2) * Rho_int_d[id * (nv + 1) + lev + 1]
                                          * KH_d[id * (nv + 1) + lev + 1]
                                          / ((Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                                             * (Altitude_d[lev + 1] - Altitude_d[lev]));
                cpr_tmp[id * nv + lev + 1] =
                    ctmp[id * nv + lev + 1]
                    / (btmp[id * nv + lev + 1] - atmp[id * nv + lev + 1] * cpr_tmp[id * nv + lev]);
                dtmp[id * nv + lev + 1] = -Rho_d[id * nv + lev] * pt_d[id * nv + lev] / time_step;
                dpr_tmp[id * nv + lev + 1] =
                    (dtmp[id * nv + lev + 1] - atmp[id * nv + lev + 1] * dpr_tmp[id * nv + (lev)])
                    / (btmp[id * nv + lev + 1] - atmp[id * nv + lev + 1] * cpr_tmp[id * nv + lev]);
            }
#ifdef DIAG_CHECK_BL_THOMAS_DIAG_DOM
            //reset values
            diagnostics_data[id * nv + lev + 1].flag = 0;
            diagnostics_data[id * nv + lev + 1].data = make_double4(0.0, 0.0, 0.0, 0.0);
            //check if matrix diagonally dominant
            if (lev < bl_top_lev_d[id] + 1) {
                if (!(fabs(btmp[id * nv + lev + 1])
                      > THOMAS_DIAG_DOM_FACTOR
                            * (fabs(atmp[id * nv + lev + 1]) + fabs(ctmp[id * nv + lev + 1])))) {
                    atomicOr(diagnostics_flag, BL_THOMAS_NOT_DD);
                    diagnostics_data[id * nv + lev + 1].flag   = BL_THOMAS_NOT_DD;
                    diagnostics_data[id * nv + lev + 1].data.x = atmp[id * nv + lev + 1];
                    diagnostics_data[id * nv + lev + 1].data.y = btmp[id * nv + lev + 1];
                    diagnostics_data[id * nv + lev + 1].data.z = ctmp[id * nv + lev + 1];
                }
            }
#endif
        }

        for (lev = bl_top_lev_d[id]; lev >= 0; lev--) {
            //backward sweep
            if (lev == bl_top_lev_d[id]) {
                pt_d[id * nv + lev] = dpr_tmp[id * nv + lev + 1];
            }
            else {
                pt_d[id * nv + lev] = (dpr_tmp[id * nv + lev + 1]
                                       - cpr_tmp[id * nv + lev + 1] * pt_d[id * nv + (lev + 1)]);
            }
            // convert back to pressure
            // pressure_d[id * nv + lev] =
            // P_Ref * pow(Rho_d[id * nv + lev] * Rd * pt_d[id * nv + lev] / P_Ref, Cp / Cv);
            pressure_d[id * nv + lev] =
                P_Ref * pow(Rd * Rho_d[id * nv + lev] * pt_d[id * nv + lev] / P_Ref, Cp / Cv);
        }

        pt_surf_d[id] = (dpr_tmp[id * nv + 0] - cpr_tmp[id * nv + 0] * pt_d[id * nv + 0]);
        //Tsurface_d[id] =
        //pow(Rho_int_d[id * (nv + 1) + 0] * Rd / P_Ref, Rd / Cv) * pow(pt_surf_d[id], Cp / Cv);
        // p_surf_d[id]   = Rho_surf_d[id] * Rd * Tsurface_d[id];
        Tsurface_d[id] = pt_surf_d[id] * pow(p_surf_d[id] / P_Ref, kappa);
        // Rho_surf_d[id] = p_surf_d[id] / (Rd * Tsurface_d[id]);
        //Tsurface_d[id] -= Fsen / Csurf * time_step;
        if (id == 10) {
            printf("\nTs = %f\n", Tsurface_d[id]);
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
                                        double *      atmp,
                                        double *      btmp,
                                        double *      ctmp,
                                        double *      cpr_tmp,
                                        double *      dtmp,
                                        double *      dpr_tmp,
                                        double *      KH_d,
                                        double *      Rho_int_d,
                                        double *      pt_surf_d,
                                        double *      p_int_d,
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
                atmp[id * (nv + 1) + lev + 1] = 0;
                btmp[id * (nv + 1) + lev + 1] =
                    -(Cp / Csurf * Rho_int_d[id * (nv + 1)] * KH_d[id * (nv + 1) + lev + 1]
                          / Altitude_d[lev + 1]
                      + 1.0 / time_step);
                ctmp[id * (nv + 1) + lev + 1] = Cp / Csurf * Rho_int_d[id * (nv + 1)]
                                                * KH_d[id * (nv + 1) + lev + 1]
                                                / Altitude_d[lev + 1];
                cpr_tmp[id * (nv + 1) + lev + 1] =
                    ctmp[id * (nv + 1) + lev + 1] / btmp[id * (nv + 1) + lev + 1];
                dtmp[id * (nv + 1) + lev + 1] = -pt_surf_d[id] / time_step;
                dpr_tmp[id * (nv + 1) + lev + 1] =
                    dtmp[id * (nv + 1) + lev + 1] / btmp[id * (nv + 1) + lev + 1];
            }
            else if (lev == 0) { //lowest layer, v at lowest boundary = 0, dz0 = Altitude0
                atmp[id * (nv + 1) + lev + 1] =
                    Cp / Cv * pow(rlow, 2) * Rho_int_d[id * (nv + 1) + lev]
                    * KH_d[id * (nv + 1) + lev] / (Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                    / Altitude_d[lev] * pow(p_int_d[id * (nv + 1) + lev] / P_Ref, kappa);
                ctmp[id * (nv + 1) + lev + 1] =
                    Cp / Cv * pow(rup, 2) * Rho_int_d[id * (nv + 1) + lev + 1]
                    * KH_d[id * (nv + 1) + lev + 1]
                    / ((Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                       * (Altitude_d[lev + 1] - Altitude_d[lev]))
                    * pow(p_int_d[id * (nv + 1) + lev + 1] / P_Ref, kappa);
                btmp[id * (nv + 1) + lev + 1] =
                    -(atmp[id * (nv + 1) + lev + 1] + ctmp[id * (nv + 1) + lev + 1]
                      + Rho_d[id * nv + lev] / time_step
                            * pow(pressure_d[id * nv + lev] / P_Ref, kappa));
                cpr_tmp[id * (nv + 1) + lev + 1] =
                    ctmp[id * (nv + 1) + lev + 1]
                    / (btmp[id * (nv + 1) + lev + 1]
                       - atmp[id * (nv + 1) + lev + 1] * cpr_tmp[id * (nv + 1) + lev]);
                dtmp[id * (nv + 1) + lev + 1] = -Rho_d[id * nv + lev] * pt_d[id * nv + lev]
                                                / time_step
                                                * pow(pressure_d[id * nv + lev] / P_Ref, kappa);
                dpr_tmp[id * (nv + 1) + lev + 1] =
                    (dtmp[id * (nv + 1) + lev + 1]
                     - atmp[id * (nv + 1) + lev + 1] * dpr_tmp[id * (nv + 1) + lev])
                    / (btmp[id * (nv + 1) + lev + 1]
                       - atmp[id * (nv + 1) + lev + 1] * cpr_tmp[id * (nv + 1) + lev]);
            }
            else if (lev == bl_top_lev_d[id]) {
                atmp[id * (nv + 1) + lev + 1] = Cp / Cv * pow(rlow, 2)
                                                * Rho_int_d[id * (nv + 1) + lev]
                                                * KH_d[id * (nv + 1) + lev]
                                                / ((Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                                                   * (Altitude_d[lev] - Altitude_d[lev - 1]))
                                                * pow(p_int_d[id * (nv + 1) + lev] / P_Ref, kappa);
                btmp[id * (nv + 1) + lev + 1] =
                    -(Cp / Cv / (Altitudeh_d[lev + 1] - Altitudeh_d[lev]) * pow(rup, 2)
                          * Rho_int_d[id * (nv + 1) + lev + 1] * KH_d[id * (nv + 1) + lev + 1]
                          / (Altitude_d[lev + 1] - Altitude_d[lev])
                          * pow(p_int_d[id * (nv + 1) + lev + 1] / P_Ref, kappa)
                      + atmp[id * (nv + 1) + lev + 1]
                      + Rho_d[id * nv + lev] / time_step
                            * pow(pressure_d[id * nv + lev] / P_Ref, kappa));
                ctmp[id * (nv + 1) + lev + 1]    = 0;
                cpr_tmp[id * (nv + 1) + lev + 1] = 0; //not used, i think
                dtmp[id * (nv + 1) + lev + 1] =
                    -Rho_d[id * nv + lev] * pt_d[id * nv + lev] / time_step
                        * pow(pressure_d[id * nv + lev] / P_Ref, kappa)
                    - Cp / Cv * pow(rup, 2) * Rho_int_d[id * (nv + 1) + lev + 1]
                          * KH_d[id * (nv + 1) + lev + 1]
                          / (Altitudeh_d[lev + 1] - Altitudeh_d[lev]) * pt_d[id * nv + (lev + 1)]
                          / (Altitude_d[lev + 1] - Altitude_d[lev])
                          * pow(p_int_d[id * (nv + 1) + lev + 1] / P_Ref, kappa);
                dpr_tmp[id * (nv + 1) + lev + 1] =
                    (dtmp[id * (nv + 1) + lev + 1]
                     - atmp[id * (nv + 1) + lev + 1] * dpr_tmp[id * (nv + 1) + (lev)])
                    / (btmp[id * (nv + 1) + lev + 1]
                       - atmp[id * (nv + 1) + lev + 1] * cpr_tmp[id * (nv + 1) + lev]);
            }
            else {
                atmp[id * (nv + 1) + lev + 1] = Cp / Cv * pow(rlow, 2)
                                                * Rho_int_d[id * (nv + 1) + lev]
                                                * KH_d[id * (nv + 1) + lev]
                                                / ((Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                                                   * (Altitude_d[lev] - Altitude_d[lev - 1]))
                                                * pow(p_int_d[id * (nv + 1) + lev] / P_Ref, kappa);
                ctmp[id * (nv + 1) + lev + 1] =
                    Cp / Cv * pow(rup, 2) * Rho_int_d[id * (nv + 1) + lev + 1]
                    * KH_d[id * (nv + 1) + lev + 1]
                    / ((Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                       * (Altitude_d[lev + 1] - Altitude_d[lev]))
                    * pow(p_int_d[id * (nv + 1) + lev + 1] / P_Ref, kappa);
                btmp[id * (nv + 1) + lev + 1] =
                    -(atmp[id * (nv + 1) + lev + 1] + ctmp[id * (nv + 1) + lev + 1]
                      + Rho_d[id * nv + lev] / time_step
                            * pow(pressure_d[id * nv + lev] / P_Ref, kappa));
                cpr_tmp[id * (nv + 1) + lev + 1] =
                    ctmp[id * (nv + 1) + lev + 1]
                    / (btmp[id * (nv + 1) + lev + 1]
                       - atmp[id * (nv + 1) + lev + 1] * cpr_tmp[id * (nv + 1) + lev]);
                dtmp[id * (nv + 1) + lev + 1] = -Rho_d[id * nv + lev] * pt_d[id * nv + lev]
                                                / time_step
                                                * pow(pressure_d[id * nv + lev] / P_Ref, kappa);
                dpr_tmp[id * (nv + 1) + lev + 1] =
                    (dtmp[id * (nv + 1) + lev + 1]
                     - atmp[id * (nv + 1) + lev + 1] * dpr_tmp[id * (nv + 1) + (lev)])
                    / (btmp[id * (nv + 1) + lev + 1]
                       - atmp[id * (nv + 1) + lev + 1] * cpr_tmp[id * (nv + 1) + lev]);
            }
#ifdef DIAG_CHECK_BL_THOMAS_DIAG_DOM
            //reset values
            diagnostics_data[id * nv + lev + 1].flag = 0;
            diagnostics_data[id * nv + lev + 1].data = make_double4(0.0, 0.0, 0.0, 0.0);
            //check if matrix diagonally dominant
            if (lev < bl_top_lev_d[id] + 1) {
                if (!(fabs(btmp[id * (nv + 1) + lev + 1])
                      > THOMAS_DIAG_DOM_FACTOR
                            * (fabs(atmp[id * (nv + 1) + lev + 1])
                               + fabs(ctmp[id * (nv + 1) + lev + 1])))) {
                    atomicOr(diagnostics_flag, BL_THOMAS_NOT_DD);
                    diagnostics_data[id * nv + lev + 1].flag   = BL_THOMAS_NOT_DD;
                    diagnostics_data[id * nv + lev + 1].data.x = atmp[id * (nv + 1) + lev + 1];
                    diagnostics_data[id * nv + lev + 1].data.y = btmp[id * (nv + 1) + lev + 1];
                    diagnostics_data[id * nv + lev + 1].data.z = ctmp[id * (nv + 1) + lev + 1];
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
            // convert back to pressure
            // pressure_d[id * nv + lev] =
            // P_Ref * pow(Rho_d[id * nv + lev] * Rd * pt_d[id * nv + lev] / P_Ref, Cp / Cv);
            //pressure_d[id * nv + lev] =
            //    P_Ref * pow(Rd * Rho_d[id * nv + lev] * pt_d[id * nv + lev] / P_Ref, Cp / Cv);
            //need to first update temp using previous pressure (to conserve energy)
            temperature_d[id * nv + lev] =
                pt_d[id * nv + lev] * pow(pressure_d[id * nv + lev] / P_Ref, kappa);
            //now we can calculate new pressure
            pressure_d[id * nv + lev] = temperature_d[id * nv + lev] * (Rd * Rho_d[id * nv + lev]);
        }

        pt_surf_d[id] =
            (dpr_tmp[id * (nv + 1) + 0] - cpr_tmp[id * (nv + 1) + 0] * pt_d[id * nv + 0]);
        //Tsurface_d[id] =
        //pow(Rho_int_d[id * (nv + 1) + 0] * Rd / P_Ref, Rd / Cv) * pow(pt_surf_d[id], Cp / Cv);
        // p_surf_d[id]   = Rho_surf_d[id] * Rd * Tsurface_d[id];
        Tsurface_d[id] = pt_surf_d[id] * pow(p_int_d[id * (nv + 1) + 0] / P_Ref, kappa);
        // Rho_surf_d[id] = p_surf_d[id] / (Rd * Tsurface_d[id]);
        //Tsurface_d[id] -= Fsen / Csurf * time_step;
    }
}


__global__ void Heat_Diff_Expl(double *pt_d,
                               double *Rho_d,
                               double *Altitude_d,
                               double *Altitudeh_d,
                               double *Tsurface_d,
                               double *dTsurf_dt_d,
                               double *KH_d,
                               double *Rho_int_d,
                               double *pt_surf_d,
                               double *Rho_surf_d,
                               double *p_surf_d,
                               double  time_step,
                               double  Rd,
                               double  Cp,
                               double  P_Ref,
                               double  Csurf,
                               double *profx_Qheat_d,
                               int     num,
                               int     nv,
                               int *   bl_top_lev_d) {

    int id = blockIdx.x * blockDim.x + threadIdx.x;
    int lev;

    if (id < num) {
        // double Cv    = Cp - Rd;
        // double kappa = Rd / Cp;

        for (lev = 0; lev < bl_top_lev_d[id] + 1; lev++) {
            if (lev == 0) {
                profx_Qheat_d[id * nv + lev] +=
                    Cp / (Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                    * (Rho_int_d[id * (nv + 1) + lev + 1] * KH_d[id * (nv + 1) + lev + 1]
                           * (pt_d[id * nv + lev + 1] - pt_d[id * nv + lev])
                           / (Altitude_d[lev + 1] - Altitude_d[lev])
                       - Rho_int_d[id * (nv + 1) + lev] * KH_d[id * (nv + 1) + lev]
                             * (pt_d[id * nv + lev] - pt_surf_d[id]) / (Altitude_d[lev]));
            }
            // else if (lev == bl_top_lev_d[id]) {
            //
            //
            // }
            else {
                profx_Qheat_d[id * nv + lev] +=
                    Cp / (Altitudeh_d[lev + 1] - Altitudeh_d[lev])
                    * (Rho_int_d[id * (nv + 1) + lev + 1] * KH_d[id * (nv + 1) + lev + 1]
                           * (pt_d[id * nv + lev + 1] - pt_d[id * nv + lev])
                           / (Altitude_d[lev + 1] - Altitude_d[lev])
                       - Rho_int_d[id * (nv + 1) + lev] * KH_d[id * (nv + 1) + lev]
                             * (pt_d[id * nv + lev] - pt_d[id * nv + lev - 1])
                             / (Altitude_d[lev] - Altitude_d[lev - 1]));
            }
        }
        //update surface temperature
        // pt_surf_d[id] += 1.0 / Csurf * Cp * KH_d[id * (nv + 1) + 0] / Altitude_d[0]
        //                  * (Rho_d[id * nv + 0] * pt_d[id * nv + 0] - Rho_surf_d[id] * pt_surf_d[id])
        //                  * time_step;
        // if (pt_surf_d[id] < 0)
        //     pt_surf_d[id] = 0;
        // // Tsurface_d[id] = pow(Rho_surf_d[id] * Rd / P_Ref, Rd / Cv) * pow(pt_surf_d[id], Cp / Cv);
        // Tsurface_d[id] = pt_surf_d[id] * pow(p_surf_d[id] / P_Ref, kappa);
        // double T0 = pow(Rho_d[id * nv + 0] * Rd / P_Ref, Rd / Cv) * pow(pt_d[id * nv + 0], Cp / Cv);

        dTsurf_dt_d[id] +=
            1.0 / Csurf * Cp * KH_d[id * (nv + 1) + 0] / Altitude_d[0]
            * (Rho_d[id * nv + 0] * pt_d[id * nv + 0] - Rho_surf_d[id] * pt_surf_d[id]);

        // printf("Tsurf= %f, T0 = %f, dTsurf = %e\n", Tsurface_d[id], T0, dTsurf_dt_d[id]);
        // Rho_surf_d[id] = p_surf_d[id] / (Rd * Tsurface_d[id]);
        // if (id == 2561) {
        //     printf("T = %f, p = %f, rho = %f, pt = %f\n",
        //            Tsurface_d[id],
        //            p_surf_d[id],
        //            Rho_surf_d[id],
        //            pt_surf_d[id]);
        // }
    }
}

__global__ void Calc_MOlength_Cdrag_BLdepth(double *pressure_d,
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
                                            double  z_therm,
                                            double *RiB_d,
                                            int *   bl_top_lev_d,
                                            double *bl_top_height_d,
                                            double  f_surf_layer,
                                            double *pt_surf_d,
                                            double *p_int_d,
                                            double *CD_d,
                                            double *CH_d,
                                            double *L_MO_d,
                                            double *F_sens_d,
                                            double *vh_lowest_d,
                                            double *Rho_int_d,
                                            int     num,
                                            int     nv) {

    //first get bulk Ri number of first layer
    //convert to zeta and compute L (MO length) for column
    //then zeta for layers above = z/L
    //from zeta, compute drag coefficients, BL depth, add'l quantities

    int id = blockIdx.x * blockDim.x + threadIdx.x;
    int lev;

    if (id < num) {
        double kappa = Rd / Cp;
        double extrap_surf, pt_interface;
        double vh_layer, vh_layer_below, vh_interface;
        double zeta_surf, RiB;

        //extrapolate some quantities to surface
        extrap_surf = -Altitude_d[1] / (Altitude_d[0] - Altitude_d[1]);
        p_int_d[id * (nv + 1) + 0] =
            pressure_d[id * nv + 1] + extrap_surf * (pressure_d[id * nv] - pressure_d[id * nv + 1]);
        // p_int_d[id * (nv + 1)] = p_surf_d[id]; //may want to eliminate p_surf_d
        pt_surf_d[id] = Tsurface_d[id] * pow(p_int_d[id * (nv + 1) + 0] / P_Ref, -kappa);
        Rho_int_d[id * (nv + 1)] =
            Rho_d[id * nv + 1] + extrap_surf * (Rho_d[id * nv] - Rho_d[id * nv + 1]);

        //things we need a center of lowest layer
        pt_d[id * nv] =
            pow(P_Ref, kappa) * pow(pressure_d[id * nv], 1.0 - kappa) / (Rho_d[id * nv] * Rd);
        vh_lowest_d[id] = sqrt((pow(Mh_d[id * nv * 3 + 0], 2) + pow(Mh_d[id * nv * 3 + 1], 2)
                                + pow(Mh_d[id * nv * 3 + 2], 2)))
                          / Rho_d[id * nv];
        vh_layer = vh_lowest_d[id];

        if (pow(vh_lowest_d[id], 2) == 0) { //zero velocity, RiB = large +number
            RiB_d[id] = LARGERiB;
        }
        else { // bulk Richardson number, wrt to surface
            RiB_d[id] = Gravit * Altitude_d[0] * (pt_d[id * nv] - pt_surf_d[id])
                        / (pt_surf_d[id] * pow(vh_lowest_d[id], 2));
        }
        zeta_surf =
            RiB_2_zeta(RiB_d[id], Ri_crit, Altitude_d[0] / z_rough, Altitude_d[0] / z_therm);

        if (RiB_d[id] > Ri_crit) { //super-critical, set L to zero
            L_MO_d[id] = 0.0;
        }
        else {
            L_MO_d[id] = Altitude_d[0] / zeta_surf; //monin-obukhov length for this column
        }

        //surface drag coefficient
        if (RiB_d[id] < 0) { //unstable (temporarily model as neutral)
            CD_d[id] = pow(KVONKARMAN, 2)
                       * pow(int_phi_M_unstable(zeta_surf, Altitude_d[0] / z_rough), -2);
            CH_d[id] = pow(KVONKARMAN, 2)
                       * pow(int_phi_M_unstable(zeta_surf, Altitude_d[0] / z_rough), -1)
                       * pow(int_phi_H_unstable(zeta_surf, Altitude_d[0] / z_therm), -1);
        }
        else if (RiB_d[id] == 0) { //neutral
            CD_d[id] = pow(KVONKARMAN, 2) * pow(log(Altitude_d[0] / z_rough), -2);
            CH_d[id] = pow(KVONKARMAN, 2) * pow(log(Altitude_d[0] / z_rough), -1)
                       * pow(log(Altitude_d[0] / z_therm), -1);
        }
        else if (RiB_d[id] <= Ri_crit) { //stable
            CD_d[id] =
                pow(KVONKARMAN, 2) * pow(log(Altitude_d[0] / z_rough) + zeta_surf / Ri_crit, -2);
            CH_d[id] = pow(KVONKARMAN, 2)
                       * pow(log(Altitude_d[0] / z_rough) + zeta_surf / Ri_crit, -1)
                       * pow(log(Altitude_d[0] / z_therm) + zeta_surf / Ri_crit, -1);
        }
        else { //super-critical
            CD_d[id] = 0;
            CH_d[id] = 0;
        }

        if (isnan(L_MO_d[id])) {
            printf("halt!");
        }

        //sensible heat flux from surface to lowest layer of atmos
        F_sens_d[id] = CH_d[id] * vh_lowest_d[id] * (pt_surf_d[id] - pt_d[id * nv]);

        //loop over layers above to determine height
        RiB              = RiB_d[id]; //temporarily store RiB of layer below
        bl_top_lev_d[id] = 0;
        for (lev = 1; lev <= nv; lev++) {
            if (RiB <= Ri_crit) { //stop if we've crossed Ri_crit at the top of previous layer
                //compute vh and pt of top interface
                vh_layer_below = vh_layer;
                vh_layer       = sqrt((pow(Mh_d[id * nv * 3 + lev * 3 + 0], 2)
                                 + pow(Mh_d[id * nv * 3 + lev * 3 + 1], 2)
                                 + pow(Mh_d[id * nv * 3 + lev * 3 + 2], 2)))
                           / Rho_d[id * nv + lev];
                vh_interface = vh_layer_below
                               + (vh_layer - vh_layer_below)
                                     * (Altitudeh_d[lev] - Altitude_d[lev - 1])
                                     / (Altitude_d[lev] - Altitude_d[lev - 1]);
                pt_d[id * nv + lev] = pow(P_Ref, kappa)
                                      * pow(pressure_d[id * nv + lev], 1.0 - kappa)
                                      / (Rho_d[id * nv + lev] * Rd);
                pt_interface = pt_d[id * nv + lev - 1]
                               + (pt_d[id * nv + lev] - pt_d[id * nv + lev - 1])
                                     * (Altitudeh_d[lev] - Altitude_d[lev - 1])
                                     / (Altitude_d[lev] - Altitude_d[lev - 1]);

                Rho_int_d[id * (nv + 1) + lev] = Rho_d[id * nv + lev - 1]
                                                 + (Rho_d[id * nv + lev] - Rho_d[id * nv + lev - 1])
                                                       * (Altitudeh_d[lev] - Altitude_d[lev - 1])
                                                       / (Altitude_d[lev] - Altitude_d[lev - 1]);
                p_int_d[id * (nv + 1) + lev] =
                    pressure_d[id * nv + lev - 1]
                    + (pressure_d[id * nv + lev] - pressure_d[id * nv + lev - 1])
                          * (Altitudeh_d[lev] - Altitude_d[lev - 1])
                          / (Altitude_d[lev] - Altitude_d[lev - 1]);
                //calc RiB at top interface of layer (lev-1)
                if (pow(vh_interface, 2) == 0) { //zero velocity, set RiB to a big +number
                    RiB = LARGERiB;
                }
                else {
                    if (RiB_d[id] < 0) { //unstable part
                        // later! include effect of thermals here
                        RiB = Gravit * Altitudeh_d[lev] * (pt_interface - pt_d[id * nv + 0])
                              / (pt_d[id * nv + 0] * pow(vh_interface, 2));
                    }
                    else { //stable & neutral
                        RiB = Gravit * Altitudeh_d[lev] * (pt_interface - pt_d[id * nv + 0])
                              / (pt_d[id * nv + 0] * pow(vh_interface, 2));
                    }
                }
                bl_top_lev_d[id] = lev; //layer above this interface is in BL
            }
        }
        bl_top_height_d[id] = Altitudeh_d[bl_top_lev_d[id] + 1];
    }
}

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
                           double  z_therm,
                           double *RiGrad_d,
                           double *pt_surf_d,
                           double *p_int_d,
                           double *CD_d,
                           double *CH_d,
                           double *Rho_int_d,
                           double *KM_d,
                           double *KH_d,
                           int     num,
                           int     nv) {

    //gradient richardson number at intefaces
    //need to iterate on level or not?

    int id  = blockIdx.x * blockDim.x + threadIdx.x;
    int nvi = gridDim.y;
    int lev = blockIdx.y;

    if (id < num) {
        double kappa = Rd / Cp;
        double pt_below, pt_int;
        double shear2, mix_length, asym_len_scale, e_mix;

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
                CD_d[id] =
                    CN * stability_fM_u(CN, RiGrad_d[id * nvi + lev], Altitude_d[lev], z_rough);
                CH_d[id] =
                    CN * stability_fH_u(CN, RiGrad_d[id * nvi + lev], Altitude_d[lev], z_rough);
            }
            else {
                CD_d[id] = CN * stability_f_s(RiGrad_d[id * nvi + lev]);
                CH_d[id] = CN * stability_f_s(RiGrad_d[id * nvi + lev]);
            }
            KM_d[id * nvi + lev] = CD_d[id] * sqrt(shear2) * pow(Altitude_d[lev], 2);
            KH_d[id * nvi + lev] = CH_d[id] * sqrt(shear2) * pow(Altitude_d[lev], 2);
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
            if (Altitudeh_d[lev] <= TRANSITION_HEIGHT) {
                asym_len_scale = ABL_ASYM_LEN;
            }
            else {
                asym_len_scale = FREE_ASYM_LEN
                                 + (ABL_ASYM_LEN - FREE_ASYM_LEN)
                                       * exp(1 - Altitudeh_d[lev] / TRANSITION_HEIGHT);
            }
            mix_length = pow(1.0 / KVONKARMAN / Altitudeh_d[lev] + 1.0 / asym_len_scale, -1.0);
            if (RiGrad_d[id * nvi + lev] < Ri_crit) { // prevent nan value if RiGrad > Ri_crit
                e_mix = pow(mix_length, 2) * shear2 * (1 - RiGrad_d[id * nvi + lev] / Ri_crit);
            }
            else {
                e_mix = 0.0;
            }
            if (e_mix < E_MIN_MIX)
                e_mix = E_MIN_MIX;
            //for now, diffusivities are the same. could modify later for convective conditions
            KH_d[id * nvi + lev] = mix_length * sqrt(e_mix);
            KM_d[id * nvi + lev] = mix_length * sqrt(e_mix);
        }
    }
}


__global__ void CalcRiB(double *pressure_d,
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
                        double  z_therm,
                        double *RiB_d,
                        int *   bl_top_lev_d,
                        double *bl_top_height_d,
                        int *   sl_top_lev_d,
                        double  f_surf_layer,
                        double *pt_surf_d,
                        double *p_surf_d,
                        double *Rho_surf_d,
                        double *CD_d,
                        double *CH_d,
                        double *zeta_d,
                        double *vh_lowest_d,
                        double *Rho_int_d,
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

                //  p_surf       = pressure_d[id * nv + lev];
                p_surf_d[id] = p_surf;

                pt_surf       = Tsurface_d[id] * pow(p_surf / P_Ref, -kappa);
                pt_surf_d[id] = pt_surf;

                Rho_surf_d[id] = Rho_d[id * nv + lev + 1]
                                 + extrap_surf * (Rho_d[id * nv + lev] - Rho_d[id * nv + lev + 1]);

                Rho_int_d[id * (nv + 1) + lev] = Rho_surf_d[id];
                // if (id == 2561) {
                //     printf(
                //         "p = %f, pt = %f, rho = %f\n", p_surf_d[id], pt_surf_d[id], Rho_surf_d[id]);
                // }
                // calculate pt and horizontal velocity of layer
                pt_layer = pow(P_Ref, kappa) * pow(pressure_d[id * nv + lev], 1.0 - kappa)
                           / (Rho_d[id * nv + lev] * Rd);
                pt_lowest           = pt_layer; //will need this later
                pt_d[id * nv + lev] = pt_layer;
                vh_layer            = sqrt((pow(Mh_d[id * nv * 3 + lev * 3 + 0], 2)
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
                    CH_d[id] = pow(KVONKARMAN, 2) * pow(log(Altitude_d[lev] / z_rough), -1)
                               * pow(log(Altitude_d[lev] / z_therm), -1);
                }
                else if (RiB_d[id * (nv + 1) + lev] == 0) { //neutral
                    CD_d[id] = pow(KVONKARMAN, 2) * pow(log(Altitude_d[lev] / z_rough), -2);
                    CH_d[id] = pow(KVONKARMAN, 2) * pow(log(Altitude_d[lev] / z_rough), -1)
                               * pow(log(Altitude_d[lev] / z_therm), -1);
                }
                else if (RiB_d[id * (nv + 1) + lev] <= Ri_crit) { //stable
                    CD_d[id] = pow(KVONKARMAN, 2)
                               * pow(log(Altitude_d[lev] / z_rough) + zeta / Ri_crit, -2);
                    CH_d[id] = pow(KVONKARMAN, 2)
                               * pow(log(Altitude_d[lev] / z_rough) + zeta / Ri_crit, -1)
                               * pow(log(Altitude_d[lev] / z_therm) + zeta / Ri_crit, -1);
                }
                else { //super-critical
                    CD_d[id] = 0;
                    CH_d[id] = 0;
                }
            }
            else if (lev == nv) {
                //what should I do at the top level??
                RiB_d[id * (nv + 1) + lev] = LARGERiB; //top level can't be incorporated into BL?
                // convert to zeta (m-o stability param)
                zeta_d[id * (nv + 1) + lev] = RiB_2_zeta(RiB_d[id * (nv + 1) + lev],
                                                         Ri_crit,
                                                         Altitudeh_d[lev] / z_rough,
                                                         Altitudeh_d[lev] / z_therm);
            }
            else {
                //potential temperatures for this layer, layer below, and interface b/w
                pt_layer_below = pt_layer;
                pt_layer       = pow(P_Ref, kappa) * pow(pressure_d[id * nv + lev], 1.0 - kappa)
                           / (Rho_d[id * nv + lev] * Rd);
                pt_d[id * nv + lev] = pt_layer;
                pt_interface        = pt_layer_below
                               + (pt_layer - pt_layer_below)
                                     * (Altitudeh_d[lev] - Altitude_d[lev - 1])
                                     / (Altitude_d[lev] - Altitude_d[lev - 1]);

                Rho_int_d[id * (nv + 1) + lev] = Rho_d[id * nv + lev - 1]
                                                 + (Rho_d[id * nv + lev] - Rho_d[id * nv + lev - 1])
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
                                                         Altitudeh_d[lev] / z_rough,
                                                         Altitudeh_d[lev] / z_therm);
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
        sl_top_lev_d[id] = 0;
        for (lev = 0; lev <= bl_top_lev_d[id]; lev++) {
            if (Altitudeh_d[lev + 1] <= f_surf_layer * bl_top_height_d[id]) {
                sl_top_lev_d[id] = lev;
            }
        }
    }
}

__device__ double phi_M_unstable(double zeta) {
    //phi_M for unstable conditions
    double gamma = GAMMA_M;
    return pow(1 - gamma * zeta, -1. / 3);
}

__device__ double int_phi_M_unstable(double zeta, double z_z0) {
    //integral of phi_M/zeta for unstable conditions
    double gamma = GAMMA_M;
    double x     = pow(1 - gamma * zeta, 1. / 3);
    double x0    = pow(1 - gamma * zeta / z_z0, 1. / 3);
    return sqrt(3.0) * atan((1 + 2 * x) / sqrt(3.0)) - sqrt(3.0) * atan((1 + 2 * x0) / sqrt(3.0))
           + log((1 - x) / (1 - x0)) - 0.5 * log((1 + x + x * x) / (1 + x0 + x0 * x0));
}

__device__ double phi_H_unstable(double zeta) {
    //phi_H for unstable conditions
    double gamma = GAMMA_H;
    return pow(1 - gamma * zeta, -1. / 2);
}

__device__ double int_phi_H_unstable(double zeta, double z_zT) {
    //integral of phi_H/zeta for unstable conditions
    double gamma = GAMMA_H;
    return -2 * asinh(1. / sqrt(-gamma * zeta)) + 2 * asinh(1. / sqrt(-gamma * zeta / z_zT));
}

__device__ double f_newton_RiB_zeta(double RiB, double zeta, double z_z0, double z_zT) {
    //function relating RiB and zeta. we're finding the root (zeta value) that makes this = 0
    return RiB - zeta * int_phi_H_unstable(zeta, z_zT) / pow(int_phi_M_unstable(zeta, z_z0), 2.0);
}

__device__ double fpr_newton_RiB_zeta(double zeta, double z_z0, double z_zT) {
    //derivative in zeta of the above function
    return -int_phi_H_unstable(zeta, z_zT) / pow(int_phi_M_unstable(zeta, z_z0), 2)
           - phi_H_unstable(zeta) / pow(int_phi_M_unstable(zeta, z_z0), 2)
           + 2 * int_phi_H_unstable(zeta, z_zT) / pow(int_phi_M_unstable(zeta, z_z0), 3)
                 * phi_M_unstable(zeta);
}

__device__ double RiB_2_zeta(double RiB, double Ri_crit, double z_z0, double z_zT) {
    if (RiB < 0) {
        double eps  = 1e-8 * fabs(RiB); //convergence for newton scheme
        double zeta = RiB, zetaf = zeta + 2 * eps;
        // double f, fp;
        int iter = 0;
        while (fabs(zetaf - zeta) > eps && iter < 1000) {
            zeta = zetaf;
            zetaf =
                zeta
                - f_newton_RiB_zeta(RiB, zeta, z_z0, z_zT) / fpr_newton_RiB_zeta(zeta, z_z0, z_zT);
            // f  = f_newton_RiB_zeta(RiB, zeta, z_z0, z_zT);
            // fp = fpr_newton_RiB_zeta(zeta, z_z0, z_zT);
            // printf("RiB = %f, zeta = %f, f = %e, fp = %f\n", RiB, zeta, f, fp);
            // printf("zeta = %f, phiH = %f, phiM = %f, intphiH = %f, intphiM = %f\n",
            //        zeta,
            //        phi_H_unstable(zeta),
            //        phi_M_unstable(zeta),
            //        int_phi_H_unstable(zeta, z_zT),
            //        int_phi_M_unstable(zeta, z_z0));
            iter++;
        }
        if (iter == 1000) {
            printf("exceeded limit in newton iteration, RiB = %e, zeta = %e\n", RiB, zetaf);
        }
        return zetaf;
    }
    else if (RiB == 0) {
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

__global__ void CalcKM_KH(double *RiB_d,
                          double *L_MO_d,
                          double *CD_d,
                          double *CH_d,
                          double *bl_top_height_d,
                          int *   bl_top_lev_d,
                          double *F_sens_d,
                          double *counter_grad_d,
                          double *vh_lowest_d,
                          double *Altitude_d,
                          double *Altitudeh_d,
                          double  Ri_crit,
                          double  z_rough,
                          double  z_therm,
                          double  f_surf_layer,
                          double *pt_surf_d,
                          double *KM_d,
                          double *KH_d,
                          double  Gravit,
                          int     num) {
    int id  = blockIdx.x * blockDim.x + threadIdx.x;
    int nvi = gridDim.y;
    int lev = blockIdx.y;

    if (id < num) {
        if (lev == 0) { // coefficient at surface
            KM_d[id * (nvi) + lev] = CD_d[id] * vh_lowest_d[id] * Altitude_d[lev];
            // this is the only place CH is used?
            KH_d[id * (nvi) + lev]         = CH_d[id] * vh_lowest_d[id] * Altitude_d[lev];
            counter_grad_d[id * nvi + lev] = 0.0;
        }
        else {
            //factor common to all coeffs above first layer
            double kvk_z_scale = KVONKARMAN * Altitudeh_d[lev]
                                 * pow(1.0 - Altitudeh_d[lev] / bl_top_height_d[id], 2);
            if (Altitudeh_d[lev] <= f_surf_layer * bl_top_height_d[id]) { //inner layer
                //condition on surface layer stability as in Frierson
                if (RiB_d[id] < 0) { //unstable
                    KM_d[id * (nvi) + lev] = kvk_z_scale * sqrt(CD_d[id]) * vh_lowest_d[id]
                                             / phi_M_unstable(Altitudeh_d[lev] / L_MO_d[id]);
                    KH_d[id * (nvi) + lev] = kvk_z_scale * sqrt(CD_d[id]) * vh_lowest_d[id]
                                             / phi_H_unstable(Altitudeh_d[lev] / L_MO_d[id]);
                }
                else if (RiB_d[id] == 0) { //neutral
                    KM_d[id * (nvi) + lev] = kvk_z_scale * sqrt(CD_d[id]) * vh_lowest_d[id];
                    KH_d[id * (nvi) + lev] =
                        kvk_z_scale * sqrt(CD_d[id]) * vh_lowest_d[id]; //scales with CD, not CH
                }
                else { //stable
                    if (L_MO_d[id] == 0) {
                        KM_d[id * nvi + lev] = 0.0;
                        KH_d[id * nvi + lev] = 0.0;
                    }
                    else {
                        KM_d[id * (nvi) + lev] = kvk_z_scale * sqrt(CD_d[id]) * vh_lowest_d[id]
                                                 / (1 + Altitudeh_d[lev] / L_MO_d[id] / Ri_crit);
                        KH_d[id * (nvi) + lev] =
                            kvk_z_scale * sqrt(CD_d[id]) * vh_lowest_d[id]
                            / (1
                               + Altitudeh_d[lev] / L_MO_d[id] / Ri_crit); //scales with CD, not CH
                    }
                }
                counter_grad_d[id * nvi + lev] = 0.0;
            }
            else if (Altitudeh_d[lev] <= bl_top_height_d[id]) { //outer layer
                if (RiB_d[id] < 0) {                            //unstable
                    //trickiest part of the whole damn thing
                    double wstar, wm, wt, Pr;
                    //convective velocity scale above inner layer
                    wstar =
                        pow(Gravit / pt_surf_d[id] * F_sens_d[id] * bl_top_height_d[id], 1. / 3);
                    //total velocity scale including surface friction velocity
                    wm = pow(pow(sqrt(CD_d[id]) * vh_lowest_d[id], 3)
                                 + GAMMA_H * f_surf_layer * KVONKARMAN * pow(wstar, 3),
                             1. / 3);
                    //Prandtl number at top of inner layer
                    Pr = phi_H_unstable(f_surf_layer * bl_top_height_d[id] / L_MO_d[id])
                             / phi_M_unstable(f_surf_layer * bl_top_height_d[id] / L_MO_d[id])
                         + A_THERM * KVONKARMAN * f_surf_layer * wstar / wm;
                    //velocity scale for heat mixing
                    wt = wm / Pr;

                    counter_grad_d[id * nvi + lev] =
                        A_THERM * wstar * F_sens_d[id] / (pow(wm, 2) * bl_top_height_d[id]);

                    KM_d[id * nvi + lev] = kvk_z_scale * wm;
                    KH_d[id * nvi + lev] = kvk_z_scale * wt;
                    if (isnan(KH_d[id * nvi + lev]) || isnan(KM_d[id * nvi + lev])) {
                        printf("%f %f %f %f", wstar, wm, Pr, wt);
                    }
                    if (KH_d[id * nvi + lev] > 14000.) {
                        printf("%f %f %f %f %f", wstar, wm, Pr, wt, kvk_z_scale);
                    }
                }
                else if (RiB_d[id] == 0) { //neutral (same as inner layer)
                    KM_d[id * (nvi) + lev] = kvk_z_scale * sqrt(CD_d[id]) * vh_lowest_d[id];
                    KH_d[id * (nvi) + lev] =
                        kvk_z_scale * sqrt(CD_d[id]) * vh_lowest_d[id]; //scales with CD, not CH
                    counter_grad_d[id * nvi + lev] = 0.0;
                }
                else { //stable (same as inner layer)
                    if (L_MO_d[id] == 0) {
                        KM_d[id * nvi + lev] = 0.0;
                        KH_d[id * nvi + lev] = 0.0;
                    }
                    else {
                        KM_d[id * (nvi) + lev] = kvk_z_scale * sqrt(CD_d[id]) * vh_lowest_d[id]
                                                 / (1 + Altitudeh_d[lev] / L_MO_d[id] / Ri_crit);
                        KH_d[id * (nvi) + lev] =
                            kvk_z_scale * sqrt(CD_d[id]) * vh_lowest_d[id]
                            / (1
                               + Altitudeh_d[lev] / L_MO_d[id] / Ri_crit); //scales with CD, not CH
                    }
                    counter_grad_d[id * nvi + lev] = 0.0;
                }
            }
            else {
                KM_d[id * (nvi) + lev]         = 0.0;
                KH_d[id * (nvi) + lev]         = 0.0;
                counter_grad_d[id * nvi + lev] = 0.0;
            }
        }
        if (KH_d[id * nvi + lev] > 14000.) {
            printf("somethings up here");
        }
        if (isnan(KH_d[id * nvi + lev]) || isnan(KM_d[id * nvi + lev])) {
            printf("stopppp");
        }
    }
}

__global__ void CalcKM_KH_old(double *RiB_d,
                              double *zeta_d,
                              double *CD_d,
                              double *CH_d,
                              double *bl_top_height_d,
                              int *   bl_top_lev_d,
                              int *   sl_top_lev_d,
                              double *vh_lowest_d,
                              double *Altitude_d,
                              double *Altitudeh_d,
                              double  Ri_crit,
                              double  z_rough,
                              double  z_therm,
                              double  f_surf_layer,
                              double *KM_d,
                              double *KH_d,
                              int     num) {
    int id  = blockIdx.x * blockDim.x + threadIdx.x;
    int nvi = gridDim.y;
    int lev = blockIdx.y;

    if (id < num) {
        if (lev == 0) { // coefficient at surface
            KM_d[id * (nvi) + lev] = CD_d[id] * vh_lowest_d[id] * Altitude_d[lev];
            // this is the only place CH is used?
            KH_d[id * (nvi) + lev] = CH_d[id] * vh_lowest_d[id] * Altitude_d[lev];
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
                    KH_d[id * (nvi) + lev] = KVONKARMAN * Altitudeh_d[lev] * sqrt(CD_d[id])
                                             * vh_lowest_d[id]; //scales with CD, not CH
                }
                else if (RiB_d[id * (nvi) + 0] == 0) { //neutral
                    KM_d[id * (nvi) + lev] =
                        KVONKARMAN * Altitudeh_d[lev] * sqrt(CD_d[id]) * vh_lowest_d[id];
                    KH_d[id * (nvi) + lev] = KVONKARMAN * Altitudeh_d[lev] * sqrt(CD_d[id])
                                             * vh_lowest_d[id]; //scales with CD, not CH
                }
                else { //stable
                    KM_d[id * (nvi) + lev] = KVONKARMAN * Altitudeh_d[lev] * sqrt(CD_d[id])
                                             * vh_lowest_d[id]
                                             / (1 + zeta_d[id * (nvi) + lev] / Ri_crit);
                    KH_d[id * (nvi) + lev] =
                        KVONKARMAN * Altitudeh_d[lev] * sqrt(CD_d[id]) * vh_lowest_d[id]
                        / (1 + zeta_d[id * (nvi) + lev] / Ri_crit); //scales with CD, not CH
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
                    KH_d[id * (nvi) + lev] = KVONKARMAN * f_surf_layer * bl_top_height_d[id]
                                             * sqrt(CD_d[id]) * vh_lowest_d[id]
                                             * z_scale; //scales with CD, not CH
                }
                else if (RiB_d[id * (nvi) + 0] == 0) { //neutral
                    KM_d[id * (nvi) + lev] = KVONKARMAN * f_surf_layer * bl_top_height_d[id]
                                             * sqrt(CD_d[id]) * vh_lowest_d[id] * z_scale;
                    KH_d[id * (nvi) + lev] = KVONKARMAN * f_surf_layer * bl_top_height_d[id]
                                             * sqrt(CD_d[id]) * vh_lowest_d[id]
                                             * z_scale; //scales with CD, not CH
                }
                else { //stable  (not quite consistent: z = f*h but zeta is at interface? above f*h)
                    KM_d[id * (nvi) + lev] = KVONKARMAN * f_surf_layer * bl_top_height_d[id]
                                             * sqrt(CD_d[id]) * vh_lowest_d[id] * z_scale
                                             / (1
                                                + zeta_d[id * (nvi) + sl_top_lev_d[id] + 1]
                                                      / Ri_crit); ///bl_top_lev_d is wrong!!
                    KH_d[id * (nvi) + lev] = KVONKARMAN * f_surf_layer * bl_top_height_d[id]
                                             * sqrt(CD_d[id]) * vh_lowest_d[id] * z_scale
                                             / (1
                                                + zeta_d[id * (nvi) + sl_top_lev_d[id] + 1]
                                                      / Ri_crit); //scales with CD, not CH
                }
                if (KM_d[id * (nvi) + lev] < 0) {
                    printf("Stop here! %f \n", z_scale);
                }
            }

            else {
                KM_d[id * (nvi) + lev] = 0.0;
                KH_d[id * (nvi) + lev] = 0.0;
            }
        }

        // hack
        // KM_d[id * nvi + lev] = 5e-3;
        // KH_d[id * nvi + lev] = 5e-3;
    }
}
