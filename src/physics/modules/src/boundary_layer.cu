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

    return true;
}


bool boundary_layer::free_memory() {

    return true;
}

bool boundary_layer::initial_conditions(const ESP &esp, const SimulationSetup &sim) {
    bool config_OK = true;

    bl_type = RAYLEIGHHS;
    if (bl_type_str == "RayleighHS") {
        bl_type = RAYLEIGHHS;
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

    BLSetup(bl_type, surf_drag_config, bl_sigma_config);

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

    rayleighHS<<<NB, NTH>>>(esp.Mh_d,
                            esp.pressure_d,
                            esp.Rho_d,
                            esp.Altitude_d,
                            surf_drag,
                            bl_sigma,
                            sim.Gravit,
                            time_step,
                            esp.point_num);

    return true;
}

bool boundary_layer::configure(config_file &config_reader) {

    config_reader.append_config_var("bl_type", bl_type_str, string(bl_type_default)); //

    // coefficient of drag strength
    config_reader.append_config_var("surf_drag", surf_drag_config, surf_drag_config);

    // percent of surface pressure where bl starts
    config_reader.append_config_var("bl_sigma", bl_sigma_config, bl_sigma_config);

    return true;
}

bool boundary_layer::store(const ESP &esp, storage &s) {

    return true;
}

bool boundary_layer::store_init(storage &s) {
    s.append_value(surf_drag, "/surf_drag", "1/s", "surface drag coefficient");
    s.append_value(bl_sigma, "/bl_sigma", " ", "boundary layer sigma coordinate");

    return true;
}

void boundary_layer::BLSetup(int bl_type, double surf_drag_, double bl_sigma_) {
    if (bl_type == RAYLEIGHHS) {
        surf_drag = surf_drag_;
        bl_sigma  = bl_sigma_;
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
    }
}
