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
// Method: Chemistry physics module
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
//
// History:
// Version Date       Comment
// ======= ====       =======
//
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////

#pragma once

#include "phy_module_base.h"

class chemistry : public phy_module_base
{
public:
    chemistry();
    ~chemistry();

    bool initialise_memory(const ESP &esp, device_RK_array_manager &phy_modules_core_arrays);

    bool initial_conditions(const ESP &esp, const SimulationSetup &sim, storage *s);

    bool dyn_core_loop_init(const ESP &esp);

    bool dyn_core_loop_slow_modes(const ESP &            esp,
                                  const SimulationSetup &sim,
                                  int                    nstep, // Step number
                                  double                 time_step);            // Time-step [s]

    bool dyn_core_loop_fast_modes(const ESP &            esp,
                                  const SimulationSetup &sim,
                                  int                    nstep, // Step number
                                  double                 times);                // Time-step [s]

    bool dyn_core_loop_end(const ESP &esp);

    bool phy_loop(ESP &                  esp,
                  const SimulationSetup &sim,
                  kernel_diagnostics &   diag,
                  int                    nstep, // Step number
                  double                 time_step);            // Time-step [s]

    bool store(const ESP &esp, storage &s);

    bool store_init(storage &s);

    bool configure(config_file &config_reader);

    virtual bool free_memory();

    void print_config();

private:
    int ntr = 5;

    // host array
    double *tauch4_h;
    double *tauco_h;
    double *tauh2o_h;
    double *tauco2_h;
    double *taunh3_h;

    double *ch4eq_h;
    double *coeq_h;
    double *h2oeq_h;
    double *co2eq_h;
    double *nh3eq_h;

    double *P_che_h;
    double *T_che_h;

    double *tracer_h;

    // device arrays
    double *tauch4_d;
    double *tauco_d;
    double *tauh2o_d;
    double *tauco2_d;
    double *taunh3_d;

    double *ch4eq_d;
    double *coeq_d;
    double *h2oeq_d;
    double *co2eq_d;
    double *nh3eq_d;

    double *tracer_d;
    double *tracers_d;
    double *tracerk_d;

    double *P_che_d;
    double *T_che_d;

    double *difftr_d;

    // filename parameters
    std::string chem_time_filename;
    std::string fEQ_filename;
};
