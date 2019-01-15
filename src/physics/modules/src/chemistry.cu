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
// Method: Chemical kinetics physics module
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
#include "chemistry.h"

#include "thor_chemistry.h"

#include "chemistry_device.h" // Simple chemistry.
#include "chemistry_host.h"

#include "binary_test.h"
#include "debug.h"
#include "directories.h"

chemistry::chemistry() {
}

chemistry::~chemistry() {
}

void chemistry::print_config() {
    log::printf("  Chemistry module\n");
    log::printf("    chem_time_filename: %s.\n", chem_time_filename.c_str());
    log::printf("    fEQ_filename: %s.\n", fEQ_filename.c_str());    
}

bool chemistry::initialise_memory(const ESP &              esp,
                                  device_RK_array_manager &phy_modules_core_arrays) {
    const int num_trace_pts = 7425;
    // TODO: find name
    const int num_a = 135;
    const int num_b = 55;

    // alloc host memory
    coeq_h  = (double *)malloc(num_trace_pts * sizeof(double));
    co2eq_h = (double *)malloc(num_trace_pts * sizeof(double));
    ch4eq_h = (double *)malloc(num_trace_pts * sizeof(double));
    h2oeq_h = (double *)malloc(num_trace_pts * sizeof(double));
    nh3eq_h = (double *)malloc(num_trace_pts * sizeof(double));

    tauco_h  = (double *)malloc(num_trace_pts * sizeof(double));
    tauco2_h = (double *)malloc(num_trace_pts * sizeof(double));
    tauch4_h = (double *)malloc(num_trace_pts * sizeof(double));
    tauh2o_h = (double *)malloc(num_trace_pts * sizeof(double));
    taunh3_h = (double *)malloc(num_trace_pts * sizeof(double));

    P_che_h = (double *)malloc(num_a * sizeof(double));
    T_che_h = (double *)malloc(num_b * sizeof(double));

    tracer_h = (double *)malloc(esp.nv * esp.point_num * ntr * sizeof(double));

    // alloc device memory
    cudaMalloc((void **)&tracer_d, esp.nv * esp.point_num * ntr * sizeof(double));
    cudaMalloc((void **)&tracers_d, esp.nv * esp.point_num * ntr * sizeof(double));
    cudaMalloc((void **)&tracerk_d, esp.nv * esp.point_num * ntr * sizeof(double));

    cudaMalloc((void **)&coeq_d, num_trace_pts * sizeof(double));
    cudaMalloc((void **)&co2eq_d, num_trace_pts * sizeof(double));
    cudaMalloc((void **)&ch4eq_d, num_trace_pts * sizeof(double));
    cudaMalloc((void **)&h2oeq_d, num_trace_pts * sizeof(double));
    cudaMalloc((void **)&nh3eq_d, num_trace_pts * sizeof(double));

    cudaMalloc((void **)&tauco_d, num_trace_pts * sizeof(double));
    cudaMalloc((void **)&tauco2_d, num_trace_pts * sizeof(double));
    cudaMalloc((void **)&tauch4_d, num_trace_pts * sizeof(double));
    cudaMalloc((void **)&tauh2o_d, num_trace_pts * sizeof(double));
    cudaMalloc((void **)&taunh3_d, num_trace_pts * sizeof(double));

    cudaMalloc((void **)&P_che_d, num_a * sizeof(double));
    cudaMalloc((void **)&T_che_d, num_b * sizeof(double));

    cudaMalloc((void **)&difftr_d, esp.nv * esp.point_num * ntr * sizeof(double));

    // register the arrays that need to be updated by RK kernel
    phy_modules_core_arrays.register_array(tracers_d, tracerk_d, tracer_d, ntr);

#ifdef BENCHMARKING
    std::map<std::string, output_def> defs =
        {
            {"tracer_d", {tracer_d, esp.nv * esp.point_num * ntr, "RK tracer", "ti", true}},
            {"tracers_d", {tracers_d, esp.nv * esp.point_num * ntr, "RK tracers", "ts", true}},
            {"tracerk_d", {tracerk_d, esp.nv * esp.point_num * ntr, "RK tracerk", "tk", true}},
            {"difftr_d", {difftr_d, esp.nv * esp.point_num * ntr, "Diffusion Tr", "dftr", true}}
        };

    std::vector<std::string> output_vars = {"tracer_d", "tracers_d", "tracerk_d", "difftr_d"};


    binary_test::get_instance().register_phy_modules_variables(defs, std::vector<std::string>(), output_vars);

#endif // BENCHMARING
    return true;
}


bool chemistry::free_memory() {
    free(tauch4_h);
    free(tauco_h);
    free(tauh2o_h);
    free(tauco2_h);
    free(taunh3_h);

    free(ch4eq_h);
    free(coeq_h);
    free(h2oeq_h);
    free(co2eq_h);
    free(nh3eq_h);

    free(P_che_h);
    free(T_che_h);

    free(tracer_h);

    cudaFree(ch4eq_d);
    cudaFree(coeq_d);
    cudaFree(h2oeq_d);
    cudaFree(co2eq_d);
    cudaFree(nh3eq_d);

    cudaFree(tauch4_d);
    cudaFree(tauco_d);
    cudaFree(tauh2o_d);
    cudaFree(tauco2_d);
    cudaFree(taunh3_d);

    cudaFree(tracer_d);
    cudaFree(tracers_d);
    cudaFree(tracerk_d);

    cudaFree(P_che_d);
    cudaFree(T_che_d);

    cudaFree(difftr_d);

    return true;
}

bool chemistry::initial_conditions(const ESP &            esp,
                                   const SimulationSetup &sim) {
    // Input for chemistry
    FILE * infile1;
    int    NT = 55;
    int    NP = 135;
    double dummy;



    
    log::printf("Chemistry: Loading fEQ file: %s.\n", fEQ_filename.c_str());
    if (!path_exists(fEQ_filename))
    {
        log::printf("\nfEQ input file %s does not exist.\n", fEQ_filename.c_str());
        exit(EXIT_FAILURE);
    }
    
    infile1 = fopen(fEQ_filename.c_str(), "r");
    if (infile1 == NULL) {
        log::printf("\nUnable to open fEQ input file %s.\n", fEQ_filename.c_str());
        return false;
    }
    
    for (int i = 0; i < NT; i++) {
        for (int j = 0; j < NP; j++) {
            if (fscanf(infile1,
                       "%lf %lf %lf %lf %lf %lf %lf",
                       &T_che_h[i],
                       &P_che_h[j],
                       &ch4eq_h[j * NT + i],
                       &coeq_h[j * NT + i],
                       &h2oeq_h[j * NT + i],
                       &co2eq_h[j * NT + i],
                       &nh3eq_h[j * NT + i])
                != 7) {
                log::printf("error parsing chem time file: %s.\n", fEQ_filename.c_str());
                fclose(infile1);
                return false;
            }
        }
    }
    
    fclose(infile1);

    log::printf("Chemistry: Loading chem time file: %s.\n", chem_time_filename.c_str());
    if (!path_exists(chem_time_filename))
    {
        log::printf("\nchem time input file %s does not exist.\n", chem_time_filename.c_str());
        exit(EXIT_FAILURE);
    }
    
    infile1 = fopen(chem_time_filename.c_str(), "r");
    if (infile1 == NULL) {
        log::printf("\nUnable to open chem time input file %s.\n", chem_time_filename.c_str());
        exit(EXIT_FAILURE);
    }
    
    for (int i = 0; i < NT; i++) {
        for (int j = 0; j < NP; j++) {
            if (fscanf(infile1,
                       "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                       &T_che_h[i],
                       &P_che_h[j],
                       &tauch4_h[j * NT + i],
                       &tauco_h[j * NT + i],
                       &dummy,
                       &dummy,
                       &tauh2o_h[j * NT + i],
                       &tauco2_h[j * NT + i],
                       &taunh3_h[j * NT + i],
                       &dummy)
                != 10) {
                log::printf("error parsing fEQ file: %s.\n", chem_time_filename.c_str());
                fclose(infile1);
                return false;
            }
        }
    }

    for (int j = 0; j < NP; j++) P_che_h[j] = log(P_che_h[j]);
    fclose(infile1);


    for (int lev = 0; lev < esp.nv; lev++) {
        for (int i = 0; i < esp.point_num; i++) {
            // CH4
            tracer_h[i * esp.nv * ntr + lev * ntr + 0] = Compute_tracer_host(
                                                             ch4eq_h,
                                                             P_che_h,
                                                             T_che_h,
                                                             esp.temperature_h[i * esp.nv + lev],
                                                             esp.pressure_h[i * esp.nv + lev])
                                                         * esp.Rho_h[i * esp.nv + lev];
            // CO
            tracer_h[i * esp.nv * ntr + lev * ntr + 1] = Compute_tracer_host(
                                                             coeq_h,
                                                             P_che_h,
                                                             T_che_h,
                                                             esp.temperature_h[i * esp.nv + lev],
                                                             esp.pressure_h[i * esp.nv + lev])
                                                         * esp.Rho_h[i * esp.nv + lev];
            // H2O
            tracer_h[i * esp.nv * ntr + lev * ntr + 2] = Compute_tracer_host(h2oeq_h,
                                                                             P_che_h,
                                                                             T_che_h,
                                                                             esp.temperature_h[i * esp.nv + lev],
                                                                             esp.pressure_h[i * esp.nv + lev])
                                                         * esp.Rho_h[i * esp.nv + lev];
            // CO2
            tracer_h[i * esp.nv * ntr + lev * ntr + 3] = Compute_tracer_host(co2eq_h,
                                                                             P_che_h,
                                                                             T_che_h,
                                                                             esp.temperature_h[i * esp.nv + lev],
                                                                             esp.pressure_h[i * esp.nv + lev])
                                                         * esp.Rho_h[i * esp.nv + lev];
            // NH3
            tracer_h[i * esp.nv * ntr + lev * ntr + 4] = Compute_tracer_host(nh3eq_h,
                                                                             P_che_h,
                                                                             T_che_h,
                                                                             esp.temperature_h[i * esp.nv + lev],
                                                                             esp.pressure_h[i * esp.nv + lev])
                                                         * esp.Rho_h[i * esp.nv + lev];
        }
    }


    cudaMemcpy(coeq_d, coeq_h, 7425 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(ch4eq_d, ch4eq_h, 7425 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(h2oeq_d, h2oeq_h, 7425 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(co2eq_d, co2eq_h, 7425 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(nh3eq_d, nh3eq_h, 7425 * sizeof(double), cudaMemcpyHostToDevice);

    cudaMemcpy(tauco_d, tauco_h, 7425 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(tauch4_d, tauch4_h, 7425 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(tauh2o_d, tauh2o_h, 7425 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(tauco2_d, tauco2_h, 7425 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(taunh3_d, taunh3_h, 7425 * sizeof(double), cudaMemcpyHostToDevice);

    cudaMemcpy(P_che_d, P_che_h, 135 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(T_che_d, T_che_h, 55 * sizeof(double), cudaMemcpyHostToDevice);

    cudaMemcpy(tracer_d, tracer_h, esp.point_num * esp.nv * ntr * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemset(tracers_d, 0, sizeof(double) * esp.nv * esp.point_num * ntr);
    cudaMemset(tracerk_d, 0, sizeof(double) * esp.nv * esp.point_num * ntr);


    return true;
}

bool chemistry::dyn_core_loop_init(const ESP &esp) {
    cudaMemcpy(tracerk_d, tracer_d, esp.point_num * esp.nv * ntr * sizeof(double), cudaMemcpyDeviceToDevice);

    cudaMemset(tracers_d, 0, sizeof(double) * esp.point_num * esp.nv * ntr);

    return true;
}

bool chemistry::dyn_core_loop_slow_modes(const ESP &            esp,
                                         const SimulationSetup &sim,
                                         int                    nstep, // Step number
                                         double                 time_step) {           // Time-step [s]
    const int LN = 16;                                                 // Size of the inner region side.
    dim3      NT(esp.nl_region, esp.nl_region, 1);                     // Number of threads in a block.
    dim3      NBPT(2, 1, ntr);                                         // Number of blocks. (POLES)
    dim3      NBTR(esp.nr, esp.nv, ntr);                               // Number of blocks in the diffusion routine for tracers.
    dim3      NBTRP(2, esp.nv, ntr);                                   // Number of blocks in the diffusion routine for tracers. (POLES)

    if (sim.HyDiff) {
        // Tracers
        // TODO: check: where should this be set to 0 ?
        cudaMemset(esp.diff_d, 0, sizeof(double) * 6 * esp.point_num * esp.nv);
        cudaDeviceSynchronize();
        Tracer_Eq_Diffusion<LN, LN><<<NBTR, NT>>>(difftr_d,
                                                  esp.diff_d,
                                                  tracerk_d,
                                                  esp.Rhok_d,
                                                  esp.areasTr_d,
                                                  esp.nvecoa_d,
                                                  esp.nvecti_d,
                                                  esp.nvecte_d,
                                                  esp.Kdh4_d,
                                                  esp.Altitude_d,
                                                  sim.A,
                                                  esp.maps_d,
                                                  ntr, //
                                                  esp.nl_region,
                                                  0,
                                                  sim.DeepModel);

        Tracer_Eq_Diffusion_Poles<5><<<NBTRP, 1>>>(difftr_d,
                                                   esp.diff_d,
                                                   tracerk_d,
                                                   esp.Rhok_d,
                                                   esp.areasTr_d,
                                                   esp.nvecoa_d,
                                                   esp.nvecti_d,
                                                   esp.nvecte_d,
                                                   esp.Kdh4_d,
                                                   esp.Altitude_d,
                                                   esp.Altitudeh_d,
                                                   sim.A,
                                                   esp.point_local_d,
                                                   ntr,
                                                   esp.point_num,
                                                   0,
                                                   sim.DeepModel);
        cudaDeviceSynchronize();
        Tracer_Eq_Diffusion<LN, LN><<<NBTR, NT>>>(difftr_d,
                                                  esp.diff_d,
                                                  tracerk_d,
                                                  esp.Rhok_d,
                                                  esp.areasTr_d,
                                                  esp.nvecoa_d,
                                                  esp.nvecti_d,
                                                  esp.nvecte_d,
                                                  esp.Kdh4_d,
                                                  esp.Altitude_d,
                                                  sim.A,
                                                  esp.maps_d,
                                                  ntr, //
                                                  esp.nl_region,
                                                  1,
                                                  sim.DeepModel);

        Tracer_Eq_Diffusion_Poles<5><<<NBTRP, 1>>>(difftr_d,
                                                   esp.diff_d,
                                                   tracerk_d,
                                                   esp.Rhok_d,
                                                   esp.areasTr_d,
                                                   esp.nvecoa_d,
                                                   esp.nvecti_d,
                                                   esp.nvecte_d,
                                                   esp.Kdh4_d,
                                                   esp.Altitude_d,
                                                   esp.Altitudeh_d,
                                                   sim.A,
                                                   esp.point_local_d,
                                                   ntr,
                                                   esp.point_num,
                                                   1,
                                                   sim.DeepModel);
    }


    return true;
}

bool chemistry::dyn_core_loop_fast_modes(const ESP &            esp,
                                         const SimulationSetup &sim,
                                         int                    nstep, // Step number
                                         double                 times) {               // Time-step [s]
    const int LN = 16;                                                 // Size of the inner region side.
    dim3      NT(esp.nl_region, esp.nl_region, 1);                     // Number of threads in a block.
    dim3      NBPT(2, 1, ntr);                                         // Number of blocks. (POLES)
    dim3      NBTR(esp.nr, esp.nv, ntr);                               // Number of blocks in the diffusion routine for tracers.
    dim3      NBTRP(2, esp.nv, ntr);                                   // Number of blocks in the diffusion routine for tracers. (POLES)


    //
    // Tracer equation.
    cudaDeviceSynchronize();
    Tracer_Eq<LN, LN><<<NBTR, NT>>>(tracers_d,
                                    tracerk_d,
                                    esp.Rhos_d,
                                    esp.Rhok_d,
                                    esp.Mhs_d,
                                    esp.Mhk_d,
                                    esp.Whs_d,
                                    esp.Whk_d,
                                    difftr_d,
                                    esp.div_d,
                                    esp.Altitude_d,
                                    esp.Altitudeh_d,
                                    sim.A,
                                    times,
                                    esp.maps_d,
                                    ntr,
                                    esp.nl_region,
                                    sim.DeepModel);

    Tracer_Eq_Poles<6><<<NBPT, 1>>>(tracers_d,
                                    tracerk_d,
                                    esp.Rhos_d,
                                    esp.Rhok_d,
                                    esp.Mhs_d,
                                    esp.Mhk_d,
                                    esp.Whs_d,
                                    esp.Whk_d,
                                    difftr_d,
                                    esp.div_d,
                                    esp.Altitude_d,
                                    esp.Altitudeh_d,
                                    sim.A,
                                    times,
                                    esp.point_local_d,
                                    ntr,
                                    esp.point_num,
                                    esp.nv,
                                    sim.DeepModel);

    return true;
}


bool chemistry::dyn_core_loop_end(const ESP &esp) {
    cudaMemcpy(tracer_d, tracerk_d, esp.point_num * esp.nv * ntr * sizeof(double), cudaMemcpyDeviceToDevice);


    return true;
}


bool chemistry::phy_loop(ESP &                  esp,
                         const SimulationSetup &sim,
                         int                    nstep, // Step number
                         double                 time_step) {           // Time-step [s]
    const int NTH = 256;
    dim3      NBTR((esp.point_num / NTH) + 1, esp.nv, ntr);

    //
    ////////////////////////
    // Simple chemistry
    cudaDeviceSynchronize();
    Tracers_relax_chemistry_co2<<<NBTR, NTH>>>(tracer_d,
                                               tauch4_d,
                                               tauco_d,
                                               tauh2o_d,
                                               tauco2_d,
                                               taunh3_d,
                                               ch4eq_d,
                                               coeq_d,
                                               h2oeq_d,
                                               co2eq_d,
                                               nh3eq_d,
                                               P_che_d,
                                               T_che_d,
                                               esp.temperature_d,
                                               esp.pressure_d,
                                               esp.Rho_d,
                                               time_step,
                                               ntr,
                                               esp.point_num);
    cudaDeviceSynchronize();
    Tracers_relax_chemistry<<<NBTR, NTH>>>(tracer_d,
                                           tauch4_d,
                                           tauco_d,
                                           tauh2o_d,
                                           tauco2_d,
                                           taunh3_d,
                                           ch4eq_d,
                                           coeq_d,
                                           h2oeq_d,
                                           co2eq_d,
                                           nh3eq_d,
                                           P_che_d,
                                           T_che_d,
                                           esp.temperature_d,
                                           esp.pressure_d,
                                           esp.Rho_d,
                                           time_step,
                                           ntr,
                                           esp.point_num);


    return true;
}

bool chemistry::configure(config_file &config_reader) {

    config_reader.append_config_var("chem_time_file", chem_time_filename, std::string("NoFilename.txt"));
    config_reader.append_config_var("chem_fEQ_file", fEQ_filename, std::string("NoFilename.txt"));
    return true;
}

bool chemistry::store(const ESP &esp,
                      storage &  s) {
    cudaMemcpy(tracer_h, tracer_d, esp.point_num * esp.nv * ntr * sizeof(double), cudaMemcpyDeviceToHost);
    s.append_table(tracer_h,
                   esp.nv * esp.point_num * ntr,
                   "/tracer",
                   " ",
                   "Volume mixing ratio");

    return true;
}

bool chemistry::store_init(storage &s) {
    return true;
}
