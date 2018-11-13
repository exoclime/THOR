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
// Build the class ESP (Exoclimes Simulation Platform)
//
//
// Description:
//   Declare and initialize variables in the model
//
// Method: -
//
//
// Known limitations: None.
//
//
// Known issues: None.
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

#include "../headers/phy/chemistry_host.h"
#include "../headers/phy/valkyrie_conservation.h"
#include "../headers/phy/valkyrie_jet_steadystate.h"
#include "directories.h"
#include "esp.h"
#include "hdf5.h"
#include "storage.h"
#include <map>
#include <stdio.h>

// physical modules
#include "phy_modules.h"

__host__ ESP::ESP(int *           point_local_,
                  int *           maps_,
                  double *        lonlat_,
                  double *        Altitude_,
                  double *        Altitudeh_,
                  double *        nvecoa_,
                  double *        nvecti_,
                  double *        nvecte_,
                  double *        areasT_,
                  double *        areasTr_,
                  double *        div_,
                  double *        grad_,
                  double *        func_r_,
                  int             nl_region_,
                  int             nr_,
                  int             nv_,
                  int             nvi_,
                  int             glevel_,
                  bool            spring_dynamics_,
                  double          spring_beta_,
                  int             nlat_,
                  int             ntr_,
                  int *           zonal_mean_tab,
                  double          Rv_sponge_,
                  double          ns_sponge_,
                  double          t_shrink_,
                  int             point_num_,
                  bool            conservation,
                  benchmark_types core_benchmark_,
                  log_writer &    logwriter_) :
    nl_region(nl_region_),
    nr(nr_),
    point_num(point_num_),
    nv(nv_),
    nvi(nvi_),
    nlat(nlat_),
    ntr(ntr_),
    glevel(glevel_),
    spring_dynamics(spring_dynamics_),
    spring_beta(spring_beta_),
    logwriter(logwriter_),
    core_benchmark(core_benchmark_) {

    point_local_h = point_local_;
    maps_h        = maps_;

    lonlat_h = lonlat_;

    Altitude_h  = Altitude_;
    Altitudeh_h = Altitudeh_;

    nvecoa_h  = nvecoa_;
    nvecti_h  = nvecti_;
    nvecte_h  = nvecte_;
    areasTr_h = areasTr_;
    areasT_h  = areasT_;

    div_h  = div_;
    grad_h = grad_;

    func_r_h = func_r_;

    zonal_mean_tab_h = zonal_mean_tab;

    Rv_sponge = Rv_sponge_;
    ns_sponge = ns_sponge_;
    t_shrink  = t_shrink_;

    //
    //  Allocate Data
    if (core_benchmark != NO_BENCHMARK)
        alloc_data(conservation);
}

__host__ void ESP::alloc_data(bool conservation) {


    //
    //  Description:
    //
    //  Allocate data on host and device.
    //
    //  Allocate data in host
    //  Diagnostics
    Rho_h         = (double *)malloc(nv * point_num * sizeof(double));
    pressure_h    = (double *)malloc(nv * point_num * sizeof(double));
    temperature_h = (double *)malloc(nv * point_num * sizeof(double));
    Mh_h          = (double *)malloc(nv * point_num * 3 * sizeof(double));
    W_h           = (double *)malloc(nv * point_num * sizeof(double));
    Wh_h          = (double *)malloc(nvi * point_num * sizeof(double));

    if (conservation == true) {
        Etotal_h  = (double *)malloc(nv * point_num * sizeof(double));
        Mass_h    = (double *)malloc(nv * point_num * sizeof(double));
        AngMomx_h = (double *)malloc(nv * point_num * sizeof(double));
        AngMomy_h = (double *)malloc(nv * point_num * sizeof(double));
        AngMomz_h = (double *)malloc(nv * point_num * sizeof(double));
    }

    coeq_h  = (double *)malloc(7425 * sizeof(double));
    co2eq_h = (double *)malloc(7425 * sizeof(double));
    ch4eq_h = (double *)malloc(7425 * sizeof(double));
    h2oeq_h = (double *)malloc(7425 * sizeof(double));
    nh3eq_h = (double *)malloc(7425 * sizeof(double));

    tauco_h  = (double *)malloc(7425 * sizeof(double));
    tauco2_h = (double *)malloc(7425 * sizeof(double));
    tauch4_h = (double *)malloc(7425 * sizeof(double));
    tauh2o_h = (double *)malloc(7425 * sizeof(double));
    taunh3_h = (double *)malloc(7425 * sizeof(double));

    P_che_h = (double *)malloc(135 * sizeof(double));
    T_che_h = (double *)malloc(55 * sizeof(double));

    tracer_h = (double *)malloc(nv * point_num * ntr * sizeof(double));

    //  Allocate data in device
    //  Grid
    cudaMalloc((void **)&point_local_d, 6 * point_num * sizeof(int));
    cudaMalloc((void **)&maps_d, (nl_region + 2) * (nl_region + 2) * nr * sizeof(int));

    //  Operators
    cudaMalloc((void **)&nvecoa_d, 6 * 3 * point_num * sizeof(double));
    cudaMalloc((void **)&nvecti_d, 6 * 3 * point_num * sizeof(double));
    cudaMalloc((void **)&nvecte_d, 6 * 3 * point_num * sizeof(double));
    cudaMalloc((void **)&areasT_d, point_num * sizeof(double));
    cudaMalloc((void **)&areasTr_d, 6 * point_num * sizeof(double));
    cudaMalloc((void **)&func_r_d, 3 * point_num * sizeof(double));
    cudaMalloc((void **)&div_d, 7 * 3 * point_num * sizeof(double));
    cudaMalloc((void **)&grad_d, 7 * 3 * point_num * sizeof(double));

    //  Altitude (grid)
    cudaMalloc((void **)&Altitude_d, nv * sizeof(double));
    cudaMalloc((void **)&Altitudeh_d, nvi * sizeof(double));

    //  Longitude-latitude
    cudaMalloc((void **)&lonlat_d, 2 * point_num * sizeof(double));

    //  Diagnostics
    cudaMalloc((void **)&Mh_d, nv * point_num * 3 * sizeof(double));
    cudaMalloc((void **)&W_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&Wh_d, nvi * point_num * sizeof(double));
    cudaMalloc((void **)&Rho_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&pressure_d, nv * point_num * sizeof(double));

    cudaMalloc((void **)&tracer_d, nv * point_num * ntr * sizeof(double));
    cudaMalloc((void **)&tracers_d, nv * point_num * ntr * sizeof(double));
    cudaMalloc((void **)&tracerk_d, nv * point_num * ntr * sizeof(double));

    cudaMalloc((void **)&coeq_d, 7425 * sizeof(double));
    cudaMalloc((void **)&co2eq_d, 7425 * sizeof(double));
    cudaMalloc((void **)&ch4eq_d, 7425 * sizeof(double));
    cudaMalloc((void **)&h2oeq_d, 7425 * sizeof(double));
    cudaMalloc((void **)&nh3eq_d, 7425 * sizeof(double));

    cudaMalloc((void **)&tauco_d, 7425 * sizeof(double));
    cudaMalloc((void **)&tauco2_d, 7425 * sizeof(double));
    cudaMalloc((void **)&tauch4_d, 7425 * sizeof(double));
    cudaMalloc((void **)&tauh2o_d, 7425 * sizeof(double));
    cudaMalloc((void **)&taunh3_d, 7425 * sizeof(double));

    cudaMalloc((void **)&P_che_d, 135 * sizeof(double));
    cudaMalloc((void **)&T_che_d, 55 * sizeof(double));

    //  Temperature
    cudaMalloc((void **)&temperature_d, nv * point_num * sizeof(double));

    //  Potential temperature
    cudaMalloc((void **)&pt_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&pth_d, nvi * point_num * sizeof(double));

    //  Entalphy
    cudaMalloc((void **)&h_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&hh_d, nvi * point_num * sizeof(double));

    //  Advection
    cudaMalloc((void **)&Adv_d, nv * point_num * 3 * sizeof(double));

    //  3D vector
    cudaMalloc((void **)&v_d, nv * point_num * 3 * sizeof(double));

    //  Effective gravity
    cudaMalloc((void **)&gtil_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&gtilh_d, nvi * point_num * sizeof(double));

    //  Slow modes
    cudaMalloc((void **)&SlowMh_d, nv * point_num * 3 * sizeof(double));
    cudaMalloc((void **)&SlowWh_d, nvi * point_num * sizeof(double));
    cudaMalloc((void **)&SlowRho_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&Slowpressure_d, nv * point_num * sizeof(double));


    //  Deviations
    cudaMalloc((void **)&pressures_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&Rhos_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&Mhs_d, nv * point_num * 3 * sizeof(double));
    cudaMalloc((void **)&Ws_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&Whs_d, nvi * point_num * sizeof(double));


    //  RK-Method
    cudaMalloc((void **)&pressurek_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&Rhok_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&Mhk_d, nv * point_num * 3 * sizeof(double));
    cudaMalloc((void **)&Wk_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&Whk_d, nvi * point_num * sizeof(double));

    //  Vertical integration
    cudaMalloc((void **)&Sp_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&Sd_d, nv * point_num * sizeof(double));

    //  Diffusion
    cudaMalloc((void **)&Kdhz_d, nv * sizeof(double));
    cudaMalloc((void **)&Kdh4_d, nv * sizeof(double));
    cudaMalloc((void **)&DivM_d, nv * point_num * 3 * sizeof(double));
    cudaMalloc((void **)&diffpr_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&diffmh_d, 3 * nv * point_num * sizeof(double));
    cudaMalloc((void **)&diffw_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&diffrh_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&diff_d, 6 * nv * point_num * sizeof(double));
    cudaMalloc((void **)&divg_Mh_d, 3 * nv * point_num * sizeof(double));
    cudaMalloc((void **)&difftr_d, nv * point_num * ntr * sizeof(double));

    //  Extras-nan
    cudaMalloc((void **)&check_d, sizeof(bool));

    cudaMalloc((void **)&vbar_d, 3 * nv * point_num * sizeof(double));
    cudaMalloc((void **)&zonal_mean_tab_d, 2 * point_num * sizeof(int));

    if (conservation == true) {
        //  Conservation quantities
        cudaMalloc((void **)&Etotal_d, nv * point_num * sizeof(double));
        cudaMalloc((void **)&Mass_d, nv * point_num * sizeof(double));
        cudaMalloc((void **)&AngMomx_d, nv * point_num * sizeof(double));
        cudaMalloc((void **)&AngMomy_d, nv * point_num * sizeof(double));
        cudaMalloc((void **)&AngMomz_d, nv * point_num * sizeof(double));
        cudaMalloc((void **)&GlobalE_d, 1 * sizeof(double));
        cudaMalloc((void **)&GlobalMass_d, 1 * sizeof(double));
        cudaMalloc((void **)&GlobalAMx_d, 1 * sizeof(double));
        cudaMalloc((void **)&GlobalAMy_d, 1 * sizeof(double));
        cudaMalloc((void **)&GlobalAMz_d, 1 * sizeof(double));
    }
    // PHY modules
    if (core_benchmark != NO_BENCHMARK)
        phy_modules_init_mem(*this);
}

__host__ bool ESP::initial_values(bool               rest,
                                  const std::string &initial_conditions_filename,
                                  const bool &       continue_sim,
                                  double             timestep_dyn,
                                  XPlanet &          xplanet,
                                  double             kb,
                                  double             mu,
                                  bool               sponge,
                                  bool               DeepModel,
                                  int                TPprof,
                                  int                chemistry,
                                  int &              nstep,
                                  double &           simulation_start_time,
                                  int &              output_file_idx,
                                  bool               conservation) {

    output_file_idx = 0;
    nstep           = 0;

    // Store some general configs
    planet         = xplanet;

    //  Set initial conditions.
    //
    //
    //  Initial atmospheric conditions
    if (rest) {
        double Ha = planet.Rd * planet.Tmean / planet.Gravit;
        for (int i = 0; i < point_num; i++) {
            //
            //          Initial conditions for an isothermal Atmosphere
            //

            for (int lev = 0; lev < nv; lev++) {
                pressure_h[i * nv + lev] = planet.P_Ref * exp(-Altitude_h[lev] / Ha);
                if (TPprof == 0) {
                    temperature_h[i * nv + lev] = planet.Tmean;
                }
                else if (TPprof == 1) {
                    double tau                  = pressure_h[i * nv + lev] / (1e4); //tau = 1 at 0.1 bar
                    double gamma                = 0.6;                              // ratio of sw to lw opacity
                    double f                    = 0.25;
                    temperature_h[i * nv + lev] = pow(3 * planet.Tmean * planet.Tmean * planet.Tmean * planet.Tmean * f * (2 / 3 + 1 / (gamma * sqrt(3)) + (gamma / sqrt(3) - 1 / (gamma * sqrt(3))) * exp(-gamma * tau * sqrt(3))), 0.25);
                }
                if (core_benchmark == HS_DEEP_HOT_JUPITER) {
                    double Ptil = 0.0;
                    if (pressure_h[i * nv + lev] >= 1e5) {
                        Ptil = log10(pressure_h[i * nv + lev] / 100000);
                    }
                    temperature_h[i * nv + lev] = 1696.6986 + 132.2318 * Ptil - 174.30459 * Ptil * Ptil
                                                  + 12.579612 * Ptil * Ptil * Ptil + 59.513639 * Ptil * Ptil * Ptil * Ptil
                                                  + 9.6706522 * Ptil * Ptil * Ptil * Ptil * Ptil
                                                  - 4.1136048 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil
                                                  - 1.0632301 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil
                                                  + 0.064400203 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil
                                                  + 0.035974396 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil
                                                  + 0.0025740066 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil;
                }
            }

            for (int lev = 0; lev < nv; lev++) {
                //              Density [kg/m3]
                Rho_h[i * nv + lev] = pressure_h[i * nv + lev] / (temperature_h[i * nv + lev] * planet.Rd);

                //              Momentum [kg/m3 m/s]
                Mh_h[i * 3 * nv + 3 * lev + 0] = 0.0;
                Mh_h[i * 3 * nv + 3 * lev + 1] = 0.0;
                Mh_h[i * 3 * nv + 3 * lev + 2] = 0.0;

                //              Vertical momentum [kg/m3 m/s]
                W_h[i * nv + lev]        = 0.0; // Center of the layer.
                Wh_h[i * (nv + 1) + lev] = 0.0; // Layers interface.
            }
            Wh_h[i * (nv + 1) + nv] = 0.0;
        }
        if (core_benchmark == JET_STEADY) {
            //  Number of threads per block.
            const int NTH = 256;

            //  Specify the block sizes.
            dim3 NB((point_num / NTH) + 1, nv, 1);

            cudaMemcpy(Altitude_d, Altitude_h, nv * sizeof(double), cudaMemcpyHostToDevice);
            cudaMemcpy(pressure_d, pressure_h, point_num * nv * sizeof(double), cudaMemcpyHostToDevice);
            cudaMemcpy(Mh_d, Mh_h, 3 * point_num * nv * sizeof(double), cudaMemcpyHostToDevice);
            cudaMemcpy(Rho_d, Rho_h, point_num * nv * sizeof(double), cudaMemcpyHostToDevice);
            cudaMemcpy(temperature_d, temperature_h, point_num * nv * sizeof(double), cudaMemcpyHostToDevice);
            cudaMemcpy(lonlat_d, lonlat_h, 2 * point_num * sizeof(double), cudaMemcpyHostToDevice);
            setup_jet<<<NB, NTH>>>(Mh_d,
                                   // setup_jet <<< 1, 1 >>>  (Mh_d,
                                   pressure_d,
                                   Rho_d,
                                   temperature_d,
                                   planet.Cp,
                                   planet.Rd,
                                   planet.Omega,
                                   planet.A,
                                   Altitude_d,
                                   lonlat_d,
                                   point_num);

            cudaMemcpy(Mh_h, Mh_d, 3 * point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
            cudaMemcpy(temperature_h, temperature_d, point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
            cudaMemcpy(pressure_h, pressure_d, point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
            cudaMemcpy(Rho_h, Rho_d, point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
        }

        simulation_start_time = 0.0;
    }
    else {
        bool load_OK = true;
        // build planet filename
        string planet_filename;

        path   p(initial_conditions_filename);
        int    file_number = 0;
        string basename    = "";

        string parent_path = p.parent();

        if (continue_sim) {
            if (!match_output_file_numbering_scheme(initial_conditions_filename,
                                                    basename,
                                                    file_number)) {
                printf("Loading initial conditions: "
                       "Could not recognise file numbering scheme "
                       "for input %s: (found base: %s, num: %d) \n",
                       initial_conditions_filename.c_str(),
                       basename.c_str(),
                       file_number);
                return false;
            }

            output_file_idx = file_number;

            planet_filename = p.parent() + "/esp_output_planet_" + basename + ".h5";
        }
        else {
            planet_filename = p.parent() + "/" + p.stem() + "_planet.h5";
        }

        // check existence of files
        if (!path_exists(initial_conditions_filename)) {
            printf("initial condition file %s not found.\n", initial_conditions_filename.c_str());
            return false;
        }

        if (!path_exists(planet_filename)) {
            printf("planet_file %s not found.\n", planet_filename.c_str());
            return false;
        }


        printf("Loading planet from: %s\n", planet_filename.c_str());
        printf("Loading initial conditions from: %s\n", initial_conditions_filename.c_str());

        // Check planet data
        {
            // values to check agains variable
            map<string, double> mapValues;

            mapValues["/A"]            = planet.A;
            mapValues["/Top_altitude"] = planet.Top_altitude;
            mapValues["/glevel"]       = glevel;
            mapValues["/vlevel"]       = nv;

            hid_t file_id;
            file_id = H5Fopen(planet_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

            bool values_match = true;

            for (const std::pair<std::string, double> &element : mapValues) {
                double value = 0.0;
                load_OK &= load_double_value_from_h5file(file_id, element.first, value);

                if (value != element.second) {
                    printf("mismatch for %s value between config value: %f and initial condition value %f.\n",
                           element.first.c_str(),
                           element.second,
                           value);
                    values_match = false;
                }
            }

            H5Fclose(file_id);

            if (load_OK == false || values_match == false) {
                printf("Could not reload full configuration.\n");

                return false;
            }
        }


        //      Restart from an existing simulation.
        {

            // Load atmospheric data
            hid_t file_id;
            file_id = H5Fopen(initial_conditions_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
            // Step number
            load_OK &= load_int_value_from_h5file(file_id, "/nstep", nstep);
            //      Density
            load_OK &= load_double_table_from_h5file(file_id, "/Rho", Rho_h, point_num * nv);

            //      Pressure
            load_OK &= load_double_table_from_h5file(file_id, "/Pressure", pressure_h, point_num * nv);

            //      Horizontal momentum
            load_OK &= load_double_table_from_h5file(file_id, "/Mh", Mh_h, point_num * nv * 3);
            //      Vertical momentum
            load_OK &= load_double_table_from_h5file(file_id, "/Wh", Wh_h, point_num * nvi);

            //      Simulation start time
            load_OK &= load_double_value_from_h5file(file_id, "/simulation_time", simulation_start_time);
            H5Fclose(file_id);
        }


        if (!load_OK)
            return false;

        for (int i = 0; i < point_num; i++)
            for (int lev = 0; lev < nv; lev++)
                temperature_h[i * nv + lev] = pressure_h[i * nv + lev] / (planet.Rd * Rho_h[i * nv + lev]);

        for (int i = 0; i < point_num; i++) {
            for (int lev = 0; lev < nv; lev++) {
                double xi   = Altitude_h[lev];
                double xim1 = Altitudeh_h[lev];
                double xip1 = Altitudeh_h[lev + 1];

                double a = (xi - xip1) / (xim1 - xip1);
                double b = (xi - xim1) / (xip1 - xim1);

                W_h[i * nv + lev] = Wh_h[i * (nv + 1) + lev] * a + Wh_h[i * (nv + 1) + lev + 1] * b;
            }
        }
    }
#ifdef BENCHMARKING
    // recompute temperature from pressure and density, to have correct rounding for binary comparison
    for (int i = 0; i < point_num; i++)
        for (int lev = 0; lev < nv; lev++)
            temperature_h[i * nv + lev] = pressure_h[i * nv + lev] / (planet.Rd * Rho_h[i * nv + lev]);
#endif // BENCHMARKING

    //  Diffusion
    //  Horizontal
    double *Kdhz_h, *Kdh4_h;
    Kdhz_h = new double[nv];
    Kdh4_h = new double[nv];
    for (int lev = 0; lev < nv; lev++) {
        //      Diffusion constant.
        double dbar = sqrt(2 * M_PI / 5) * planet.A / (pow(2, glevel));
        Kdh4_h[lev] = planet.Diffc * pow(dbar, 4.) / timestep_dyn;
        Kdhz_h[lev] = planet.Diffc * pow(dbar, 4.) / timestep_dyn;
    }

    // Input for chemistry
    FILE * infile1;
    int    NT = 55;
    int    NP = 135;
    double dummy;
    if (chemistry == 1) {
        infile1 = fopen("ifile/solar_fEQ_THOR.txt", "r");
        if (infile1 == NULL) {
            printf("\nUnable to open input file.\n");
            exit(EXIT_FAILURE);
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
                    printf("error parsing ifile/solar_fEQ_THOR.txt\n");
                    fclose(infile1);
                    return false;
                }
            }
        }


        fclose(infile1);

        infile1 = fopen("ifile/solar_chem_time.txt", "r");
        if (infile1 == NULL) {
            printf("\nUnable to open input file.\n");
            return false;
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
                    printf("error parsing ifile/solar_chem_time.txt\n");
                    fclose(infile1);
                    return false;
                }
            }
        }

        for (int j = 0; j < NP; j++) P_che_h[j] = log(P_che_h[j]);
        fclose(infile1);

        // CH4
        for (int lev = 0; lev < nv; lev++) {
            for (int i = 0; i < point_num; i++) {
                tracer_h[i * nv * ntr + lev * ntr + 0] = Compute_tracer_host(ch4eq_h,
                                                                             P_che_h,
                                                                             T_che_h,
                                                                             temperature_h[i * nv + lev],
                                                                             pressure_h[i * nv + lev])
                                                         * Rho_h[i * nv + lev];
            }
        }
        // CO
        for (int lev = 0; lev < nv; lev++) {
            for (int i = 0; i < point_num; i++) {
                tracer_h[i * nv * ntr + lev * ntr + 1] = Compute_tracer_host(coeq_h,
                                                                             P_che_h,
                                                                             T_che_h,
                                                                             temperature_h[i * nv + lev],
                                                                             pressure_h[i * nv + lev])
                                                         * Rho_h[i * nv + lev];
            }
        }
        // H2O
        for (int lev = 0; lev < nv; lev++) {
            for (int i = 0; i < point_num; i++) {
                tracer_h[i * nv * ntr + lev * ntr + 2] = Compute_tracer_host(h2oeq_h,
                                                                             P_che_h,
                                                                             T_che_h,
                                                                             temperature_h[i * nv + lev],
                                                                             pressure_h[i * nv + lev])
                                                         * Rho_h[i * nv + lev];
            }
        }
        // CO2
        for (int lev = 0; lev < nv; lev++) {
            for (int i = 0; i < point_num; i++) {
                tracer_h[i * nv * ntr + lev * ntr + 3] = Compute_tracer_host(co2eq_h,
                                                                             P_che_h,
                                                                             T_che_h,
                                                                             temperature_h[i * nv + lev],
                                                                             pressure_h[i * nv + lev])
                                                         * Rho_h[i * nv + lev];
            }
        }
        // NH3
        for (int lev = 0; lev < nv; lev++) {
            for (int i = 0; i < point_num; i++) {
                tracer_h[i * nv * ntr + lev * ntr + 4] = Compute_tracer_host(nh3eq_h,
                                                                             P_che_h,
                                                                             T_che_h,
                                                                             temperature_h[i * nv + lev],
                                                                             pressure_h[i * nv + lev])
                                                         * Rho_h[i * nv + lev];
            }
        }
    }

    //  Copy memory to the devide
    cudaMemcpy(point_local_d, point_local_h, 6 * point_num * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(maps_d, maps_h, (nl_region + 2) * (nl_region + 2) * nr * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(Altitude_d, Altitude_h, nv * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(Altitudeh_d, Altitudeh_h, nvi * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(nvecoa_d, nvecoa_h, 6 * 3 * point_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(nvecti_d, nvecti_h, 6 * 3 * point_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(nvecte_d, nvecte_h, 6 * 3 * point_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(areasTr_d, areasTr_h, 6 * point_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(areasT_d, areasT_h, point_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(lonlat_d, lonlat_h, 2 * point_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(func_r_d, func_r_h, 3 * point_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(temperature_d, temperature_h, point_num * nv * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(Mh_d, Mh_h, point_num * nv * 3 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(W_d, W_h, point_num * nv * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(Wh_d, Wh_h, point_num * nvi * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(Rho_d, Rho_h, point_num * nv * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(pressure_d, pressure_h, point_num * nv * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(div_d, div_h, 7 * 3 * point_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(grad_d, grad_h, 7 * 3 * point_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(Kdhz_d, Kdhz_h, nv * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(Kdh4_d, Kdh4_h, nv * sizeof(double), cudaMemcpyHostToDevice);

    if (sponge == true)
        cudaMemcpy(zonal_mean_tab_d, zonal_mean_tab_h, 2 * point_num * sizeof(int), cudaMemcpyHostToDevice);

    if (chemistry == 1) {
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

        cudaMemcpy(tracer_d, tracer_h, point_num * nv * ntr * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemset(tracers_d, 0, sizeof(double) * nv * point_num * ntr);
        cudaMemset(tracerk_d, 0, sizeof(double) * nv * point_num * ntr);
    }

    //  Initialize arrays
    cudaMemset(Adv_d, 0, sizeof(double) * 3 * point_num * nv);
    cudaMemset(v_d, 0, sizeof(double) * nv * point_num * 3);
    cudaMemset(pt_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(pth_d, 0, sizeof(double) * nvi * point_num);
    cudaMemset(SlowMh_d, 0, sizeof(double) * nv * point_num * 3);
    cudaMemset(SlowWh_d, 0, sizeof(double) * nvi * point_num);
    cudaMemset(SlowRho_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(Slowpressure_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(h_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(hh_d, 0, sizeof(double) * nvi * point_num);
    cudaMemset(Rhos_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(pressures_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(Mhs_d, 0, sizeof(double) * nv * point_num * 3);
    cudaMemset(Ws_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(Whs_d, 0, sizeof(double) * nvi * point_num);
    cudaMemset(gtil_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(gtilh_d, 0, sizeof(double) * nvi * point_num);
    cudaMemset(Rhok_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(pressurek_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(Mhk_d, 0, sizeof(double) * nv * point_num * 3);
    cudaMemset(Wk_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(Whk_d, 0, sizeof(double) * nvi * point_num);
    cudaMemset(Sp_d, 0, sizeof(double) * point_num * nv);
    cudaMemset(Sd_d, 0, sizeof(double) * point_num * nv);
    cudaMemset(DivM_d, 0, sizeof(double) * point_num * 3 * nv);
    cudaMemset(diffpr_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(diffmh_d, 0, sizeof(double) * 3 * nv * point_num);
    cudaMemset(diffw_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(diffrh_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(diff_d, 0, sizeof(double) * 6 * nv * point_num);
    cudaMemset(divg_Mh_d, 0, sizeof(double) * 3 * nv * point_num);
    cudaMemset(difftr_d, 0, sizeof(double) * nv * point_num * ntr);

    delete[] Kdh4_h;
    delete[] Kdhz_h;

    return true;
}

__host__ ESP::~ESP() {

    //
    //  Description: Frees the memory space.
    //
    //  Host
    free(point_local_h);
    free(maps_h);
    free(lonlat_h);
    free(Altitude_h);
    free(Altitudeh_h);
    free(nvecoa_h);
    free(nvecti_h);
    free(nvecte_h);
    free(areasTr_h);
    free(div_h);
    free(grad_h);
    free(func_r_h);
    free(Rho_h);
    free(pressure_h);
    free(temperature_h);
    free(Mh_h);
    free(W_h);
    free(Wh_h);

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

    //  Device
    cudaFree(point_local_d);
    cudaFree(maps_d);
    cudaFree(Altitude_d);
    cudaFree(Altitudeh_d);
    cudaFree(nvecoa_d);
    cudaFree(nvecti_d);
    cudaFree(nvecte_d);
    cudaFree(areasT_d);
    cudaFree(areasTr_d);
    cudaFree(lonlat_d);
    cudaFree(div_d);
    cudaFree(grad_d);
    cudaFree(func_r_d);
    cudaFree(Rho_d);
    cudaFree(pressure_d);
    cudaFree(temperature_d);
    cudaFree(W_d);
    cudaFree(Wh_d);
    cudaFree(h_d);
    cudaFree(hh_d);
    cudaFree(Adv_d);
    cudaFree(gtil_d);
    cudaFree(gtilh_d);
    cudaFree(v_d);
    cudaFree(pt_d);
    cudaFree(pth_d);
    cudaFree(SlowMh_d);
    cudaFree(SlowWh_d);
    cudaFree(SlowRho_d);
    cudaFree(Slowpressure_d);
    cudaFree(Rhok_d);
    cudaFree(pressurek_d);
    cudaFree(Mhk_d);
    cudaFree(Whk_d);
    cudaFree(Wk_d);
    cudaFree(Rhos_d);
    cudaFree(pressures_d);
    cudaFree(Mhs_d);
    cudaFree(Whs_d);
    cudaFree(Ws_d);

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

    cudaFree(Sd_d);
    cudaFree(Sp_d);
    cudaFree(Kdhz_d);
    cudaFree(Kdh4_d);
    cudaFree(DivM_d);
    cudaFree(diffpr_d);
    cudaFree(diffmh_d);
    cudaFree(diffw_d);
    cudaFree(diffrh_d);
    cudaFree(diff_d);
    cudaFree(difftr_d);
    cudaFree(divg_Mh_d);

    //  Conservation quantities
    cudaFree(Etotal_d);
    cudaFree(Mass_d);
    cudaFree(AngMomx_d);
    cudaFree(AngMomy_d);
    cudaFree(AngMomz_d);
    cudaFree(GlobalE_d);
    cudaFree(GlobalMass_d);
    cudaFree(GlobalAMx_d);
    cudaFree(GlobalAMy_d);
    cudaFree(GlobalAMz_d);

    if (core_benchmark != NO_BENCHMARK)
        phy_modules_free_mem();


    printf("\n\n Free memory!\n\n");
}
