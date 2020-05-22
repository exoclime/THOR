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
// Description: binary correctness test of output, enabled with compile time switches
//
//
//
// Method: [1] - Dumps output to binary file on a flag
//         [2] - Reads data from binary files on a flag and compare to
//               dumped output
//
//
//
// Known limitations: None.
//
// Known issues: None.
//
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

#include "binary_test.h"

#ifdef BENCHMARKING
#    include "debug_helpers.h"
#    include "directories.h"
#    include "esp.h"
#    include "grid.h"
#    include "log_writer.h"
#    include <functional>
#    include <iomanip>
#    include <iostream>
#    include <sstream>
#    include <vector>

using namespace std;


// define all the variables we can check in the debugging functions
map<string, output_def> build_definitions(ESP& esp, Icogrid& grid) {

    map<string, output_def> out = {
        // {"map name, {variable pointer, table size, name, short name, on device, function call back for loc}}
        {"Rho_d",
         {esp.Rho_d,
          esp.nv * esp.point_num,
          "Density",
          "rho",
          true,
          std::bind(
              &ESP::index_to_location_scalar, &esp, std::placeholders::_1, std::placeholders::_2)}},

        {"pressure_d",
         {esp.pressure_d,
          esp.nv * esp.point_num,
          "Pressure",
          "P",
          true,
          std::bind(
              &ESP::index_to_location_scalar, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"Mh_d",
         {esp.Mh_d,
          esp.nv * esp.point_num * 3,
          "Horizontal Momentum",
          "Mh",
          true,
          std::bind(
              &ESP::index_to_location_vector, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"Wh_d",
         {esp.Wh_d,
          esp.nvi * esp.point_num,
          "Vertical Momentum",
          "Wh",
          true,
          std::bind(&ESP::index_to_location_scalar_mid,
                    &esp,
                    std::placeholders::_1,
                    std::placeholders::_2)}},
        {"temperature_d",
         {esp.temperature_d,
          esp.nv * esp.point_num,
          "Temperature",
          "T",
          true,
          std::bind(
              &ESP::index_to_location_scalar, &esp, std::placeholders::_1, std::placeholders::_2)}},

        {"W_d",
         {esp.W_d,
          esp.nv * esp.point_num,
          "W vert momentum",
          "W",
          true,
          std::bind(
              &ESP::index_to_location_scalar, &esp, std::placeholders::_1, std::placeholders::_2)}},

        {"pt_d",
         {esp.pt_d,
          esp.nv * esp.point_num,
          "Potential Temp pt",
          "pt",
          true,
          std::bind(
              &ESP::index_to_location_scalar, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"pth_d",
         {esp.pth_d,
          esp.nvi * esp.point_num,
          "Potential Temp pth",
          "pth",
          true,
          std::bind(&ESP::index_to_location_scalar_mid,
                    &esp,
                    std::placeholders::_1,
                    std::placeholders::_2)}},


        {"h_d",
         {esp.h_d,
          esp.nv * esp.point_num,
          "Enthalpy h",
          "h",
          true,
          std::bind(
              &ESP::index_to_location_scalar, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"hh_d",
         {esp.hh_d,
          esp.nvi * esp.point_num,
          "Enthalpy hh",
          "hh",
          true,
          std::bind(&ESP::index_to_location_scalar_mid,
                    &esp,
                    std::placeholders::_1,
                    std::placeholders::_2)}},

        {"pressures_d",
         {esp.pressures_d,
          esp.nv * esp.point_num,
          "RK pressures",
          "ps",
          true,
          std::bind(
              &ESP::index_to_location_scalar, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"pressurek_d",
         {esp.pressurek_d,
          esp.nv * esp.point_num,
          "RK pressurek",
          "pk",
          true,
          std::bind(
              &ESP::index_to_location_scalar, &esp, std::placeholders::_1, std::placeholders::_2)}},


        // slow modes
        {"SlowRho_d",
         {esp.SlowRho_d,
          esp.nv * esp.point_num,
          "Slow Density",
          "Srho",
          true,
          std::bind(
              &ESP::index_to_location_scalar, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"SlowMh_d",
         {esp.SlowMh_d,
          esp.nv * esp.point_num * 3,
          "Slow Mh",
          "SMh",
          true,
          std::bind(
              &ESP::index_to_location_vector, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"SlowWh_d",
         {esp.SlowWh_d,
          esp.nvi * esp.point_num,
          "Slow Wh",
          "SWh",
          true,
          std::bind(&ESP::index_to_location_scalar_mid,
                    &esp,
                    std::placeholders::_1,
                    std::placeholders::_2)}},
        {"Slowpressure_d",
         {esp.Slowpressure_d,
          esp.nv * esp.point_num,
          "Slow Pressure",
          "SP",
          true,
          std::bind(
              &ESP::index_to_location_scalar, &esp, std::placeholders::_1, std::placeholders::_2)}},

        // Diffusion
        {"diffrh_d",
         {esp.diffrh_d,
          esp.nv * esp.point_num,
          "Diff Rho",
          "dfrh",
          true,
          std::bind(
              &ESP::index_to_location_scalar, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"diffpr_d",
         {esp.diffpr_d,
          esp.nv * esp.point_num,
          "Diff Press",
          "dfpr",
          true,
          std::bind(
              &ESP::index_to_location_scalar, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"diffmh_d",
         {esp.diffmh_d,
          esp.nv * esp.point_num * 3,
          "Diff Mom",
          "dfmh",
          true,
          std::bind(
              &ESP::index_to_location_vector, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"diffw_d",
         {esp.diffw_d,
          esp.nv * esp.point_num,
          "Diff VMom",
          "dfw",
          true,
          std::bind(
              &ESP::index_to_location_scalar, &esp, std::placeholders::_1, std::placeholders::_2)}},

        {"diff_d",
         {esp.diff_d,
          6 * esp.nv * esp.point_num,
          "Diff",
          "dif",
          true,
          std::bind(&ESP::dummy, &esp, std::placeholders::_1, std::placeholders::_2)}},

        {"diffrv_d",
         {esp.diffrv_d,
          esp.nv * esp.point_num,
          "Diff Rho Vert",
          "dfrv",
          true,
          std::bind(
              &ESP::index_to_location_scalar, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"diffprv_d",
         {esp.diffprv_d,
          esp.nv * esp.point_num,
          "Diff Press Vert",
          "dfprv",
          true,
          std::bind(
              &ESP::index_to_location_scalar, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"diffmv_d",
         {esp.diffmv_d,
          esp.nv * esp.point_num * 3,
          "Diff Mom Vert",
          "dfmv",
          true,
          std::bind(
              &ESP::index_to_location_vector, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"diffwv_d",
         {esp.diffwv_d,
          esp.nv * esp.point_num,
          "Diff VMom Vert",
          "dfwv",
          true,
          std::bind(
              &ESP::index_to_location_scalar, &esp, std::placeholders::_1, std::placeholders::_2)}},

        {"diffv_d",
         {esp.diffv_d2,
          6 * esp.nv * esp.point_num,
          "Diff Vert",
          "difv",
          true,
          std::bind(&ESP::dummy, &esp, std::placeholders::_1, std::placeholders::_2)}},


        {"DivM_d",
         {esp.DivM_d,
          esp.nv * esp.point_num * 3,
          "DivM",
          "dM",
          true,
          std::bind(
              &ESP::index_to_location_vector, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"divg_Mh_d",
         {esp.divg_Mh_d,
          esp.nv * esp.point_num * 3,
          "divg_Mh_d",
          "dgMh",
          true,
          std::bind(
              &ESP::index_to_location_vector, &esp, std::placeholders::_1, std::placeholders::_2)}},

        {"gtil_d",
         {esp.gtil_d,
          esp.nv * esp.point_num,
          "gtil_d",
          "gt",
          true,
          std::bind(
              &ESP::index_to_location_scalar, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"gtilh_d",
         {esp.gtilh_d,
          esp.nvi * esp.point_num,
          "gtilh_d",
          "gth",
          true,
          std::bind(&ESP::index_to_location_scalar_mid,
                    &esp,
                    std::placeholders::_1,
                    std::placeholders::_2)}},

        {"Rhos_d",
         {esp.Rhos_d,
          esp.nv * esp.point_num,
          "RK Rhos",
          "rhos",
          true,
          std::bind(
              &ESP::index_to_location_scalar, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"Mhs_d",
         {esp.Mhs_d,
          esp.nv * esp.point_num * 3,
          "RK Mhs",
          "Mhs",
          true,
          std::bind(
              &ESP::index_to_location_vector, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"Ws_d",
         {esp.Ws_d,
          esp.nv * esp.point_num,
          "RK Ws",
          "Ws",
          true,
          std::bind(
              &ESP::index_to_location_scalar, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"Whs_d",
         {esp.Whs_d,
          esp.nvi * esp.point_num,
          "RK Whs",
          "Whs",
          true,
          std::bind(&ESP::index_to_location_scalar_mid,
                    &esp,
                    std::placeholders::_1,
                    std::placeholders::_2)}},

        {"Rhok_d",
         {esp.Rhok_d,
          esp.nv * esp.point_num,
          "RK Rhok",
          "Rhok",
          true,
          std::bind(
              &ESP::index_to_location_scalar, &esp, std::placeholders::_1, std::placeholders::_2)}},

        {"Mhk_d",
         {esp.Mhk_d,
          esp.nv * esp.point_num * 3,
          "RK Mhk",
          "Mhk",
          true,
          std::bind(
              &ESP::index_to_location_vector, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"Wk_d",
         {esp.Wk_d,
          esp.nv * esp.point_num,
          "RK Wk",
          "Wk",
          true,
          std::bind(
              &ESP::index_to_location_scalar, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"Whk_d",
         {esp.Whk_d,
          esp.nvi * esp.point_num,
          "RK Whk",
          "Whk",
          true,
          std::bind(&ESP::index_to_location_scalar_mid,
                    &esp,
                    std::placeholders::_1,
                    std::placeholders::_2)}},
        // local variables

        {"Adv_d",
         {esp.Adv_d,
          esp.nv * esp.point_num * 3,
          "Advection",
          "Adv",
          true,
          std::bind(
              &ESP::index_to_location_vector, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"v_d",
         {esp.v_d,
          esp.nv * esp.point_num * 3,
          "Velocity",
          "v",
          true,
          std::bind(
              &ESP::index_to_location_vector, &esp, std::placeholders::_1, std::placeholders::_2)}},

        // grid
        {"point_xyz",
         {grid.point_xyz,
          grid.point_num * 3,
          "xyz",
          "x",
          false,
          std::bind(
              &ESP::index_to_location_3xn, &esp, std::placeholders::_1, std::placeholders::_2)}},

        {"point_xyzq",
         {grid.point_xyzq,
          grid.point_num * 3 * 6,
          "xyzq",
          "xq",
          false,
          std::bind(
              &ESP::index_to_location_3x6xn, &esp, std::placeholders::_1, std::placeholders::_2)}},

        /*
              disable int tables for now
              {"pent_ind",      { grid.pent_ind, 12, "pent_ind", "pi", false}},

              {"point_local",   { grid.point_local, 6*grid.point_num, "point_local", "pl", false}},
              {"halo",          { grid.halo, grid.nh, "halo", "halo", false}},
              {"maps",          { grid.maps, (grid.nl_region+2)*(grid.nl_region+2)*grid.nr, "maps", "m", false}},
            */

        {"func_r",
         {grid.func_r,
          3 * grid.point_num,
          "func_r",
          "f",
          false,
          std::bind(
              &ESP::index_to_location_3xn, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"areas",
         {grid.areas,
          6 * 3 * grid.point_num,
          "areas",
          "a",
          false,
          std::bind(
              &ESP::index_to_location_3x6xn, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"areasTr",
         {grid.areasTr,
          6 * grid.point_num,
          "areasTr",
          "aTr",
          false,
          std::bind(
              &ESP::index_to_location_6xn, &esp, std::placeholders::_1, std::placeholders::_2)}},

        {"areasT",
         {grid.areasT,
          grid.point_num,
          "areasT",
          "aT",
          false,
          std::bind(
              &ESP::index_to_location_1xn, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"nvec",
         {grid.nvec,
          6 * 3 * grid.point_num,
          "nvec",
          "nc",
          false,
          std::bind(
              &ESP::index_to_location_3x6xn, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"nvecoa",
         {grid.nvecoa,
          6 * 3 * grid.point_num,
          "nvecoa",
          "na",
          false,
          std::bind(
              &ESP::index_to_location_3x6xn, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"nvecti",
         {grid.nvecti,
          6 * 3 * grid.point_num,
          "nvecti",
          "nti",
          false,
          std::bind(
              &ESP::index_to_location_3x6xn, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"nvecte",
         {grid.nvecte,
          6 * 3 * grid.point_num,
          "nvecte",
          "nte",
          false,
          std::bind(
              &ESP::index_to_location_3x6xn, &esp, std::placeholders::_1, std::placeholders::_2)}},

        {"Altitude",
         {grid.Altitude,
          grid.nv,
          "Altitude",
          "Alt",
          false,
          std::bind(&ESP::dummy, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"Altitudeh",
         {grid.Altitudeh,
          grid.nvi,
          "Altitudeh",
          "Alth",
          false,
          std::bind(&ESP::dummy, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"lonlat",
         {grid.lonlat,
          2 * grid.point_num,
          "lonlat",
          "ll",
          false,
          std::bind(
              &ESP::index_to_location_2xn, &esp, std::placeholders::_1, std::placeholders::_2)}},

        {"Sp_d",
         {esp.Sp_d,
          esp.nv * esp.point_num,
          "Sp_d",
          "Sp",
          true,
          std::bind(
              &ESP::index_to_location_scalar, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"Sd_d",
         {esp.Sd_d,
          esp.nv * esp.point_num,
          "Sd_d",
          "Sd",
          true,
          std::bind(
              &ESP::index_to_location_scalar, &esp, std::placeholders::_1, std::placeholders::_2)}},

        {"div",
         {grid.div,
          7 * 3 * grid.point_num,
          "div",
          "d",
          false,
          std::bind(
              &ESP::index_to_location_3x7xn, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"grad",
         {grid.grad,
          7 * 3 * grid.point_num,
          "grad",
          "g",
          false,
          std::bind(
              &ESP::index_to_location_3x7xn, &esp, std::placeholders::_1, std::placeholders::_2)}},
        {"Qheat",
         {esp.profx_Qheat_d,
          esp.nv * esp.point_num,
          "Qheat",
          "qh",
          true,
          std::bind(&ESP::index_to_location_scalar,
                    &esp,
                    std::placeholders::_1,
                    std::placeholders::_2)}}};

    return out;
}


binary_test::binary_test(string output_dir_, string output_base_name_) :
    output_dir(output_dir_),
    output_base_name(output_base_name_),
    nan_check_d(nullptr) {
    create_output_dir(output_dir);
}
binary_test::~binary_test() {
    if (nan_check_d != nullptr)
        deinit_device_mem_check(nan_check_d);
}
// generic data test function
void binary_test::check_data(const string&         iteration,
                             const string&         ref_name,
                             const vector<string>& input_vars,
                             const vector<string>& output_vars,
                             const bool&           trace_phy_modules) {
#    ifdef BENCH_CHECK_LAST_CUDA_ERROR
    check_last_cuda_error(ref_name);
#    endif // BENCH_CHECK_LAST_CUDA_ERROR

    // load definitions for variables
    vector<output_def> data_output;
    for (auto& name : output_vars) {
        auto&& it = output_definitions.find(name);
        if (it != output_definitions.end()) {
            data_output.push_back(it->second);
        }
    }

    if (trace_phy_modules) {
        for (auto& name : phy_modules_data_out) {
            auto&& it = output_definitions.find(name);
            if (it != output_definitions.end()) {
                data_output.push_back(it->second);
            }
        }
    }


// Specific debugging functions
#    ifdef BENCH_POINT_WRITE
    output_reference(iteration, ref_name, data_output);
#    endif // BENCH_POINT_WRITE

#    ifdef BENCH_POINT_COMPARE
    compare_to_reference(iteration, ref_name, data_output);

#    endif // BENCH_POINT_COMPARE

#    ifdef BENCH_NAN_CHECK
    if (nan_check_d == nullptr) {
        nan_check_d = init_device_mem_check(nan_check_d);
    }

    bool out;
    out = check_nan(iteration, ref_name, data_output);
    if (!out) {
        exit(-1);
    }
#    endif // BENCH_NAN_CHECK
}
// Check for NaNs
// Returns true if data is ok (contains no NaNs!)
bool binary_test::check_nan(const string&             iteration,
                            const string&             ref_name,
                            const vector<output_def>& data_output) {
    bool               out = true;
    std::ostringstream oss;
    oss << std::setw(6) << iteration << " NaNChk: " << std::setw(30) << ref_name;

    for (auto& def : data_output) {
        bool nan = false;
        // copy data to host if needed
        // and write it to the output file
        nan = check_array_for_nan(def.data, def.size, def.device_ptr, nan_check_d);

        oss << " " << def.short_name << ": " << nan;
        out &= !nan; // need to use binary neg, nan is true if data has NaNs, out should be true if
                     // data is OK
    }

#    ifndef BENCH_PRINT_DEBUG // if unset, print only failures
    if (!out)
#    endif // BENCH_PRINT_DEBUG
    {
        oss << "\n";

        log::printf(oss.str().c_str());


        for (auto& def : data_output) {
            if (!path_exists(output_crash))
                create_output_dir(output_crash);
            crash_report(def, output_crash, iteration);
        }
        printf("Crash report output in %s\n", output_crash.c_str());
    }

    return out;
}


void binary_test::output_reference(const string&             iteration,
                                   const string&             ref_name,
                                   const vector<output_def>& data_output) {
    // open file
    string output_name = output_dir + "/" + output_base_name + ref_name + "_" + iteration + ".h5";

    if (!path_exists(output_dir))
        create_output_dir(output_dir);


    storage s(output_name);

    for (auto& def : data_output) {
        // copy data to host if needed
        // and write it to the output file
        if (def.device_ptr) {
            getDeviceData(def.data, mem_buf.get(), def.size * sizeof(double));
            check_last_cuda_error(string("output ") + ref_name + string(" ") + def.name);


            s.append_table(mem_buf.get(), def.size, def.short_name, "-", def.short_name);
        }
        else
            s.append_table(def.data, def.size, def.short_name, "-", def.short_name);
    }
}

binary_test& binary_test::get_instance() {
    static binary_test bt(BENCHMARK_DUMP_REF_PATH, BENCHMARK_DUMP_BASENAME);


    return bt;
}


bool binary_test::compare_to_reference(const string&             iteration,
                                       const string&             ref_name,
                                       const vector<output_def>& data_output) {

    string output_name = output_dir + "/" + output_base_name + ref_name + "_" + iteration + ".h5";

    if (!path_exists(output_name)) {
        log::printf("No compare file: %s\n.", output_name.c_str());

        return true;
    }

    storage s(output_name, true);

    bool               out = true;
    std::ostringstream oss;
    oss << std::left << std::setw(50) << output_name;
    oss << std::setw(8) << iteration << " ref: " << std::setw(30) << ref_name;

    std::map<std::string, compare_statistics> stats_table;

    for (auto& def : data_output) {
        bool comp = false;

        // copy data to host if needed
        // and write it to the output file
        if (!s.has_table(def.short_name)) {

            // no table in input
            oss << " " << def.short_name << ": "
                << "N/A";

            continue;
        }

        compare_statistics stats;

        if (def.device_ptr) {
            getDeviceData(def.data, mem_buf.get(), def.size * sizeof(double));
            comp = compare_to_saved_data(s, def.short_name, mem_buf.get(), def.size, stats);
        }
        else {
            comp = compare_to_saved_data(s, def.short_name, def.data, def.size, stats);
        }
#    ifdef BENCH_COMPARE_PRINT_STATISTICS
        if (comp == false) {
            stats_table[def.short_name] = stats;
        }
#    endif // BENCH_COMPARE_PRINT_STATISTICS

        oss << " " << def.short_name << ": " << comp;
        out &= comp;
    }

#    ifndef BENCH_PRINT_DEBUG // if unset, print only failures
    if (!out)
#    endif
    {
        oss << "\n";
        log::printf(oss.str().c_str());
    }


#    ifdef BENCH_COMPARE_PRINT_STATISTICS
    for (auto const& v : stats_table) {
        auto const& key   = v.first;
        auto const& value = v.second;
        log::printf("  %5s - num (fail/tot/fst idx): %8d/%8d - %5d Δabs(mx:%11g,mn:%11g) - "
                    "Δrel(mx:%11g,mn:%11g) -ref(mx:%11g,mn:%11g) -val(mx:%11g,mn:%11g)\n",
                    key.c_str(),
                    value.num_failures,
                    value.num_values,
                    value.first_failure_idx,
                    value.max_abs_delta,
                    value.mean_abs_delta,
                    value.max_rel_delta,
                    value.mean_rel_delta,
                    value.max_ref,
                    value.mean_ref,
                    value.max_val,
                    value.mean_val);
    }

#    endif // BENCH_COMPARE_PRINT_STATISTICS

    return out;
}

void binary_test::append_definitions(const map<string, output_def>& defs) {
    output_definitions.insert(defs.begin(), defs.end());
    int memsize = 0;
    for (auto& d : output_definitions) {
        if (d.second.size > memsize)
            memsize = d.second.size;
    }

    mem_buf = std::unique_ptr<double[]>(new double[memsize], std::default_delete<double[]>());
}

void binary_test::register_phy_modules_variables(const std::map<string, output_def>& defs,
                                                 const vector<string>& phy_modules_data_in_,
                                                 const vector<string>& phy_modules_data_out_) {
    append_definitions(defs);
    phy_modules_data_in.insert(
        phy_modules_data_in.end(), phy_modules_data_in_.begin(), phy_modules_data_in_.end());
    phy_modules_data_out.insert(
        phy_modules_data_out.end(), phy_modules_data_out_.begin(), phy_modules_data_out_.end());
}


#endif // BENCHMARKING
