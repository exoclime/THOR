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
// Description: test for grid setup
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
#include "debug.h"

#include "grid.h"

#include "log_writer.h"
#include <iostream>
#include <map>
#include <vector>

using namespace std;


int main() {
    bool output  = false;
    bool compare = true;

    for (int l = 4; l < 7; l++) {
        int max_count = 0;

        Icogrid grid(true, // Spring dynamics option
                     1.15, // Parameter beta for spring dynamics
                     l,    // Horizontal resolution level
                     32,   // Number of vertical layers
                     5,
                     6371000.0, // Planet radius
                     36000.0,   // Top of the model's domain
                     false,     // sponge layer
                     &max_count);

        std::map<string, output_def> output_definitions = {
            // {"map name, {variable pointer, table size, name, short name, on device}}
            // grid
            {"point_xyz", {grid.point_xyz, grid.point_num * 3, "xyz", "x", false}},

            {"point_xyzq", {grid.point_xyzq, grid.point_num * 3 * 6, "xyzq", "xq", false}},

            /*
              disable int tables for now
              {"pent_ind",      { grid.pent_ind, 12, "pent_ind", "pi", false}},

              {"point_local",   { grid.point_local, 6*grid.point_num, "point_local", "pl", false}},
              {"halo",          { grid.halo, grid.nh, "halo", "halo", false}},
              {"maps",          { grid.maps, (grid.nl_region+2)*(grid.nl_region+2)*grid.nr, "maps", "m", false}},
            */

            {"func_r", {grid.func_r, 3 * grid.point_num, "func_r", "f", false}},
            {"areas", {grid.areas, 6 * 3 * grid.point_num, "areas", "a", false}},
            {"areasTr", {grid.areasTr, 6 * grid.point_num, "areasTr", "aTr", false}},

            {"areasT", {grid.areasT, grid.point_num, "areasT", "aT", false}},
            {"nvec", {grid.nvec, 6 * 3 * grid.point_num, "nvec", "nc", false}},
            {"nvecoa", {grid.nvecoa, 6 * 3 * grid.point_num, "nvecoa", "na", false}},
            {"nvecti", {grid.nvecti, 6 * 3 * grid.point_num, "nvecti", "nti", false}},
            {"nvecte", {grid.nvecte, 6 * 3 * grid.point_num, "nvecte", "nte", false}},

            {"Altitude", {grid.Altitude, grid.nv, "Altitude", "Alt", false}},
            {"Altitudeh", {grid.Altitudeh, grid.nvi, "Altitudeh", "Alth", false}},
            {"lonlat", {grid.lonlat, 2 * grid.point_num, "lonlat", "ll", false}},

            {"div", {grid.div, 7 * 3 * grid.point_num, "div", "d", false}},
            {"grad", {grid.grad, 7 * 3 * grid.point_num, "grad", "g", false}},
        };

        std::vector<string> output_vars({"point_xyz",
                                         "point_xyzq",
                                         "func_r",
                                         "areasTr",
                                         "areasT",
                                         "nvec",
                                         "nvecoa",
                                         "nvecti",
                                         "nvecte",
                                         "Altitude",
                                         "Altitudeh",
                                         "lonlat",
                                         "div",
                                         "grad"});

        // load definitions for variables
        std::vector<output_def> data_output;
        for (auto& name : output_vars) {
            auto&& it = output_definitions.find(name);
            if (it != output_definitions.end()) {
                data_output.push_back(it->second);
            }
        }

        binary_test& btester = binary_test::get_instance();
        btester.set_output("grid_test_" + std::to_string(l), "./results/grid_test/");
        btester.append_definitions(output_definitions);

        if (output)
            btester.output_reference(std::to_string(l), "test", data_output);


        if (compare) {
            cout << "start compare " << l << endl;

            bool result = btester.compare_to_reference(std::to_string(l), "test", data_output);
            if (result)
                cout << "grid compare " << l << " SUCCESS" << endl;
            else
                cout << "grid compare " << l << " FAIL" << endl;
        }
    }


    exit(0);
}
