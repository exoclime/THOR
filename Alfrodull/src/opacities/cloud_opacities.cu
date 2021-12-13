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
// Cloud opacity loading functions. For version 1.0, this simply loads on value per frequency bin
// could be extended to load lookup tables and do interpolation.
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


#include "cloud_opacities.h"

#include "physics_constants.h"
#include <cstdio>
#include <stdexcept>
#include <tuple>

#include "opacities_helpers.h"

cloud_opacity_table::cloud_opacity_table() {
}

bool cloud_opacity_table::load(const string& filename) {

    log::printf("Loading cloud opacity tables from %s\n", filename.c_str());
    storage s(filename, true);

    int num_asym = read_table_to_device<double>(s, "/asymmetry", dev_asymmetry);

    int num_scat = read_table_to_device<double>(s, "/scattering", dev_scat_cross_sections);

    int num_abs = read_table_to_device<double>(s, "/absorption", dev_abs_cross_sections);

    int num_wave = read_table_to_device<double>(s, "/wavelength", dev_wavelength);

    log::printf("Loaded %d asym, %d scat and %d abs values.\n", num_asym, num_scat, num_abs);

    cudaError_t err = cudaGetLastError();

    // Check device query
    if (err != cudaSuccess) {
        log::printf("[%s:%d] CUDA error check reports error: %s\n",
                    __FILE__,
                    __LINE__,
                    cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    return true;
}
