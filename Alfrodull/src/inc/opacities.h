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

#pragma once

#include <memory>
#include <string>


#include "cuda_device_memory.h"
#include "storage.h"

#include "opacities_helpers.h"

// class to handle opacities and meanmolmass from k_table
// loads data table
// initialise variables for ref table management
// prepares GPU data
// loads data to GPU
// interpolates data
class opacity_table
{
public:
    opacity_table();

    bool load_opacity_table(const string& filename, bool is_CGS);

    std::unique_ptr<double[]> data_opac_wave = nullptr;
    //private:
    string opacity_filename;

    double experimental_opacities_offset = 0.0;
    // device tables
    cuda_device_memory<double> dev_kpoints;

    cuda_device_memory<double> dev_temperatures;
    int                        n_temperatures = 0;
    cuda_device_memory<double> dev_pressures;
    int                        n_pressures = 0;

    // wieghted Rayeigh c
    cuda_device_memory<double> dev_scat_cross_sections;

    // Mean molecular mass
    cuda_device_memory<double> dev_meanmolmass;

    cuda_device_memory<double> dev_opac_wave;
    int                        nbin = 0;

    cuda_device_memory<double> dev_opac_y;
    int                        ny = 0;

    cuda_device_memory<double> dev_opac_interwave;

    cuda_device_memory<double> dev_opac_deltawave;
};
