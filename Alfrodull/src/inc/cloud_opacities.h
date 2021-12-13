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

// class to handle cloud scattering cross scection, absorption cross section and asymmetry
// loads data table
// initialise variables for ref table management
// prepares GPU data
// loads data to GPU
// interpolates data
class cloud_opacity_table
{
public:
    cloud_opacity_table();

    bool load(const string& filename);

    // device tables
    // scattering cross section
    cuda_device_memory<double> dev_scat_cross_sections;

    // absorption cross section
    cuda_device_memory<double> dev_abs_cross_sections;

    // asymmetry
    cuda_device_memory<double> dev_asymmetry;

    // wavelength grid
    cuda_device_memory<double> dev_wavelength;
};
