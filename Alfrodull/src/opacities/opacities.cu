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
// Opacity tables loading and interpolation. This loads opacity tables from HELIOS'
// mixted opacity file generator.
//
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

#include "opacities.h"

#include "physics_constants.h"
#include <cstdio>
#include <stdexcept>
#include <tuple>

#include "opacities_helpers.h"

opacity_table::opacity_table() {
}

bool opacity_table::load_opacity_table(const string& filename, bool is_CGS) {
#ifdef CGS_UNITS
#    warning "Compiling with CGS units"
    // For usage in Helios without unit conversion
    const double temperatures_unit_conv = 1.0;
    const double pressures_unit_conv    = 1.0;
    const double wavelength_unit_conv   = 1.0;
    const double opacity_unit_conv      = 1.0;
    const double scat_cross_unit_conv   = 1.0;
#else  // SI units
    printf("Loading opacity tables from %s as %s units\n", filename.c_str(), is_CGS ? "CGS" : "SI");
    double temperatures_unit_conv = 1.0;
    double pressures_unit_conv    = 1.0;
    double wavelength_unit_conv   = 1.0;
    double opacity_unit_conv      = 1.0;
    double scat_cross_unit_conv   = 1.0;

    if (is_CGS) {
        // for SI units, convert CGS to SI
        temperatures_unit_conv = 1.0;
        pressures_unit_conv    = 1.0e-1;
        wavelength_unit_conv   = 1.0e-2;
        opacity_unit_conv      = 1.0e-1;
        scat_cross_unit_conv   = 1.0e-4;
    }
    else {
        // Paththrough
        temperatures_unit_conv = 1.0;
        pressures_unit_conv    = 1.0;
        wavelength_unit_conv   = 1.0;
        opacity_unit_conv      = 1.0;
        scat_cross_unit_conv   = 1.0;
    }
#endif // CGS_UNIT

    storage s(filename, true);

    read_table_to_device<double>(
        s, "/kpoints", dev_kpoints, opacity_unit_conv, experimental_opacities_offset);
    n_temperatures =
        read_table_to_device<double>(s, "/temperatures", dev_temperatures, temperatures_unit_conv);
    n_pressures = read_table_to_device<double>(s, "/pressures", dev_pressures, pressures_unit_conv);
    read_table_to_device<double>(
        s, "/weighted Rayleigh cross-sections", dev_scat_cross_sections, scat_cross_unit_conv);
    {
        std::unique_ptr<double[]> data = nullptr;
        int                       size = 0;

        tie(data, size) = read_table_to_host<double>(s, "/meanmolmass");
        for (int i = 0; i < size; i++)
            data[i] *= AMU;
        push_table_to_device<double>(data, size, dev_meanmolmass);
    }

    if (s.has_table("/center wavelengths")) {
        tie(data_opac_wave, nbin) = read_table_to_host<double>(s, "/center wavelengths");
    }
    else {
        tie(data_opac_wave, nbin) = read_table_to_host<double>(s, "/wavelengths");
    }
    for (int i = 0; i < nbin; i++)
        data_opac_wave[i] *= wavelength_unit_conv;
    push_table_to_device<double>(data_opac_wave, nbin, dev_opac_wave);

    if (s.has_table("/ypoints")) {
        ny = read_table_to_device<double>(s, "/ypoints", dev_opac_y);
    }
    else {
        std::unique_ptr<double[]> data(new double[1]);
        data[0]  = 0;
        int size = 1;
        push_table_to_device<double>(data, size, dev_opac_y);
        ny = 1;
    }


    std::unique_ptr<double[]> data_opac_interwave(new double[nbin + 1]);
    if (s.has_table("/interface wavelengths")) {
        read_table_to_device<double>(
            s, "/interface wavelengths", dev_opac_interwave, wavelength_unit_conv);
    }
    else {
        // quick and dirty way to get the lamda interface values
        data_opac_interwave[0] = data_opac_wave[0] - (data_opac_wave[1] - data_opac_wave[0]) / 2.0;
        for (int i = 0; i < nbin - 1; i++)
            data_opac_interwave[i + 1] = (data_opac_wave[i + 1] + data_opac_wave[i]) / 2.0;
        data_opac_interwave[nbin] =
            data_opac_wave[nbin - 1] + (data_opac_wave[nbin - 1] - data_opac_wave[nbin - 2]) / 2.0;
        push_table_to_device<double>(data_opac_interwave, nbin + 1, dev_opac_interwave);
    }
    if (s.has_table("/wavelength width of bins")) {
        read_table_to_device<double>(
            s, "/wavelength width of bins", dev_opac_deltawave, wavelength_unit_conv);
    }
    else {
        if (nbin == 1) {
            std::unique_ptr<double[]> data_opac_deltawave(new double[1]);
            data_opac_deltawave[0] = 0.0;
            push_table_to_device<double>(data_opac_deltawave, 1, dev_opac_deltawave);
        }
        else {

            std::unique_ptr<double[]> data_opac_deltawave(new double[nbin]);
            // unit already converted with interwave
            for (int i = 0; i < nbin; i++)
                data_opac_deltawave[i] = data_opac_interwave[i + 1] - data_opac_interwave[i];
            push_table_to_device<double>(data_opac_deltawave, nbin, dev_opac_deltawave);
        }
    }

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
