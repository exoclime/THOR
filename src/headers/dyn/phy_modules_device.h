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
// datatype to store values to update on device in THOR loop
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


#pragma once
#include <vector>

#define NUM_PHY_MODULES_DYN_CORE_ARRAYS 10

struct device_RK_array {
    double* array_d;
    double* arrayk_d;
    double* arrayi_d;
    int     dimensions;
    device_RK_array(){};

    device_RK_array(double* array_d_, double* arrayk_d_, double* arrayi_d_, int dimensions_) :
        array_d(array_d_),
        arrayk_d(arrayk_d_),
        arrayi_d(arrayi_d_),
        dimensions(dimensions_){};
};

class device_RK_array_manager
{
public:
    device_RK_array_manager();

    bool register_array(double* array_d, double* arrayk_d, double* arrayi, int dimensions);

    void allocate_device_array();

private:
    std::vector<device_RK_array> data;
};
