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
//
//
// Description: Defines planet parameters.
//
//
// Method: -
//
//
// Known limitations: None
//
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

// Physical Constants
#define kb_constant 1.38e-23  // Boltzmann constant [J/K]
#define mu_constant 1.660e-27 // Atomic mass unit   [kg]


class XPlanet
{

public:
    char simulation_ID[160];

    //////////////
    // BULK     //
    //////////////

    double A;
    double Omega;
    double Gravit;

    ////////////////
    // ATMOSPHERE //
    ////////////////

    double Rd;
    double Cp;
    double Tmean;
    double P_Ref;
    double Top_altitude;
    double Diffc;

    XPlanet();
};
