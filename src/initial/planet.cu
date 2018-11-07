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
// Defines Planet's properties
//
//
// Description: Planet parameters.
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

#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "../headers/planet.h"

XPlanet::XPlanet() {
    //
    //  Earth
    // ID
    sprintf(simulation_ID, "%s", "Earth");
    //////////////
    // BULK     //
    //////////////
    A      = 72427000.0; // Radius [m]
    Omega  = 9.09E-5;    // Rotation rate [s-1]
    Gravit = 47.0;       // Gravitational acceleration [m/s^2]
    ////////////////
    // ATMOSPHERE //
    ////////////////
    Rd           = 3714;       // Gas constant [J/(Kg K)]
    Cp           = 13000;      // Specific heat capacities [J/(Kg K)]
    Tmean        = 1400;       // Mean atmospheric temperature [K]
    P_Ref        = 10000000.0; // Reference surface pressure [Pa]
    Top_altitude = 1235376.0;  // Altitude of the top of the model domain [m]
    Diffc        = 0.009973;   // Strength of diffusion
}
