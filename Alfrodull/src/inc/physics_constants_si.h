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
///
// physical constants
const double PI              = 3.141592653589793;
const double HCONST          = 6.62607015e-34;         // J s
const double CSPEED          = 299792458.0;            // m / s
const double KBOLTZMANN      = 1.38064852e-23;         // J / K
const double STEFANBOLTZMANN = 5.6703744191844314e-08; // W / (K^4 m^2)


// constants in CGS
const double C        = 299792458.0;            // speed of light in m / s
const double K_B      = 1.38064852e-23;         // Boltzmann constant in J / K
const double H        = 6.62607015e-34;         // Planck constant in J s
const double R_UNIV   = 8.31446261815324;       // universal gas constant in J / ( K mol )
const double SIGMA_SB = 5.6703744191844314e-08; // Stefan-Boltzmann constant in W / (K^4 m^2)
const double AU       = 149597870700.0;         // astronomical unit in m
const double AMU      = 1.6605390666e-27;       // atomic mass unit in kg

const double Mass_E = 9.1093837015e-31; // mass of electron in kg

const double Q_E     = 1.602176634e-19;        // charge of electron in Coulomb
const double R_SUN   = 695700000.0;            // solar radius in m
const double M_SUN   = 1.988409870698051e+30;  // solar mass in kg
const double R_JUP   = 71492000.0;             // radius Jupiter in m
const double M_JUP   = 1.8981245973360505e+27; // mass Jupiter in kg
const double R_EARTH = 6378100.0;              // radius Earth in m
const double M_EARTH = 5.972167867791379e+24;  // mass Earth in kg
const double G       = 6.6743e-11;             // gravitational constant m^3 / ( kg s^2)
const double GAMMA   = 0.5772156649;           // Euler-Mascheroni constant
