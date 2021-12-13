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

#ifdef CGS_UNITS
// To check what units are used. Is pretty noisy
// #warning "Compiling with CGS units"
// physical constants
const double PI         = 3.141592653589793;
const double HCONST     = 6.62607004e-27;
const double CSPEED     = 29979245800.0;
const double KBOLTZMANN = 1.38064852e-16;
const double STEFANBOLTZMANN =
    5.6703669999999995e-5; // yes, it needs to have this exact value to be consistent with astropy

// constants in CGS
const double C        = 29979245800.0;          // speed of light in cm / s
const double K_B      = 1.380649e-16;           // Boltzmann constant in erg / K
const double H        = 6.62607015e-27;         // Planck constant in erg s
const double R_UNIV   = 83144626.1815324;       // universal gas constant in erg / mol / K
const double SIGMA_SB = 5.6703744191844314e-05; // Stefan-Boltzmann constant in erg / cm2 / K
const double AU       = 14959787070000.0;       // astronomical unit in cm
const double AMU      = 1.6605390666e-24;       // atomic mass unit in g

const double Mass_E = 9.1093837015e-28; // mass of electron in g

const double Q_E     = 4.803204712570263e-10;  // charge of electron in Franklin, cgs-unit of charge
const double R_SUN   = 69570000000.0;          // solar radius in cm
const double M_SUN   = 1.988409870698051e+33;  // solar mass
const double R_JUP   = 7149200000.0;           // radius Jupiter in cm, old value -- 6.9911e9
const double M_JUP   = 1.8981245973360504e+30; // mass Jupiter
const double R_EARTH = 637810000.0;            // radius Earth in cm, old value
const double M_EARTH = 5.972167867791379e+27;  // mass Earth
const double G       = 6.674299999999999e-08;  // gravitational constant (cgs)
const double GAMMA   = 0.5772156649;           // Euler-Mascheroni constant
#else
// To check what units are used. Is pretty noisy
// #warning "Compiling with SI units"
#    include "physics_constants_si.h"
#endif // CGS_constants
