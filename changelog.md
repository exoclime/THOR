# changes since v2.0.0

 * argument parser fails and complains when wrong number of argument for keys
 * added test script to run multiple sims
 * fix conservation output (was not copied back before output)
 * clean up RT module startup output
 * added new routines to mjolnir/hamarr (KE, SR)
 * fixed issues in regrid
 * changed many of the plot functions to use regrid data (some still need work)
 * added entropy to global & grid conservation outputs
 * fixed pot temp update error in Density_Pressure_Eqs_Poles *major*
 * started generalizing plot functions in hamarr to make adding new quantities easier
 * added acoustic wave experiment from Tomita & Satoh 2004
 * added gravity wave experiment from Tomita & Satoh 2004
 * added python scripts to allow editing of initial h5 files (preliminary)
 * removed TPprof option, which didn't work as desired anyway: users should use the python tools to edit initial h5 files if non-isothermal state is needed
