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
 * removed TPprof option, which didn't work as desired anyway: users should use the python tools to edit initial h5 files if non-isothermal initial state is needed
 * removal of spurious vertical component of horizontal momentum diffusion has been moved
 to a separate function "Correct_Horizontal" to avoid potential issues related to the order
 in which threads are called
 * incorrect diffusion quantity (diffmh_d) was being passed to Momentum_Eq in the dynamical core, changed now to correct value (DivM_d)
 * added crash_report tool that dumps location of nans when the model crashes
 * further updates to python plotting tools to make them more flexible
 * added "custom_example.py" to mjolnir to demonstrate how to make multipanel plots with the python tools
 * grid output now contains differential operators grad, div, and vertical component of curl (the latter is not used by the model, but may be useful in post-processing)
 * added new features to double-grey RT: latitude variations, power law scaling of optical
 depth, and surface heating
 * added boundary layer module which calculates drag against the lower surface (rayleigh drag only for now)
 * changed name of 'diff_fac' input parameter to 'diff_ang' to be more consistent ('diff_ang' = 1/diffusivity factor)
 * fixed incompatibility with recent versions (> 1.10.1) of hdf5 libraries
 * updated behavior of NonHydro = false and DeepModel = false model options to be consistent with White+ 2005
