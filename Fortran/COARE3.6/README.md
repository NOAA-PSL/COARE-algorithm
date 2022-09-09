This repository provides a fortran implementation of the COARE3.6 algorithm for the computation of surface fluxes of momentum, heat and moisture over oceans.  Within the file module_coare36_model.f90, the algorithm is implemented in two subroutines:
 - coare36flux_coolskin, which assumes that an ocean bulk temperature is input in ts.  This subroutine will compute a skin temperature that accounts for the heating/cooling of the surface by radiative, sensible and latent heat fluxes.
 - coare36flux, assumes that a skin temperature is input in ts and only computes the surface energy fluxes.

This repository includes several files that have been extracted from a preliminary implementation in NOAA's Unified Forecast System (UFS), and the two model_coare36*f90 files are intended for use in that model.

The current versions of these files have been verified against a matlab implementation of COARE3.6 in coare36vn_zrf.m.  Sample output from the matlab code is provided in test_36_output_matlab_072821.txt.  Do note that acceleration due to gravity was fixed at 9.81 m/s2 in that test.  The default algorithm allows this to change with latitude.

The fortran program test_coare_36.f90 is a driver for the test suite in the data file test_35_data.txt.  This program may be compiled and run (using the intel compiler) with the script run_test.bash.  This program will produce an output file test_36_output_f90coolskin.txt, which may be compared against test_36_output_matlab_072821.txt.  An additional matlab script PlotCOARE36Errors.m produces plots of absolute and relative error for many quantities output by the scheme and, optionally, writes them to files in the Figures/ directory.  The default version of the algorithm agrees to the matlab version within 1% for most variables.  Better agreement can be achieved by decreasing dT_skin_min in module_coare36_parameters.f90 to 1.e-8, which tightens the tolerance for stopping iterations within the (iterative) flux computation.

Peter Blossey, July 2021
