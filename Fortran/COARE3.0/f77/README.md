Brief notes about the algorithm and the test data set appears at the head of the fortran file. Here we provide some comments on the structure of the code:

1. The input "read" statement is set up for the test data file test3_0.txt . This consists of four days of Moana Wave COARE data, 26-29 Nov 1992, prepared from Chris Fairall's hourly data file wavhr2_5.asc dated 31/10/96. A full description of the Moana Wave operations, instruments and data set is given at http://www.ncdc.noaa.gov/coare/catalog/data/air_sea_fluxes/moana_flux.html

2. A list of input variables is given, with units etc..  Only the first 11, the critical environmental and position variables, appear in the test data.  If time series of p and/or zi are available to the user, they may be added to the “700 read…” input string, and the default values disabled.  The remaining input parameters relate to the instrument set-up and are expected to remain fixed for the duration of the observations.  They are therefore set in the main program, as are the switches for cool skin, warm layer and wave state options.

3. Properly, wind speed should be relative to the water surface.  “u” should be the vector sum of measured wind speed and surface current if available, calculated by the user outside the program.

4. Two sea temperature measurements are given in the test data, one at only 0.05m depth, the other at 6m from Mike Gregg's Advanced Microstructure Profiler which operated from the Moana Wave during leg 1. The data was kindly provided in suitable form by Hemantha Wijesekera (COAS, OSU).  It allows a better demonstration of the calculation of skin temperature from the bulk via the warm layer option.  “ts_depth” should be set to correspond with whichever of “ts” or “hwt” is selected.

5. If skin temperature (sst) is measured directly with an infra-red radiometer, the warm layer and cool skin codes should be bypassed by setting jwarm and jcool to zero

6. jwave selects the method of calculating zo at the beginning of the main iteration loop in subroutine(ASL), i.e. Smith (1988) or one of the two wind/wave parameterizations.  The latter require values for the significant wave height (hwave) and dominant wave period (twave), which are calculated in the code from formulas given by Taylor and Yelland (2001) for fully developed seas.  If measurements of hwave and twave are available, they should be added to the input data string and the default values disabled.

7. Structure of the fortran code.
The main program (fluxes) opens the input and output files, reads the data and sets fixed values, defaults and options.  It adds a data line number (index) and calls subroutine bulk_flux, passing input data in COMMON.
Bulk_flux defines most physical constants and coefficients and, after determining the proper conditions, calculates and integrates the diurnal warming of the ocean surface, using fluxes and net longwave radiation from the previous time-step and the solar absorption profile.  The fraction of warming above the temperature sensor is added to the measurement, and subroutine ASL is called for the flux and boundary layer calculations.
ASL is a descendant of the original LKB code, but almost all operations and parameterizations are changed.  After a series of first guesses and operations to characterize the atmospheric surface layer within the framework of Monin-Obhukov similarity theory, the core of the subroutine is an iteration loop.  This iterates three times over the fluxes (in the form u*, t*, q*), the roughness parameters (zo, zot, zoq), the M-O stability parameter and profile phi functions, and also calculates gustiness and the cool skin within the loop.  Final values are returned to bulk_flux in COMMON.
Finally, bulk_flux calculates the surface fluxes (Wm-2), skin temperature (sst), heat and momentum fluxes due to rainfall, neutral transfer coefficients, values of state variables at standard height, etc., and saves the fluxes for the next timestep warm layer integrals.  Output files are written before returning to the main program for the next line of input data.

8. Outputs available from bulk_flux are listed at the head of the program.  The outputs in tst3_0bf.out are given below. To illustrate the warm layer and cool skin, we output the respective delta-temperatures and layer thicknesses. Note that dt_warm is the temperature across the entire warm layer.  Only if tk_pwp is less than the sensor depth (ts_depth = 6m in the test case) will tsw=ts +dt_warm. Otherwise, a linear profile is assumed, and the warming above the bulk sensor calculated.  The measurement of “ts” at 0.05m depth will generally include most of the warm layer but not the cool skin effect.

Programs:
*	cor3_0bf.for: this version is set up to use the Fairall near-surface temperature sensor for Ts bulk.
*	cor3_0bh.for: this version is setup to use the MSP 6m-depth temperature sensor for Ts bulk.
 
Input data:
*	test3_0.txt (Test data)

1	Date: YYYYMMDDHHmmss.ss, YYYY=year, MM=month, DD=day, HH=hour, mm=minute,ss.ss=sec
2	U:     true wind speed at 15-m height  m/s corrected for surface currents
3	Tsea:  sea surface temp (at about 0.05m depth)  deg.C
4	Tair:  Vaisala air temperature (about 15 m)  deg.C
5	qair:  Vaisala air specific humidity (about 15 m)  g/kg
6	Rs:    solar irradiance  W/m2
7	Rl:    downwelling longwave irradiance  W/m2
8	Rain:  precipitation mm/hr
9	Lat:   Latitude (N=+)
10	Lon:   Longitude (E=+)
11	MSP:   MSP temperature at 6m depth  deg.C

Output files:
*	tst3_0bf.out and tst3_0bh.out  	Fortran output files from test data

1	index:	data line number
2	xtime:	YYYYMMDDHHmmss, date and time as read in (without dec. sec.)
3	hf:	sensible heat flux   W/m2
4	ef:	latent heat flux    W/m2
5	sst:	sea skin temperature   deg.C
6	tau:	surface stress   N/m2
7	Wbar:	mean Webb vertical velocity m/s
8	rf:	sensible heat flux due to precipitation   W/m2
9	dter:	cool skin effect  deg.C
10	dt_wrm: warming across entire warm layer deg.C
11	tk_pwp:warm layer thickness  m
12	tkt*1000:tkt=cool skin thickness
13	Wg:	gustiness velocity m/s
