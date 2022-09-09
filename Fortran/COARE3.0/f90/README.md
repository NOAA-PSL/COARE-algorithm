Programs:
*	cor3_0af.F90: this version is set up to use the Fairall near-surface temperature sensor for Ts bulk.
*	cor3_0ah.F90: this version is setup to use the MSP 6m-depth temperature sensor for Ts bulk.

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
*	tst3_0af_out.txt and tst3_0ah_out.txt				Matlab output file from test data 


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

