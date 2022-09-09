"""
Functions for COARE model version 3.6 bulk flux calculations.

Translated from MATLAB scripts written by Jim Edson and Chris Fairall, 
further edited by Elizabeth Thompson. Main Matalb script source is coare36vn_zrf_et.m

Execute '%run coare36vn_zrf_et' from the iPython command line for test run with
'test_36_data.txt' input data file. 
Code can also be imported as a module so that subfunctions can be used independently if desired.
List of functions in this code are:
    ['RHcalc',
    'albedo_vector',
    'bucksat',
    'coare36vn_zrf_et',
    'grv',
    'psit_26',
    'psiu_26',
    'psiu_40',
    'qsat26air',
    'qsat26sea']

ludovic Bariteau, CU/CIRES, NOAA/ESRL/PSL
v1: August 2022
"""
import numpy as np
import os
    
def coare36vn_zrf_et(u, zu , t, zt, rh, zq, P, ts, sw_dn, lw_dn, lat, lon,jd, zi,rain, Ss, cp=None, sigH=None, zrf_u=10.0, zrf_t=10.0, zrf_q=10.0):   
#**************************************************************************
# VERSION INFO:
    
# Vectorized version of COARE 3.6 code (Fairall et al, 2003) with
# modification based on the CLIMODE, MBL and CBLAST experiments
# (Edson et al., 2012). The cool skin and surface wave options are included.
# A separate warm layer function can be used to call this function.
    
# This 3.6 version include parameterizations using wave height and wave
# slope using cp and sigH.  If these are set to NaN, then the wind
# speed dependent formulation is used.  The parameterizations are based
# on fits to the Banner-Norison wave model and the Fairall-Edson flux
# database.  This version also allows salinity as a input.
# Open ocean example Ss=35; Great Lakes Ss=0;
    
#**************************************************************************
# COOL SKIN:
    
# An important component of this code is whether the inputed ts
# represents the ocean skin temperature or a subsurface temperature.
# How this variable is treated is determined by the jcool parameter:
#   set jcool=1 if ts is subsurface or bulk ocean temperature (default);
#   set jcool=0 if ts is skin temperature or to not run cool skin model.
# The code updates the cool-skin temperature depression dT_skin and
# thickness dz_skin during iteration loop for consistency. The number of
# iterations set to nits = 6.
    
    jcoolx = 1
#**************************************************************************
### INPUTS:
    
# Notes on input default values, missing values, vectors vs. single values:
#   - the code assumes u,t,rh,ts,P,sw_dn,lw_dn,rain,Ss,cp,sigH are vectors;
#   - sensor heights (zu,zt,zl) latitude lat, longitude lon, julian date jd,
#       and PBL height zi may be constants;
#   - air pressure P and radiation sw_dn, lw_dn may be vectors or constants.
#   - input NaNs as vectors or single values to indicate no data.
#   - assign a default value to P, lw_dn, sw_dn, lat, zi if unknown, single
#       values of these inputs are okay.
# Notes about signs and units:
#   - radiation signs: positive warms the ocean
#   - signs and units change throughout the program for ease of calculations.
#   - the signs and units noted here are for the inputs.
    
    
#  u = water-relative wind speed magnitude (m/s) at height zu (m)
#             i.e. mean wind speed accounting for the ocean current vector.
#             i.e. the magnitude of the difference between the wind vector
#             (at height zu) and ocean surface current vector.
#             If not available, use true wind speed to compute fluxes in
#             earth-coordinates only which will be ignoring the stress
#             contribution from the ocean current to all fluxes
#  t = air temperature (degC) at height zt (m)
#  rh = relative humidity (#) at height zq (m)
#  P = sea level air pressure (mb)
#  ts = seawater temperature (degC), see jcool below for cool skin
#             calculation, and separate warm layer code for specifying
#             sensor depth and whether warm layer is computed
#  sw_dn = downward (positive) shortwave radiation (W/m^2)
#  lw_dn = downward (positive) longwave radiation (W/m^2)
#  lat = latitude defined positive to north
#  lon = longitude defined positive to east, if using other version,
#             adjust the eorw string input to albedo_vector function
#  jd = year day or julian day, where day Jan 1 00:00 UTC = 0
#  zi = PBL height (m) (default or typical value = 600m)
#  rain = rain rate (mm/hr)
#  Ss = sea surface salinity (PSU)
#  cp = phase speed of dominant waves (m/s) computed from peak period
#  sigH = significant wave height (m)
#  zu, zt, zq heights of the observations (m)
#  zrf_u, zrf_t, zrf_q  reference height for profile.  Use this to compare observations at different heights
    
#**************************************************************************
#### OUTPUTS: the user controls the output array A at the end of the code.
    
# Note about signs and units:
#   - radiation signs: positive warms the ocean
#   - sensible, rain, and latent flux signs: positive cools the ocean
#   - signs and units change throughout the program for ease of calculations.
#   - the signs and units noted here are for the final outputs.
    
#    usr = friction velocity that includes gustiness (m/s), u*
#    tau = wind stress that includes gustiness (N/m^2)
#    hsb = sensible heat flux (W/m^2) ... positive for Tair < Tskin
#    hlb = latent heat flux (W/m^2) ... positive for qair < qs
#    hbb = atmospheric buoyany flux (W/m^2)... positive when hlb and hsb heat the atmosphere
#    hsbb = atmospheric buoyancy flux from sonic ... as above, computed with sonic anemometer T
#    hlwebb = webb factor to be added to hl covariance and ID latent heat fluxes
#    tsr = temperature scaling parameter (K), t*
#    qsr = specific humidity scaling parameter (kg/kg), q*
#    zo = momentum roughness length (m)
#    zot = thermal roughness length (m)
#    zoq = moisture roughness length (m)
#    Cd = wind stress transfer (drag) coefficient at height zu (unitless)
#    Ch = sensible heat transfer coefficient (Stanton number) at height zu (unitless)
#    Ce = latent heat transfer coefficient (Dalton number) at height zu (unitless)
#    L = Monin-Obukhov length scale (m)
#    zeta = Monin-Obukhov stability parameter zu/L (dimensionless)
#    dT_skin = cool-skin temperature depression (degC), pos value means skin is cooler than subskin
#    dq_skin = cool-skin humidity depression (g/kg)
#    dz_skin = cool-skin thickness (m)
#    Urf = wind speed at reference height (user can select height at input)
#    Trf = air temperature at reference height
#    Qrf = air specific humidity at reference height
#    RHrf = air relative humidity at reference height
#    UrfN = neutral value of wind speed at reference height
#    TrfN = neutral value of air temp at reference height
#    qarfN = neutral value of air specific humidity at reference height
#    lw_net = Net IR radiation computed by COARE (W/m2)... positive heating ocean
#    sw_net = Net solar radiation computed by COARE (W/m2)... positive heating ocean
#    Le = latent heat of vaporization (J/K)
#    rhoa = density of air at input parameter height zt, typically same as zq (kg/m3)
#    UN = neutral value of wind speed at zu (m/s)
#    U10 = wind speed adjusted to 10 m (m/s)
#    UN10 = neutral value of wind speed at 10m (m/s)
#    Cdn_10 = neutral value of drag coefficient at 10m (unitless)
#    Chn_10 = neutral value of Stanton number at 10m (unitless)
#    Cen_10 = neutral value of Dalton number at 10m (unitless)
#    hrain = rain heat flux (W/m^2)... positive cooling ocean
#    Qs = sea surface specific humidity, i.e. assuming saturation (g/kg)
#    Evap = evaporation rate (mm/h)
#    T10 = air temperature at 10m (deg C)
#    Q10 = air specific humidity at 10m (g/kg)
#    RH10 = air relative humidity at 10m (#)
#    P10 = air pressure at 10m (mb)
#    rhoa10 = air density at 10m (kg/m3)
#    gust = gustiness velocity (m/s)
#    wc_frac = whitecap fraction (ratio)
#    Edis = energy dissipated by wave breaking (W/m^2)
    
#**************************************************************************
#### ADDITONAL CALCULATIONS:
    
#   using COARE output, one can easily calculate the following using the
#   sign conventions and names herein:
    
#     ### Skin sea surface temperature or interface temperature; neglect
#     ### dT_warm_to_skin if warm layer code is not used as the driver
#     #### program for this program
#           Tskin = ts + dT_warm_to_skin - dT_skin;
    
#     ### Upwelling radiative fluxes: positive heating the ocean
#           lw_up = lw_net - lw_dn;
#           sw_up = sw_net - sw_dn;
    
#     ### Net heat flux: positive heating ocean
#     ### note that hs, hl, hrain are defined when positive cooling
#     ### ocean by COARE, so their signs are flipped here:
#           hnet = sw_net + lw_net - hs - hl - hrain;
    
#**************************************************************************
#### REFERENCES:
    
    #  Fairall, C. W., E. F. Bradley, J. S. Godfrey, G. A. Wick, J. B. Edson,
#  and G. S. Young, 1996a: Cool-skin and warm-layer effects on sea surface
#  temperature. J. Geophys. Res., 101, 1295?1308.
    
    #  Fairall, C. W., E. F. Bradley, D. P. Rogers, J. B. Edson, and G. S. Young,
#  1996b: Bulk parameterization of air-sea fluxes for Tropical Ocean- Global
#  Atmosphere Coupled- Ocean Atmosphere Response Experiment. J. Geophys. Res.,
#  101, 3747?3764.
    
    #  Fairall, C. W., A. B. White, J. B. Edson, and J. E. Hare, 1997: Integrated
#  shipboard measurements of the marine boundary layer. Journal of Atmospheric
#  and Oceanic Technology, 14, 338?359
    
    #  Fairall, C.W., E.F. Bradley, J.E. Hare, A.A. Grachev, and J.B. Edson (2003),
#  Bulk parameterization of air sea fluxes: updates and verification for the
#  COARE algorithm, J. Climate, 16, 571-590.
    
    #  Edson, J.B., J. V. S. Raju, R.A. Weller, S. Bigorre, A. Plueddemann, C.W.
#  Fairall, S. Miller, L. Mahrt, Dean Vickers, and Hans Hersbach, 2013: On
#  the Exchange of momentum over the open ocean. J. Phys. Oceanogr., 43,
#  15891610. doi: http://dx.doi.org/10.1175/JPO-D-12-0173.1
    
#**************************************************************************
# CODE HISTORY:
    
# 1. 12/14/05 - created based on scalar version coare26sn.m with input
#    on vectorization from C. Moffat.
# 2. 12/21/05 - sign error in psiu_26 corrected, and code added to use variable
#    values from the first pass through the iteration loop for the stable case
#    with very thin M-O length relative to zu (zetau>50) (as is done in the
#    scalar coare26sn and COARE3 codes).
# 3. 7/26/11 - S = dT was corrected to read S = ut.
# 4. 7/28/11 - modification to roughness length parameterizations based
#    on the CLIMODE, MBL, Gasex and CBLAST experiments are incorporated
# 5. 9/20/2017 - New wave parameterization added based on fits to wave model
# 6. 9/2020 - tested and updated to give consistent readme info and units,
#    and so that no external functions are required. They are all included at
#    end of this program now. Changed names for a few things... including skin
#    dter -> dT_skin; dt -> dT; dqer -> dq_skin; tkt -> dz_skin
#    and others to avoid ambiguity:
#    Rnl -> lw_net; Rns -> sw_net; Rl -> lw_dn; Rs -> sw_dn;
#    SST -> Tskin; Also corrected heights at which q and P are
#    computed to be more accurate, changed units of qstar to kg/kg, removed
#    extra 1000 on neutral 10 m transfer coefficients;
# 7. 10/2021 - implemented zenith angle dependent sw_up and sw_net;
#    changed buoyancy flux calculation to follow Stull
#    textbook version of tv* and tv_sonic*; reformatted preamble of program for
#    consistent formatting and to reduce redundancy; resolved issues of
#    nomenclature around T adjusted to heights vs theta potential
#    temperature when computing dT for the purpose of sensible heat flux.
#-----------------------------------------------------------------------

#***********  prep input data *********************************************
    
### Make sure INPUTS are consistent in size. 
# Best to avoid NaNs as inputs as well. Will prevent weird results

    # be sure array inputs are ndarray floats for single value function
    # if inputs are already ndarray float this does nothing
    # otherwise copies are created in the local namespace
    # .flatten() return a 1D version in case single value input is already an array (array([[]]) vs array([]))
    if u.size ==1 and t.size ==1: 
        u = np.copy(np.asarray([u], dtype=float)).flatten()
        zu = np.copy(np.asarray([zu], dtype=float)).flatten()
        t = np.copy(np.asarray([t], dtype=float)).flatten()
        zt = np.copy(np.asarray([zt], dtype=float)).flatten()
        rh = np.copy(np.asarray([rh], dtype=float)).flatten()
        zq = np.copy(np.asarray([zq], dtype=float)).flatten()
        P = np.copy(np.asarray([P], dtype=float)).flatten()
        ts = np.copy(np.asarray([ts], dtype=float)).flatten()
        sw_dn = np.copy(np.asarray([sw_dn], dtype=float)).flatten()
        lw_dn = np.copy(np.asarray([lw_dn], dtype=float)).flatten()
        lat = np.copy(np.asarray([lat], dtype=float)).flatten()
        lon = np.copy(np.asarray([lon], dtype=float)).flatten()
        jd = np.copy(np.asarray([jd], dtype=float)).flatten()
        zi = np.copy(np.asarray([zi], dtype=float)).flatten()
        rain = np.copy(np.asarray([rain], dtype=float)).flatten()
        Ss = np.copy(np.asarray([Ss], dtype=float)).flatten()
        zrf_u = np.copy(np.asarray([zrf_u], dtype=float)).flatten()
        zrf_t = np.copy(np.asarray([zrf_t], dtype=float)).flatten()
        zrf_q = np.copy(np.asarray([zrf_q], dtype=float)).flatten()
    
    N = np.size(u)
    jcool = jcoolx * np.ones(N)
    
    if cp is not None and cp.size==1:
        cp = np.copy(np.asarray([cp], dtype=float)).flatten()
    elif cp is None:
        cp = np.nan * np.ones(N)

    if sigH is not None and sigH.size==1:
        sigH = np.copy(np.asarray([sigH], dtype=float)).flatten()
    elif sigH is None:
        sigH = np.nan * np.ones(N)
     
# Option to set local variables to default values if input is NaN... can do
# single value or fill each individual. Warning... this will fill arrays
# with the dummy values and produce results where no input data are valid
# ii=find(isnan(P)); P(ii)=1013;    # pressure
# ii=find(isnan(sw_dn)); sw_dn(ii)=200;   # incident shortwave radiation
# ii=find(isnan(lat)); lat(ii)=45;  # latitude
# ii=find(isnan(lw_dn)); lw_dn(ii)=400-1.6*abs(lat(ii)); # incident longwave radiation
# ii=find(isnan(zi)); zi(ii)=600;   # PBL height
# ii=find(isnan(Ss)); Ss(ii)=35;    # Salinity
    
# find missing input data
# iip = np.where(np.isnan(P))
# iirs = np.where(np.isnan(sw_dn))
# iilat = np.where(np.isnan(lat))
# iirl = np.where(np.isnan(lw_dn))
# iizi = np.where(np.isnan(zi))
# iiSs = np.where(np.isnan(Ss))
# Input variable u is assumed to be wind speed corrected for surface current
# (magnitude of difference between wind and surface current vectors). To
# follow orginal Fairall code, set surface current speed us=0. If us surface
# current data are available, construct u prior to using this code and
# input us = 0*u here;
    us = 0 * u
    # convert rh to specific humidity after accounting for salt effect on freezing
    # point of water
    Tf = - 0.0575 * Ss + 0.00171052 * Ss ** 1.5 - np.multiply(0.0002154996 * Ss,Ss)
    Qs = qsat26sea(ts,P,Ss,Tf) / 1000
    P_tq = P - (0.125 * zt)
    Q,Pv = qsat26air(t,P_tq,rh)
    
    # Assumes rh relative to ice T<0
    # Pv is the partial pressure due to wate vapor in mb
    Q = Q / 1000

    ice = np.zeros(N)
    iice = np.array(np.where(ts < Tf))
    ice[iice] = 1
    jcool[iice] = 0
    zos = 0.0005
    #***********  set constants ***********************************************
    zref = 10
    Beta = 1.2
    von = 0.4
    fdg = 1.0
    T2K = 273.16
    grav = grv(lat)
    #***********  air constants ***********************************************
    Rgas = 287.1
    Le = (2.501 - 0.00237 * ts) * 1000000.0
    cpa = 1004.67
    cpv = cpa * (1 + 0.84 * Q)
    rhoa = P_tq * 100.0 / (np.multiply(Rgas * (t + T2K),(1 + 0.61 * Q)))
    # Pv is the partial pressure due to wate vapor in mb
    rhodry = (P_tq - Pv) * 100.0 / (Rgas * (t + T2K))
    visa = 1.326e-05 * (1 + np.multiply(0.006542,t) + 8.301e-06 * t ** 2 - 4.84e-09 * t ** 3)
    lapse = grav / cpa
    
    #***********  cool skin constants  ***************************************
    ### includes salinity dependent thermal expansion coeff for water
    tsw = ts
    ii = np.array(np.where(ts < Tf))
    if np.size(ii) != 0:
        tsw[ii] = Tf[ii]
    Al35 = 2.1e-05 * (tsw + 3.2) ** 0.79
    # Al0 = (2.2 * real((tsw - 1) ** 0.82) - 5) * 1e-05
    Al0_i=(tsw - 1) ** 0.82
    Al0 = (2.2 * Al0_i.real - 5) * 1e-05
    Al = Al0 + np.multiply((Al35 - Al0),Ss) / 35
    ###################
    bets = 0.00075
    be = bets * Ss
    ####  see "Computing the seater expansion coefficients directly from the
    ####  1980 equation of state".  J. Lillibridge, J.Atmos.Oceanic.Tech, 1980.
    cpw = 4000
    rhow = 1022
    visw = 1e-06
    tcw = 0.6
    bigc = 16 * grav * cpw * (rhow * visw) ** 3.0 / (tcw ** 2 * rhoa ** 2)
    wetc = np.multiply(0.622 * Le,Qs) / (Rgas * (ts + T2K) ** 2)
    #***********  net solar and IR radiation fluxes ***************************
    ### net solar flux, aka sw, aka shortwave
    
    # *** for time-varying, i.e. zenith angle varying albedo using Payne 1972:
    # insert 'E' for input to albedo function if longitude is defined positive
    # to E (normal), in this case lon sign will be flipped for the calculation.
    # Otherwise specify 'W' and the sign will not be changed in the function.
    # Check: albedo should usually peak at sunrise not at sunset, though it may
    # vary based on sw_dn.
    alb,T_sw,solarmax_sw,psi_sw = albedo_vector(sw_dn,jd,lon,lat,eorw='E')
    sw_net = np.multiply((1 - alb),sw_dn)
    
    # *** for constant albedo:
    # sw_net = 0.945.*sw_dn; # constant albedo correction, positive heating ocean
    
    ### net longwave aka IR aka infrared
    # initial value here is positive for cooling ocean in the calculations
    # below. However it is returned at end of program as -lw_net in final output so
    # that it is positive heating ocean like the other input radiation values.
    lw_net = 0.97 * (5.67e-08 * (ts - 0.3 * jcool + T2K) ** 4 - lw_dn)
    #***********  begin bulk loop *********************************************
    
    #***********  first guess *************************************************
    
    # wind speed minus current speed
    du = u - us
    # air sea temperature difference for the purpose of sensible heat flux
    dT = ts - t - np.multiply(lapse,zt)
    # air-sea T diff must account for lapse rate between surface and instrument height
    # t is air temperature in C, ts is surface water temperature in C. dT is
    # an approximation that is equivalent to  dtheta where theta is the
    # potential temperature, and the pressure at sea level and instrument level
    # are used. They are equivalent (max difference = 0.0022 K). This way
    # elimniates the need to involve the pressures at different heights.
    # Using or assuming dry adiabatic lapse rate between the two heights
    # doesn't matter because if real pressures are used the result is the
    # unchanged. The dT need not include conversion to K either. Here's an example:
    # grav = grv(lat);
    # lapse=grav/cpa;
    # P_at_tq_height=(psealevel - (0.125*zt)); # P at tq measurement height (mb)
    # note psealevel is adjusted using same expression from pa height
    # Ta is originally in C and C2K = 273.15 to convert from C to K
    # theta = (b10.Ta+C2K).*(1000./P_tq).^(Rgas/cpa);
    # TadjK = (b10.Ta+C2K) + lapse*zt;
    # Tadj = b10.Ta + lapse*zt;
    # theta_sfc = (b10.Tskin+C2K).*(1000./b10.psealevel).^(Rgas/cpa);
    # TadjK_sfc = b10.Tskin+C2K;
    # Tadj_sfc = b10.Tskin;
        
    ### the adj versions are only 0.0022 K smaller than theta versions)
    # dtheta = theta_sfc - theta;
    # dTadjK = TadjK_sfc - TadjK;
    # dTadj = Tadj_sfc - Tadj; # so dT = Tskin - (Ta + lapse*zt) = Tskin - Ta - lapse*zt
        
    # put things into different units and expressions for more calculations,
    # including first guesses that get redone later
    dq = Qs - Q
    ta = t + T2K
    tv = np.multiply(ta,(1 + 0.61 * Q))
    gust = 0.5
    dT_skin = 0.3
    ut = np.sqrt(du ** 2 + gust ** 2)
    u10 = np.multiply(ut,np.log(10 / 0.0001)) / np.log(zu / 0.0001)
    usr = 0.035 * u10
    zo10 = 0.011 * usr ** 2.0 / grav + 0.11 * visa / usr
    Cd10 = (von / np.log(10.0 / zo10)) ** 2
    Ch10 = 0.00115
    Ct10 = Ch10 / np.sqrt(Cd10)
    zot10 = 10.0 / np.exp(von / Ct10)
    Cd = (von / np.log(zu / zo10)) ** 2
    Ct = von / np.log(zt / zot10)
    CC = von * Ct / Cd
    Ribcu = - zu / zi / 0.004 / Beta ** 3
    Ribu = np.multiply(np.multiply(- grav,zu) / ta,((dT - np.multiply(dT_skin,jcool)) + np.multiply(0.61 * ta,dq))) / ut ** 2
    zetau = np.multiply(np.multiply(CC,Ribu),(1 + 27 / 9 * Ribu / CC))
    k50 = np.array(np.where(zetau > 50))
    
    k = np.array(np.where(Ribu < 0))
    if np.size(Ribcu) == 1:
        zetau[k] = np.multiply(CC[k],Ribu[k]) / (1 + Ribu[k] / Ribcu)
        del k
    else:
        zetau[k] = np.multiply(CC[k],Ribu[k]) / (1 + Ribu[k] / Ribcu[k])
        del k
    
    L10 = zu / zetau
    gf = ut / du
    usr = np.multiply(ut,von) / (np.log(zu / zo10) - psiu_40(zu / L10))
    tsr = np.multiply(- (dT - np.multiply(dT_skin,jcool)),von) * fdg / (np.log(zt / zot10) - psit_26(zt / L10))
    qsr = - (dq - np.multiply(np.multiply(wetc,dT_skin),jcool)) * von * fdg / (np.log(zq / zot10) - psit_26(zq / L10))
    dz_skin = 0.001 * np.ones(N)
    #**********************************************************
    #  The following gives the new formulation for the
    #  Charnock variable
    #**********************************************************
    #############   COARE 3.5 wind speed dependent charnock
    charnC = 0.011 * np.ones(N)
    umax = 19
    a1 = 0.0017
    a2 = - 0.005
    # charnC = a1 * u10 + a2
    charnC=np.copy(np.asarray(a1 * u10 + a2, dtype=float))
    k = np.array(np.where(u10 > umax))
    if k.size!=0:
        charnC[k] = a1 * umax + a2
    #########   if wave age is given but not wave height, use parameterized
    #########   wave height based on wind speed
    hsig = np.multiply((0.02 * (cp / u10) ** 1.1 - 0.0025),u10 ** 2)
    hsig = np.maximum(hsig,0.25)
    ii = np.array(np.where(np.logical_and(np.logical_not(np.isnan(cp)) ,np.isnan(sigH))))
    if ii.size!=0:
        sigH[ii] = hsig[ii]
    Ad = 0.2
    #Ad=0.73./sqrt(u10);
    Bd = 2.2
    zoS = np.multiply(np.multiply(sigH,Ad),(usr / cp) ** Bd)
    charnS = np.multiply(zoS,grav) / usr / usr
    nits = 10
    
    charn = np.copy(charnC)  # creates a deep copy of charnC - if shallow copy (= only) charnC may change too below!
    ii = np.array(np.where(np.logical_not(np.isnan(cp))))
    charn[ii] = charnS[ii]
    #**************  bulk loop ************************************************
    
    for i in np.arange(1,nits+1).reshape(-1):
        zeta = np.multiply(np.multiply(np.multiply(von,grav),zu) / ta,(tsr + np.multiply(0.61 * ta,qsr))) / (usr ** 2)
        L = zu / zeta
        zo = np.multiply(charn,usr ** 2.0) / grav + 0.11 * visa / usr
        zo[iice] = zos
        rr = np.multiply(zo,usr) / visa
        rt = np.zeros(u.size)
        rq = np.zeros(u.size)
        # This thermal roughness length Stanton number is close to COARE 3.0 value
        zoq = np.minimum(0.00016,5.8e-05 / rr ** 0.72)
        # Andreas 1987 for snow/ice
        ik = np.array(np.where(rr[iice] <= 0.135))
        rt[iice[ik]] = rr[iice[ik]] * np.exp(1.25)
        rq[iice[ik]] = rr[iice[ik]] * np.exp(1.61)
        ik = np.array(np.where(rr[iice] > np.logical_and(0.135,rr[iice]) <= 2.5))
        rt[iice[ik]] = np.multiply(rr[iice[ik]],np.exp(0.149 - 0.55 * np.log(rr[iice[ik]])))
        rq[iice[ik]] = np.multiply(rr[iice[ik]],np.exp(0.351 - 0.628 * np.log(rr[iice[ik]])))
        ik = np.array(np.where(rr[iice] > np.logical_and(2.5,rr[iice]) <= 1000))
        rt[iice[ik]] = np.multiply(rr[iice[ik]],np.exp(0.317 - 0.565 * np.log(rr[iice[ik]]) - np.multiply(0.183 * np.log(rr[iice[ik]]),np.log(rr[iice[ik]]))))
        rq[iice[ik]] = np.multiply(rr[iice[ik]],np.exp(0.396 - 0.512 * np.log(rr[iice[ik]]) - np.multiply(0.18 * np.log(rr[iice[ik]]),np.log(rr[iice[ik]]))))
        # Dalton number is close to COARE 3.0 value
        zot = zoq
        cdhf = von / (np.log(zu / zo) - psiu_26(zu / L))
        cqhf = np.multiply(von,fdg) / (np.log(zq / zoq) - psit_26(zq / L))
        cthf = np.multiply(von,fdg) / (np.log(zt / zot) - psit_26(zt / L))
        usr = np.multiply(ut,cdhf)
        qsr = np.multiply(- (dq - np.multiply(np.multiply(wetc,dT_skin),jcool)),cqhf)
        tsr = np.multiply(- (dT - np.multiply(dT_skin,jcool)),cthf)
        # original COARE version buoyancy flux
        tvsr1 = tsr + np.multiply(0.61 * ta,qsr)
        tssr1 = tsr + np.multiply(0.51 * ta,qsr)
        # new COARE version buoyancy flux from Stull (1988) page 146
        # tsr here uses dT with the lapse rate adjustment (see code above). The
        # Q and ta values should be at measurement height, not adjusted heights
        tvsr = np.multiply(tsr,(1 + np.multiply(0.61,Q))) + np.multiply(0.61 * ta,qsr)
        tssr = np.multiply(tsr,(1 + np.multiply(0.51,Q))) + np.multiply(0.51 * ta,qsr)
        Bf = np.multiply(np.multiply(- grav / ta,usr),tvsr)
        gust = 0.2 * np.ones(N)
        k = np.array(np.where(Bf > 0))
        ### gustiness in this way is from the original code. Notes:
        # we measured the actual gustiness by measuring the variance of the
        # wind speed and empirically derived the the scaling. It's empirical
        # but it seems appropriate... the longer the time average then the larger
        # the gustiness factor should be, to account for the gustiness averaged
        # or smoothed out by the averaging. wstar is the convective velocity.
        # gustiness is beta times wstar. gustiness is different between mean of
        # velocity and square of the mean of the velocity vector components.
        # The actual wind (mean + fluctuations) is still the most relavent
        # for the flux. The models do u v w, and then compute vector avg to get
        # speed, so we've done the same thing. coare alg input is the magnitude
        # of the mean vector wind relative to water.
        if np.size(zi) == 1:
            gust[k] = Beta * (np.multiply(Bf[k],zi)) ** 0.333
            del k
        else:
            gust[k] = Beta * (np.multiply(Bf[k],zi[k])) ** 0.333
            del k
        ut = np.sqrt(du ** 2 + gust ** 2)
        gf = ut / du
        hsb = np.multiply(np.multiply(- rhoa * cpa,usr),tsr)
        hlb = np.multiply(np.multiply(np.multiply(- rhoa,Le),usr),qsr)
        qout = lw_net + hsb + hlb
        ### rain heat flux is not included in qout because we don't fully
        # understand the evolution or gradient of the cool skin layer in the
        # presence of rain, and the sea snake subsurface measurement input
        # value will capture some of the rain-cooled water already. TBD.
        ### solar absorption:
        # The absorption function below is from a Soloviev paper, appears as
        # Eq 17 Fairall et al. 1996 and updated/tested by Wick et al. 2005. The
        # coefficient was changed from 1.37 to 0.065 ~ about halved.
        # Most of the time this adjustment makes no difference. But then there
        # are times when the wind is weak, insolation is high, and it matters a
        # lot. Using the original 1.37 coefficient resulted in many unwarranted
        # warm-skins that didn't seem realistic. See Wick et al. 2005 for details.
        # That's the last time the cool-skin routine was updated. The
        # absorption is not from Paulson & Simpson because that was derived in a lab.
        # It absorbed too much and produced too many warm layers. It likely
        # approximated too much near-IR (longerwavelength solar) absorption
        # which probably doesn't make it to the ocean since it was probably absorbed
        # somewhere in the atmosphere first. The below expression could
        # likely use 2 exponentials if you had a shallow mixed layer...
        # but we find better results with 3 exponentials. That's the best so
        # far we've found that covers the possible depths.
        dels = np.multiply(sw_net,(0.065 + 11 * dz_skin - np.multiply(6.6e-05 / dz_skin,(1 - np.exp(- dz_skin / 0.0008)))))
        qcol = qout - dels
        # only needs stress, water temp, sum of sensible, latent, ir, solar,
        # and latent individually.
        alq = np.multiply(Al,qcol) + np.multiply(np.multiply(be,hlb),cpw) / Le
        xlamx = 6.0 * np.ones(N)
        #     the other is the salinity part caused by latent heat flux (evap) leaving behind salt.
        dz_skin = np.minimum(0.01,np.multiply(xlamx,visw) / (np.multiply(np.sqrt(rhoa / rhow),usr)))
        k = np.array(np.where(alq > 0))
        xlamx[k] = 6.0 / (1 + (np.multiply(bigc[k],alq[k]) / usr[k] ** 4) ** 0.75) ** 0.333
        dz_skin[k] = np.multiply(xlamx[k],visw) / (np.multiply(np.sqrt(rhoa[k] / rhow),usr[k]))
        del k
        dT_skin = np.multiply(qcol,dz_skin) / tcw
        dq_skin = np.multiply(wetc,dT_skin)
        lw_net = 0.97 * (5.67e-08 * (ts - np.multiply(dT_skin,jcool) + T2K) ** 4 - lw_dn)
        if i == 1:
            usr50 = usr[k50]
            tsr50 = tsr[k50]
            qsr50 = qsr[k50]
            L50 = L[k50]
            zeta50 = zeta[k50]
            dT_skin50 = dT_skin[k50]
            dq_skin50 = dq_skin[k50]
            tkt50 = dz_skin[k50]
        u10N = np.multiply(usr / von / gf,np.log(10.0 / zo))
        charnC = a1 * u10N + a2
        k = u10N > umax
        charnC[k] = a1 * umax + a2
        charn = charnC
        zoS = np.multiply(np.multiply(sigH,Ad),(usr / cp) ** Bd)
        charnS = np.multiply(zoS,grav) / usr / usr
        ii = np.array(np.where(np.logical_not(np.isnan(cp))))
        charn[ii] = charnS[ii]
    
    # end bulk loop
    
    # insert first iteration solution for case with zetau>50
    usr[k50] = usr50
    tsr[k50] = tsr50
    qsr[k50] = qsr50
    L[k50] = L50
    zeta[k50] = zeta50
    dT_skin[k50] = dT_skin50
    dq_skin[k50] = dq_skin50
    dz_skin[k50] = tkt50
    #****************  compute fluxes  ****************************************
    tau = np.multiply(np.multiply(rhoa,usr),usr) / gf
    
    hsb = np.multiply(np.multiply(np.multiply(- rhoa,cpa),usr),tsr)
    
    hlb = np.multiply(np.multiply(np.multiply(- rhoa,Le),usr),qsr)
    
    hbb = np.multiply(np.multiply(np.multiply(- rhoa,cpa),usr),tvsr)
    
    hbb1 = np.multiply(np.multiply(np.multiply(- rhoa,cpa),usr),tvsr1)
    
    hsbb = np.multiply(np.multiply(np.multiply(- rhoa,cpa),usr),tssr)
    
    hsbb1 = np.multiply(np.multiply(np.multiply(- rhoa,cpa),usr),tssr1)
    
    wbar = 1.61 * hlb / Le / (1 + 1.61 * Q) / rhoa + hsb / rhoa / cpa / ta
    
    hlwebb = np.multiply(np.multiply(np.multiply(rhoa,wbar),Q),Le)
    
    Evap = 1000 * hlb / Le / 1000 * 3600
    
    #*****  compute transfer coeffs relative to ut @ meas. ht  ****************
    Cd = tau / rhoa / ut / np.maximum(0.1,du)
    Ch = np.multiply(- usr,tsr) / ut / (dT - np.multiply(dT_skin,jcool))
    Ce = np.multiply(- usr,qsr) / (dq - np.multiply(dq_skin,jcool)) / ut
    #***##  compute 10-m neutral coeff relative to ut *************************
    Cdn_10 = von ** 2.0 / np.log(10.0 / zo) ** 2
    Chn_10 = von ** 2.0 * fdg / np.log(10.0 / zo) / np.log(10.0 / zot)
    Cen_10 = von ** 2.0 * fdg / np.log(10.0 / zo) / np.log(10.0 / zoq)
    #***##  compute 10-m neutral coeff relative to ut *************************
    
    # Find the stability functions for computing values at user defined
    # reference heights and 10 m
    psi = psiu_26(zu / L)
    psi10 = psiu_26(10.0 / L)
    psirf = psiu_26(zrf_u / L)
    psiT = psit_26(zt / L)
    psi10T = psit_26(10.0 / L)
    psirfT = psit_26(zrf_t / L)
    psirfQ = psit_26(zrf_q / L)
    gf = ut / du
    #*********************************************************
    #  Determine the wind speeds relative to ocean surface at different heights
    #  Note that usr is the friction velocity that includes
    #  gustiness usr = sqrt(Cd) S, which is equation (18) in
    #  Fairall et al. (1996)
    #*********************************************************
    S = ut
    U = du
    S10 = S + np.multiply(usr / von,(np.log(10.0 / zu) - psi10 + psi))
    U10 = S10 / gf
    # or U10 = U + usr./von./gf.*(log(10/zu)-psi10+psi);
    Urf = U + np.multiply(usr / von / gf,(np.log(zrf_u / zu) - psirf + psi))
    UN = U + np.multiply(psi,usr) / von / gf
    U10N = U10 + np.multiply(psi10,usr) / von / gf
    
    UrfN = Urf + np.multiply(psirf,usr) / von / gf
    UN2 = np.multiply(usr / von / gf,np.log(zu / zo))
    U10N2 = np.multiply(usr / von / gf,np.log(10.0 / zo))
    UrfN2 = np.multiply(usr / von / gf,np.log(zrf_u / zo))
    #******** rain heat flux *****************************
    dwat = 2.11e-05 * ((t + T2K) / T2K) ** 1.94
    dtmp = np.multiply((1.0 + 0.003309 * t - np.multiply(np.multiply(1.44e-06,t),t)),0.02411) / (np.multiply(rhoa,cpa))
    dqs_dt = np.multiply(Q,Le) / (np.multiply(Rgas,(t + T2K) ** 2))
    alfac = 1.0 / (1 + 0.622 * (np.multiply(np.multiply(dqs_dt,Le),dwat)) / (np.multiply(cpa,dtmp)))
    hrain = np.multiply(np.multiply(np.multiply(rain,alfac),cpw),((ts - t - np.multiply(dT_skin,jcool)) + np.multiply((Qs - Q - np.multiply(dq_skin,jcool)),Le) / cpa)) / 3600
    
    Tskin = ts - np.multiply(dT_skin,jcool)
    
    # P is sea level pressure, so use subtraction through hydrostatic equation
    # to get P10 and P at reference height
    P10 = P - (0.125 * 10)
    Prf = P - (0.125 * zref)
    T10 = t + np.multiply(tsr / von,(np.log(10.0 / zt) - psi10T + psiT)) + np.multiply(lapse,(zt - 10))
    Trf = t + np.multiply(tsr / von,(np.log(zrf_t / zt) - psirfT + psiT)) + np.multiply(lapse,(zt - zrf_t))
    TN = t + np.multiply(psiT,tsr) / von
    T10N = T10 + np.multiply(psi10T,tsr) / von
    TrfN = Trf + np.multiply(psirfT,tsr) / von
    # unused... these are here to make sure you gets the same answer whether
    # you used the thermal calculated roughness lengths or the values at the
    # measurement height. So at this point they are just illustrative and can
    # be removed or ignored if you want.
    TN2 = Tskin + np.multiply(tsr / von,np.log(zt / zot)) - np.multiply(lapse,zt)
    T10N2 = Tskin + np.multiply(tsr / von,np.log(10.0 / zot)) - np.multiply(lapse,10)
    TrfN2 = Tskin + np.multiply(tsr / von,np.log(zrf_t / zot)) - np.multiply(lapse,zrf_t)
    dq_skin = np.multiply(np.multiply(wetc,dT_skin),jcool)
    Qs = Qs - dq_skin
    dq_skin = dq_skin * 1000
    Qs = Qs * 1000
    Q = Q * 1000
    Q10 = Q + np.multiply(np.multiply(1000.0,qsr) / von,(np.log(10.0 / zq) - psi10T + psiT))
    Qrf = Q + np.multiply(np.multiply(1000.0,qsr) / von,(np.log(zrf_q / zq) - psirfQ + psiT))
    QN = Q + np.multiply(np.multiply(psiT,1000.0),qsr) / von / np.sqrt(gf)
    Q10N = Q10 + np.multiply(np.multiply(psi10T,1000.0),qsr) / von
    QrfN = Qrf + np.multiply(np.multiply(psirfQ,1000.0),qsr) / von
    # unused... these are here to make sure you gets the same answer whether
    # you used the thermal calculated roughness lengths or the values at the
    # measurement height. So at this point they are just illustrative and can
    # be removed or ignored if you want.
    QN2 = Qs + np.multiply(np.multiply(1000.0,qsr) / von,np.log(zq / zoq))
    Q10N2 = Qs + np.multiply(np.multiply(1000.0,qsr) / von,np.log(10.0 / zoq))
    QrfN2 = Qs + np.multiply(np.multiply(1000.0,qsr) / von,np.log(zrf_q / zoq))
    RHrf = RHcalc(Trf,Prf,Qrf / 1000,Tf)
    RH10 = RHcalc(T10,P10,Q10 / 1000,Tf)
    # recompute rhoa10 with 10-m values of everything else.
    rhoa10 = P10 * 100.0 / (np.multiply(Rgas * (T10 + T2K),(1 + 0.61 * (Q10 / 1000))))
    ############  Other wave breaking statistics from Banner-Morison wave model
    wc_frac = 0.00073 * (U10N - 2) ** 1.43
    wc_frac[U10 < 2.1] = 1e-05
    
    kk = np.array(np.where(np.isfinite(cp) == 1))
    wc_frac[kk] = 0.0016 * U10N[kk] ** 1.1 / np.sqrt(cp[kk] / U10N[kk])
    
    Edis = np.multiply(np.multiply(0.095 * rhoa,U10N),usr ** 2)
    wc_frac[iice] = 0
    Edis[iice] = 0
    #****************  output  ****************************************************
    # only return values if jcool = 1; if cool skin model was intended to be run
    dT_skinx = np.multiply(dT_skin,jcool)
    dq_skinx = np.multiply(dq_skin,jcool)
    # get rid of filled values where nans are present in input data
    bad_input = np.array(np.where(np.isnan(u) == 1))
    gust[bad_input] = np.nan
    dz_skin[bad_input] = np.nan
    zot[bad_input] = np.nan
    zoq[bad_input] = np.nan
    # flip lw_net sign for standard radiation sign convention: positive heating ocean
    lw_net = - lw_net
    # this sign flip means lw_net, net long wave flux, is equivalent to:
    # lw_net = 0.97*(lw_dn_best - 5.67e-8*(Tskin+C2K).^4);
    
    # adjust A output as desired:
    out = np.array([usr,tau,hsb,hlb,hbb,hsbb,hlwebb,tsr,qsr,zo,zot,zoq,Cd,Ch,Ce,L,zeta,dT_skinx,dq_skinx,dz_skin,Urf,Trf,Qrf,RHrf,UrfN,TrfN,QrfN,lw_net,sw_net,Le,rhoa,UN,U10,U10N,Cdn_10,Chn_10,Cen_10,hrain,Qs,Evap,T10,T10N,Q10,Q10N,RH10,P10,rhoa10,gust,wc_frac,Edis])
    #                1   2   3   4   5   6    7      8   9  10  11  12 13 14 15 16  17   18       19        20    21  22  23  24   25   26   27     28    29   30  31  32 33   34    35     36    37      38  39  40   41  42   43   44   45  46   47    48     49    50
    A = np.column_stack(out)
    return A
    
#------------------------------------------------------------------------------
    
def psit_26(zeta = None): 
    # computes temperature structure function
    dzeta = np.minimum(50,0.35 * zeta)
    psi = - ((1 + 0.6667 * zeta) ** 1.5 + np.multiply(0.6667 * (zeta - 14.28),np.exp(- dzeta)) + 8.525)
    k = np.array(np.where(zeta < 0))
    x = (1 - 15 * zeta[k]) ** 0.5
    psik = 2 * np.log((1 + x) / 2)
    x = (1 - 34.15 * zeta[k]) ** 0.3333
    psic = 1.5 * np.log((1 + x + x ** 2) / 3) - np.sqrt(3) * np.arctan((1 + 2 * x) / np.sqrt(3)) + 4 * np.arctan(1) / np.sqrt(3)
    f = zeta[k] ** 2.0 / (1 + zeta[k] ** 2)
    psi[k] = np.multiply((1 - f),psik) + np.multiply(f,psic)
    return psi
    
#------------------------------------------------------------------------------
    
def psiu_26(zeta = None): 
    # computes velocity structure function
    dzeta = np.minimum(50,0.35 * zeta)
    a = 0.7
    b = 3 / 4
    c = 5
    d = 0.35
    psi = - (a * zeta + np.multiply(b * (zeta - c / d),np.exp(- dzeta)) + b * c / d)
    k = np.array(np.where(zeta < 0))
    x = (1 - 15 * zeta[k]) ** 0.25
    psik = 2 * np.log((1 + x) / 2) + np.log((1 + np.multiply(x,x)) / 2) - 2 * np.arctan(x) + 2 * np.arctan(1)
    x = (1 - 10.15 * zeta[k]) ** 0.3333
    psic = 1.5 * np.log((1 + x + x ** 2) / 3) - np.sqrt(3) * np.arctan((1 + 2 * x) / np.sqrt(3)) + 4 * np.arctan(1) / np.sqrt(3)
    f = zeta[k] ** 2.0 / (1 + zeta[k] ** 2)
    psi[k] = np.multiply((1 - f),psik) + np.multiply(f,psic)
    return psi
    
#------------------------------------------------------------------------------
    
def psiu_40(zeta = None): 
    # computes velocity structure function
    dzeta = np.minimum(50,0.35 * zeta)
    a = 1
    b = 3 / 4
    c = 5
    d = 0.35
    psi = - (a * zeta + np.multiply(b * (zeta - c / d),np.exp(- dzeta)) + b * c / d)
    k = np.array(np.where(zeta < 0))
    x = (1 - 18 * zeta[k]) ** 0.25
    psik = 2 * np.log((1 + x) / 2) + np.log((1 + np.multiply(x,x)) / 2) - 2 * np.arctan(x) + 2 * np.arctan(1)
    x = (1 - 10 * zeta[k]) ** 0.3333
    psic = 1.5 * np.log((1 + x + x ** 2) / 3) - np.sqrt(3) * np.arctan((1 + 2 * x) / np.sqrt(3)) + 4 * np.arctan(1) / np.sqrt(3)
    f = zeta[k] ** 2.0 / (1 + zeta[k] ** 2)
    psi[k] = np.multiply((1 - f),psik) + np.multiply(f,psic)
    return psi
    
#------------------------------------------------------------------------------
    
def bucksat(T = None,P = None,Tf = None): 
    # computes saturation vapor pressure [mb]
    # given T [degC] and P [mb] Tf is freezing pt
    exx = np.multiply(np.multiply(6.1121,np.exp(np.multiply(17.502,T) / (T + 240.97))),(1.0007 + np.multiply(3.46e-06,P)))
    ii = np.array(np.where(T < Tf))
    if np.size(ii) != 0:
        exx[ii] = np.multiply(np.multiply((1.0003 + 4.18e-06 * P[ii]),6.1115),np.exp(np.multiply(22.452,T[ii]) / (T[ii] + 272.55)))
    
    return exx
    
#------------------------------------------------------------------------------
    
def qsat26sea(T = None,P = None,Ss = None,Tf = None): 
    # computes surface saturation specific humidity [g/kg]
    # given T [degC] and P [mb]
    ex = bucksat(T,P,Tf)
    fs = 1 - 0.02 * Ss / 35
    es = np.multiply(fs,ex)
    qs = 622 * es / (P - 0.378 * es)
    return qs
    
#------------------------------------------------------------------------------
    
def qsat26air(T = None,P = None,rh = None): 
    # computes saturation specific humidity [g/kg]
    # given T [degC] and P [mb]
    Tf = 0
    es = bucksat(T,P,Tf)
    em = np.multiply(0.01 * rh,es)
    q = 622 * em / (P - 0.378 * em)
    return q,em
    
#------------------------------------------------------------------------------
    
def grv(lat = None): 
    # computes g [m/sec^2] given lat in deg
    gamma = 9.7803267715
    c1 = 0.0052790414
    c2 = 2.32718e-05
    c3 = 1.262e-07
    c4 = 7e-10
    phi = lat * np.pi / 180
    x = np.sin(phi)
    g = gamma * (1 + c1 * x ** 2 + c2 * x ** 4 + c3 * x ** 6 + c4 * x ** 8)
    return g
    
#------------------------------------------------------------------------------
    
def RHcalc(T = None,P = None,Q = None,Tf = None): 
    # computes relative humidity given T,P, & Q
    es = np.multiply(np.multiply(6.1121,np.exp(np.multiply(17.502,T) / (T + 240.97))),(1.0007 + np.multiply(3.46e-06,P)))
    ii = np.array(np.where(T < Tf))
    if np.size(ii) != 0:
        es[ii] = np.multiply(np.multiply(6.1115,np.exp(np.multiply(22.452,T[ii]) / (T[ii] + 272.55))),(1.0003 + 4.18e-06 * P[ii]))
    em = np.multiply(Q,P) / (np.multiply(0.378,Q) + 0.622)
    RHrf = 100 * em / es
    return RHrf
    
#------------------------------------------------------------------------------
    
def albedo_vector(sw_dn = None,jd = None,lon = None,lat = None,eorw = None): 
    #  Computes transmission and albedo from downwelling sw_dn using
    #  lat   : latitude in degrees (positive to the north)
    #  lon   : longitude in degrees (positive to the west)
    #  jd    : yearday
    #  sw_dn : downwelling solar radiation measured at surface
    #  eorw  : 'E' if longitude is positive to the east, or 'W' if otherwise
        
    # updates:
    #   20-10-2021: ET vectorized function
    
    if eorw == 'E':
        #     disp('lon is positive to east so negate for albedo calculation');
        lon = - lon
    elif eorw == 'W':
        #     disp('lon is already positive to west so go ahead with albedo calculation');
        pass
    else:
        print('please provide sign information on whether lon is deg E or deg W')
        
    
    alb = np.full([np.size(sw_dn)],np.nan)
    lat = lat * np.pi / 180
    lon = lon * np.pi / 180
    SC = 1380
    utc = (jd - np.fix(jd)) * 24
    h = np.pi * utc / 12 - lon
    declination = 23.45 * np.cos(2 * np.pi * (jd - 173) / 365.25)
    solarzenithnoon = (lat * 180 / np.pi - declination)
    solaraltitudenoon = 90 - solarzenithnoon
    sd = declination * np.pi / 180
    gamma = 1
    gamma2 = gamma * gamma
    
    sinpsi = np.multiply(np.sin(lat),np.sin(sd)) - np.multiply(np.multiply(np.cos(lat),np.cos(sd)),np.cos(h))
    psi = np.multiply(np.arcsin(sinpsi),180) / np.pi
    solarmax = np.multiply(SC,sinpsi) / gamma2
    #solarmax=1380*sinpsi*(0.61+0.20*sinpsi);
    
    T = np.minimum(2,sw_dn / solarmax)
    
    Ts = np.arange(0,1+0.05,0.05)
    As = np.arange(0,90+2,2)
    
    #  Look up table from Payne (1972)  Only adjustment is to T=0.95 Alt=10 value
    #       0     2     4     6     8     10    12    14    16    18    20   22     24    26    28    30    32    34    36    38    40    42    44    46    48    50    52    54    56    58    60    62    64    66    68    70    72    74    76    78    80    82    84    86    88    90
    a = np.array([[0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061],[0.062,0.062,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061],[0.072,0.07,0.068,0.065,0.065,0.063,0.062,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.06,0.061,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06],[0.087,0.083,0.079,0.073,0.07,0.068,0.066,0.065,0.064,0.063,0.062,0.061,0.061,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06],[0.115,0.108,0.098,0.086,0.082,0.077,0.072,0.071,0.067,0.067,0.065,0.063,0.062,0.061,0.061,0.06,0.06,0.06,0.06,0.061,0.061,0.061,0.061,0.06,0.059,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.059,0.059,0.059],[0.163,0.145,0.13,0.11,0.101,0.092,0.084,0.079,0.072,0.072,0.068,0.067,0.064,0.063,0.062,0.061,0.061,0.061,0.06,0.06,0.06,0.06,0.06,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.058],[0.235,0.198,0.174,0.15,0.131,0.114,0.103,0.094,0.083,0.08,0.074,0.074,0.07,0.067,0.065,0.064,0.063,0.062,0.061,0.06,0.06,0.06,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.058,0.058,0.058],[0.318,0.263,0.228,0.192,0.168,0.143,0.127,0.113,0.099,0.092,0.084,0.082,0.076,0.072,0.07,0.067,0.065,0.064,0.062,0.062,0.06,0.06,0.06,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.058,0.058,0.058,0.058,0.058,0.058,0.058,0.058,0.057,0.058,0.058,0.058,0.058,0.057,0.057,0.057],[0.395,0.336,0.29,0.248,0.208,0.176,0.151,0.134,0.117,0.107,0.097,0.091,0.085,0.079,0.075,0.071,0.068,0.067,0.065,0.063,0.062,0.061,0.06,0.06,0.06,0.059,0.059,0.058,0.058,0.058,0.057,0.057,0.057,0.057,0.057,0.057,0.057,0.056,0.056,0.056,0.056,0.056,0.056,0.056,0.056,0.055],[0.472,0.415,0.357,0.306,0.252,0.21,0.176,0.154,0.135,0.125,0.111,0.102,0.094,0.086,0.081,0.076,0.072,0.071,0.068,0.066,0.065,0.063,0.062,0.061,0.06,0.059,0.058,0.057,0.057,0.057,0.056,0.055,0.055,0.055,0.055,0.055,0.055,0.054,0.053,0.054,0.053,0.053,0.054,0.054,0.053,0.053],[0.542,0.487,0.424,0.36,0.295,0.242,0.198,0.173,0.15,0.136,0.121,0.11,0.101,0.093,0.086,0.081,0.076,0.073,0.069,0.067,0.065,0.064,0.062,0.06,0.059,0.058,0.057,0.056,0.055,0.055,0.054,0.053,0.053,0.052,0.052,0.052,0.051,0.051,0.05,0.05,0.05,0.05,0.051,0.05,0.05,0.05],[0.604,0.547,0.498,0.407,0.331,0.272,0.219,0.185,0.16,0.141,0.127,0.116,0.105,0.097,0.089,0.083,0.077,0.074,0.069,0.066,0.063,0.061,0.059,0.057,0.056,0.055,0.054,0.053,0.053,0.052,0.051,0.05,0.05,0.049,0.049,0.049,0.048,0.047,0.047,0.047,0.046,0.046,0.047,0.047,0.046,0.046],[0.655,0.595,0.556,0.444,0.358,0.288,0.236,0.19,0.164,0.145,0.13,0.119,0.107,0.098,0.09,0.084,0.076,0.073,0.068,0.064,0.06,0.058,0.056,0.054,0.053,0.051,0.05,0.049,0.048,0.048,0.047,0.046,0.046,0.045,0.045,0.045,0.044,0.043,0.043,0.043,0.042,0.042,0.043,0.042,0.042,0.042],[0.693,0.631,0.588,0.469,0.375,0.296,0.245,0.193,0.165,0.145,0.131,0.118,0.106,0.097,0.088,0.081,0.074,0.069,0.065,0.061,0.057,0.055,0.052,0.05,0.049,0.047,0.046,0.046,0.044,0.044,0.043,0.042,0.042,0.041,0.041,0.04,0.04,0.039,0.039,0.039,0.038,0.038,0.038,0.038,0.038,0.038],[0.719,0.656,0.603,0.48,0.385,0.3,0.25,0.193,0.164,0.145,0.131,0.116,0.103,0.092,0.084,0.076,0.071,0.065,0.061,0.057,0.054,0.051,0.049,0.047,0.045,0.043,0.043,0.042,0.041,0.04,0.039,0.039,0.038,0.038,0.037,0.036,0.036,0.035,0.035,0.034,0.034,0.034,0.034,0.034,0.034,0.034],[0.732,0.67,0.592,0.474,0.377,0.291,0.246,0.19,0.162,0.144,0.13,0.114,0.1,0.088,0.08,0.072,0.067,0.062,0.058,0.054,0.05,0.047,0.045,0.043,0.041,0.039,0.039,0.038,0.037,0.036,0.036,0.035,0.035,0.034,0.033,0.032,0.032,0.032,0.031,0.031,0.031,0.03,0.03,0.03,0.03,0.03],[0.73,0.652,0.556,0.444,0.356,0.273,0.235,0.188,0.16,0.143,0.129,0.113,0.097,0.086,0.077,0.069,0.064,0.06,0.055,0.051,0.047,0.044,0.042,0.039,0.037,0.035,0.035,0.035,0.034,0.033,0.033,0.032,0.032,0.032,0.029,0.029,0.029,0.029,0.028,0.028,0.028,0.028,0.027,0.027,0.028,0.028],[0.681,0.602,0.488,0.386,0.32,0.252,0.222,0.185,0.159,0.142,0.127,0.111,0.096,0.084,0.075,0.067,0.062,0.058,0.054,0.05,0.046,0.042,0.04,0.036,0.035,0.033,0.032,0.032,0.031,0.03,0.03,0.03,0.03,0.029,0.027,0.027,0.027,0.027,0.026,0.026,0.026,0.026,0.026,0.026,0.026,0.026],[0.581,0.494,0.393,0.333,0.288,0.237,0.211,0.182,0.158,0.141,0.126,0.11,0.095,0.083,0.074,0.066,0.061,0.057,0.053,0.049,0.045,0.041,0.039,0.034,0.033,0.032,0.031,0.03,0.029,0.028,0.028,0.028,0.028,0.027,0.026,0.026,0.026,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025],[0.453,0.398,0.342,0.301,0.266,0.226,0.205,0.18,0.157,0.14,0.125,0.109,0.095,0.083,0.074,0.065,0.061,0.057,0.052,0.048,0.044,0.04,0.038,0.033,0.032,0.031,0.03,0.029,0.028,0.027,0.027,0.026,0.026,0.026,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025],[0.425,0.37,0.325,0.29,0.255,0.22,0.2,0.178,0.157,0.14,0.122,0.108,0.095,0.083,0.074,0.065,0.061,0.056,0.052,0.048,0.044,0.04,0.038,0.033,0.032,0.031,0.03,0.029,0.028,0.027,0.026,0.026,0.026,0.026,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025]])

    if T.size==1:   ### for single value function
        Tchk = np.abs(Ts - T)
        i = np.array(np.where(Tchk == Tchk.min()))
        if psi < 0:
            alb = np.array([0])
            solarmax = np.array([0])
            T = np.array([0])
            j = np.array([0])
            psi = np.array([0])
        else:
            Achk = np.abs(As - psi)
            j = np.array(np.where(Achk == Achk.min()))
            szj = j.shape
            if szj[0] > 0:
                alb = a[i,j].flatten()
            else:
                #       print('no j found, not assigning alb to anything');
                pass
    else:  ### for vectorized function
        for k in np.arange(0,np.size(sinpsi)).reshape(-1):
            Tchk = np.abs(Ts - T[k])
            i = np.array(np.where(Tchk == Tchk.min()))
            if psi[k] < 0:
                alb[k] = 0
                solarmax[k] = 0
                T[k] = 0
                j = 0
                psi[k] = 0
            else:
                Achk = np.abs(As - psi[k])
                j = np.array(np.where(Achk == Achk.min()))
                szj = j.shape
                if szj[0] > 0:
                    alb[k] = a[i,j]
                else:
                    #       disp('no j found, not assigning alb to anything');
                    pass
    
    #disp([num2str(jd) '  ' num2str(sw_dn) '  ' num2str(alb) '  ' num2str(T) '  ' num2str(i) '  ' num2str(j)])
    return alb,T,solarmax,psi

#------------------------------------------------------------------------------

# This code executes if 'run coare36vn_zrf_et.py' is executed from iPython cmd line
# Edit line 959 to indicate path to test data file
if __name__ == '__main__':
    # import numpy as np
    # import os
    # import util
    # import matplotlib.pyplot as plt
    
    path = '/Users/ludo/Documents/Work/COARE/conversion2python_tests/'
    fil = 'test_36_data.txt'   
    data = np.genfromtxt(path+fil, skip_header=1)
    u = data[:,1]
    t = data[:,3]
    rh = data[:,5]
    P = data[:,7]
    ts = data[:,8]
    sw_dn = data[:,9]
    lw_dn = data[:,10]
    lat = data[:,11]
    lon = data[:,12]
    zi = data[:,13]
    rain = data[:,14]
    zu= data[:,2]
    zt= data[:,4]
    zq= data[:,6]
    zrf_u=10.0;
    zrf_t=10.0;
    zrf_q=10.0;
    Ss = data[:,15]
    jd = data[:,0]
    cp = data[:,16]
    sigH = data[:,17]
    
    A=coare36vn_zrf_et(u, zu , t, zt, rh, zq, P, ts, sw_dn, lw_dn, lat, lon,jd, zi,rain, Ss, cp , sigH, zrf_u, zrf_t, zrf_q)
    fnameA = os.path.join(path,'test_36_output_py_082022_withwavesinput.txt')
    # A=coare36vn_zrf_et(u, zu , t, zt, rh, zq, P, ts, sw_dn, lw_dn, lat, lon,jd, zi,rain, Ss, None , None, zrf_u, zrf_t, zrf_q)
    # fnameA = os.path.join(path,'test_36_output_py_082022_withnowavesinput.txt')
    A_hdr = 'usr\ttau\thsb\thlb\thbb\thlwebb\ttsr\tqsr\tzo\tzot\tzoq\tCd\t'
    A_hdr += 'Ch\tCe\tL\tzeta\tdT_skinx\tdq_skinx\tdz_skin\tUrf\tTrf\tQrf\t'
    A_hdr += 'RHrf\tUrfN\tTrfN\tQrfN\tlw_net\tsw_net\tLe\trhoa\tUN\tU10\tU10N\t'
    A_hdr += 'Cdn_10\tChn_10\tCen_10\thrain\tQs\tEvap\tT10\tT10N\tQ10\tQ10N\tRH10\t'
    A_hdr += 'P10\trhoa10\tgust\twc_frac\tEdis'
    np.savetxt(fnameA,A,fmt='%.18e',delimiter='\t',header=A_hdr)
    
   # test on signle value
   # A=coare36vn_zrf_et(u[0], zu[0], t[0], zt[0], rh[0], zq[0], P[0], ts[0], sw_dn[0], lw_dn[0], lat[0], lon[0],jd[0], zi[0],rain[0], Ss[0], cp[0] , sigH[0], zrf_u, zrf_t, zrf_q)


    