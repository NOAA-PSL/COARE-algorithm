"""
Functions for COARE model bulk flux calculations.

Translated and vectorized from J Edson/ C Fairall MATLAB scripts.

Execute '%run bulk.py' from the iPython command line for test run with
'test_35_data.txt' input data file.

Byron Blomquist, CU/CIRES, NOAA/ESRL/PSD3
v1: May 2015
"""

def coare35vn(u, t, rh, ts, P=1015, Rs=150, Rl=370, zu=18, zt=18, zq=18, lat=45,
             zi=600, rain=None, cp=None, sigH=None, jcool=1):
    """
    usage: A = coare35vn(u, t, rh, ts)  -  include other kwargs as desired

    Vectorized version of COARE 3 code (Fairall et al, 2003) with modification
    based on the CLIMODE, MBL and CBLAST experiments (Edson et al., 2013).
    The cool skin option is retained but warm layer and surface wave options
    have been removed.

    This version includes parameterizations of wave height and wave slope using
    cp and sigH.  Unless these are provided the wind speed dependent formulation
    is used.

    AN IMPORTANT COMPONENT OF THIS CODE IS WHETHER INPUT 'ts' REPRESENTS
    THE SKIN TEMPERATURE OR A NEAR SURFACE TEMPERATURE.  How this variable is
    treated is determined by the jcool parameter:  jcool=1 if Ts is bulk
    ocean temperature (default), jcool=0 if Ts is ocean skin temperature.

    The code assumes u, t, rh, and ts are vectors; rain, if given, is a vector;
    P, Rs, Rl, lat, zi, cp and sigH may be passed as vectors or constants;
    sensor heights (zu, zt, zq) are only constants.  All vectors must be of
    equal length.

    Default values are assigned for all variables except u,t,rh,ts.  Input
    arrays may contain NaNs to indicate missing values.  Defaults should be set
    to representative regional values if possible.

    Input definitions:

    u = ocean surface relative wind speed (m/s) at height zu(m)
    t = bulk air temperature (degC) at height zt(m)
    rh = relative humidity (%) at height zq(m)
    ts = sea water temperature (degC) - see jcool below
    P = surface air pressure (mb) (default = 1015)
    Rs = downward shortwave radiation (W/m^2) (default = 150)
    Rl = downward longwave radiation (W/m^2) (default = 370)
    zu = wind sensor height (m) (default = 18m)
    zt = bulk temperature sensor height (m) (default = 18m)
    zq = RH sensor height (m) (default = 18m)
    lat = latitude (default = 45 N)
    zi = PBL height (m) (default = 600m)
    rain = rain rate (mm/hr)
    cp = phase speed of dominant waves (m/s)
    sigH =  significant wave height (m)
    jcool = cool skin option (default = 1 for bulk SST)

    Output is a 2-D ndarray with the following variables as 37 columns.
    Other quantities may be added to output by editing lines 536/537.

    col    var     description
    -------------------------------------------------------------------------
    0      usr     friction velocity that includes gustiness (m/s)
    1      tau     wind stress (N/m^2)
    2      hsb     sensible heat flux into ocean (W/m^2)
    3      hlb     latent heat flux into ocean (W/m^2)
    4      hbb     buoyancy flux into ocean (W/m^2)
    5      hsbb    "sonic" buoyancy flux measured directly by sonic anemometer
    6      hlwebb  Webb correction for latent heat flux, add this to directly
                   measured eddy covariance latent heat flux from water vapor
                   mass concentration sensors (e.g. Licor 7500).
    7      tsr     temperature scaling parameter (K)
    8      qsr     specific humidity scaling parameter (g/Kg)
    9      zot     thermal roughness length (m)
    10     zoq     moisture roughness length (m)
    11     Cd      wind stress transfer (drag) coefficient at height zu
    12     Ch      sensible heat transfer coefficient (Stanton number) at ht zu
    13     Ce      latent heat transfer coefficient (Dalton number) at ht zq
    14     L       Obukhov length scale (m)
    15     zet     Monin-Obukhov stability parameter zu/L
    16     dter    cool-skin temperature depression (degC)
    17     dqer    cool-skin humidity depression (degC)
    18     tkt     cool-skin thickness (m)
    19     Urf     wind speed at reference height (user can select height below)
    20     Trf     temperature at reference height
    21     Qrf     specific humidity at reference height
    22     RHrf    relative humidity at reference height
    23     UrfN    neutral value of wind speed at reference height
    24     Rnl     Upwelling IR radiation computed by COARE
    25     Le      latent heat of vaporization
    26     rhoa    density of air
    27     UN      neutral value of wind speed at zu
    28     U10     wind speed adjusted to 10 m
    29     U10N    neutral value of wind speed at 10m
    30     Cdn_10  neutral value of drag coefficient at 10m
    31     Chn_10  neutral value of Stanton number at 10m
    32     Cen_10  neutral value of Dalton number at 10m
    33     RF      rain heat flux (W/m2)
    34     Evap    evaporation (mm/hr)
    35     Qs      sea surface specific humidity (g/kg)
    36     Q10     specific humidity at 10m (g/kg)
    37     RH10    RH at 10m (%)

    Notes:
    1) u is the ocean-relative wind speed, i.e., the magnitude of the
       difference between the wind (at zu) and ocean surface current
       vectors.
    2) Set jcool=0 if ts is true surface skin temperature,
       otherwise ts is assumed the bulk temperature and jcool=1.
    3) The code to compute the heat flux caused by precipitation is
       included if rain data is available (default is no rain).
    4) Code updates the cool-skin temperature depression dter and thickness
       tkt during iteration loop for consistency.
    5) Number of iterations set to nits = 6.
    6) The warm layer is not implemented in this version.

    Reference:

    Fairall, C.W., E.F. Bradley, J.E. Hare, A.A. Grachev, and J.B. Edson (2003),
    Bulk parameterization of air sea fluxes: updates and verification for the
    COARE algorithm, J. Climate, 16, 571-590.

    Code history:

    1) 12/14/05 - created based on scalar version coare26sn.m with input
       on vectorization from C. Moffat.
    2) 12/21/05 - sign error in psiu_26 corrected, and code added to use
       variable values from the first pass through the iteration loop for the
       stable case with very thin M-O length relative to zu (zetu>50) (as is
       done in the scalar coare26sn and COARE3 codes).
    3) 7/26/11 - S = dt was corrected to read S = ut.
    4) 7/28/11 - modification to roughness length parameterizations based
       on the CLIMODE, MBL, Gasex and CBLAST experiments are incorporated
    5) Python translation by BWB, Oct 2014.  Modified to allow user specified
       vectors for lat and zi.  Defaults added for zu, zt, zq.
    """

    import numpy as np
    import meteo
    import util

    # be sure array inputs are ndarray floats
    # if inputs are already ndarray float this does nothing
    # otherwise copies are created in the local namespace
    u = np.copy(np.asarray(u, dtype=float))
    t = np.copy(np.asarray(t, dtype=float))
    rh = np.copy(np.asarray(rh, dtype=float))
    ts = np.copy(np.asarray(ts, dtype=float))
    # these default to 1 element arrays
    P = np.copy(np.asarray(P, dtype=float))
    Rs = np.copy(np.asarray(Rs, dtype=float))
    Rl = np.copy(np.asarray(Rl, dtype=float))
    zi = np.copy(np.asarray(zi, dtype=float))
    lat = np.copy(np.asarray(lat, dtype=float))

    # check for mandatory input variable consistency
    len = u.size
    if not np.all([t.size==len, rh.size==len, ts.size==len]):
        raise ValueError, 'coare35vn: u, t, rh, ts arrays of different length'

    # format optional array inputs
    if P.size != len and P.size != 1:
        raise ValueError, 'coare35vn: P array of different length'
    elif P.size == 1:
        P = P * np.ones(len)

    if Rl.size != len and Rl.size != 1:
        raise ValueError, 'coare35vn: Rl array of different length'
    elif Rl.size == 1:
        Rl = Rl * np.ones(len)

    if Rs.size != len and Rs.size != 1:
        raise ValueError, 'coare35vn: Rs array of different length'
    elif Rs.size == 1:
        Rs = Rs * np.ones(len)

    if zi.size != len and zi.size != 1:
        raise ValueError, 'coare35vn: zi array of different length'
    elif zi.size == 1:
        zi = zi * np.ones(len)

    if lat.size != len and lat.size != 1:
        raise ValueError, 'coare35vn: lat array of different length'
    elif lat.size == 1:
        lat = lat * np.ones(len)

    if rain is not None:
        rain = np.asarray(rain, dtype=float)
        if rain.size != len:
            raise ValueError, 'coare35vn: rain array of different length'

    if cp is not None:
        waveage_flag = True
        cp = np.copy(np.asarray(cp, dtype=float))
        if cp.size != len:
            raise ValueError, 'coare35vn: cp array of different length'
        elif cp.size == 1:
            cp = cp * np.ones(len)
    else:
        waveage_flag = False
        cp = np.nan * np.ones(len)

    if sigH is not None:
        seastate_flag = True
        sigH = np.copy(np.asarray(sigH, dtype=float))
        if sigH.size != len:
            raise ValueError, 'coare35vn: sigH array of different length'
        elif sigH.size == 1:
            sigH = sigH * np.ones(len)
    else:
        seastate_flag = False
        sigH = np.nan * np.ones(len)

    if waveage_flag and seastate_flag:
        print 'Using seastate dependent parameterization'

    if waveage_flag and not seastate_flag:
        print 'Using waveage dependent parameterization'

    # check jcool
    if jcool != 0:
        jcool = 1   # all input other than 0 defaults to jcool=1

    # check sensor heights
    test = [type(zu) is int or type(zu) is float]
    test.append(type(zt) is int or type(zt) is float)
    test.append(type(zq) is int or type(zq) is float)
    if not np.all(test):
        raise ValueError, 'coare35vn: zu, zt, zq, should be constants'
    zu = zu * np.ones(len)
    zt = zt * np.ones(len)
    zq = zq * np.ones(len)

    # input variable u is surface relative wind speed (magnitude of difference
    # between wind and surface current vectors). To follow orginal Fairall
    # code, we set surface current speed us=0. If us data are available
    # construct u prior to using this code.
    us = np.zeros(len)

    # convert rh to specific humidity
    Qs = meteo.qsea(ts,P)/1000.0  # surface water specific humidity (kg/kg)
    Q, Pv = meteo.qair(t,P,rh)    # specific hum. and partial Pv (mb)
    Q /= 1000.0                   # Q (kg/kg)

    # set constants
    zref = 10.          # ref height, m (edit as desired)
    Beta = 1.2
    von  = 0.4          # von Karman const
    fdg  = 1.00         # Turbulent Prandtl number
    tdk  = 273.16
    grav = meteo.grv(lat)

    # air constants
    Rgas = 287.1
    Le   = (2.501 - 0.00237*ts) * 1e6
    cpa  = 1004.67
    cpv  = cpa * (1 + 0.84*Q)
    rhoa = P*100. / (Rgas * (t + tdk) * (1 + 0.61*Q))
    rhodry = (P - Pv)*100. / (Rgas * (t + tdk))
    visa = 1.326e-5 * (1 + 6.542e-3*t + 8.301e-6*t**2 - 4.84e-9*t**3)

    # cool skin constants
    Al   = 2.1e-5 * (ts + 3.2)**0.79
    be   = 0.026
    cpw  = 4000.
    rhow = 1022.
    visw = 1.e-6
    tcw  = 0.6
    bigc = 16. * grav * cpw * (rhow * visw)**3 / (tcw**2 * rhoa**2)
    wetc = 0.622 * Le * Qs / (Rgas * (ts + tdk)**2)

    # net radiation fluxes
    Rns = 0.945 * Rs        #albedo correction
    #    IRup = eps * sigma*T**4 + (1 - eps)*IR
    #    Rnl = IRup - IR
    #    Rnl = eps * sigma*T**4 - eps*IR  as below
    Rnl = 0.97 * (5.67e-8 * (ts - 0.3*jcool + tdk)**4 - Rl) # initial value
    #    IRup = Rnl + IR

    #####     BEGIN BULK LOOP

    # first guess
    du = u - us
    dt = ts - t - 0.0098*zt
    dq = Qs - Q
    ta = t + tdk
    ug = 0.5
    dter = 0.3
    ut = np.sqrt(du**2 + ug**2)
    u10 = ut * np.log(10/1e-4) / np.log(zu/1e-4)
    usr = 0.035 * u10
    zo10 = 0.011 * usr**2 / grav + 0.11*visa / usr
    Cd10 = (von / np.log(10/zo10))**2
    Ch10 = 0.00115
    Ct10 = Ch10 / np.sqrt(Cd10)
    zot10 = 10 / np.exp(von/Ct10)
    Cd = (von / np.log(zu/zo10))**2
    Ct = von / np.log(zt/zot10)
    CC = von * Ct/Cd
    Ribcu = -zu / zi / 0.004 / Beta**3
    Ribu = -grav * zu/ta * ((dt - dter*jcool) + 0.61*ta*dq) / ut**2
    zetu = CC * Ribu * (1 + 27/9  * Ribu/CC)

    k50 = util.find(zetu > 50)   # stable with thin M-O length relative to zu

    k = util.find(Ribu < 0)
    if Ribcu.size == 1:
        zetu[k] = CC[k] * Ribu[k] / (1 + Ribu[k] / Ribcu)
    else:
        zetu[k] = CC[k] * Ribu[k] / (1 + Ribu[k] / Ribcu[k])

    L10 = zu / zetu
    gf = ut / du
    usr = ut * von / (np.log(zu/zo10) - meteo.psiu_40(zu/L10))
    tsr = -(dt - dter*jcool)*von*fdg / (np.log(zt/zot10) -
            meteo.psit_26(zt/L10))
    qsr = -(dq - wetc*dter*jcool)*von*fdg / (np.log(zq/zot10) -
            meteo.psit_26(zq/L10))
    tkt = 0.001 * np.ones(len)

    # The following gives the new formulation for the Charnock variable

    charnC = 0.011 * np.ones(len)
    umax = 19
    a1 = 0.0017
    a2 = -0.0050

    charnC = a1 * u10 + a2
    j = util.find(u10 > umax)
    charnC[j] = a1 * umax + a2

    A = 0.114   # wave-age dependent coefficients
    B = 0.622

    Ad = 0.091  # Sea-state/wave-age dependent coefficients
    Bd = 2.0

    charnW = A * (usr/cp)**B
    zoS = sigH * Ad * (usr/cp)**Bd
    charnS = zoS * grav / usr / usr

    charn = 0.011 * np.ones(len)
    k = util.find(ut > 10)
    charn[k] = 0.011 + (ut[k] - 10) / (18 - 10)*(0.018 - 0.011)
    k = util.find(ut > 18)
    charn[k] = 0.018

    # begin bulk loop
    nits = 10   # number of iterations
    for i in range(nits):
        zet = von*grav*zu / ta*(tsr + 0.61*ta*qsr) / (usr**2)
        if waveage_flag:
            if seastate_flag:
                charn = charnS
            else:
                charn = charnW
        else:
            charn = charnC

        L = zu / zet
        zo = charn*usr**2/grav + 0.11*visa/usr  # surface roughness
        rr = zo*usr/visa
        # thermal roughness lengths give Stanton and Dalton numbers that
        # closely approximate COARE 3.0
        zoq = np.minimum(1.6e-4, 5.8e-5/rr**0.72)
        zot = zoq
        cdhf = von / (np.log(zu/zo) - meteo.psiu_26(zu/L))
        cqhf = von*fdg / (np.log(zq/zoq) - meteo.psit_26(zq/L))
        cthf = von*fdg / (np.log(zt/zot) - meteo.psit_26(zt/L))
        usr = ut*cdhf
        qsr = -(dq - wetc*dter*jcool)*cqhf
        tsr = -(dt - dter*jcool)*cthf
        tvsr = tsr + 0.61*ta*qsr
        tssr = tsr + 0.51*ta*qsr
        Bf = -grav / ta*usr*tvsr
        ug = 0.2 * np.ones(len)

        k = util.find(Bf > 0)
        if zi.size == 1:
            ug[k] = Beta*(Bf[k]*zi)**0.333
        else:
            ug[k] = Beta*(Bf[k]*zi[k])**0.333

        ut = np.sqrt(du**2  + ug**2)
        gf = ut/du
        hsb = -rhoa*cpa*usr*tsr
        hlb = -rhoa*Le*usr*qsr
        qout = Rnl + hsb + hlb
        dels = Rns * (0.065 + 11*tkt - 6.6e-5/tkt*(1 - np.exp(-tkt/8.0e-4)))
        qcol = qout - dels
        alq = Al*qcol + be*hlb*cpw/Le
        xlamx = 6.0 * np.ones(len)
        tkt = np.minimum(0.01, xlamx*visw/(np.sqrt(rhoa/rhow)*usr))
        k = util.find(alq > 0)
        xlamx[k] = 6/(1 + (bigc[k]*alq[k]/usr[k]**4)**0.75)**0.333
        tkt[k] = xlamx[k]*visw / (np.sqrt(rhoa[k]/rhow)*usr[k])
        dter = qcol*tkt/tcw
        dqer = wetc*dter
        Rnl = 0.97*(5.67e-8*(ts - dter*jcool + tdk)**4 - Rl)   # update dter

        # save first iteration solution for case of zetu>50
        if i == 0:
            usr50 = usr[k50]
            tsr50 = tsr[k50]
            qsr50 = qsr[k50]
            L50 = L[k50]
            zet50 = zet[k50]
            dter50 = dter[k50]
            dqer50 = dqer[k50]
            tkt50 = tkt[k50]

        u10N = usr/von/gf*np.log(10/zo)
        charnC = a1*u10N + a2
        k = util.find(u10N > umax)
        charnC[k] = a1*umax + a2
        charnW = A*(usr/cp)**B
        zoS = sigH*Ad*(usr/cp)**Bd - 0.11*visa/usr
        charnS = zoS*grav/usr/usr

    # end bulk loop

    # insert first iteration solution for case with zetu > 50
    usr[k50] = usr50
    tsr[k50] = tsr50
    qsr[k50] = qsr50
    L[k50] = L50
    zet[k50] = zet50
    dter[k50] = dter50
    dqer[k50] = dqer50
    tkt[k50] = tkt50

    # compute fluxes
    tau = rhoa*usr*usr/gf           # wind stress
    hsb = -rhoa*cpa*usr*tsr         # sensible heat flux
    hlb = -rhoa*Le*usr*qsr          # latent heat flux
    hbb = -rhoa*cpa*usr*tvsr        # buoyancy flux
    hsbb = -rhoa*cpa*usr*tssr       # sonic buoyancy flux
    wbar = 1.61*hlb/Le/(1+1.61*Q)/rhoa + hsb/rhoa/cpa/ta
    hlwebb = rhoa*wbar*Q*Le
    Evap = 1000*hlb/Le/1000*3600 # mm/hour

    # compute transfer coeffs relative to ut @ meas. ht
    Cd = tau/rhoa/ut/np.maximum(0.1,du)
    Ch = -usr*tsr/ut/(dt - dter*jcool)
    Ce = -usr*qsr/(dq - dqer*jcool)/ut

    # compute 10-m neutral coeff relative to ut
    Cdn_10 = 1000*von**2 / np.log(10/zo)**2
    Chn_10 = 1000*von**2 * fdg/np.log(10/zo) / np.log(10/zot)
    Cen_10 = 1000*von**2 * fdg/np.log(10/zo) / np.log(10/zoq)

    # compute the stability functions
    zrf_u = 10      # User defined reference heights
    zrf_t = 10
    zrf_q = 10
    psi = meteo.psiu_26(zu/L)
    psi10 = meteo.psiu_26(10/L)
    psirf = meteo.psiu_26(zrf_u/L)
    psiT = meteo.psit_26(zt/L)
    psi10T = meteo.psit_26(10/L)
    psirfT = meteo.psit_26(zrf_t/L)
    psirfQ = meteo.psit_26(zrf_q/L)
    gf = ut/du

    # Determine the wind speeds relative to ocean surface
    # Note that usr is the friction velocity that includes
    # gustiness usr = sqrt(Cd) S, which is equation (18) in
    # Fairall et al. (1996)
    S = ut
    U = du
    S10 = S + usr/von*(np.log(10/zu) - psi10 + psi)
    U10 = S10/gf
    # or U10 = U + usr/von/gf*(np.log(10/zu) - psi10 + psi)
    Urf = U + usr/von/gf*(np.log(zrf_u/zu) - psirf + psi)
    UN = U + psi*usr/von/gf
    U10N = U10 + psi10*usr/von/gf
    UrfN = Urf + psirf*usr/von/gf

    UN2 = usr/von/gf * np.log(zu/zo)
    U10N2 = usr/von/gf * np.log(10/zo)
    UrfN2 = usr/von/gf * np.log(zrf_u/zo)

    # rain heat flux after Gosnell et al., JGR, 1995
    if rain is None:
        RF = np.zeros(usr.size)
    else:
        # water vapour diffusivity
        dwat = 2.11e-5*((t + tdk)/tdk)**1.94
        # heat diffusivity
        dtmp = (1 + 3.309e-3*t - 1.44e-6*t**2) * 0.02411/(rhoa*cpa)
        # Clausius-Clapeyron
        dqs_dt = Q*Le / (Rgas*(t + tdk)**2)
        # wet bulb factor
        alfac = 1/(1 + 0.622*(dqs_dt*Le*dwat)/(cpa*dtmp))
        RF = rain*alfac*cpw*((ts-t-dter*jcool)+(Qs-Q-dqer*jcool)*Le/cpa)/3600

    lapse = grav/cpa
    SST = ts - dter*jcool

    T = t
    T10 = T + tsr/von*(np.log(10/zt) - psi10T + psiT) + lapse*(zt - 10)
    Trf = T + tsr/von*(np.log(zrf_t/zt) - psirfT + psiT) + lapse*(zt - zrf_t)
    TN = T + psiT*tsr/von
    T10N = T10 + psi10T*tsr/von
    TrfN = Trf + psirfT*tsr/von

    TN2 = SST + tsr/von * np.log(zt/zot) - lapse*zt
    T10N2 = SST + tsr/von * np.log(10/zot) - lapse*10;
    TrfN2 = SST + tsr/von * np.log(zrf_t/zot) - lapse*zrf_t

    dqer = wetc*dter*jcool
    SSQ = Qs - dqer
    SSQ = SSQ*1000
    Q = Q*1000
    qsr = qsr*1000
    Q10 = Q + qsr/von*(np.log(10/zq) - psi10T + psiT)
    Qrf = Q + qsr/von*(np.log(zrf_q/zq) - psirfQ + psiT)
    QN = Q + psiT*qsr/von/np.sqrt(gf)
    Q10N = Q10 + psi10T*qsr/von
    QrfN = Qrf + psirfQ*qsr/von

    QN2 = SSQ + qsr/von * np.log(zq/zoq)
    Q10N2 = SSQ + qsr/von * np.log(10/zoq)
    QrfN2 = SSQ + qsr/von * np.log(zrf_q/zoq)
    RHrf = meteo.rhcalc(Trf, P, Qrf/1000)
    RH10 = meteo.rhcalc(T10, P, Q10/1000)

    # output: edit these lists to add other values as desired
#     list1 = [usr,tau,hsb,hlb,hbb,hsbb,hlwebb,tsr,qsr,zot,zoq,Cd,Ch,Ce,L,zet]
#     list2 = [dter,dqer,tkt,Urf,Trf,Qrf,RHrf,UrfN,Rnl,Le,rhoa,UN,U10,U10N]
#     list3 = [Cdn_10,Chn_10,Cen_10,RF,Evap,Qs,Q10,RH10]
#     out = tuple(list1 + list2 + list3)
    # basic default output...
    list1 = [usr,tau,hsb,hlb,hlwebb,tsr,qsr,zot,zoq,Cd,Ch,Ce,L,zet]
    list2 = [dter,dqer,tkt,RF,Cdn_10,Chn_10,Cen_10]
    out = tuple(list1 + list2)
    A = np.column_stack(out)
    return A





# This code executes if 'run bulk.py' is executed from iPython cmd line
# Edit line 533 to indicate path to test data file
if __name__ == '__main__':
    import numpy as np
    import os
    import util
    import matplotlib.pyplot as plt

    path = '/Volumes/MyPassport/pyCOARE_NOAA'
    fil = 'test_35_data.txt'
    cols = 15
    data, varNames = util.load_txt_file(path,fil,cols)
    u = data[:,0]
    ta = data[:,2]
    rh = data[:,4]
    Pa = data[:,6]
    ts = data[:,7]
    rs = data[:,8]
    rl = data[:,9]
    Lat = data[:,10]
    ZI = data[:,11]
    Rain = data[:,12]

    A = coare35vn(u, ta, rh, ts, P=Pa, Rs=rs, Rl=rl, zu=16, zt=16, zq=16,
                lat=Lat, zi=ZI, rain=Rain, jcool=1)
    fnameA = os.path.join(path,'test_35_output_py_082020.txt')
    A_hdr = 'usr\ttau\thsb\thlb\thlwebb\ttsr\tqsr\tzot\tzoq\tCd\t'
    A_hdr += 'Ch\tCe\tL\tzet\tdter\tdqer\ttkt\tRF\tCdn_10\tChn_10\tCen_10'
    np.savetxt(fnameA,A,fmt='%.18e',delimiter='\t',header=A_hdr)

