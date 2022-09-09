"""
Functions for COARE model version 3.6 bulk flux calculations with warm layer computations.

Translated from MATLAB scripts written by Jim Edson and Chris Fairall, 
further edited by Elizabeth Thompson. Main Matalb script source is coare36vnWarm_et.m

Execute '%run coare36vnWarm_et' from the iPython command line for test run with
'test_36_data.txt' input data file. 
List of functions in this code are:
    ['ComputeEsat',
    'scalarv',
    'coare36vnWarm_et']
The rest of the functions are imported via the module 'coare36vn_zrf_et'

ludovic Bariteau, CU/CIRES, NOAA/ESRL/PSL
v1: September 2022
"""

import numpy as np
import coare36vn_zrf_et as c36
import os

def coare36vnWarm_et(Jd, U, Zu, Tair, Zt, RH, Zq, P, Tsea, SW_dn, LW_dn, Lat, Lon, Zi, Rainrate, Ts_depth, Ss, cp=None, sigH=None,zrf_u = 10.0,zrf_t = 10.0,zrf_q = 10.0): 
    print('WarmCoolLayer')
#***********   input data **************
#       Jd = day-of-year or julian day
#	    U = wind speed magnitude (m/s) corrected for currents, i.e. relative to water at height zu
#	    Zu = height (m) of wind measurement
#	    Tair = air temp (degC) at height zt
#	    Zt = height (m) of air temperature measurement
#	    RH = relative humidity (#) at height zq
#	    Zq = height (m) of air humidity measurement
#	    P = air pressure at sea level (mb)
#       Tsea = surface sea temp (degC) at ts_depth
#	    SW_dn = downward solar flux (w/m^2) defined positive down
#	    LW_dn = downward IR flux (w/m^2) defined positive down
#	    Lat = latitude (deg N=+)
#	    Lon = longitude (deg E=+) # If using other version, see
#               usage of lon in albedo_vector function and adjust 'E' or
#               'W' string input
#       Zi = inversion height (m)
#       Rainrate = rain rate (mm/hr)
#       Ts_depth = depth (m) of water temperature measurement, positive for below
#                   surface
#       Ss = sea surface salinity (PSU)
#       cp = phase speed of dominant waves (m/s)
#       sigH = significant wave height (m)
#       zu, zt, zq = heights of the observations (m)
#       zrf_u, zrf_t, zrf_q = reference height for profile.
#           Use this to compare observations at different heights
    
    
#********** output data  ***************
# Outputs
# Adds onto output from coare36vn_zrf_et
# .... see that function for updated output. It can change. This function adds 4 variables onto it:
# previously named dt_wrm, tk_pwp, dsea, du_wrm... renamed to the following
# for clarity:
#         dT_warm = dT from base of warm layer to skin, i.e. warming across entire warm layer depth (deg C)
#         dz_warm = warm layer thickness (m)
#         dT_warm_to_skin = dT from measurement depth to skin due to warm layer,
#                       such that Tskin = tsea + dT_warm_to_skin - dT_skin
#         du_warm = total current accumulation in warm layer (m/s ?...  unsure of units but likely m/s)
    
#********** history ********************
# updated 09/2020 for consistency with units, readme info, and coare 3.6 main function
#    Changed names for a few things...
#    dt_wrm -> dT_warm; tk_pwp -> dz_warm; dsea -> dT_warm_to_skin;
#    skin: dter -> dT_skin
#    and fluxes to avoid ambiguity: RF -> hrain
#    Rnl -> lw_net; Rns -> sw_net; Rl -> lw_dn; Rs -> sw_dn;
#  updated 05/2022 - fix a glitch in code when using wave input - LB
#   added     cpi=cp[ibg];   # air density
#             sigHi=sigH[ibg];   # air density
# and edited line 258 to use cpi and sigHi
# Bx=coare36vn_zrf_et(u,zu,t,zt,rh,zq,p,ts,sw_dn,lw_dn,lat,lon,jd,zi,rain,ss,cpi,sigHi,zrf_u,zrf_t,zrf_q);
    
#********** Set cool skin options ******************
    jcool = 1
    icount = 1
#********** Set wave ******************
### ... not sure if this is necessary ...
# if no wave info is provided, fill array with nan.
# if length(cp) == 1 && isnan(cp) == 1
#     cp = jd*nan;
# end
# if length(sigH) == 1 && isnan(sigH) == 1
#     sigH = jd*nan;
# end
    
    # be sure array inputs are ndarray floats for single value function
    # if inputs are already ndarray float this does nothing
    # otherwise copies are created in the local namespace
    if U.size ==1 and Tair.size ==1: 
        U = np.copy(np.asarray([U], dtype=float))
        Zu = np.copy(np.asarray([Zu], dtype=float))
        Tair = np.copy(np.asarray([Tair], dtype=float))
        Zt = np.copy(np.asarray([Zt], dtype=float))
        RH = np.copy(np.asarray([RH], dtype=float))
        Zq = np.copy(np.asarray([Zq], dtype=float))
        P = np.copy(np.asarray([P], dtype=float))
        Tsea = np.copy(np.asarray([Tsea], dtype=float))
        SW_dn = np.copy(np.asarray([SW_dn], dtype=float))
        LW_dn = np.copy(np.asarray([LW_dn], dtype=float))
        Lat = np.copy(np.asarray([Lat], dtype=float))
        Lon = np.copy(np.asarray([Lon], dtype=float))
        Jd = np.copy(np.asarray([Jd], dtype=float))
        Zi = np.copy(np.asarray([Zi], dtype=float))
        Rainrate = np.copy(np.asarray([Rainrate], dtype=float))
        Ss = np.copy(np.asarray([Ss], dtype=float))
        zrf_u = np.copy(np.asarray([zrf_u], dtype=float))
        zrf_t = np.copy(np.asarray([zrf_t], dtype=float))
        zrf_q = np.copy(np.asarray([zrf_q], dtype=float))
        Ts_depth = np.copy(np.asarray([Ts_depth], dtype=float))

    N = np.size(U)
 
    if cp is not None and cp.size==1:
        cp = np.copy(np.asarray([cp], dtype=float))
    elif cp is None:
        cp = np.nan * np.ones(N)
    
    if sigH is not None and sigH.size==1:
        sigH = np.copy(np.asarray([sigH], dtype=float))
    elif sigH is None:
        sigH = np.nan * np.ones(N)

#*********************  housekeep variables  ********
# Call coare36vn to get initial flux values
    Bx = c36.coare36vn_zrf_et(U[0],Zu[0],Tair[0],Zt[0],RH[0],Zq[0],P[0],Tsea[0],SW_dn[0],LW_dn[0],Lat[0],Lon[0],Jd[0],Zi[0],Rainrate[0],Ss[0],cp[0],sigH[0],zrf_u,zrf_t,zrf_q)

    ### check these indices for you latest version of coare!
    tau_old = Bx[0,1]
    hs_old = Bx[0,2]
    hl_old = Bx[0,3]
    dT_skin_old = Bx[0,17]
    hrain_old = Bx[0,37]
    
    qcol_ac = 0.0
    tau_ac = 0.0
    dT_warm = 0.0
    du_warm = 0.0
    max_pwp = 19.0
    dz_warm = max_pwp
    dT_warm_to_skin = 0.0
    q_pwp = 0.0
    fxp = 0.5
    
    rich = 0.65
    
    jtime = 0
    jamset = 0
    jump = 1
    
    #*******************  set constants  ****************
    T2K = 273.16
    Rgas = 287.1
    cpa = 1004.67
    cpw = 4000.0
    rhow = 1022.0
    visw = 1e-06
    
    be = 0.026
    tcw = 0.6
    #**********************************************************
    #******************  setup read data loop  ****************
    P_tq = P - (0.125 * Zt)
    Press,tseak,tairk,Qsatsea,Qsat,Qair,Rhoair,Rhodry = scalarv(P,Tsea,Tair,RH,Zt)

    nx = np.size(Jd)
    
    # this is an empty array for saving warm layer code output values. Will be
    # added to coare output at the end.
    warm_output = np.full([nx,4],np.nan)
    for ibg in np.arange(0,nx).reshape(-1):  #np.arange(0,53).reshape(-1):
        jd = Jd[ibg]
        p = P[ibg]
        u = U[ibg]
        tsea = Tsea[ibg]
        t = Tair[ibg]
        ss = Ss[ibg]
        qs = Qsatsea[ibg]
        q = Qair[ibg]
        rh = RH[ibg]
        sw_dn = SW_dn[ibg]
        lw_dn = LW_dn[ibg]
        rain = Rainrate[ibg]
        grav = c36.grv(Lat[ibg])
        lat = Lat[ibg]
        lon = Lon[ibg]
        rhoa = Rhoair[ibg]
        cpi = cp[ibg]
        sigHi = sigH[ibg]
        zi=Zi[ibg]
        zu=Zu[ibg]
        zt=Zt[ibg]
        zq=Zq[ibg]
        ts_depth=Ts_depth[ibg]
        
        #*****  variables for warm layer  ***
        ### for constant albedo
        # sw_net=.945*sw_dn;     #Net Solar: positive warming ocean, constant albedo
        ### for albedo that is time-varying, i.e. zenith angle varying
        # insert 'E' for input to albedo function if longitude is defined positive
        # to E, so that lon can be flipped. The function ideally works with
        # longitude positive to the west. Check: albedo should peak at sunrise not
        # sunset.
        alb,T_sw,solarmax_sw,psi_sw = c36.albedo_vector(sw_dn,jd,lon,lat,'E')
        sw_net = np.multiply((1 - alb[0]),sw_dn)
        lw_net = 0.97 * (5.67e-08 * (tsea - dT_skin_old * jcool + T2K) ** 4 - lw_dn)
        cpv = cpa * (1 + 0.84 * q / 1000)
        visa = 1.326e-05 * (1 + 0.006542 * t + 8.301e-06 * t * t - 4.84e-09 * t * t * t)
        Al = 2.1e-05 * (tsea + 3.2) ** 0.79
        ctd1 = np.sqrt(2 * rich * cpw / (Al * grav * rhow))
        ctd2 = np.sqrt(2 * Al * grav / (rich * rhow)) / (cpw ** 1.5)
        
        #********************************************************
        #****  Compute apply warm layer  correction *************
        #********************************************************
        intime = jd - np.fix(jd)
        loc = (lon + 7.5) / 15
        chktime = loc + intime * 24
        if chktime > 24:
            chktime = chktime - 24
        newtime = (chktime - 24 * np.fix(chktime / 24)) * 3600
        if icount > 1:
            if newtime <= 21600 or jump == 0:
                jump = 0
                if newtime < jtime:
                    jamset = 0
                    fxp = 0.5
                    dz_warm = max_pwp
                    tau_ac = 0.0
                    qcol_ac = 0.0
                    dT_warm = 0.0
                    du_warm = 0.0
                else:
                    #************************************
                    #****   set warm layer constants  ***
                    #************************************
                    dtime = newtime - jtime
                    qr_out = lw_net + hs_old + hl_old + hrain_old
                    q_pwp = fxp * sw_net - qr_out
                    # qqrx[ibg] = hs_old
                    if q_pwp >= 50 or jamset == 1:
                        jamset = 1
                        tau_ac = tau_ac + np.maximum(0.002,tau_old) * dtime
                        if qcol_ac + q_pwp * dtime > 0:
                            #******************************************
                            # Compute the absorption profile
                            #******************************************
                            for i in np.arange(1,5+1).reshape(-1):
                                #### The original version since Fairall et al. 1996:
                                fxp = 1 - (0.28 * 0.014 * (1 - np.exp(- dz_warm / 0.014)) + 0.27 * 0.357 * (1 - np.exp(- dz_warm / 0.357)) + 0.45 * 12.82 * (1 - np.exp(- dz_warm / 12.82))) / dz_warm
                                # the above integrated flux formulation is wrong for the warm layer,
                                # but it has been used in this scheme since 1996 without
                                # making bad predictions.
                                # Simon recognized that fxp should be the
                                # fraction absorbed in the layer of the form
                                # 1-sum(ai*exp(-tk_pwp/gammai)) with sum(ai)=1
                                # One can choose different exponential
                                # absorption parameters, but it has to be of the right form.
                                # Simon idealized coefficients from profiles of
                                # absorption in DYNAMO by C. Ohlmann.
                                # see /Users/sdeszoek/Data/cruises/DYNAMO_2011/solar_absorption/test_absorption_fcns.m
                                # Correct form of Soloviev 3-band absorption:
                                # --using original F96 absorption bands:
                                #  fxp=1-(0.32*exp(-tk_pwp/22.0) + 0.34*exp(-tk_pwp/1.2) + 0.34*exp(-tk_pwp/0.014)); # NOT TESTED!!
                                # --using DYNAMO absorption bands (F, invgamma defined above):
                                #### NOT TESTED!! Correction of fxp from Simon ***
                                #fxp=1-sum(F.*(exp(-tk_pwp*invgamma)),2);
                                qjoule = (fxp * sw_net - qr_out) * dtime
                                if qcol_ac + qjoule > 0:
                                    dz_warm = np.minimum(max_pwp,ctd1 * tau_ac / np.sqrt(qcol_ac + qjoule))
                        else:
                            fxp = 0.75
                            dz_warm = max_pwp
                            qjoule = (fxp * sw_net - qr_out) * dtime
                        qcol_ac = qcol_ac + qjoule
                        #*******  compute dt_warm  ******
                        if qcol_ac > 0:
                            dT_warm = ctd2 * (qcol_ac) ** 1.5 / tau_ac
                            du_warm = 2 * tau_ac / (dz_warm * rhow)
                        else:
                            dT_warm = 0
                            du_warm = 0
                # Compute warm layer dT between input measurement and skin layer
                if dz_warm < ts_depth:
                    dT_warm_to_skin = dT_warm
                else:
                    dT_warm_to_skin = dT_warm * ts_depth / dz_warm
        jtime = newtime
        
        #************* output from routine  *****************************
        # Adjust tsea for warm layer above measurement. Even if tsea is at 5 cm for the sea snake,
        # there could be warming present between the interface and the subsurface
        # temperature. dT_warm_to_skin estimates that warming between the levels.
        ts = tsea + dT_warm_to_skin
        # Rerun COARE with the warm-layer corrected tsea temperature. COARE will
        # apply a cool skin to this, completing all calculations needed for Tskin and fluxes.
        # Using COARE ouput from this function, Tskin = Tsnake - dT_skin + dT_warm_to_skin
        # note: in prior COARE lingo/code: dT_warm_to_skin used to be dsea and dT_skin used to be dter
        Bx = c36.coare36vn_zrf_et(u,zu,t,zt,rh,zq,p,ts,sw_dn,lw_dn,lat,lon,jd,zi,rain,ss,cpi,sigHi,zrf_u,zrf_t,zrf_q)
        # save values from this time step to be used in next time step, this is how
        # the integral is computed
        tau_old = Bx[0,1]
        hs_old = Bx[0,2]
        hl_old = Bx[0,3]
        dT_skin_old = Bx[0,17]
        hrain_old = Bx[0,37]
        warm_output[ibg,0] = dT_warm
        warm_output[ibg,1] = dz_warm
        warm_output[ibg,2] = dT_warm_to_skin
        warm_output[ibg,3] = du_warm
        
        icount = icount + 1
    #end of data line loop 
    
    # get rid of filled values where nans are present in input data
    bad_input = np.where(np.isnan(sw_dn) == 1)
    # disp(['bad solar values = ' sprintf('#i',length(bad_input))]);
    warm_output[bad_input,:] = np.nan
    
    #**************************************************
    # Recompute entire time series of fluxes with seawater T adjusted for warm layer
    #**************************************************
    del Bx
    Tsea = Tsea + warm_output[:,2]
    Bx = c36.coare36vn_zrf_et(U,Zu,Tair,Zt,RH,Zq,P,Tsea,SW_dn,LW_dn,Lat,Lon,Jd,Zi,Rainrate,Ss,cp,sigH,zrf_u,zrf_t,zrf_q)
    B = np.hstack((Bx,warm_output))
    
    #************* output from routine  *****************************
    ### adds to coarevn_zrf output the following 4 vars:
    # B = [ <<< outputs from main coare function >>> .... dT_warm  dz_warm  dT_warm_to_skin  du_warm ]
    return B
    
#------------------------------------------------------------------------------  
def scalarv(P0 = None,tsea = None,tair = None,rh = None,zt = None): 
# Compute the require scalar variables for the bulk code
# Vectorized when needed
# Inputs:
# P0    Air pressure (mb)
# tsea  sea temp (C)
# tair  Air temperature (C)
# rh    Relative humidity (#)
    
    Press = P0 * 100
    P_tq = 100 * (P0 - (0.125 * zt))
    tseak = tsea + 273.15
    tairk = tair + 273.15
    
    if np.size(tsea) > 1:
        #**********************COMPUTES QSEA***************************
        Exx = 6.1121 * np.exp(17.502 * tsea / (240.97 + tsea))
        Exx = np.multiply(Exx,(1.0007 + P0 * 3.46e-06))
        Esatsea = Exx * 0.98
        Qsatsea = 0.622 * Esatsea / (P0 - 0.378 * Esatsea) * 1000
        #**********************COMPUTES QAIR***************************
        Exx = 6.1121 * np.exp(17.502 * tair / (240.97 + tair))
        Exx = np.multiply(Exx,(1.0007 + P_tq * 3.46e-06))
        Qsatms = 0.622 * Exx / (P_tq - 0.378 * Exx) * 1000
        Ems = np.multiply(Exx,rh) / 100
        Qms = 0.622 * Ems / (P_tq - 0.378 * Ems) * 1000
        E = Ems * 100
        #******************COMPUTES AIR DENSITY*******************
        Rhoair = P_tq / (np.multiply(tairk,(1 + 0.61 * Qms / 1000)) * 287.05)
        Rhodry = (P_tq - E) / (tairk * 287.05)
    else:
        #**********************COMPUTES QSEA***************************
        Ex = ComputeEsat(tsea,P0)
        Esatsea = Ex * 0.98
        Qsatsea = 0.622 * Esatsea / (P0 - 0.378 * Esatsea) * 1000
        #**********************COMPUTES QAIR***************************
        Esatms = ComputeEsat(tair,P_tq)
        Qsatms = 0.622 * Esatms / (P_tq - 0.378 * Esatms) * 1000
        Ems = Esatms * rh / 100
        Qms = 0.622 * Ems / (P_tq - 0.378 * Ems) * 1000
        E = Ems * 100
        #******************COMPUTES AIR DENSITY*******************
        Rhoair = P_tq / (tairk * (1 + 0.61 * Qms / 1000) * 287.05)
        Rhodry = (P_tq - E) / (tairk * 287.05)
    
    return Press,tseak,tairk,Qsatsea,Qsatms,Qms,Rhoair,Rhodry
    
#------------------------------------------------------------------------------
def ComputeEsat(T = None,P = None): 
#   Given temperature (C) and pressure (mb), returns
#   saturation vapor pressure (mb).
    Exx = 6.1121 * np.exp(17.502 * T / (240.97 + T))
    Exx = np.multiply(Exx,(1.0007 + P * 3.46e-06))
    return Exx
    
#------------------------------------------------------------------------------

# This code executes if 'run coare36vnWarm_et.py' is executed from iPython cmd line
# Edit line 403 to indicate path to test data file
if __name__ == '__main__':
    # import numpy as np
    # import os
    # import util
    # import matplotlib.pyplot as plt
    
    path = '/Users/ludo/Documents/Work/COARE/conversion2python_tests/'
    fil = 'test_36_data.txt'   
    data = np.genfromtxt(path+fil, skip_header=1)
    U = data[:,1]
    Tair = data[:,3]
    RH = data[:,5]
    P = data[:,7]
    Tsea = data[:,8]
    SW_dn = data[:,9]
    LW_dn = data[:,10]
    Lat = data[:,11]
    Lon = data[:,12]
    Zi = data[:,13]
    Rainrate = data[:,14]
    Zu= data[:,2]
    Zt= data[:,4]
    Zq= data[:,6]
    zrf_u=10.0;
    zrf_t=10.0;
    zrf_q=10.0;
    Ss = data[:,15]
    Jd = data[:,0]
    cp = data[:,16]
    sigH = data[:,17]
    Tsg = data[:,18]
    Ts_depth = data[:,19]
  
    # A=coare36vnWarm_et(Jd, U, Zu, Tair, Zt, RH, Zq, P, Tsg, SW_dn, LW_dn, Lat, Lon, Zi, Rainrate, Ts_depth, Ss, cp = None, sigH = None, zrf_u, zrf_t, zrf_q)
    fnameA = os.path.join(path,'test_36_output_py_082022_withwavesinput_withwarmlayer.txt')
    A=coare36vnWarm_et(Jd, U, Zu, Tair, Zt, RH, Zq, P, Tsg, SW_dn, LW_dn, Lat, Lon, Zi, Rainrate, Ts_depth, Ss, None, None, zrf_u, zrf_t, zrf_q)
    fnameA = os.path.join(path,'test_36_output_py_082022_withnowavesinput_withwarmlayer.txt')
    A_hdr = 'usr\ttau\thsb\thlb\thbb\thlwebb\ttsr\tqsr\tzo\tzot\tzoq\tCd\t'
    A_hdr += 'Ch\tCe\tL\tzeta\tdT_skinx\tdq_skinx\tdz_skin\tUrf\tTrf\tQrf\t'
    A_hdr += 'RHrf\tUrfN\tTrfN\tQrfN\tlw_net\tsw_net\tLe\trhoa\tUN\tU10\tU10N\t'
    A_hdr += 'Cdn_10\tChn_10\tCen_10\thrain\tQs\tEvap\tT10\tT10N\tQ10\tQ10N\tRH10\t'
    A_hdr += 'P10\trhoa10\tgust\twc_frac\tEdis\tdT_warm\tdz_warm\tdT_warm_to_skin\tdu_warm'
    np.savetxt(fnameA,A,fmt='%.18e',delimiter='\t',header=A_hdr)
    
    # test on signle value
    # A=coare36vnWarm_et(Jd[0], U[0], Zu[0], Tair[0], Zt[0], RH[0], Zq[0], P[0], Tsg[0], SW_dn[0], LW_dn[0], Lat[0], Lon[0], Zi[0], Rainrate[0], Ts_depth[0], Ss[0], cp[0], sigH[0], zrf_u, zrf_t, zrf_q)



