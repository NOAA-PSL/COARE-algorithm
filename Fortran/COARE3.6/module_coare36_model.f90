module coare36_module

  use machine , only : kind_phys
  use module_coare36_parameters

  implicit none

contains

  subroutine coare36flux_coolskin( &
    ! input arguments
    u,zu,t,zt,rh,zq,P,ts,sw_dn,lw_dn,grav,zi,rain,Ss,cp,sigH, &
    ! output arguments
    tau,hsb,hlb,Le,sw_net,lw_dn_abs,lw_up,dT_skinx,dz_skin, &
    output_vec)

! -------------------------------------------------------------------------- !
!   This Fortran code is based on the COARE 3.6 Matlab code (Fairall et al,
! 2003; Elizabeth Thompson, personal communication, Feb & Mar 2021). The
! latter is a modified version based on the CLIMODE, MBL, CBLAST, ATOMIC
! experiments (Edson et al. 2012; Elizabeth Thompson et al., personal
! ocmmunication, Feb & Mar 2021)
!   The COARE 3.6 algorithm in the public-released UFS (uncoupled; basically
! GFSv16beta) for calculating latent heat flux (hlb), sensible heat flux
! (hsb), wind stress (tau), cool-skin temperature correction (dT_skinx),
! cool-skin thickness (dz_skin), and latent heat of vaporation (Le; as a
! by-product).
!   The cool-skin schemne of the COARE 3.6 is implemented here because the
! cool-skin scheme in the GFSv16beta Near-Surface Sea Temperature (NSST)
! scheme appears to be older.
!   The input arguments Ss (sea surface salinity), cp (dominant wave speed),
! and sigH (signiciant wave height) are missing in the uncoupled UFS, but
! they are kept here for future full coupling work. They are currently
! assigned constant values inside this subroutine.
! The output dT_skinx
!   Many variables in the COARE 3.6 Matlabe code are discarded here, because
! the UFS or GFS does not need them and it is good practice to save time &
! space in a GCM run. The dicarded variables, however, may be needed for
! comparison with observations. If so, they need to be added in later.
! -------------------------------------------------------------------------- !

! ********************************************************************
!   An important component of this code is whether the inputed ts 
! represents the skin temperature of a near surface temperature. In the
! COARE 3.6 Matlabe code, there is a flag called "jcool" to determine
! whether ts is true skin temperature (jcool=0) or ocean bulk
! temperature (jcool=1). Here we assume the SST in the UFS (or GFS)
! is always ocean bulk temperature, so the jcool is equal to 1 and
! dropped in the associated calculations.
!   If one does not wish to run cool-skin scheme (only energy fluxes are
! wanted), then coare36flux, instead of coare36flux_coolskin should be
! called.
! ********************************************************************

    implicit none

    ! input arguments
    real(kind=kind_phys), intent(in) :: u, zu, t, zt, rh, zq, P, ts,    &
      sw_dn, lw_dn, grav, zi, rain, Ss, cp
    ! input/output arguments
    real(kind_phys), intent(inout) :: sigH
    ! output arguments
    real(kind=kind_phys), intent(out) :: tau, hsb, hlb, Le, &
      sw_net, lw_dn_abs, lw_up, dT_skinx, dz_skin

    !bloss: add output_vec for full output tests against the python code
    real(kind=kind_phys), dimension(21), intent(out), optional :: output_vec

    ! local variables
    real :: jcool, Tf, us, Qs, P_tq, Q, Pv, zos1, cpv, rhoa, &
      visa, tsw, Al35, Al0, Al, be, bigc, &
      lw_net, du, dT, dq, ta, gust, dT_skin, dT_skin_pre, ut, u10, usr, &
      zo10, Cd10, &
      Ch10, Ct10, zot10, Cd, Ct, CC, Ribcu, Ribu, zetu, L10, usr_no_gust, &
      tsr, qsr, wetc, charnC, charn, charnS, dq_skin, umax, &
      a1, a2, Ad, Bd, zoS2, nits, zet, L, zo, rr, zoq, zot, cdhf, &
      cqhf, cthf, tvsr, tssr, Bf, qout, dwat, dtmp, dqs_dt, alfac, &
      dels, qcol, alq, xlamx, usr50, tsr50, qsr50, &
      L50, zet50, dT_skin50, dq_skin50, tkt50, u10N
    integer :: ice, i 
    integer :: k50 = 0 !bloss: Initial value

    integer :: ii !bloss: output_vec index

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Input:  
!
!     u = water-relative mean wind speed magnitude (m/s) at height zu (m)
!         ... if not available, use true wind speed to compute fluxes in
!         earth-coordinates only
!         which will be ignoring the stress contribution by ocean
!         current, which affects all other fluxes too.
!     t = air temperature (degC) at height zt (m)
!     q = specific humidity at height zq (m)
!     P = surface sea level air pressure (mb) (typical value = 1013)
!    ts = seawater temperature (degC)
! sw_dn = downward (positive) shortwave radiation (W/m^2) (typical value = 150) 
! lw_dn = downward (positive) longwave radiation (W/m^2) (typical value = 370)
!  grav = gravitational force, varying with latitude
!    zi = PBL height (m) (default or typical value = 600m)
!  rain = rain rate (mm/hr)
!    Ss = sea surface salinity (PSU)
!    cp = phase speed of dominant waves (m/s)  
!  sigH = significant wave height (m)
!  zu, zt, zq heights of the u, t, rh (m)
!
! The user controls the output at the end of the code.

! Output:  
!
!    tau  = wind stress that includes gustiness (N/m^2)
!    hsb  = sensible heat flux (W/m^2) ... positive for Tair < Tskin
!    hlb  = latent heat flux (W/m^2) ... positive for qair < qs
!     Le  = latent heat of vaporation (J/K)
! dT_skin = cool-skin temperature depression (degC), pos value means skin is
!           cooler than subskin
! dz_skin = cool-skin thickness (m)

! Others:
!
!% radiation signs: positive warms the ocean
!% sensible, rain, and latent flux signs: positive cools the ocean
!% NOTE: signs change throughout the program for ease of calculations. The
!% units noted here are for the final outputs.
!
!%    usr = friction velocity that includes gustiness (m/s)
!%    tau = wind stress that includes gustiness (N/m^2)
!%    hsb = sensible heat flux (W/m^2) ... positive for Tair < Tskin
!%    hlb = latent heat flux (W/m^2) ... positive for qair < qs
!%    hbb = atmospheric buoyany flux (W/m^2)... positive when hlb and hsb heat
!the atmosphere
!%   hsbb = atmospheric buoyancy flux from sonic ... as above, computed with
!sonic anemometer T
!% hlwebb = webb factor to be added to hl covariance and ID latent heat fluxes
!%    tsr = temperature scaling parameter (K)
!%    qsr = specific humidity scaling parameter (g/kg)
!%     zo = momentum roughness length (m) 
!%    zot = thermal roughness length (m) 
!%    zoq = moisture roughness length (m)
!%     Cd = wind stress transfer (drag) coefficient at height zu   
!%     Ch = sensible heat transfer coefficient (Stanton number) at height zu   
!%     Ce = latent heat transfer coefficient (Dalton number) at height zu
!%      L = Monin-Obukhov length scale (m) 
!%    zet = Monin-Obukhov stability parameter zu/L (dimensionless)
!%dq_skin = cool-skin humidity depression (g/kg)
!% lw_net = Net IR radiation computed by COARE (W/m2)... positive heating ocean
!% sw_net = Net solar radiation computed by COARE (W/m2)... positive heating
!ocean
!%   rhoa = density of air at input parameter height zt, typically same as zq
!(kg/m3)
!%   gust = gustiness velocity (m/s)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! JWWang, 2021/3/10
!   The following variables exist in COARE 3.6 Matlabe code but discarded
! in the current UFS or GFS.
!%    Urf = wind speed at reference height (user can select height at input)
!%    Trf = air temperature at reference height
!%    Qrf = air specific humidity at reference height
!%   RHrf = air relative humidity at reference height
!%   UrfN = neutral value of wind speed at reference height
!%     UN = neutral value of wind speed at zu (m/s)
!%    U10 = wind speed adjusted to 10 m (m/s)
!%   UN10 = neutral value of wind speed at 10m (m/s)
!% Cdn_10 = neutral value of drag coefficient at 10m (unitless)
!% Chn_10 = neutral value of Stanton number at 10m (unitless)
!% Cen_10 = neutral value of Dalton number at 10m (unitless)
!%  hrain = rain heat flux (W/m^2)... positive cooling ocean
!%     Qs = sea surface specific humidity, i.e. assuming saturation (g/kg)
!%   Evap = evaporation rate (mm/h)
!%    T10 = air temperature at 10m (deg C)
!%    Q10 = air specific humidity at 10m (g/kg)
!%   RH10 = air relative humidity at 10m (%)
!%    P10 = air pressure at 10m (mb)
!% rhoa10 = air density at 10m (kg/m3)
!%wc_frac = whitecap fraction (ratio)
!%   Edis = energy dissipated by wave breaking (W/m^2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!% Notes: 1) u is the surface-relative wind speed, i.e., the magnitude of the
!%           difference between the wind (at zu) and ocean surface current 
!%           vectors.
!%        3) Assign a default value to P, lw_dn, sw_dn, zi if unknown.
!%        4) Code updates the cool-skin temperature depression dT_skin and
!%           thickness dz_skin during iteration loop for consistency.
!%        5) Number of iterations set to nits = 10 or escape from the loop
!            when abs(dT_skin) is too small.
!
!% Reference:
!%
!%  Fairall, C.W., E.F. Bradley, J.E. Hare, A.A. Grachev, and J.B. Edson (2003),
!%  Bulk parameterization of air sea fluxes: updates and verification for the 
!%  COARE algorithm, J. Climate, 16, 571-590.
!%
!%  Edson, J.B., J. V. S. Raju, R.A. Weller, S. Bigorre, A. Plueddemann, C.W.
!Fairall, 
!%  S. Miller, L. Mahrt, Dean Vickers, and Hans Hersbach, 2013: On the Exchange
!of momentum
!%  over the open ocean. J. Phys. Oceanogr., 43, 1589–1610. doi:
!http://dx.doi.org/10.1175/JPO-D-12-0173.1 
!
!% Code history:
!% 
!% 1. 12/14/05 - created based on scalar version coare26sn.m with input
!%    on vectorization from C. Moffat.  
!% 2. 12/21/05 - sign error in psiu_26 corrected, and code added to use variable
!%    values from the first pass through the iteration loop for the stable case
!%    with very thin M-O length relative to zu (zetu>50) (as is done in the 
!%    scalar coare26sn and COARE3 codes).
!% 3. 7/26/11 - S = dT was corrected to read S = ut.
!% 4. 7/28/11 - modification to roughness length parameterizations based 
!%    on the CLIMODE, MBL, Gasex and CBLAST experiments are incorporated
!% 5. New wave parameterization added 9/20/2017  based on fits to wave model
!% 6. tested and updated in 9/2020 to give consistent readme info and units,
!%    and so that no external functions are required. They are all included at
!%    end of this program now. Changed names for a few things... including skin
!%    dter -> dT_skin; dt -> dT; dqer -> dq_skin; tkt -> dz_skin
!%    and others to avoid ambiguity:
!%    Rnl -> lw_net; Rns -> sw_net; Rl -> lw_dn; Rs -> sw_dn;
!%    SST -> Tskin
!% 7. 3/10/21 - Adapted for UFS or GFS fortran code 
!%-----------------------------------------------------------------------

    ! skip sanity check
    ! if(P .lt. 900.) P=101300.

!% input variable u is assumed to be wind speed corrected for surface current
!% (magnitude of difference between wind and surface current vectors). to 
!% follow orginal Fairall code, set surface current speed us=0. if us data 
!% are available, construct u prior to using this code.

    us = 0. * u

!% convert rh to specific humidity after accounting for salt effect on freezing
!% point of water
    !%freezing point of seawater
    Tf = -0.0575*Ss + 1.71052E-3*Ss**1.5 - 2.154996E-4*Ss*Ss
    !% surface water specific humidity (g/kg)
    Qs = qsat26sea(ts,P,Ss,Tf)/1000.
    !% P at tq measurement height
!bloss(Test without height correction)    P_tq = P - (0.125*zt)
    P_tq = P 
    call qsat26air(t,P_tq,rh,Q,Pv) !% specific humidity of air (g/kg).  
    !% Assumes rh relative to ice T<0
    !% Pv is the partial pressure due to wate vapor in mb
    Q = Q/1000. !% change back to g/g

    ! JWWang, 2021/2/26
    if(ts < Tf) then
      ice = 1
      jcool = 0.
    else
      ice = 0
      jcool = 1.
    endif

    zos1 = 5.E-4

    !%***********  air constants **********************************************
    Le   = (2.501-.00237*ts)*1.e6
    cpv  = cpa*(1.+0.84*Q)
    rhoa = P_tq*100./(Rgas*(t+T2K)*(1.+0.61*Q))
    !% Pv is the partial pressure due to wate vapor in mb
    visa = 1.326e-5*(1.+6.542e-3*t+8.301e-6*t*t-4.84e-9*t*t*t)

    !%***********  cool skin constants  ***************************************
    !%%% includes salinity dependent thermal expansion coeff for water
    tsw = ts
    if(ts < Tf) tsw=Tf
    Al35 = 2.1e-5*(tsw+3.2)**0.79
    Al0 = (2.2*(tsw-1.)**0.82-5.)*1.e-5
    Al = Al0 + (Al35-Al0)*Ss/35.
    !%%%%%%%%%%%%%%%%%%%
    be = bets*Ss !% be is beta*Salinity
    !%%%%  see "Computing the seater expansion coefficients directly from the
    !%%%%  1980 equation of state".  J. Lillibridge, J.Atmos.Oceanic.Tech, 1980.
    bigc = 16.*grav*cpw*(rhow*visw)**3/(tcw*tcw*rhoa*rhoa)
    wetc = eps*Le*Qs/(Rgas*(ts+T2K)**2)

    !%***********  net solar and IR radiation fluxes ***************************
    !% net solar aka sw aka shortwave
    sw_net = sw_onema*sw_dn !% albedo correction, positive heating ocean

    !% lw_up = eps*sigma*Tskin^4 + (1-eps)*lw_dn
    !% lw_net = lw_up + lw_dn
    !% lw_net = eps(sigma*Tskin^4 - lw-dn)  
    !%     as below for lw_dn>0 heating ocean and lw_net>0 cooling ocean
    !%     Tskin is in Kelvin: ts-0.3*jcool+T2K... approximates a cool-skin of
    !%     amount 0.3 K as initial guess. This gets updated as this code
    !%     progresses.

    !% net longwave aka IR aka infrared
    !% initial value here is positive for cooling ocean in the calculations
    !% below. However it is returned at end of program as -lw_net in final
    !% output so that it is positive heating ocean like the other input
    !% radiation values.

    ! JWWang, 2021/3/23
    !bloss, 2021-07-30: Use physcons.F90 value for Stefan-Boltzman constant
    lw_dn_abs = lw_emiss*lw_dn
    lw_up = lw_emiss*(sigma_r*(ts-0.3*jcool+T2K)**4) + (1.-lw_emiss)*lw_dn
    lw_net = lw_up - lw_dn

!%****************  begin bulk loop ********************************************

!%***********  first guess ************************************************
    du = u - us
    dT = ts - t - .0098*zt
    dq = Qs - Q
    ta = t + T2K
    gust = 0.5
    dT_skin  = 0.3
    ut    = sqrt(du*du+gust*gust)
    u10   = ut*log(10./1e-4)/log(zu/1.e-4)
    usr   = 0.035*u10
    zo10  = 0.011*usr*usr/grav + 0.11*visa/usr
    Cd10  = (von/log(10./zo10))**2
    Ch10  = 0.00115
    Ct10  = Ch10/sqrt(Cd10)
    zot10 = 10./exp(von/Ct10)
    Cd    = (von/log(zu/zo10))**2
    Ct    = von/log(zt/zot10)
    CC    = von*Ct/Cd
    Ribcu = -zu/zi/.004/Beta**3
    Ribu  = -grav*zu/ta*((dT-dT_skin*jcool)+.61*ta*dq)/(ut*ut)
!    zetu  = CC*Ribu*(1.+27./9.*Ribu/CC)
    zetu  = CC*Ribu*(1.+3.*Ribu/CC)
    if(zetu > 50.) k50=1 !% stable with very thin M-O length relative to zu
    if(Ribu < 0.) then
      zetu = CC*Ribu/(1.+Ribu/Ribcu) ! clear k;
    endif
    L10 = zu/zetu

    usr = ut*von/(log(zu/zo10)-psiu_40(zu/L10))
    tsr = -(dT-dT_skin*jcool)*von*fdg/(log(zt/zot10)-psit_26(zt/L10))
    qsr = -(dq-wetc*dT_skin*jcool)*von*fdg/(log(zq/zot10)-psit_26(zq/L10))
    tvsr=tsr+0.61*ta*qsr; !bloss(2021-07-29)
    dz_skin = 0.001

    !%**********************************************************
    !%  The following gives the new formulation for the
    !%  Charnock variable
    !%**********************************************************
    !%%%%%%%%%%%%%   COARE 3.5 wind speed dependent charnock
    umax = 19.
    a1 = 0.0017
    a2 = -0.0050
    charnC = a1*min(u10,umax) + a2

!%%%%%%%%%   if wave age is given but not wave height, use parameterized
!%%%%%%%%%   wave height
    if(cp > 0. .and. sigH < 0.) then
      sigH = (0.02*(cp/u10)**1.1-0.0025)*u10*u10
      sigH = max(sigH,.25)
    elseif(cp <= 0. .and. sigH < 0.) then
      sigH = 0.25
    endif
   
    Ad = 0.2  !%Sea-state/wave-age dependent coefficients from wave model
    !%Ad=0.73./sqrt(u10);
    Bd = 2.2

    charn = charnC
    if(cp > 0.) then
      zoS2 = sigH*Ad*(usr/cp)**Bd
      charnS = zoS2*grav/usr/usr
      charn = charnS
    endif

    nits = 10 !% number of iterations

!%**************  bulk loop **************************************************

    dT_skin_pre = 0.
    do i = 1, nits
!bloss(2021-07-29)      zet = von*grav*zu/ta*(tsr+.61*ta*qsr)/(usr*usr)
      zet = von*grav*zu*(tvsr/ta)/(usr*usr)
      L = zu/zet
      zo = charn*usr*usr/grav + 0.11*visa/usr !% surface roughness
      if(ice == 1) zo=zos1
      rr = zo*usr/visa

      !% This thermal roughness length Stanton number is close to COARE 3.0
      !  value
      zoq = min(1.6e-4,5.8e-5/rr**.72)  

      !% Dalton number is close to COARE 3.0 value
      !% why is it the same as
      zot = zoq                               
      cdhf = von/(log(zu/zo)-psiu_26(zu/L))
      cqhf = von*fdg/(log(zq/zoq)-psit_26(zq/L))
      cthf = von*fdg/(log(zt/zot)-psit_26(zt/L))
      usr_no_gust = du*cdhf ! ustar without gustiness -- avoids division by wind speed which could be zero.
      usr = ut*cdhf ! ustar with gustiness
      qsr = -(dq-wetc*dT_skin*jcool)*cqhf
      tsr = -(dT-dT_skin*jcool)*cthf
      tvsr=tsr+0.61*ta*qsr; !bloss(2021-07-29)
      Bf = -grav/ta*usr*tvsr
      gust = 0.2
      if(Bf > 0.) then
        gust = Beta*(Bf*zi)**.333 ! clear k;
      endif
      ut = sqrt(du*du+gust*gust)
!bloss: Define usr_no_gust above instead of using gf
!bloss      gf = ut/du
      hsb = -rhoa*cpa*usr*tsr
      hlb = -rhoa*Le*usr*qsr
      ! JWWang, 2021/3/3
      ! From Elizabeth Thompson (2021/3/3):
      !   Technically we actually advise leaving out hrain from qout for the
      ! cool skin model part of the code you are asking about. We don't
      ! actually know how rain affects the cool skin, throughout the 0-1 m
      ! depth of water the T gradient totally changes, and there is a cool 
      ! skin on top of that but we don't have enough measurements in rain to
      ! understand it yet. It'd be good to analyze soon with some new research
      ! data, but let's not put it in GFS yet. Technically, the "bulk" water
      ! temperature used as input to COARE (typically collected at 1 m or 5
      ! cm) should contain the rain-cooled effect already, particularly if it
      ! is from as shallow as 5 cm. So this is another reason for us, who
      ! typically input a 5 cm T from our floating sea snake, to not include
      ! hrain in qout.
      !   We do include hrain in calculating hnet at the end outside of the
      !code. 
      qout = lw_net + hsb + hlb
 
      !% The absorption function below is from a Soloviev paper, appears as 
      !% Eq 17 Fairall et al. 1996 and updated/tested by Wick et al. 2005. The
      !% coefficient was changed from 1.37 to 0.065 ~ about halved.
      !% Most of the time this adjustment makes no difference. But then there
      !% are times when the wind is weak, insolation is high, and it matters a
      !% lot. Using the original 1.37 coefficient resulted in many unwarranted
      !% warm-skins that didn't seem realistic. See Wick et al. 2005 for details.
      !% That's the last time the cool-skin routine was updated. The absorption
      !% is not from Paulson & Simpson because that was derived in a lab.
      !% It absorbed to much and produced too many warm layers. It likely had 
      !% too much near-IR (longerwavelength solar) absorption in ocean
      !% which probably doesn't make it to the ocean, was probably absorbed
      !% somewhere in the atmosphere first. The below expression could 
      !% likely use 2 exponentials if you had a shallow mixed layer... 
      !% but we find better results with 3 exponentials. That's the best so 
      !% far we've found that covers the possible depths. 
      dels = sw_net * &
        (0.065+11.*dz_skin-6.6e-5/dz_skin*(1.-exp(-dz_skin/8.0e-4)))  
      qcol = qout - dels
      alq = Al*qcol + be*hlb*cpw/Le
      xlamx = 6.0
      dz_skin = min(0.01,xlamx*visw/(sqrt(rhoa/rhow)*usr))
      if(alq > 0.) then
        xlamx = 6./(1.+(bigc*alq/usr**4)**0.75)**0.333
        dz_skin = xlamx*visw/(sqrt(rhoa/rhow)*usr) ! clear k;
      endif
      dT_skin = qcol*dz_skin/tcw
      dq_skin = wetc*dT_skin
!!$993      format('n/lw_dn/lw_net/sw_dn/sw_net/hlb/alq = ',i4,7e12.5)
!!$      write(*,993) i,lw_dn,lw_net,sw_dn,sw_net,hlb,alq
!!$994      format('n/qout/qcol/dels/alq/dz_skin/dT_skin/dq_skin = ',i4,7e12.5)
!!$      write(*,994) i,qout,qcol,dels,alq,dz_skin,dT_skin,dq_skin
      ! JWWang, 2021/3/23
      !bloss, 2021-07-30: Use physcons.F90 value for Stefan-Boltzman constant
      lw_dn_abs = lw_emiss*lw_dn
      lw_up = lw_emiss*(sigma_r*(ts-dT_skin*jcool+T2K)**4) + &
        (1.-lw_emiss)*lw_dn
      lw_net = lw_up - lw_dn


      if(i == 1 .and. k50 == 1) then
        !% save first iteration solution for case of zetu>50;
        usr50 = usr
        tsr50 = tsr
        qsr50 = qsr
        dT_skin50 = dT_skin
        tkt50 = dz_skin
      endif
                                                             !  dT_skin
      ! compute u10N using ustar without gustiness enhancement
      u10N = usr_no_gust/von*log(10./zo)
      charnC = a1*u10N + a2
      if(u10N > umax) then 
        charnC = a1*umax + a2
      endif

      charn = charnC
      if(cp > 0.) then
        zoS2 = sigH*Ad*(usr/cp)**Bd !%-0.11*visa./usr;
        charnS = zoS2*grav/usr/usr
        charn = charnS
      endif

      if(abs(dT_skin-dT_skin_pre) .lt. dT_skin_min) then
        exit
      else
        dT_skin_pre = dT_skin
      endif
    enddo

    !% insert first iteration solution for case with zetu>50
    if(k50 == 1) then
      usr = usr50
      tsr = tsr50
      qsr = qsr50
      dT_skin = dT_skin50
      dz_skin = tkt50
    endif

    !%****************  compute fluxes  ****************************************
    tau = rhoa*usr*usr_no_gust    !% wind stress
    hsb = -rhoa*cpa*usr*tsr       !% sensible heat flux
    hlb = -rhoa*Le*usr*qsr        !% latent heat flux

    dT_skinx = dT_skin*jcool

    if(PRESENT(output_vec)) then
      ii = 0
      ii = ii + 1; output_vec(ii) = usr
      ii = ii + 1; output_vec(ii) = tau
      ii = ii + 1; output_vec(ii) = hsb
      ii = ii + 1; output_vec(ii) = hlb
      ii = ii + 1; output_vec(ii) = -999999. !bloss hlwebb

      ii = ii + 1; output_vec(ii) = tsr
      ii = ii + 1; output_vec(ii) = qsr
      ii = ii + 1; output_vec(ii) = zot
      ii = ii + 1; output_vec(ii) = zoq
      ii = ii + 1; output_vec(ii) = tau/rhoa/ut/max(0.1,du)

      ii = ii + 1; output_vec(ii) = -usr*tsr/ut/(dt-dT_skin*jcool)
      ii = ii + 1; output_vec(ii) = -usr*qsr/ut/(dq-dq_skin*jcool)
      ii = ii + 1; output_vec(ii) = L
      ii = ii + 1; output_vec(ii) = zet
      ii = ii + 1; output_vec(ii) = dT_skin  !dter

      ii = ii + 1; output_vec(ii) = dq_skin !dqer
      ii = ii + 1; output_vec(ii) = dz_skin !tkt
      ii = ii + 1; output_vec(ii) = -999999. !RF
      ii = ii + 1; output_vec(ii) = 1000.*(von/log(10./zo))*(von/log(10./zo)) !Cdn_10
      ii = ii + 1; output_vec(ii) = 1000.*fdg*(von/log(10./zo))*(von/log(10./zot)) !Chn_10

      ii = ii + 1; output_vec(ii) = 1000.*fdg*(von/log(10./zo))*(von/log(10./zoq)) !Cen_10
    end if

  end subroutine coare36flux_coolskin


  subroutine coare36flux( &
    ! input arguments
    u,zu,t,zt,rh,zq,P,ts,sw_dn,lw_dn,grav,zi,rain,Ss,cp,sigH, &
    ! output arguments
    tau,hsb,hlb,Le,sw_net,lw_dn_abs,lw_up, &
    output_vec)

! -------------------------------------------------------------------------- !
!   This Fortran code is based on the COARE 3.6 Matlab code (Fairall et al,
! 2003; Elizabeth Thompson, personal communication, Feb & Mar 2021). The
! latter is a modified version based on the CLIMODE, MBL, CBLAST, ATOMIC
! experiments (Edson et al. 2012; Elizabeth Thompson et al., personal
! ocmmunication, Feb & Mar 2021)
!   The COARE 3.6 algorithm in the public-released UFS (uncoupled; basically
! GFSv16beta) for calculating latent heat flux (hlb), sensible heat flux
! (hsb), wind stress (tau), cool-skin temperature correction (dT_skinx),
! cool-skin thickness (dz_skin), and latent heat of vaporation (Le; as a
! by-product).
!   The cool-skin schemne of the COARE 3.6 is implemented here because the
! cool-skin scheme in the GFSv16beta Near-Surface Sea Temperature (NSST)
! scheme appears to be older.
!   The input arguments Ss (sea surface salinity), cp (dominant wave speed),
! and sigH (signiciant wave height) are missing in the uncoupled UFS, but
! they are kept here for future full coupling work. They are currently
! assigned constant values inside this subroutine.
! The output dT_skinx
!   Many variables in the COARE 3.6 Matlabe code are discarded here, because
! the UFS or GFS does not need them and it is good practice to save time &
! space in a GCM run. The dicarded variables, however, may be needed for
! comparison with observations. If so, they need to be added in later.
! -------------------------------------------------------------------------- !

! ********************************************************************
!   subroutine coare36flux is a simplified version of subroutine
! coare36flux_coolskin. subroutine coare36flux is not inteneded to perform
! cool-skin correction. Instead, only energy fluxes are wanted as output.
! ********************************************************************

    implicit none

    ! input arguments
    real(kind=kind_phys), intent(in) :: u, zu, t, zt, rh, zq, P, ts,    &
      sw_dn, lw_dn, grav, zi, rain, Ss, cp
    ! input/output arguments
    real(kind_phys), intent(inout) :: sigH
    ! output arguments
    real(kind=kind_phys), intent(out) :: tau, hsb, hlb, Le, &
      sw_net, lw_dn_abs, lw_up

    !bloss: add output_vec for full output tests against the python code
    real(kind=kind_phys), dimension(21), intent(out), optional :: output_vec

    ! local variables
    real ::  Tf, us, Qs, P_tq, Q, Pv, zos1, cpv, rhoa, &
      visa, tsw, Al35, Al0, Al, rhow, visw, &
      lw_net, du, dT, dq, ta, gust, ut, u10, usr, zo10, Cd10, &
      Ch10, Ct10, zot10, Cd, Ct, CC, Ribcu, Ribu, zetu, L10, usr_no_gust, &
      tsr, qsr, charnC, charn, charnS, dq_skin, umax, &
      a1, a2, Ad, Bd, zoS2, nits, zet, L, zo, rr, zoq, zot, cdhf, &
      cqhf, cthf, tvsr, tssr, Bf, qout, dwat, dtmp, dqs_dt, alfac, &
      dels, qcol, alq, xlamx, u10N
    integer :: ice, i 

    integer :: ii !bloss: output_vec index

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Input:  
!
!     u = water-relative mean wind speed magnitude (m/s) at height zu (m)
!         ... if not available, use true wind speed to compute fluxes in
!         earth-coordinates only
!         which will be ignoring the stress contribution by ocean
!         current, which affects all other fluxes too.
!     t = air temperature (degC) at height zt (m)
!     q = specific humidity at height zq (m)
!     P = surface sea level air pressure (mb) (typical value = 1013)
!    ts = seawater temperature (degC)
! sw_dn = downward (positive) shortwave radiation (W/m^2) (typical value = 150) 
! lw_dn = downward (positive) longwave radiation (W/m^2) (typical value = 370)
!  grav = gravitationl force, varying with latitudes
!    zi = PBL height (m) (default or typical value = 600m)
!  rain = rain rate (mm/hr)
!    Ss = sea surface salinity (PSU)
!    cp = phase speed of dominant waves (m/s)  
!  sigH = significant wave height (m)
!  zu, zt, zq heights of the u, t, rh (m)
!
! The user controls the output at the end of the code.

! Output:  
!
!    tau  = wind stress that includes gustiness (N/m^2)
!    hsb  = sensible heat flux (W/m^2) ... positive for Tair < Tskin
!    hlb  = latent heat flux (W/m^2) ... positive for qair < qs
!     Le  = latent heat of vaporation (J/K)

! Others:
!
!% radiation signs: positive warms the ocean
!% sensible, rain, and latent flux signs: positive cools the ocean
!% NOTE: signs change throughout the program for ease of calculations. The
!% units noted here are for the final outputs.
!
!%    usr = friction velocity that includes gustiness (m/s)
!%    tau = wind stress that includes gustiness (N/m^2)
!%    hsb = sensible heat flux (W/m^2) ... positive for Tair < Tskin
!%    hlb = latent heat flux (W/m^2) ... positive for qair < qs
!%    hbb = atmospheric buoyany flux (W/m^2)... positive when hlb and hsb heat
!the atmosphere
!%   hsbb = atmospheric buoyancy flux from sonic ... as above, computed with
!sonic anemometer T
!% hlwebb = webb factor to be added to hl covariance and ID latent heat fluxes
!%    tsr = temperature scaling parameter (K)
!%    qsr = specific humidity scaling parameter (g/kg)
!%     zo = momentum roughness length (m) 
!%    zot = thermal roughness length (m) 
!%    zoq = moisture roughness length (m)
!%     Cd = wind stress transfer (drag) coefficient at height zu   
!%     Ch = sensible heat transfer coefficient (Stanton number) at height zu   
!%     Ce = latent heat transfer coefficient (Dalton number) at height zu
!%      L = Monin-Obukhov length scale (m) 
!%    zet = Monin-Obukhov stability parameter zu/L (dimensionless)
!%dq_skin = cool-skin humidity depression (g/kg)
!ocean
!%   rhoa = density of air at input parameter height zt, typically same as zq
!(kg/m3)
!%   gust = gustiness velocity (m/s)

!% Notes: 1) u is the surface-relative wind speed, i.e., the magnitude of the
!%           difference between the wind (at zu) and ocean surface current 
!%           vectors.
!%        3) Assign a default value to P, lw_dn, sw_dn, zi if unknown.
!
!% Reference:
!%
!%  Fairall, C.W., E.F. Bradley, J.E. Hare, A.A. Grachev, and J.B. Edson (2003),
!%  Bulk parameterization of air sea fluxes: updates and verification for the 
!%  COARE algorithm, J. Climate, 16, 571-590.
!%
!%  Edson, J.B., J. V. S. Raju, R.A. Weller, S. Bigorre, A. Plueddemann, C.W.
!Fairall, 
!%  S. Miller, L. Mahrt, Dean Vickers, and Hans Hersbach, 2013: On the Exchange
!of momentum
!%  over the open ocean. J. Phys. Oceanogr., 43, 1589–1610. doi:
!http://dx.doi.org/10.1175/JPO-D-12-0173.1 
!
!% Code history:
!% 
!% 1. 12/14/05 - created based on scalar version coare26sn.m with input
!%    on vectorization from C. Moffat.  
!% 2. 12/21/05 - sign error in psiu_26 corrected, and code added to use variable
!%    values from the first pass through the iteration loop for the stable case
!%    with very thin M-O length relative to zu (zetu>50) (as is done in the 
!%    scalar coare26sn and COARE3 codes).
!% 3. 7/26/11 - S = dT was corrected to read S = ut.
!% 4. 7/28/11 - modification to roughness length parameterizations based 
!%    on the CLIMODE, MBL, Gasex and CBLAST experiments are incorporated
!% 5. New wave parameterization added 9/20/2017  based on fits to wave model
!% 6. tested and updated in 9/2020 to give consistent readme info and units,
!%    and so that no external functions are required. They are all included at
!%    end of this program now.
!%    SST -> Tskin
!% 7. 3/10/21 - Adapted for UFS or GFS fortran code 
!% 8. 3/13/21 - Adapted from 3/10/21 subroutine coare36flux_coolskin to 
!%    only calculate energy fluxes
!%-----------------------------------------------------------------------

!% input variable u is assumed to be wind speed corrected for surface current
!% (magnitude of difference between wind and surface current vectors). to 
!% follow orginal Fairall code, set surface current speed us=0. if us data 
!% are available, construct u prior to using this code.
    us = 0. * u

!% convert rh to specific humidity after accounting for salt effect on freezing
!% point of water
    !%freezing point of seawater
    Tf = -0.0575*Ss + 1.71052E-3*Ss**1.5 - 2.154996E-4*Ss*Ss
    !% surface water specific humidity (g/kg)
    Qs = qsat26sea(ts,P,Ss,Tf)/1000.
    !% P at tq measurement height
    P_tq = P - (0.125*zt)
    call qsat26air(t,P_tq,rh,Q,Pv) !% specific humidity of air (g/kg).  
    !% Assumes rh relative to ice T<0
    !% Pv is the partial pressure due to wate vapor in mb
    Q = Q/1000. !% change back to g/g

    ! JWWang, 2021/2/26
    if(ts < Tf) then
      ice = 1
    else
      ice = 0
    endif

    zos1 = 5.E-4

    !%***********  air constants **********************************************
    Le   = (2.501-.00237*ts)*1.e6
    cpv  = cpa*(1.+0.84*Q)
    rhoa = P_tq*100./(Rgas*(t+T2K)*(1.+0.61*Q))
    !% Pv is the partial pressure due to wate vapor in mb
    visa = 1.326e-5*(1.+6.542e-3*t+8.301e-6*t*t-4.84e-9*t*t*t)

    !%***********  cool skin constants  ***************************************
    !%%% includes salinity dependent thermal expansion coeff for water
    tsw = ts
    if(ts < Tf) tsw=Tf
    Al35 = 2.1e-5*(tsw+3.2)**0.79
    Al0 = (2.2*(tsw-1.)**0.82-5.)*1.e-5
    Al = Al0 + (Al35-Al0)*Ss/35.

    !%***********  net solar and IR radiation fluxes ***************************
    !% net solar aka sw aka shortwave
    sw_net = sw_onema*sw_dn !% albedo correction, positive heating ocean
    !% lw_up = eps*sigma*Tskin^4 + (1-eps)*lw_dn
    !% lw_net = lw_up + lw_dn
    !% lw_net = eps(sigma*Tskin^4 - lw-dn)  
    !%     as below for lw_dn>0 heating ocean and lw_net>0 cooling ocean
    !%     Tskin is in Kelvin: ts-0.3*jcool+T2K... approximates a cool-skin of
    !%     amount 0.3 K as initial guess. This gets updated as this code
    !%     progresses.

    !% net longwave aka IR aka infrared
    !% initial value here is positive for cooling ocean in the calculations
    !% below. However it is returned at end of program as -lw_net in final
    !% output so that it is positive heating ocean like the other input
    !% radiation values.

    ! JWWang, 2021/3/23
    !bloss, 2021-07-30: Use physcons.F90 value for Stefan-Boltzman constant
    lw_dn_abs = lw_emiss*lw_dn
    lw_up = lw_emiss*(sigma_r*(ts+T2K)**4) + &
      (1.-lw_emiss)*lw_dn
    lw_net = lw_up - lw_dn

!%***********  first guess ************************************************
    du = u - us
    dT = ts - t - .0098*zt
    dq = Qs - Q
    ta = t + T2K
    gust = 0.5
    ut    = sqrt(du*du+gust*gust)
    u10   = ut*log(10./1e-4)/log(zu/1.e-4)
    usr   = 0.035*u10
    zo10  = 0.011*usr*usr/grav + 0.11*visa/usr
    Cd10  = (von/log(10./zo10))**2
    Ch10  = 0.00115
    Ct10  = Ch10/sqrt(Cd10)
    zot10 = 10./exp(von/Ct10)
    Cd    = (von/log(zu/zo10))**2
    Ct    = von/log(zt/zot10)
    CC    = von*Ct/Cd
    Ribcu = -zu/zi/.004/Beta**3
    Ribu  = -grav*zu/ta*(dT+.61*ta*dq)/(ut*ut)
!    zetu  = CC*Ribu*(1.+27./9.*Ribu/CC)
    zetu  = CC*Ribu*(1.+3.*Ribu/CC)
    if(Ribu < 0.) then
      zetu = CC*Ribu/(1.+Ribu/Ribcu) ! clear k;
    endif
    L10 = zu/zetu
    usr = ut*von/(log(zu/zo10)-psiu_40(zu/L10))
    tsr = -dT*von*fdg/(log(zt/zot10)-psit_26(zt/L10))
    qsr = -dq*von*fdg/(log(zq/zot10)-psit_26(zq/L10))

    !%**********************************************************
    !%  The following gives the new formulation for the
    !%  Charnock variable
    !%**********************************************************
    !%%%%%%%%%%%%%   COARE 3.5 wind speed dependent charnock
    umax = 19.
    a1 = 0.0017
    a2 = -0.0050
    charnC = a1*min(u10,umax) + a2

!%%%%%%%%%   if wave age is given but not wave height, use parameterized
!%%%%%%%%%   wave height
    if(cp > 0. .and. sigH < 0.) then
      sigH = (0.02*(cp/u10)**1.1-0.0025)*u10*u10
      sigH = max(sigH,.25)
    elseif(cp <= 0. .and. sigH < 0.) then
      sigH = 0.25
    endif
    
    Ad = 0.2  !%Sea-state/wave-age dependent coefficients from wave model
    !%Ad=0.73./sqrt(u10);
    Bd = 2.2

    charn = charnC
    if(cp > 0.) then
      zoS2 = sigH*Ad*(usr/cp)**Bd
      charnS = zoS2*grav/usr/usr
      charn = charnS
    endif

    zet = von*grav*zu/ta*(tsr+.61*ta*qsr)/(usr*usr)
    L = zu/zet
    zo = charn*usr*usr/grav + 0.11*visa/usr !% surface roughness
    if(ice == 1) zo=zos1
    rr = zo*usr/visa

    !% This thermal roughness length Stanton number is close to COARE 3.0
    !  value
    zoq = min(1.6e-4,5.8e-5/rr**.72)  

    !% Dalton number is close to COARE 3.0 value
    !% why is it the same as
    zot = zoq                               
    cdhf = von/(log(zu/zo)-psiu_26(zu/L))
    cqhf = von*fdg/(log(zq/zoq)-psit_26(zq/L))
    cthf = von*fdg/(log(zt/zot)-psit_26(zt/L))
    !bloss: Eliminate gf to avoid division by wind speed (which could be zero).
    usr_no_gust = du*cdhf ! ustar without gustiness
    usr = ut*cdhf
    qsr = -dq*cqhf
    tsr = -dT*cthf
    tvsr=tsr+0.61*ta*qsr; !bloss(2021-07-29)
    Bf = -grav/ta*usr*tvsr
    gust = 0.2
    if(Bf > 0.) then
      gust = Beta*(Bf*zi)**.333 ! clear k;
    endif
    ut = sqrt(du*du+gust*gust)
!bloss: Define usr_no_gust above instead of using gf
!bloss    gf = ut/du

    !%****************  compute fluxes  ****************************************
    tau = rhoa*usr*usr_no_gust    !% wind stress
    hsb = -rhoa*cpa*usr*tsr       !% sensible heat flux
    hlb = -rhoa*Le*usr*qsr        !% latent heat flux

    if(PRESENT(output_vec)) then
      !bloss: Output a vector of values that can be compared with output
      !  from the matlab or python version of COARE for verification.
      ii = 0
      ii = ii + 1; output_vec(ii) = usr
      ii = ii + 1; output_vec(ii) = tau
      ii = ii + 1; output_vec(ii) = hsb
      ii = ii + 1; output_vec(ii) = hlb
      ii = ii + 1; output_vec(ii) = -999999. !hlwebb

      ii = ii + 1; output_vec(ii) = tsr
      ii = ii + 1; output_vec(ii) = qsr
      ii = ii + 1; output_vec(ii) = zot
      ii = ii + 1; output_vec(ii) = zoq
      ii = ii + 1; output_vec(ii) = Cd

      ii = ii + 1; output_vec(ii) = -999999. !Ch
      ii = ii + 1; output_vec(ii) = -999999. !Ce
      ii = ii + 1; output_vec(ii) = L
      ii = ii + 1; output_vec(ii) = zet
      ii = ii + 1; output_vec(ii) = -999999.  !dter

      ii = ii + 1; output_vec(ii) = -999999. !dqer
      ii = ii + 1; output_vec(ii) = -999999. !tkt
      ii = ii + 1; output_vec(ii) = -999999. !RF
      ii = ii + 1; output_vec(ii) = Cd10 !Cdn_10
      ii = ii + 1; output_vec(ii) = Ch10 !Chn_10

      ii = ii + 1; output_vec(ii) = -999999. !Cen_10
    end if

  end subroutine coare36flux

!%------------------------------------------------------------------------------
  real function psit_26(zet)
!% computes temperature structure function
 
    real(kind_phys) :: zet
    real(kind_phys) :: dzet, psi, x, psik, psic, f

    dzet = min(50.,0.35*zet) !% stable
    psit_26 = -((1.+0.6667*zet)**1.5+0.6667*(zet-14.28)*exp(-dzet)+8.525)
    if(zet < 0.) then !% unstable
      x = (1.-15.*zet)**0.5
      psik = 2.*log((1.+x)/2.)
      x = (1.-34.15*zet)**0.3333
      psic = 1.5*log((1.+x+x*x)/3.) - sqrt(3.)*atan((1.+2.*x)/sqrt(3.)) + &
        4.*atan(1.)/sqrt(3.)
      f = zet*zet/(1.+zet*zet)
      psit_26 = (1.-f)*psik + f*psic
    endif
  end function psit_26

  !%----------------------------------------------------------------------------
  real function psiu_26(zet)
    !% computes velocity structure function
 
    real(kind_phys) :: zet
    real(kind_phys) :: dzet, a, b, c, d, psi, x, psik, psic, f

    dzet = min(50.,0.35*zet) !% stable
    a = 0.7
    b = 3./4.
    c = 5.
    d = 0.35
    psiu_26 = -(a*zet+b*(zet-c/d)*exp(-dzet)+b*c/d)

    if(zet < 0.) then !% unstable
      x = (1.-15.*zet)**0.25
      psik = 2.*log((1.+x)/2.)+log((1.+x*x)/2.)-2.*atan(x)+2.*atan(1.)
      x = (1.-10.15*zet)**0.3333
      psic = 1.5*log((1.+x+x*x)/3.) - sqrt(3.)*atan((1.+2.*x)/sqrt(3.)) + &
        4.*atan(1.)/sqrt(3.)
      f = zet*zet/(1.+zet*zet)
      psiu_26 = (1.-f)*psik + f*psic
    endif
  end function psiu_26

  !%----------------------------------------------------------------------------
  real function psiu_40(zet)
  !% computes velocity structure function

    real(kind_phys) :: zet
    real(kind_phys) :: dzet, a, b, c, d, psi, x, psik, psic, f

    dzet = min(50.,0.35*zet) !% stable
    a = 1.
    b = 3./4.
    c = 5.
    d = 0.35
    psiu_40 = -(a*zet+b*(zet-c/d)*exp(-dzet)+b*c/d)

    if(zet < 0.) then !% unstable
      x = (1.-18.*zet)**0.25
      psik = 2.*log((1.+x)/2.)+log((1.+x*x)/2.) - 2.*atan(x) + 2.*atan(1.)
      x = (1.-10.*zet)**0.3333
      psic = 1.5*log((1.+x+x*x)/3.) - sqrt(3.)*atan((1.+2.*x)/sqrt(3.)) + &
        4.*atan(1.)/sqrt(3.)
      f = zet*zet/(1.+zet*zet)
      psiu_40 = (1.-f)*psik+f*psic
    endif
  end function psiu_40

  !%----------------------------------------------------------------------------
  real function bucksat(T,P,Tf)
  !% computes saturation vapor pressure [mb]
  !% given T [degC] and P [mb] Tf is freezing pt 

    real(kind_phys) :: T, P, Tf

    bucksat = 6.1121*exp(17.502*T/(T+240.97))*(1.0007+3.46e-6*P)
    if(T < Tf) then
      bucksat = (1.0003+4.18e-6*P)*6.1115*exp(22.452*T/(T+272.55)) !% vapor
                                                               !  pressure ice
    endif
  end function bucksat

!%------------------------------------------------------------------------------
  real function qsat26sea(T,P,Ss,Tf)
  !% computes surface saturation specific humidity [g/kg]
  !% given T [degC] and P [mb]

    real(kind_phys) :: T, P, Ss, Tf
    real(kind_phys) :: ex, fs, es

    ex = bucksat(T,P,Tf)
    fs = 1.-0.02*Ss/35. !% reduction sea surface vapor pressure by salinity
    es = fs*ex 
    qsat26sea = 622.*es/(P-0.378*es)
  end function qsat26sea
!%------------------------------------------------------------------------------
!function [q,em]=qsat26air(T,P,rh)
  subroutine qsat26air(T,P,rh,q,em)
  !% computes saturation specific humidity [g/kg]
  !% given T [degC] and P [mb]

    real(kind_phys), intent(in) :: T, P, rh
    real(kind_phys) :: Tf, es
    real(kind_phys), intent(out) :: q, em

    Tf = 0. !%assumes relative humidity for pure water
    es = bucksat(T,P,Tf)
    em = 0.01*rh*es  !% in mb, partial pressure of water vapor
    q = 622.*em/(P-0.378*em)
  end subroutine qsat26air

  !%----------------------------------------------------------------------------
  real function RHcalc(T,P,Q,Tf)
  !% computes relative humidity given T,P, & Q
    real(kind_phys) :: T, P, Q, Tf
    real(kind_phys) :: es, em

    es = 6.1121*exp(17.502*T/(T+240.97))*(1.0007+3.46e-6*P)
    if(T < Tf) then !%ice case
      es = 6.1115*exp(22.452*T/(T+272.55))*(1.0003+4.18e-6*P)
      em = Q*P/(0.378*Q+eps)
      RHcalc = 100.*em/es
    endif
  end function RHcalc

end module coare36_module
