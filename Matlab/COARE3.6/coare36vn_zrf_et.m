function A=coare36vn_zrf_et(u,zu,t,zt,rh,zq,P,ts,sw_dn,lw_dn,lat,lon,jd,zi,rain,Ss,cp,sigH,zrf_u,zrf_t,zrf_q)
%

%**************************************************************************
% VERSION INFO:

% Vectorized version of COARE 3.6 code (Fairall et al, 2003) with 
% modification based on the CLIMODE, MBL and CBLAST experiments 
% (Edson et al., 2012). The cool skin and surface wave options are included.
% A separate warm layer function can be used to call this function.
%
% This 3.6 version include parameterizations using wave height and wave
% slope using cp and sigH.  If these are set to NaN, then the wind
% speed dependent formulation is used.  The parameterizations are based
% on fits to the Banner-Norison wave model and the Fairall-Edson flux
% database.  This version also allows salinity as a input.  
% Open ocean example Ss=35; Great Lakes Ss=0;
 
%**************************************************************************
% COOL SKIN:
%
% An important component of this code is whether the inputed ts 
% represents the ocean skin temperature or a subsurface temperature.  
% How this variable is treated is determined by the jcool parameter:
%   set jcool=1 if ts is subsurface or bulk ocean temperature (default);
%   set jcool=0 if ts is skin temperature or to not run cool skin model.
% The code updates the cool-skin temperature depression dT_skin and 
% thickness dz_skin during iteration loop for consistency. The number of 
% iterations set to nits = 6.

    jcoolx=1;
    
%**************************************************************************
%%% INPUTS:  

% Notes on input default values, missing values, vectors vs. single values:
%   - the code assumes u,t,rh,ts,P,sw_dn,lw_dn,rain,Ss,cp,sigH are vectors; 
%   - sensor heights (zu,zt,zl) latitude lat, longitude lon, julian date jd, 
%       and PBL height zi may be constants;
%   - air pressure P and radiation sw_dn, lw_dn may be vectors or constants. 
%   - input NaNs as vectors or single values to indicate no data. 
%   - assign a default value to P, lw_dn, sw_dn, lat, zi if unknown, single
%       values of these inputs are okay.
% Notes about signs and units: 
%   - radiation signs: positive warms the ocean
%   - signs and units change throughout the program for ease of calculations.
%   - the signs and units noted here are for the inputs.

%
%     u = water-relative wind speed magnitude (m/s) at height zu (m)
%             i.e. mean wind speed accounting for the ocean current vector. 
%             i.e. the magnitude of the difference between the wind vector 
%             (at height zu) and ocean surface current vector.
%             If not available, use true wind speed to compute fluxes in 
%             earth-coordinates only which will be ignoring the stress 
%             contribution from the ocean current to all fluxes
%     t = air temperature (degC) at height zt (m)
%    rh = relative humidity (%) at height zq (m)
%     P = sea level air pressure (mb) 
%    ts = seawater temperature (degC), see jcool below for cool skin 
%             calculation, and separate warm layer code for specifying 
%             sensor depth and whether warm layer is computed
% sw_dn = downward (positive) shortwave radiation (W/m^2)
% lw_dn = downward (positive) longwave radiation (W/m^2) 
%   lat = latitude defined positive to north
%   lon = longitude defined positive to east, if using other version,
%             adjust the eorw string input to albedo_vector function
%    jd = year day or julian day, where day Jan 1 00:00 UTC = 0
%    zi = PBL height (m) (default or typical value = 600m)
%  rain = rain rate (mm/hr)
%    Ss = sea surface salinity (PSU)
%    cp = phase speed of dominant waves (m/s) computed from peak period 
%  sigH = significant wave height (m)
%  zu, zt, zq heights of the observations (m)
%  zrf_u, zrf_t, zrf_q  reference height for profile.  Use this to compare observations at different heights  

%**************************************************************************
%%%% OUTPUTS: the user controls the output array A at the end of the code.

% Note about signs and units: 
%   - radiation signs: positive warms the ocean
%   - sensible, rain, and latent flux signs: positive cools the ocean
%   - signs and units change throughout the program for ease of calculations.
%   - the signs and units noted here are for the final outputs.

%    usr = friction velocity that includes gustiness (m/s), u*
%    tau = wind stress that includes gustiness (N/m^2)
%    hsb = sensible heat flux (W/m^2) ... positive for Tair < Tskin
%    hlb = latent heat flux (W/m^2) ... positive for qair < qs
%    hbb = atmospheric buoyany flux (W/m^2)... positive when hlb and hsb heat the atmosphere
%   hsbb = atmospheric buoyancy flux from sonic ... as above, computed with sonic anemometer T
% hlwebb = webb factor to be added to hl covariance and ID latent heat fluxes
%    tsr = temperature scaling parameter (K), t*
%    qsr = specific humidity scaling parameter (kg/kg), q*
%     zo = momentum roughness length (m) 
%    zot = thermal roughness length (m) 
%    zoq = moisture roughness length (m)
%     Cd = wind stress transfer (drag) coefficient at height zu (unitless)
%     Ch = sensible heat transfer coefficient (Stanton number) at height zu (unitless)
%     Ce = latent heat transfer coefficient (Dalton number) at height zu (unitless)
%      L = Monin-Obukhov length scale (m) 
%    zeta = Monin-Obukhov stability parameter zu/L (dimensionless)
%dT_skin = cool-skin temperature depression (degC), pos value means skin is cooler than subskin
%dq_skin = cool-skin humidity depression (g/kg)
%dz_skin = cool-skin thickness (m)
%    Urf = wind speed at reference height (user can select height at input)
%    Trf = air temperature at reference height
%    Qrf = air specific humidity at reference height
%   RHrf = air relative humidity at reference height
%   UrfN = neutral value of wind speed at reference height
%   TrfN = neutral value of air temp at reference height
%  qarfN = neutral value of air specific humidity at reference height
% lw_net = Net IR radiation computed by COARE (W/m2)... positive heating ocean
% sw_net = Net solar radiation computed by COARE (W/m2)... positive heating ocean
%     Le = latent heat of vaporization (J/K)
%   rhoa = density of air at input parameter height zt, typically same as zq (kg/m3)
%     UN = neutral value of wind speed at zu (m/s)
%    U10 = wind speed adjusted to 10 m (m/s)
%   UN10 = neutral value of wind speed at 10m (m/s)
% Cdn_10 = neutral value of drag coefficient at 10m (unitless)
% Chn_10 = neutral value of Stanton number at 10m (unitless)
% Cen_10 = neutral value of Dalton number at 10m (unitless)
%  hrain = rain heat flux (W/m^2)... positive cooling ocean
%     Qs = sea surface specific humidity, i.e. assuming saturation (g/kg)
%   Evap = evaporation rate (mm/h)
%    T10 = air temperature at 10m (deg C)
%    Q10 = air specific humidity at 10m (g/kg)
%   RH10 = air relative humidity at 10m (%)
%    P10 = air pressure at 10m (mb)
% rhoa10 = air density at 10m (kg/m3)
%   gust = gustiness velocity (m/s)
%wc_frac = whitecap fraction (ratio)
%   Edis = energy dissipated by wave breaking (W/m^2)

%**************************************************************************
%%%% ADDITONAL CALCULATIONS: 

%   using COARE output, one can easily calculate the following using the
%   sign conventions and names herein:

%     %%% Skin sea surface temperature or interface temperature; neglect
%     %%% dT_warm_to_skin if warm layer code is not used as the driver
%     %%%% program for this program
%           Tskin = ts + dT_warm_to_skin - dT_skin;
%     
%     %%% Upwelling radiative fluxes: positive heating the ocean
%           lw_up = lw_net - lw_dn;
%           sw_up = sw_net - sw_dn;
%  
%     %%% Net heat flux: positive heating ocean
%     %%% note that hs, hl, hrain are defined when positive cooling
%     %%% ocean by COARE, so their signs are flipped here: 
%           hnet = sw_net + lw_net - hs - hl - hrain;

%**************************************************************************
%%%% REFERENCES:
%
%  Fairall, C. W., E. F. Bradley, J. S. Godfrey, G. A. Wick, J. B. Edson, 
%  and G. S. Young, 1996a: Cool-skin and warm-layer effects on sea surface 
%  temperature. J. Geophys. Res., 101, 1295?1308.
%
%  Fairall, C. W., E. F. Bradley, D. P. Rogers, J. B. Edson, and G. S. Young,
%  1996b: Bulk parameterization of air-sea fluxes for Tropical Ocean- Global
%  Atmosphere Coupled- Ocean Atmosphere Response Experiment. J. Geophys. Res.,
%  101, 3747?3764.
%
%  Fairall, C. W., A. B. White, J. B. Edson, and J. E. Hare, 1997: Integrated
%  shipboard measurements of the marine boundary layer. Journal of Atmospheric
%  and Oceanic Technology, 14, 338?359
%
%  Fairall, C.W., E.F. Bradley, J.E. Hare, A.A. Grachev, and J.B. Edson (2003),
%  Bulk parameterization of air sea fluxes: updates and verification for the 
%  COARE algorithm, J. Climate, 16, 571-590.
%
%  Edson, J.B., J. V. S. Raju, R.A. Weller, S. Bigorre, A. Plueddemann, C.W.
%  Fairall, S. Miller, L. Mahrt, Dean Vickers, and Hans Hersbach, 2013: On 
%  the Exchange of momentum over the open ocean. J. Phys. Oceanogr., 43, 
%  1589–1610. doi: http://dx.doi.org/10.1175/JPO-D-12-0173.1 

%**************************************************************************
% CODE HISTORY:
% 
% 1. 12/14/05 - created based on scalar version coare26sn.m with input
%    on vectorization from C. Moffat.  
% 2. 12/21/05 - sign error in psiu_26 corrected, and code added to use variable
%    values from the first pass through the iteration loop for the stable case
%    with very thin M-O length relative to zu (zetau>50) (as is done in the 
%    scalar coare26sn and COARE3 codes).
% 3. 7/26/11 - S = dT was corrected to read S = ut.
% 4. 7/28/11 - modification to roughness length parameterizations based 
%    on the CLIMODE, MBL, Gasex and CBLAST experiments are incorporated
% 5. 9/20/2017 - New wave parameterization added based on fits to wave model
% 6. 9/2020 - tested and updated to give consistent readme info and units,
%    and so that no external functions are required. They are all included at
%    end of this program now. Changed names for a few things... including skin
%    dter -> dT_skin; dt -> dT; dqer -> dq_skin; tkt -> dz_skin
%    and others to avoid ambiguity:
%    Rnl -> lw_net; Rns -> sw_net; Rl -> lw_dn; Rs -> sw_dn;
%    SST -> Tskin; Also corrected heights at which q and P are
%    computed to be more accurate, changed units of qstar to kg/kg, removed
%    extra 1000 on neutral 10 m transfer coefficients;
% 7. 10/2021 - implemented zenith angle dependent sw_up and sw_net;
%    changed buoyancy flux calculation to follow Stull 
%    textbook version of tv* and tv_sonic*; reformatted preamble of program for
%    consistent formatting and to reduce redundancy; resolved issues of
%    nomenclature around T adjusted to heights vs theta potential
%    temperature when computing dT for the purpose of sensible heat flux.
%-----------------------------------------------------------------------

%***********  prep input data *********************************************

% convert input to column vectors
u=u(:);t=t(:);rh=rh(:);P=P(:);ts=ts(:);
sw_dn=sw_dn(:);lw_dn=lw_dn(:);lat=lat(:);zi=zi(:);
zu=zu(:);zt=zt(:);zq=zq(:);
zrf_u=zrf_u(:);zrf_t=zrf_t(:);zrf_q=zrf_q(:);
rain=rain(:);  
Ss=Ss(:);cp=cp(:);sigH=sigH(:);
N=length(u);
jcool=jcoolx*ones(N,1);jcool=jcool(:);

% Option to set local variables to default values if input is NaN... can do
% single value or fill each individual. Warning... this will fill arrays
% with the dummy values and produce results where no input data are valid
% ii=find(isnan(P)); P(ii)=1013;    % pressure
% ii=find(isnan(sw_dn)); sw_dn(ii)=200;   % incident shortwave radiation
% ii=find(isnan(lat)); lat(ii)=45;  % latitude
% ii=find(isnan(lw_dn)); lw_dn(ii)=400-1.6*abs(lat(ii)); % incident longwave radiation
% ii=find(isnan(zi)); zi(ii)=600;   % PBL height
% ii=find(isnan(Ss)); Ss(ii)=35;    % Salinity

% find missing input data
iip=find(isnan(P)); 
iirs=find(isnan(sw_dn)); 
iilat=find(isnan(lat)); 
iirl=find(isnan(lw_dn)); 
iizi=find(isnan(zi)); 
iiSs=find(isnan(Ss)); 

% Input variable u is assumed to be wind speed corrected for surface current
% (magnitude of difference between wind and surface current vectors). To 
% follow orginal Fairall code, set surface current speed us=0. If us surface
% current data are available, construct u prior to using this code and
% input us = 0*u here; 
us = 0*u;

% convert rh to specific humidity after accounting for salt effect on freezing
% point of water
Tf=-0.0575*Ss+1.71052E-3*Ss.^1.5-2.154996E-4*Ss.*Ss; %freezing point of seawater
Qs = qsat26sea(ts,P,Ss,Tf)./1000; % surface water specific humidity (g/kg)
P_tq = P - (0.125*zt); % P at tq measurement height (mb)
[Q,Pv]  = qsat26air(t,P_tq,rh); % specific humidity of air (g/kg).  
    % Assumes rh relative to ice T<0
    % Pv is the partial pressure due to wate vapor in mb
Q=Q./1000;  % change Q to g/g

ice=zeros(size(u));
iice=find(ts<Tf);ice(iice)=1;jcool(iice)=0;
zos=5E-4;

%***********  set constants ***********************************************
zref=10;
Beta = 1.2;
von  = 0.4;
fdg  = 1.00; % Turbulent Prandtl number
T2K  = 273.16;
grav = grv(lat);

%***********  air constants ***********************************************
Rgas = 287.1;
Le   = (2.501-.00237*ts)*1e6;
cpa  = 1004.67;
cpv  = cpa*(1+0.84*Q);
rhoa = P_tq*100./(Rgas*(t+T2K).*(1+0.61*Q));
% Pv is the partial pressure due to wate vapor in mb
rhodry = (P_tq-Pv)*100./(Rgas*(t+T2K)); % dry air density. 
visa = 1.326e-5*(1+6.542e-3.*t+8.301e-6*t.^2-4.84e-9*t.^3);
lapse=grav/cpa; % dry adiabatic lapse rate, K/km

%***********  cool skin constants  ***************************************
%%% includes salinity dependent thermal expansion coeff for water
tsw=ts;ii=find(ts<Tf);tsw(ii)=Tf(ii);
Al35   = 2.1e-5*(tsw+3.2).^0.79;
Al0   =(2.2*real((tsw-1).^0.82)-5)*1e-5;
Al=Al0+(Al35-Al0).*Ss/35;
%%%%%%%%%%%%%%%%%%%
bets = 7.5e-4; % salintity expansion coeff; assumes beta is independent of T
be   = bets*Ss; % be is beta*salinity
%%%%  see "Computing the seater expansion coefficients directly from the
%%%%  1980 equation of state".  J. Lillibridge, J.Atmos.Oceanic.Tech, 1980.
cpw  = 4000;
rhow = 1022;
visw = 1e-6;
tcw  = 0.6;
bigc = 16*grav*cpw*(rhow*visw)^3./(tcw.^2*rhoa.^2);
wetc = 0.622*Le.*Qs./(Rgas*(ts+T2K).^2);

%***********  net solar and IR radiation fluxes ***************************
%%% net solar flux, aka sw, aka shortwave

% *** for time-varying, i.e. zenith angle varying albedo using Payne 1972:
% insert 'E' for input to albedo function if longitude is defined positive
% to E (normal), in this case lon sign will be flipped for the calculation.
% Otherwise specify 'W' and the sign will not be changed in the function. 
% Check: albedo should usually peak at sunrise not at sunset, though it may
% vary based on sw_dn.
[alb,T_sw,solarmax_sw,psi_sw] = albedo_vector(sw_dn,jd,lon,lat,'E');
sw_net = (1-alb).*sw_dn; % varying albedo correction, positive heating ocean

% *** for constant albedo:
% sw_net = 0.945.*sw_dn; % constant albedo correction, positive heating ocean

%%% net longwave aka IR aka infrared
% initial value here is positive for cooling ocean in the calculations
% below. However it is returned at end of program as -lw_net in final output so 
% that it is positive heating ocean like the other input radiation values.
lw_net = 0.97*(5.67e-8*(ts-0.3*jcool+T2K).^4-lw_dn); 

%***********  begin bulk loop *********************************************

%***********  first guess *************************************************

% wind speed minus current speed
du = u-us;
% air sea temperature difference for the purpose of sensible heat flux
dT = ts-t-lapse.*zt;

% air-sea T diff must account for lapse rate between surface and instrument height
% t is air temperature in C, ts is surface water temperature in C. dT is
% an approximation that is equivalent to  dtheta where theta is the
% potential temperature, and the pressure at sea level and instrument level
% are used. They are equivalent (max difference = 0.0022 K). This way 
% elimniates the need to involve the pressures at different heights. 
% Using or assuming dry adiabatic lapse rate between the two heights
% doesn't matter because if real pressures are used the result is the 
% unchanged. The dT need not include conversion to K either. Here's an example: 
% grav = grv(lat);
% lapse=grav/cpa;
% P_at_tq_height=(psealevel - (0.125*zt)); % P at tq measurement height (mb)
% note psealevel is adjusted using same expression from pa height
% Ta is originally in C and C2K = 273.15 to convert from C to K
% theta = (b10.Ta+C2K).*(1000./P_tq).^(Rgas/cpa); 
% TadjK = (b10.Ta+C2K) + lapse*zt;
% Tadj = b10.Ta + lapse*zt;
% theta_sfc = (b10.Tskin+C2K).*(1000./b10.psealevel).^(Rgas/cpa); 
% TadjK_sfc = b10.Tskin+C2K;
% Tadj_sfc = b10.Tskin;
% 
%%% the adj versions are only 0.0022 K smaller than theta versions)
% dtheta = theta_sfc - theta;
% dTadjK = TadjK_sfc - TadjK;
% dTadj = Tadj_sfc - Tadj; % so dT = Tskin - (Ta + lapse*zt) = Tskin - Ta - lapse*zt


% put things into different units and expressions for more calculations,
% including first guesses that get redone later
dq = Qs-Q;
ta = t+T2K;
tv = ta.*(1+0.61*Q);
gust = 0.5;
dT_skin  = 0.3;
ut    = sqrt(du.^2+gust.^2);
u10   = ut.*log(10/1e-4)./log(zu/1e-4);
usr   = 0.035*u10;
zo10  = 0.011*usr.^2./grav + 0.11*visa./usr;
Cd10  = (von./log(10./zo10)).^2;
Ch10  = 0.00115;
Ct10  = Ch10./sqrt(Cd10);
zot10 = 10./exp(von./Ct10);
Cd    = (von./log(zu./zo10)).^2;
Ct    = von./log(zt./zot10);
CC    = von*Ct./Cd;
Ribcu = -zu./zi./.004/Beta^3;
Ribu  = -grav.*zu./ta.*((dT-dT_skin.*jcool)+.61*ta.*dq)./ut.^2;
zetau = CC.*Ribu.*(1+27/9*Ribu./CC);
k50=find(zetau>50); % stable with very thin M-O length relative to zu
k=find(Ribu<0); 
if length(Ribcu)==1
    zetau(k)=CC(k).*Ribu(k)./(1+Ribu(k)./Ribcu); clear k;
else
    zetau(k)=CC(k).*Ribu(k)./(1+Ribu(k)./Ribcu(k)); clear k;
end
L10 = zu./zetau;
gf=ut./du;
usr = ut.*von./(log(zu./zo10)-psiu_40(zu./L10));
tsr = -(dT-dT_skin.*jcool).*von*fdg./(log(zt./zot10)-psit_26(zt./L10));
qsr = -(dq-wetc.*dT_skin.*jcool)*von*fdg./(log(zq./zot10)-psit_26(zq./L10));
dz_skin = 0.001*ones(N,1);

%**********************************************************
%  The following gives the new formulation for the
%  Charnock variable
%**********************************************************
%%%%%%%%%%%%%   COARE 3.5 wind speed dependent charnock
charnC = 0.011*ones(N,1);
umax=19;
a1=0.0017;
a2=-0.0050;
charnC=a1*u10+a2;
k=find(u10>umax);
charnC(k)=a1*umax+a2;


%%%%%%%%%   if wave age is given but not wave height, use parameterized
%%%%%%%%%   wave height based on wind speed
    hsig=(0.02*(cp./u10).^1.1-0.0025).*u10.^2;
    hsig=max(hsig,.25);
    ii=find(~isnan(cp) & isnan(sigH));
    sigH(ii)=hsig(ii);
    
Ad=0.2;  %Sea-state/wave-age dependent coefficients from wave model
%Ad=0.73./sqrt(u10);
Bd=2.2;
zoS=sigH.*Ad.*(usr./cp).^Bd;
charnS=zoS.*grav./usr./usr;

nits=10; % number of iterations
charn=charnC;
ii=find(~isnan(cp));charn(ii)=charnS(ii);
%**************  bulk loop ************************************************

for i=1:nits
    zeta=von.*grav.*zu./ta.*(tsr +.61*ta.*qsr)./(usr.^2);
    L=zu./zeta;
    zo=charn.*usr.^2./grav+0.11*visa./usr; % surface roughness
    zo(iice)=zos;
    rr=zo.*usr./visa;
    
    % This thermal roughness length Stanton number is close to COARE 3.0 value
    zoq=min(1.6e-4,5.8e-5./rr.^.72);  
     
    % Andreas 1987 for snow/ice
    ik=find(rr(iice)<=.135);
   		rt(iice(ik))=rr(iice(ik))*exp(1.250);
     	rq(iice(ik))=rr(iice(ik))*exp(1.610);
     ik=find(rr(iice)>.135 & rr(iice)<=2.5);
        rt(iice(ik))=rr(iice(ik)).*exp(0.149-.55*log(rr(iice(ik))));
     	rq(iice(ik))=rr(iice(ik)).*exp(0.351-0.628*log(rr(iice(ik))));
     ik=find(rr(iice)>2.5 & rr(iice)<=1000);
     	rt(iice(ik))=rr(iice(ik)).*exp(0.317-0.565*log(rr(iice(ik)))-0.183*log(rr(iice(ik))).*log(rr(iice(ik))));
      	rq(iice(ik))=rr(iice(ik)).*exp(0.396-0.512*log(rr(iice(ik)))-0.180*log(rr(iice(ik))).*log(rr(iice(ik))));
 
    % Dalton number is close to COARE 3.0 value
    zot=zoq;                               
    cdhf=von./(log(zu./zo)-psiu_26(zu./L));
    cqhf=von.*fdg./(log(zq./zoq)-psit_26(zq./L));
    cthf=von.*fdg./(log(zt./zot)-psit_26(zt./L));
    usr=ut.*cdhf;
    qsr=-(dq-wetc.*dT_skin.*jcool).*cqhf;
    tsr=-(dT-dT_skin.*jcool).*cthf;
    
    % original COARE version buoyancy flux
    tvsr1=tsr+0.61*ta.*qsr;
    tssr1=tsr+0.51*ta.*qsr;
    
    % new COARE version buoyancy flux from Stull (1988) page 146
    % tsr here uses dT with the lapse rate adjustment (see code above). The
    % Q and ta values should be at measurement height, not adjusted heights
    tvsr=tsr.*(1+0.61.*Q) + 0.61*ta.*qsr;
    tssr=tsr.*(1+0.51.*Q) + 0.51*ta.*qsr;
    
    Bf=-grav./ta.*usr.*tvsr;
    gust=0.2*ones(N,1);
    k=find(Bf>0); 
    %%% gustiness in this way is from the original code. Notes: 
    % we measured the actual gustiness by measuring the variance of the
    % wind speed and empirically derived the the scaling. It's empirical
    % but it seems appropriate... the longer the time average then the larger
    % the gustiness factor should be, to account for the gustiness averaged
    % or smoothed out by the averaging. wstar is the convective velocity.
    % gustiness is beta times wstar. gustiness is different between mean of
    % velocity and square of the mean of the velocity vector components.
    % The actual wind (mean + fluctuations) is still the most relavent 
    % for the flux. The models do u v w, and then compute vector avg to get
    % speed, so we've done the same thing. coare alg input is the magnitude
    % of the mean vector wind relative to water. 

    if length(zi)==1
        gust(k)=Beta*(Bf(k).*zi).^.333; clear k;
    else
        gust(k)=Beta*(Bf(k).*zi(k)).^.333; clear k;
    end
    ut=sqrt(du.^2+gust.^2);
    gf=ut./du;
    hsb=-rhoa*cpa.*usr.*tsr;
    hlb=-rhoa.*Le.*usr.*qsr;
    qout=lw_net+hsb+hlb; 
    %%% rain heat flux is not included in qout because we don't fully
    % understand the evolution or gradient of the cool skin layer in the
    % presence of rain, and the sea snake subsurface measurement input
    % value will capture some of the rain-cooled water already. TBD.
    
    %%% solar absorption: 
    % The absorption function below is from a Soloviev paper, appears as 
    % Eq 17 Fairall et al. 1996 and updated/tested by Wick et al. 2005. The
    % coefficient was changed from 1.37 to 0.065 ~ about halved.
    % Most of the time this adjustment makes no difference. But then there
    % are times when the wind is weak, insolation is high, and it matters a
    % lot. Using the original 1.37 coefficient resulted in many unwarranted
    % warm-skins that didn't seem realistic. See Wick et al. 2005 for details.
    % That's the last time the cool-skin routine was updated. The
    % absorption is not from Paulson & Simpson because that was derived in a lab.
    % It absorbed too much and produced too many warm layers. It likely  
    % approximated too much near-IR (longerwavelength solar) absorption
    % which probably doesn't make it to the ocean since it was probably absorbed
    % somewhere in the atmosphere first. The below expression could 
    % likely use 2 exponentials if you had a shallow mixed layer... 
    % but we find better results with 3 exponentials. That's the best so 
    % far we've found that covers the possible depths. 
    dels=sw_net.*(0.065+11*dz_skin-6.6e-5./dz_skin.*(1-exp(-dz_skin/8.0e-4)));  
    qcol=qout-dels;
    % only needs stress, water temp, sum of sensible, latent, ir, solar,
    % and latent individually. 
    alq=Al.*qcol+be.*hlb.*cpw./Le; % buoyancy flux... to make it sink. Positive buoyancy fluxes generate tke. cooling at the interface creates turb in the ocean through the buoyancy term. and destroys it on the air side. 
    xlamx=6.0*ones(N,1); % contains buoyancy effect. There are two sources.Net cooling of interface by latent, sensible, ir, 
%     the other is the salinity part caused by latent heat flux (evap) leaving behind salt. 
    dz_skin=min(0.01, xlamx.*visw./(sqrt(rhoa./rhow).*usr));
    k=find(alq>0);
    xlamx(k)=6./(1+(bigc(k).*alq(k)./usr(k).^4).^0.75).^0.333;
    dz_skin(k)=xlamx(k).*visw./(sqrt(rhoa(k)./rhow).*usr(k)); clear k;
    dT_skin=qcol.*dz_skin./tcw;
    dq_skin=wetc.*dT_skin;
    lw_net=0.97*(5.67e-8*(ts-dT_skin.*jcool+T2K).^4-lw_dn); % used to update dT_skin
    if i==1 % save first iteration solution for case of zetau>50;
        usr50=usr(k50);tsr50=tsr(k50);qsr50=qsr(k50);L50=L(k50);
        zeta50=zeta(k50);dT_skin50=dT_skin(k50);dq_skin50=dq_skin(k50);tkt50=dz_skin(k50);
    end
    u10N = usr./von./gf.*log(10./zo);
    charnC=a1*u10N+a2;
    k=u10N>umax;
    charnC(k)=a1*umax+a2;
    charn=charnC;
    zoS=sigH.*Ad.*(usr./cp).^Bd;%-0.11*visa./usr;
    charnS=zoS.*grav./usr./usr;
    ii=find(~isnan(cp));charn(ii)=charnS(ii);
end

% insert first iteration solution for case with zetau>50
usr(k50)=usr50;tsr(k50)=tsr50;qsr(k50)=qsr50;L(k50)=L50;
zeta(k50)=zeta50;dT_skin(k50)=dT_skin50;dq_skin(k50)=dq_skin50;dz_skin(k50)=tkt50;

%****************  compute fluxes  ****************************************
tau=rhoa.*usr.*usr./gf;         % wind stress
hsb=-rhoa.*cpa.*usr.*tsr;       % sensible heat flux
hlb=-rhoa.*Le.*usr.*qsr;        % latent heat flux
hbb=-rhoa.*cpa.*usr.*tvsr;      % atmospheric buoyancy flux new from Jim Edson
hbb1=-rhoa.*cpa.*usr.*tvsr1;    % atmospheric buoyancy flux original
hsbb=-rhoa.*cpa.*usr.*tssr;     % atmospheric buoyancy flux new from Jim Edson
hsbb1=-rhoa.*cpa.*usr.*tssr1;   % atmospheric buoyancy flux original

wbar=1.61*hlb./Le./(1+1.61*Q)./rhoa+hsb./rhoa./cpa./ta;  % mean w
hlwebb=rhoa.*wbar.*Q.*Le;       % webb correction to add to covariance hl
Evap=1000*hlb./Le./1000*3600;   % evap rate mm/hour

%*****  compute transfer coeffs relative to ut @ meas. ht  ****************
Cd= tau./rhoa./ut./max(.1,du);
Ch=-usr.*tsr./ut./(dT-dT_skin.*jcool);
Ce=-usr.*qsr./(dq-dq_skin.*jcool)./ut;

%***%%  compute 10-m neutral coeff relative to ut *************************
Cdn_10=von.^2./log(10./zo).^2;
Chn_10=von.^2.*fdg./log(10./zo)./log(10./zot);
Cen_10=von.^2.*fdg./log(10./zo)./log(10./zoq);

%***%%  compute 10-m neutral coeff relative to ut *************************

% Find the stability functions for computing values at user defined 
% reference heights and 10 m
psi=psiu_26(zu./L);
psi10=psiu_26(10./L);
psirf=psiu_26(zrf_u./L);
psiT=psit_26(zt./L);
psi10T=psit_26(10./L);
psirfT=psit_26(zrf_t./L);
psirfQ=psit_26(zrf_q./L);
gf=ut./du;

%*********************************************************
%  Determine the wind speeds relative to ocean surface at different heights
%  Note that usr is the friction velocity that includes 
%  gustiness usr = sqrt(Cd) S, which is equation (18) in
%  Fairall et al. (1996)
%*********************************************************
S = ut;
U = du;
S10 = S + usr./von.*(log(10./zu)-psi10+psi);
U10 = S10./gf;
% or U10 = U + usr./von./gf.*(log(10/zu)-psi10+psi);
Urf = U + usr./von./gf.*(log(zrf_u./zu)-psirf+psi);
UN = U + psi.*usr/von./gf;
U10N = U10 + psi10.*usr/von./gf; % technically this removes gustiness because the average wind is the U10
UrfN = Urf + psirf.*usr/von./gf;

UN2 = usr/von./gf.*log(zu./zo);
U10N2 = usr./von./gf.*log(10./zo);
UrfN2  = usr./von./gf.*log(zrf_u./zo);

%******** rain heat flux *****************************
dwat=2.11e-5*((t+T2K)./T2K).^1.94; % water vapor diffusivity
dtmp=(1. + 3.309e-3*t - 1.44e-6.*t.*t).*0.02411./(rhoa.*cpa); % heat diffusivity
dqs_dt=Q.*Le./(Rgas.*(t+T2K).^2); % Clausius-Clapeyron
alfac= 1./(1+0.622*(dqs_dt.*Le.*dwat)./(cpa.*dtmp)); % wet bulb factor
hrain= rain.*alfac.*cpw.*((ts-t-dT_skin.*jcool)+(Qs-Q-dq_skin.*jcool).*Le./cpa)./3600; % rain heat flux

Tskin=ts-dT_skin.*jcool; % ocean skin surface temperature, or interface temperature, C

% P is sea level pressure, so use subtraction through hydrostatic equation 
% to get P10 and P at reference height
P10 = P - (0.125*10);
Prf = P - (0.125*zref);

T10 = t + tsr./von.*(log(10./zt)-psi10T+psiT) + lapse.*(zt-10);
Trf = t + tsr./von.*(log(zrf_t./zt)-psirfT+psiT) + lapse.*(zt-zrf_t);
TN = t + psiT.*tsr/von;
T10N = T10 + psi10T.*tsr/von;
TrfN = Trf + psirfT.*tsr/von;

% unused... these are here to make sure you gets the same answer whether 
% you used the thermal calculated roughness lengths or the values at the 
% measurement height. So at this point they are just illustrative and can
% be removed or ignored if you want.
TN2 = Tskin + tsr/von.*log(zt./zot)-lapse.*zt;
T10N2 = Tskin + tsr/von.*log(10./zot)-lapse.*10;
TrfN2 = Tskin + tsr/von.*log(zrf_t./zot)-lapse.*zrf_t;

dq_skin=wetc.*dT_skin.*jcool;
Qs=Qs-dq_skin;
dq_skin = dq_skin*1000;
Qs=Qs*1000;
Q=Q*1000;
Q10 = Q + 1E3.*qsr./von.*(log(10./zq)-psi10T+psiT);
Qrf = Q + 1E3.*qsr./von.*(log(zrf_q./zq)-psirfQ+psiT);
QN = Q + psiT.*1E3.*qsr/von./sqrt(gf);
Q10N = Q10 + psi10T.*1E3.*qsr/von;
QrfN = Qrf + psirfQ.*1E3.*qsr/von;

% unused... these are here to make sure you gets the same answer whether 
% you used the thermal calculated roughness lengths or the values at the 
% measurement height. So at this point they are just illustrative and can
% be removed or ignored if you want.
QN2 = Qs + 1E3.*qsr/von.*log(zq./zoq);
Q10N2 = Qs + 1E3.*qsr/von.*log(10./zoq);
QrfN2 = Qs + 1E3.*qsr/von.*log(zrf_q./zoq);

RHrf=RHcalc(Trf,Prf,Qrf/1000,Tf);
RH10=RHcalc(T10,P10,Q10/1000,Tf);

% recompute rhoa10 with 10-m values of everything else.
rhoa10 = P10*100./(Rgas*(T10+T2K).*(1+0.61*(Q10/1000)));

%%%%%%%%%%%%  Other wave breaking statistics from Banner-Morison wave model
wc_frac=7.3E-4*(U10N-2).^1.43;      % wind only
wc_frac(U10<2.1)=1e-5;

kk = find(isfinite(cp) == 1);   % wind and waves if wave info is available
wc_frac(kk)=1.6e-3*U10N(kk).^1.1./sqrt(cp(kk)./U10N(kk));  

Edis=0.095*rhoa.*U10N.*usr.^2;  %  energy dissipation rate from breaking waves W/m^2
wc_frac(iice)=0;Edis(iice)=0;

%****************  output  ****************************************************
% only return values if jcool = 1; if cool skin model was intended to be run
dT_skinx=dT_skin.*jcool;
dq_skinx=dq_skin.*jcool;

% get rid of filled values where nans are present in input data
bad_input = find(isnan(u) ==1);
gust(bad_input) = nan;
dz_skin(bad_input) = nan;
zot(bad_input) = nan;
zoq(bad_input) = nan;

% flip lw_net sign for standard radiation sign convention: positive heating ocean
lw_net = -lw_net; 
% this sign flip means lw_net, net long wave flux, is equivalent to:
% lw_net = 0.97*(lw_dn_best - 5.67e-8*(Tskin+C2K).^4);

% adjust A output as desired:
A=[usr tau hsb hlb hbb hsbb hlwebb tsr qsr zo  zot zoq Cd Ch Ce  L  zeta dT_skinx dq_skinx dz_skin Urf Trf Qrf RHrf UrfN TrfN QrfN  lw_net sw_net Le rhoa UN U10 U10N Cdn_10 Chn_10 Cen_10 hrain Qs Evap T10 T10N Q10 Q10N  RH10 P10 rhoa10 gust wc_frac Edis];
%   1   2   3   4   5   6    7      8   9  10  11  12  13 14 15  16  17  18       19        20     21  22  23   24   25   26   27     28      29  30  31  32 33   34    35     36   37      38   39  40  41  42   43   44   45    46   47    48   49      50   
end
%------------------------------------------------------------------------------
function psi=psit_26(zeta)
% computes temperature structure function
dzeta=min(50,0.35*zeta); % stable
psi=-((1+0.6667*zeta).^1.5+0.6667*(zeta-14.28).*exp(-dzeta)+8.525);
k=find(zeta<0); % unstable
x=(1-15*zeta(k)).^0.5;
psik=2*log((1+x)./2);
x=(1-34.15*zeta(k)).^0.3333;
psic=1.5*log((1+x+x.^2)./3)-sqrt(3)*atan((1+2*x)./sqrt(3))+4*atan(1)./sqrt(3);
f=zeta(k).^2./(1+zeta(k).^2);
psi(k)=(1-f).*psik+f.*psic;
end
%------------------------------------------------------------------------------
function psi=psiu_26(zeta)
% computes velocity structure function
dzeta=min(50,0.35*zeta); % stable
a=0.7;
b=3/4;
c=5;
d=0.35;
psi=-(a*zeta+b*(zeta-c/d).*exp(-dzeta)+b*c/d);
k=find(zeta<0); % unstable
x=(1-15*zeta(k)).^0.25;
psik=2*log((1+x)/2)+log((1+x.*x)/2)-2*atan(x)+2*atan(1);
x=(1-10.15*zeta(k)).^0.3333;
psic=1.5*log((1+x+x.^2)/3)-sqrt(3)*atan((1+2*x)./sqrt(3))+4*atan(1)./sqrt(3);
f=zeta(k).^2./(1+zeta(k).^2);
psi(k)=(1-f).*psik+f.*psic;
end
%------------------------------------------------------------------------------
function psi=psiu_40(zeta)
% computes velocity structure function
dzeta=min(50,0.35*zeta); % stable
a=1;
b=3/4;
c=5;
d=0.35;
psi=-(a*zeta+b*(zeta-c/d).*exp(-dzeta)+b*c/d);
k=find(zeta<0); % unstable
x=(1-18*zeta(k)).^0.25;
psik=2*log((1+x)/2)+log((1+x.*x)/2)-2*atan(x)+2*atan(1);
x=(1-10*zeta(k)).^0.3333;
psic=1.5*log((1+x+x.^2)/3)-sqrt(3)*atan((1+2*x)./sqrt(3))+4*atan(1)./sqrt(3);
f=zeta(k).^2./(1+zeta(k).^2);
psi(k)=(1-f).*psik+f.*psic;
end
%------------------------------------------------------------------------------
function exx=bucksat(T,P,Tf)
% computes saturation vapor pressure [mb]
% given T [degC] and P [mb] Tf is freezing pt 
exx=6.1121.*exp(17.502.*T./(T+240.97)).*(1.0007+3.46e-6.*P);
ii=find(T<Tf);
exx(ii)=(1.0003+4.18e-6*P(ii)).*6.1115.*exp(22.452.*T(ii)./(T(ii)+272.55));%vapor pressure ice
end
%------------------------------------------------------------------------------
function qs=qsat26sea(T,P,Ss,Tf)
% computes surface saturation specific humidity [g/kg]
% given T [degC] and P [mb]
ex=bucksat(T,P,Tf);
fs=1-0.02*Ss/35;% reduction sea surface vapor pressure by salinity
es=fs.*ex; 
qs=622*es./(P-0.378*es);
end
%------------------------------------------------------------------------------
function [q,em]=qsat26air(T,P,rh)
% computes saturation specific humidity [g/kg]
% given T [degC] and P [mb]
Tf=0;%assumes relative humidity for pure water
es=bucksat(T,P,Tf);
em=0.01*rh.*es; % in mb, partial pressure of water vapor
q=622*em./(P-0.378*em);
end
%------------------------------------------------------------------------------
function g=grv(lat)
% computes g [m/sec^2] given lat in deg
gamma=9.7803267715;
c1=0.0052790414;
c2=0.0000232718;
c3=0.0000001262;
c4=0.0000000007;
phi=lat*pi/180;
x=sin(phi);
g=gamma*(1+c1*x.^2+c2*x.^4+c3*x.^6+c4*x.^8);
end
%------------------------------------------------------------------------------
function RHrf=RHcalc(T,P,Q,Tf)
% computes relative humidity given T,P, & Q
es=6.1121.*exp(17.502.*T./(T+240.97)).*(1.0007+3.46e-6.*P);
ii=find(T<Tf);%ice case
es(ii)=6.1115.*exp(22.452.*T(ii)./(T(ii)+272.55)).*(1.0003+4.18e-6*P(ii));
em=Q.*P./(0.378.*Q+0.622);
RHrf=100*em./es;
end
%------------------------------------------------------------------------------
function [alb,T,solarmax,psi] = albedo_vector(sw_dn,jd,lon,lat,eorw)

%  Computes transmission and albedo from downwelling sw_dn using
%  lat   : latitude in degrees (positive to the north)
%  lon   : longitude in degrees (positive to the west)
%  jd    : yearday
%  sw_dn : downwelling solar radiation measured at surface
%  eorw  : 'E' if longitude is positive to the east, or 'W' if otherwise

% updates:
%   20-10-2021: ET vectorized function

if strcmp(eorw,'E') == 1
%     disp('lon is positive to east so negate for albedo calculation');
    lon = -lon;
elseif strcmp(eorw,'W') == 1
%     disp('lon is already positive to west so go ahead with albedo calculation');
else
    disp('please provide sign information on whether lon is deg E or deg W');
end

alb = nan(length(sw_dn),1); % allocate array of albedo values as vector
lat=lat*pi/180;       %Convert to radians
lon=lon*pi/180;       %Convert to radians
SC=1380;              %Solar constant W/m^2
utc=(jd-fix(jd))*24;  %UTC time (decimal hours)
h=pi*utc/12-lon;
declination=23.45*cos(2*pi*(jd-173)/365.25);       % Solar declination angle
solarzenithnoon=(lat*180/pi-declination);          % Zenith angle at noon
solaraltitudenoon=90-solarzenithnoon;              % Altitude at noon
sd=declination*pi/180;                             % Convert to radians
gamma=1;                                           % ratio of actual to mean earth-sun separation (set to 1 for now)
gamma2=gamma*gamma;
%
sinpsi = sin(lat).*sin(sd)-cos(lat).*cos(sd).*cos(h);  % Local elevation angle
psi = asin(sinpsi).*180/pi;
solarmax=SC.*sinpsi/gamma2;                            % No atmosphere
%solarmax=1380*sinpsi*(0.61+0.20*sinpsi);

T = min(2,sw_dn./solarmax);

Ts = 0:0.05:1;
As = 0:2:90;

%%% for vectorized function
for k = 1:length(sinpsi)
 Tchk = abs(Ts-T(k));
 i = find(Tchk == min(Tchk));

%%% for single value function
% Tchk = abs(Ts-T);
% i=find(Tchk==min(Tchk));


%  Look up table from Payne (1972)  Only adjustment is to T=0.95 Alt=10 value
%
%       0     2     4     6     8     10    12    14    16    18    20   22     24    26    28    30    32    34    36    38    40    42    44    46    48    50    52    54    56    58    60    62    64    66    68    70    72    74    76    78    80    82    84    86    88    90

a = [ 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061; ...
      0.062 0.062 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061; ...
      0.072 0.070 0.068 0.065 0.065 0.063 0.062 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.060 0.061 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060; ...
      0.087 0.083 0.079 0.073 0.070 0.068 0.066 0.065 0.064 0.063 0.062 0.061 0.061 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060; ...
      0.115 0.108 0.098 0.086 0.082 0.077 0.072 0.071 0.067 0.067 0.065 0.063 0.062 0.061 0.061 0.060 0.060 0.060 0.060 0.061 0.061 0.061 0.061 0.060 0.059 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.059 0.059 0.059; ...
      0.163 0.145 0.130 0.110 0.101 0.092 0.084 0.079 0.072 0.072 0.068 0.067 0.064 0.063 0.062 0.061 0.061 0.061 0.060 0.060 0.060 0.060 0.060 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.058; ...
      0.235 0.198 0.174 0.150 0.131 0.114 0.103 0.094 0.083 0.080 0.074 0.074 0.070 0.067 0.065 0.064 0.063 0.062 0.061 0.060 0.060 0.060 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.058 0.058 0.058; ...
      0.318 0.263 0.228 0.192 0.168 0.143 0.127 0.113 0.099 0.092 0.084 0.082 0.076 0.072 0.070 0.067 0.065 0.064 0.062 0.062 0.060 0.060 0.060 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.058 0.058 0.058 0.058 0.058 0.058 0.058 0.058 0.057 0.058 0.058 0.058 0.058 0.057 0.057 0.057; ...
      0.395 0.336 0.290 0.248 0.208 0.176 0.151 0.134 0.117 0.107 0.097 0.091 0.085 0.079 0.075 0.071 0.068 0.067 0.065 0.063 0.062 0.061 0.060 0.060 0.060 0.059 0.059 0.058 0.058 0.058 0.057 0.057 0.057 0.057 0.057 0.057 0.057 0.056 0.056 0.056 0.056 0.056 0.056 0.056 0.056 0.055; ...
      0.472 0.415 0.357 0.306 0.252 0.210 0.176 0.154 0.135 0.125 0.111 0.102 0.094 0.086 0.081 0.076 0.072 0.071 0.068 0.066 0.065 0.063 0.062 0.061 0.060 0.059 0.058 0.057 0.057 0.057 0.056 0.055 0.055 0.055 0.055 0.055 0.055 0.054 0.053 0.054 0.053 0.053 0.054 0.054 0.053 0.053; ...
      0.542 0.487 0.424 0.360 0.295 0.242 0.198 0.173 0.150 0.136 0.121 0.110 0.101 0.093 0.086 0.081 0.076 0.073 0.069 0.067 0.065 0.064 0.062 0.060 0.059 0.058 0.057 0.056 0.055 0.055 0.054 0.053 0.053 0.052 0.052 0.052 0.051 0.051 0.050 0.050 0.050 0.050 0.051 0.050 0.050 0.050; ...
      0.604 0.547 0.498 0.407 0.331 0.272 0.219 0.185 0.160 0.141 0.127 0.116 0.105 0.097 0.089 0.083 0.077 0.074 0.069 0.066 0.063 0.061 0.059 0.057 0.056 0.055 0.054 0.053 0.053 0.052 0.051 0.050 0.050 0.049 0.049 0.049 0.048 0.047 0.047 0.047 0.046 0.046 0.047 0.047 0.046 0.046; ...
      0.655 0.595 0.556 0.444 0.358 0.288 0.236 0.190 0.164 0.145 0.130 0.119 0.107 0.098 0.090 0.084 0.076 0.073 0.068 0.064 0.060 0.058 0.056 0.054 0.053 0.051 0.050 0.049 0.048 0.048 0.047 0.046 0.046 0.045 0.045 0.045 0.044 0.043 0.043 0.043 0.042 0.042 0.043 0.042 0.042 0.042; ... 
      0.693 0.631 0.588 0.469 0.375 0.296 0.245 0.193 0.165 0.145 0.131 0.118 0.106 0.097 0.088 0.081 0.074 0.069 0.065 0.061 0.057 0.055 0.052 0.050 0.049 0.047 0.046 0.046 0.044 0.044 0.043 0.042 0.042 0.041 0.041 0.040 0.040 0.039 0.039 0.039 0.038 0.038 0.038 0.038 0.038 0.038; ... 
      0.719 0.656 0.603 0.480 0.385 0.300 0.250 0.193 0.164 0.145 0.131 0.116 0.103 0.092 0.084 0.076 0.071 0.065 0.061 0.057 0.054 0.051 0.049 0.047 0.045 0.043 0.043 0.042 0.041 0.040 0.039 0.039 0.038 0.038 0.037 0.036 0.036 0.035 0.035 0.034 0.034 0.034 0.034 0.034 0.034 0.034; ... 
      0.732 0.670 0.592 0.474 0.377 0.291 0.246 0.190 0.162 0.144 0.130 0.114 0.100 0.088 0.080 0.072 0.067 0.062 0.058 0.054 0.050 0.047 0.045 0.043 0.041 0.039 0.039 0.038 0.037 0.036 0.036 0.035 0.035 0.034 0.033 0.032 0.032 0.032 0.031 0.031 0.031 0.030 0.030 0.030 0.030 0.030; ... 
      0.730 0.652 0.556 0.444 0.356 0.273 0.235 0.188 0.160 0.143 0.129 0.113 0.097 0.086 0.077 0.069 0.064 0.060 0.055 0.051 0.047 0.044 0.042 0.039 0.037 0.035 0.035 0.035 0.034 0.033 0.033 0.032 0.032 0.032 0.029 0.029 0.029 0.029 0.028 0.028 0.028 0.028 0.027 0.027 0.028 0.028; ... 
      0.681 0.602 0.488 0.386 0.320 0.252 0.222 0.185 0.159 0.142 0.127 0.111 0.096 0.084 0.075 0.067 0.062 0.058 0.054 0.050 0.046 0.042 0.040 0.036 0.035 0.033 0.032 0.032 0.031 0.030 0.030 0.030 0.030 0.029 0.027 0.027 0.027 0.027 0.026 0.026 0.026 0.026 0.026 0.026 0.026 0.026; ... 
      0.581 0.494 0.393 0.333 0.288 0.237 0.211 0.182 0.158 0.141 0.126 0.110 0.095 0.083 0.074 0.066 0.061 0.057 0.053 0.049 0.045 0.041 0.039 0.034 0.033 0.032 0.031 0.030 0.029 0.028 0.028 0.028 0.028 0.027 0.026 0.026 0.026 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.025; ... 
      0.453 0.398 0.342 0.301 0.266 0.226 0.205 0.180 0.157 0.140 0.125 0.109 0.095 0.083 0.074 0.065 0.061 0.057 0.052 0.048 0.044 0.040 0.038 0.033 0.032 0.031 0.030 0.029 0.028 0.027 0.027 0.026 0.026 0.026 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.025; ... 
      0.425 0.370 0.325 0.290 0.255 0.220 0.200 0.178 0.157 0.140 0.122 0.108 0.095 0.083 0.074 0.065 0.061 0.056 0.052 0.048 0.044 0.040 0.038 0.033 0.032 0.031 0.030 0.029 0.028 0.027 0.026 0.026 0.026 0.026 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.025]; 

if psi(k)<0 
   alb(k)=0;
   solarmax(k)=0;
   T(k)=0;
   j=0;
   psi(k)=0;
else
   Achk = abs(As-psi(k));
   j=find(Achk==min(Achk));
   szj = size(j);
   if szj(2) > 0
       alb(k)=a(i,j);
   else
%       disp('no j found, not assigning alb to anything');
   end
end

end % end for k list of alb array
%disp([num2str(jd) '  ' num2str(sw_dn) '  ' num2str(alb) '  ' num2str(T) '  ' num2str(i) '  ' num2str(j)])
end