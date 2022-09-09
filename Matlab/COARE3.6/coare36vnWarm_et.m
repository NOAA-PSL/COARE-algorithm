function B=coare36vnWarm_et(Jd,U,zu,Tair,zt,RH,zq,P,Tsea,SW_dn,LW_dn,Lat,Lon,zi,Rainrate,ts_depth,Ss,cp,sigH,zrf_u,zrf_t,zrf_q)

disp('WarmCoolLayer')

%***********   input data **************
% *** capitals are used for input arrays and lower case are used for single
% values used within the program 
%       Jd = day-of-year or julian day
%	     U = wind speed magnitude (m/s) corrected for currents, i.e. relative to water at height zu
%	    zu = height (m) of wind measurement
%	  Tair = air temp (degC) at height zt
%	    zt = height (m) of air temperature measurement
%	    RH = relative humidity (%) at height zq
%	    zq = height (m) of air humidity measurement
%	     P = air pressure at sea level (mb) 
%     Tsea = surface sea temp (degC) at ts_depth
%	 SW_dn = downward solar flux (w/m^2) defined positive down
%	 LW_dn = downward IR flux (w/m^2) defined positive down
%	   Lat = latitude (deg N=+)
%	   Lon = longitude (deg E=+) % If using other version, see
%               usage of lon in albedo_vector function and adjust 'E' or 
%               'W' string input
%       zi = inversion height (m)
% Rainrate = rain rate (mm/hr)
% ts_depth = depth (m) of water temperature measurement, positive for below
% surface
%       Ss = sea surface salinity (PSU)
%       cp = phase speed of dominant waves (m/s)  
%     sigH = significant wave height (m)
% zu, zt, zq = heights of the observations (m)
% zrf_u, zrf_t, zrf_q = reference height for profile.  
%   Use this to compare observations at different heights  
%
%  
%********** output data  ***************
% Outputs
% Adds onto output from coare36vn_zrf_et
% .... see that function for updated output. It can change. This function adds 4 variables onto it:
% previously named dt_wrm, tk_pwp, dsea, du_wrm... renamed to the following
% for clarity:
%         dT_warm = dT from base of warm layer to skin, i.e. warming across entire warm layer depth (deg C) 
%         dz_warm = warm layer thickness (m)
% dT_warm_to_skin = dT from measurement depth to skin due to warm layer, 
%                       such that Tskin = Tsea + dT_warm_to_skin - dT_skin
%         du_warm = total current accumulation in warm layer (m/s ?...  unsure of units but likely m/s)

%********** history ********************
% updated 09/2020 for consistency with units, readme info, and coare 3.6 main function
%    Changed names for a few things... 
%    dt_wrm -> dT_warm; tk_pwp -> dz_warm; dsea -> dT_warm_to_skin;
%    skin: dter -> dT_skin
%    and fluxes to avoid ambiguity: RF -> hrain
%    Rnl -> lw_net; Rns -> sw_net; Rl -> lw_dn; Rs -> sw_dn;
% 05/2022 - fix a glitch in code when using wave input - LB
%   added     cpi=cp(ibg);   % air density
%             sigHi=sigH(ibg);   % air density
% and edited line 258 to use cpi and sigHi
% Bx=coare36vn_zrf_et(u,zu,t,zt,rh,zq,p,ts,sw_dn,lw_dn,lat,lon,jd,zi,rain,ss,cpi,sigHi,zrf_u,zrf_t,zrf_q);



%********** Set cool skin options ******************
jcool=1;  % 0=no cool skin calc; 1=do cool skin calc
icount=1;

%********** Set wave ******************
%%% ... not sure if this is necessary ... 
% if no wave info is provided, fill array with nan. 
% if length(cp) == 1 && isnan(cp) == 1
%     cp = Jd*nan;
% end
% if length(sigH) == 1 && isnan(sigH) == 1
%     sigH = Jd*nan;
% end

%*********************  housekeep variables  ********
% Call coare35vn to get initial flux values
Bx=coare36vn_zrf_et(U(1),zu(1),Tair(1),zt(1),RH(1),zq(1),P(1),Tsea(1),SW_dn(1),LW_dn(1),Lat(1),Lon(1),Jd(1),zi(1),Rainrate(1),Ss(1),cp(1),sigH(1),zrf_u(1),zrf_t(1),zrf_q(1));
%%% check these indices for you latest version of coare!
tau_old=Bx(2);  %stress
hs_old=Bx(3);   %sensible heat flux
hl_old=Bx(4);   %latent heat flux 
dT_skin_old=Bx(18);%cool skin
hrain_old=Bx(38);  %rain heat flux

qcol_ac=0;      %accumulates heat from integral, J/s
tau_ac=0;       %accumulates stress from integral, N/m2
dT_warm=0;      %total warming (amplitude deg C) in warm layer
du_warm=0;      %total current (amplitude m/s?) in warm layer, m/s most likely
max_pwp=19;     %maximum depth of warm layer (adjustable)
dz_warm=max_pwp;%initial depth set to max value
dT_warm_to_skin=0; %dT from measurement input depth to skin, initially set to 0
q_pwp=0;        %total heat absorped in warm layer
fxp=.5;         %initial value of solar flux absorption

rich=.65;       %critical Richardson number tuned for COARE

jtime=0;
jamset=0;
jump=1;

%*******************  set constants  ****************
T2K=273.16;   %Converts to Kelvin
Rgas=287.1;   %Universal gas constant
cpa=1004.67;  %Specific heat of air at constant pressure
cpw=4000;     %Specific heat of water
rhow=1022;    %Density of water
visw=1e-6;    %Viscosity of water

be=0.026;
tcw=0.6;

%**********************************************************
%******************  setup read data loop  ****************
P_tq = P - (0.125*zt); % P at tq measurement height
[Press,Tseak,Tairk,Qsatsea,Qsat,Qair,Rhoair,Rhodry]=scalarv(P,Tsea,Tair,RH,zt);

nx=length(Jd);        %# of lines of data

% this is an empty array for saving warm layer code output values. Will be
% added to coare output at the end.
warm_output = nan(nx,4);

for ibg = 1:nx 			% major read loop discritized by time
    jd=Jd(ibg);       % yearday where 1 is Jan 1
    p=P(ibg);           % air pressure at sea level
    u=U(ibg);           % wind speed corrected for surface currents... i.e. wind speed adjused or relative to moving water surface
    tsea=Tsea(ibg);     % subsurface sea temp
    t=Tair(ibg);        % air temp
    ss=Ss(ibg);         % salinity at or near the surface
    qs=Qsatsea(ibg);    % sea surface humidity
    q=Qair(ibg);        % air specific humidity 
    rh=RH(ibg);         % air relative humidity
    sw_dn=SW_dn(ibg);      % downward solar flux (positive down)
    lw_dn=LW_dn(ibg);         % doward IR flux (positive down)
    rain=Rainrate(ibg); % rain rate
    grav=grv(Lat(ibg)); % gravity
    lat=Lat(ibg);      % latitude
    lon=Lon(ibg);      % longitude
    rhoa=Rhoair(ibg);   % air density
    cpi=cp(ibg);   % air density
    sigHi=sigH(ibg);   % air density
    
    %*****  variables for warm layer  ***
    
    %%% for constant albedo
    % sw_net=.945*sw_dn;     %Net Solar: positive warming ocean, constant albedo

    %%% for albedo that is time-varying, i.e. zenith angle varying
    % insert 'E' for input to albedo function if longitude is defined positive
    % to E, so that lon can be flipped. The function ideally works with
    % longitude positive to the west. Check: albedo should peak at sunrise not
    % sunset.
    [alb,T_sw,solarmax_sw,psi_sw] = albedo_vector(sw_dn,jd,lon,lat,'E');
    sw_net = (1-alb).*sw_dn; % constant albedo correction, positive heating ocean

    lw_net=.97*(5.67e-8*(tsea-dT_skin_old*jcool+T2K)^4-lw_dn); %Net IR: positive cooling ocean
    cpv=cpa*(1+0.84*q/1000);
    visa=1.326e-5*(1+6.542e-3*t+8.301e-6*t*t-4.84e-9*t*t*t);
    Al=2.1e-5*(tsea+3.2)^0.79;
    ctd1=sqrt(2*rich*cpw/(Al*grav*rhow));       %mess-o-constants 1
    ctd2=sqrt(2*Al*grav/(rich*rhow))/(cpw^1.5); %mess-o-constants 2
    
    %********************************************************
    %****  Compute apply warm layer  correction *************
    %********************************************************
    
    intime=jd-fix(jd);
    loc=(lon+7.5)/15;
    chktime=loc+intime*24;
    if chktime>24
        chktime=chktime-24;
    end
    newtime=(chktime-24*fix(chktime/24))*3600;
    if icount>1  %not first time thru
        if newtime<=21600 || jump==0
            jump=0;
            if newtime < jtime	%re-zero at midnight
                jamset=0;
                fxp=.5;
                dz_warm=max_pwp;
                tau_ac=0;
                qcol_ac=0;
                dT_warm=0;
                du_warm=0;
            else
                %************************************
                %****   set warm layer constants  ***
                %************************************
                dtime=newtime-jtime;                %delta time for integrals
                qr_out=lw_net+hs_old+hl_old+hrain_old; %total cooling at surface
                q_pwp=fxp*sw_net-qr_out;               %tot heat abs in warm layer
                  qqrx(ibg)=hs_old;
                if q_pwp>=50 || jamset==1           %Check for threshold
                    jamset=1;                       %indicates threshold crossed
                    tau_ac=tau_ac+max(.002,tau_old)*dtime;	%momentum integral
                    if qcol_ac+q_pwp*dtime>0	            %check threshold for warm layer existence
                        %******************************************
                        % Compute the absorption profile
                        %******************************************
                        for i=1:5                           %loop 5 times for fxp

                            %%%% The original version since Fairall et al. 1996:
                            fxp=1-(0.28*0.014*(1-exp(-dz_warm/0.014))+0.27*0.357*(1-exp(-dz_warm/0.357))+0.45*12.82*(1-exp(-dz_warm/12.82)))/dz_warm;
                            % the above integrated flux formulation is wrong for the warm layer,
                            % but it has been used in this scheme since 1996 without
                            % making bad predictions.
                            % Simon recognized that fxp should be the
                            % fraction absorbed in the layer of the form
                            % 1-sum(ai*exp(-tk_pwp/gammai)) with sum(ai)=1
                            % One can choose different exponential
                            % absorption parameters, but it has to be of the right form.
                            % Simon idealized coefficients from profiles of
                            % absorption in DYNAMO by C. Ohlmann.
                            % see /Users/sdeszoek/Data/cruises/DYNAMO_2011/solar_absorption/test_absorption_fcns.m
                            % 
                            % Correct form of Soloviev 3-band absorption:
                            % --using original F96 absorption bands:
                            %  fxp=1-(0.32*exp(-tk_pwp/22.0) + 0.34*exp(-tk_pwp/1.2) + 0.34*exp(-tk_pwp/0.014)); % NOT TESTED!!
                            % --using DYNAMO absorption bands (F, invgamma defined above):

                            %%%% NOT TESTED!! Correction of fxp from Simon ***
%                            fxp=1-sum(F.*(exp(-tk_pwp*invgamma)),2);   
                            
                            qjoule=(fxp*sw_net-qr_out)*dtime;
                            if qcol_ac+qjoule>0             %Compute warm-layer depth
                                dz_warm=min(max_pwp,ctd1*tau_ac/sqrt(qcol_ac+qjoule));
                            end
                        end
                    else   %warm layer wiped out
                        fxp=0.75;
                        dz_warm=max_pwp;
                        qjoule=(fxp*sw_net-qr_out)*dtime;
                    end
                    qcol_ac=qcol_ac+qjoule; %heat integral
                    %*******  compute dt_warm  ******
                    if qcol_ac>0
                        dT_warm=ctd2*(qcol_ac)^1.5/tau_ac;
                        du_warm=2*tau_ac/(dz_warm*rhow); % check units... should be m/s?
                    else
                        dT_warm=0;
                        du_warm=0;
                    end
                end%                    end threshold check
            end%                        end midnight reset
            % Compute warm layer dT between input measurement and skin layer
            if dz_warm<ts_depth          
                dT_warm_to_skin=dT_warm;
            else
                dT_warm_to_skin=dT_warm*ts_depth/dz_warm; % assumes a linear profile of T throughout warm layer
            end
            
        end%                            end 6am start first time thru
    end%                                end first time thru check
    jtime=newtime;
    %************* output from routine  *****************************

% Adjust tsea for warm layer above measurement. Even if tsea is at 5 cm for the sea snake,
% there could be warming present between the interface and the subsurface
% temperature. dT_warm_to_skin estimates that warming between the levels.
    ts=tsea+dT_warm_to_skin;
    
% Rerun COARE with the warm-layer corrected tsea temperature. COARE will
% apply a cool skin to this, completing all calculations needed for Tskin and fluxes. 
% Using COARE ouput from this function, Tskin = Tsnake - dT_skin + dT_warm_to_skin
% note: in prior COARE lingo/code: dT_warm_to_skin used to be dsea and dT_skin used to be dter
    Bx=coare36vn_zrf_et(u,zu,t,zt,rh,zq,p,ts,sw_dn,lw_dn,lat,lon,jd,zi,rain,ss,cpi,sigHi,zrf_u,zrf_t,zrf_q);

% save values from this time step to be used in next time step, this is how
% the integral is computed
    tau_old=Bx(2);            % stress
    hs_old=Bx(3);             % sensible heat flux
    hl_old=Bx(4);             % latent heat flux
    dT_skin_old=Bx(18);          % cool skin
    hrain_old=Bx(38);         % rain heat flux
    
    warm_output(ibg,1) = dT_warm;         % dT between base of warm layer and skin, across entire layer, deg C
    warm_output(ibg,2) = dz_warm;         % warm layer thickness, m
    warm_output(ibg,3) = dT_warm_to_skin; % dT between measurement and skin level due to warm layer, degC
    warm_output(ibg,4) = du_warm;         % current across entire layer, m/s ?
   
    icount=icount+1;
    
end %  data line loop

% get rid of filled values where nans are present in input data
bad_input = find(isnan(SW_dn) == 1);
% disp(['bad solar values = ' sprintf('%i',length(bad_input))]);
warm_output(bad_input,:) = nan;

%**************************************************
% Recompute entire time series of fluxes with seawater T adjusted for warm layer
%**************************************************
clear Bx
Tsea=Tsea+warm_output(:,3);

Bx=coare36vn_zrf_et(U,zu,Tair,zt,RH,zq,P,Tsea,SW_dn,LW_dn,Lat,Lon,Jd,zi,Rainrate,Ss,cp,sigH,zrf_u,zrf_t,zrf_q);

B=[Bx warm_output];    %Add all warm layer variables to the end

%************* output from routine  *****************************
%%% adds to coarevn_zrf output the following 4 vars:
% B = [ <<< outputs from main coare function >>> .... dT_warm  dz_warm  dT_warm_to_skin  du_warm ]
end

%------------------------------------------------------------------------------
function g=grv(latx)
% computes g [m/sec^2] given lat in deg
gamma=9.7803267715;
c1=0.0052790414;
c2=0.0000232718;
c3=0.0000001262;
c4=0.0000000007;
phi=latx*pi/180;
x=sin(phi);
g=gamma*(1+c1*x.^2+c2*x.^4+c3*x.^6+c4*x.^8);
end

%------------------------------------------------------------------------------
function [Press,Tseak,Tairk,Qsatsea,Qsatms,Qms,Rhoair,Rhodry]=scalarv(P0,Tsea,Tair,RH,zt)

% Compute the require scalar variables for the bulk code
% Vectorized when needed
% Inputs:
% P0    Air pressure (mb)
% Tsea  sea temp (C)
% Tair  Air temperature (C)
% RH    Relative humidity (%)

Press=P0*100;                   %ATMOSPHERIC SEA LEVEL PRESSURE (Pa)
P_tq=100*(P0 - (0.125*zt));   % P at tq measurement height (Pa)
Tseak=Tsea+273.15;              %SEA TEMPERATURE (K)
Tairk=Tair+273.15;              %AIR TEMPERATURE (K)

if length(Tsea)>1
    %**********************COMPUTES QSEA***************************
    Exx=6.1121*exp(17.502*Tsea./(240.97+Tsea));
    Exx=Exx.*(1.0007+P0*3.46E-6);
    Esatsea=Exx*0.98;
    Qsatsea=.622*Esatsea./(P0-.378*Esatsea)*1000;  %SPEC HUM SEA (g/kg)
    
    %**********************COMPUTES QAIR***************************
    Exx=6.1121*exp(17.502*Tair./(240.97+Tair));
    Exx=Exx.*(1.0007+P_tq*3.46E-6);
    Qsatms=.622*Exx./(P_tq-.378*Exx)*1000;
    Ems=Exx.*RH/100;
    Qms=.622*Ems./(P_tq-.378*Ems)*1000;   %SPEC HUM AIR (g/kg)
    E=Ems*100;
    %******************COMPUTES AIR DENSITY*******************
    Rhoair=P_tq./(Tairk.*(1+.61*Qms/1000)*287.05);
    Rhodry=(P_tq-E)./(Tairk*287.05);
else
    %**********************COMPUTES QSEA***************************
    Ex=ComputeEsat(Tsea,P0);                         %SATURATION VALUE
    Esatsea=Ex*.98;
    Qsatsea=.622*Esatsea/(P0-.378*Esatsea)*1000;  %SPEC HUM SEA (g/kg)
    
    %**********************COMPUTES QAIR***************************
    Esatms=ComputeEsat(Tair,P_tq);          %SATURATION VALUE
    Qsatms=.622*Esatms/(P_tq-.378*Esatms)*1000;
    Ems=Esatms*RH/100;
    Qms=.622*Ems/(P_tq-.378*Ems)*1000;   %SPEC HUM AIR (g/kg)
    E=Ems*100;
    
    %******************COMPUTES AIR DENSITY*******************
    Rhoair=P_tq/(Tairk*(1+.61*Qms/1000)*287.05);
    Rhodry=(P_tq-E)/(Tairk*287.05);
end
end

%------------------------------------------------------------------------------

function Exx=ComputeEsat(T,P)
%         Given temperature (C) and pressure (mb), returns
%         saturation vapor pressure (mb).
Exx=6.1121*exp(17.502*T./(240.97+T));
Exx=Exx.*(1.0007+P*3.46E-6);
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

  
%   whos psi
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