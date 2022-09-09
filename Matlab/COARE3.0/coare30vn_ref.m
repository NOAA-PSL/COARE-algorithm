function A=coare30vn(u,zu,t,zt,rh,zq,P,ts,Rs,Rl,lat,zi,zref_u,zref_t,zref_q,U10n)
%
% A=[usr tau hsb hlb hbb hsbb tsr qsr zot zoq Cd Ch Ce L zet dter tkt Urf Trf Qrf RHrf];    
disp('WHOI vectorized COARE 3.0a');
% Vectorized version of COARE3 code (Fairall et al, 2003) with cool skin option
% retained but warm layer and surface wave options removed. Assumes u,t,rh,ts 
% are vectors; sensor heights zu,zt,zl, latitude lat, and PBL height zi are
% constants; air pressure P and radiation Rs,Rl may be vectors or constants. 
% Default values are assigned for P,Rs,Rl,lat,and zi if these data are not 
% available.  Input NaNs to indicate no data. Defaults should be set to 
% representative regional values if possible.
%
% Input:  
%
%     u = relative wind speed (m/s) at height zu(m)
%     t = bulk air temperature (degC) at height zt(m)
%    rh = relative humidity (%) at height zq(m)
%     P = surface air pressure (mb) (default = 1015)
%    ts = water temperature (degC)
%    Rs = downward shortwave radiation (W/m^2) (default = 150) 
%    Rl = downward longwave radiation (W/m^2) (default = 370)
%   lat = latitude (default = +45 N)
%    zi = PBL height (m) (default = 600m)
%    zref_u,zref_t,zreef_q are reference heights for wind, temp, spec hum
%    values.  If you use zref_u=3, then 
%

%outputs of A
%A=[usr tau hsb hlb hbb hsbb tsr qsr zot zoq Cd Ch Ce L zet dter tkt Urf Trf Qrf RHrf];    
%   usr = friction velocity (m/s)
%   tau = wind stress (N/m^2)
%   hsb = sensible heat flux into ocean (W/m^2)
%   hlb = latent heat flux into ocean (W/m^2)
%   hbb = bouyancy flux into ocean (W/m^2)
%   hsbb = sonic bouyancy flux into ocean (W/m^2)
%   tsr = t*
%   qsr = q*
%   zot = z_o for temperature (m)
%   zoq = z_o for humidity (m)
%    Cd = wind stress transfer coefficient at height zu   
%    Ch = sensible heat transfer coefficient at height zt   
%    Ce = latent heat transfer coefficient at height zq
%     L = Obukhov length scale (m) 
%   zet = Monin-Obukhov stability parameter zu/L 
%  dter = cool-skin temperature depression (degC)
%   tkt = cool-skin thickness (m)
%   Urf = windspeed at reference height, zref_u (m/s)
%   Trf = temperture at reference height, zref_t (C)
%   Qrf = specific humidity at reference height, zref_q (g/kg)
%   RHrf = relative humidity at reference height (%)

% Notes: 1) u is the relative wind speed, i.e., the magnitude of the
%           difference between the wind (at zu) and ocean surface current 
%           vectors.
%        2) Set jcool=0 in code if ts is true surface skin temperature,
%           otherwise ts is assumed the bulk temperature and jcool=1.
%        3) Set P=NaN to assign default value if no air pressure data 
%           available. 
%        4) Set Rs=NaN, Rl=NaN if no radiation data available.  This assigns 
%           default values to Rs, Rl so that cool skin option can be applied. 
%        5) Set lat=NaN and/or zi=NaN to assign default values if latitude
%           and/or PBL height not given. 
%        6) The code to compute the heat flux caused by precipitation is 
%           included if rain data is available (default is no rain).
%        7) Code updates the cool-skin temperature depression dter and thickness
%           tkt during iteration loop for consistency.
%        8) Number of iterations set to nits = 6.

% Reference:
%
%  Fairall, C.W., E.F. Bradley, J.E. Hare, A.A. Grachev, and J.B. Edson (2003),
%  Bulk parameterization of air sea fluxes: updates and verification for the 
%  COARE algorithm, J. Climate, 16, 571-590.

% Code history:
% 
% 1. 12/14/05 - created based on scalar version coare26sn.m with input
%    on vectorization from C. Moffat.  
% 2. 12/21/05 - sign error in psiu_26 corrected, and code added to use variable
%    values from the first pass through the iteration loop for the stable case
%    with very thin M-O length relative to zu (zetu>50) (as is done in the 
%    scalar coare26sn and COARE3 codes).


%-----------------------------------------------------------------------

% set jcool=1 if Ts is bulk ocean temperature (default),
%     jcool=0 if Ts is true ocean skin temperature. 
      jcool=0;

% convert input to column vectors
     u=u(:);t=t(:);
     rh=rh(:);P=P(:);ts=ts(:);
     Rs=Rs(:);Rl=Rl(:);lat=lat(:);zi=zi(:);     
     N=length(u);zt=zt(:);zu=zu(:);zq=zq(:);
      
% set local variables to default values if input is NaN
     if isnan(P); P=1015*ones(N,1); end;      % pressure
     if isnan(Rs); Rs=370*ones(N,1); end;     % incident shortwave radiation 
     if isnan(Rl); Rl=150*ones(N,1); end;     % incident longwave radiation
     if isnan(lat); lat=45; end;              % latitude
     if isnan(zi); zi=600; end;               % PBL height      

% input variable u is assumed relative wind speed (magnitude of difference
% between wind and surface current vectors). to follow orginal Fairall code, set
% surface current speed us=0. if us data are available, construct u prior to 
% using this code.
     us = 0*u;    

% convert rh to specific humidity
     Qs = qsat26sea(ts,P)./1000;    % surface water specific humidity (g/kg)
     Q  = qsat26air(t,P,rh)./1000;  % specific humidity of air (g/kg)

% set rain to zero 
     rain = 0*u; % rain rate (mm/hr) - keep as option
   
     %***********  set constants **********************************************
     %zref=10.0;
     Beta = 1.2;
     von  = 0.4;
     fdg  = 1.00;
     tdk  = 273.16;
     grav = grv(lat); 
     
     %***********  air constants **********************************************
     Rgas = 287.1;
     Le   = (2.501-.00237*ts)*1e6;
     cpa  = 1004.67;
     cpv  = cpa*(1+0.84*Q);
     rhoa = P*100./(Rgas*(t+tdk).*(1+0.61*Q));
     visa = 1.326e-5*(1+6.542e-3.*t+8.301e-6*t.^2-4.84e-9*t.^3);
     
     %***********  cool skin constants  ***************************************
     Al   = 2.1e-5*(ts+3.2).^0.79; 
     be   = 0.026;
     cpw  = 4000;
     rhow = 1022;
     visw = 1e-6;
     tcw  = 0.6;
     bigc = 16*grav*cpw*(rhow*visw)^3./(tcw.^2*rhoa.^2);
     wetc = 0.622*Le.*Qs./(Rgas*(ts+tdk).^2);                
     
     %***********  net radiation fluxes ***************************************
     Rns = 0.945.*Rs; % albedo correction
     Rnl = 0.97*(5.67e-8*(ts-0.3*jcool+tdk).^4-Rl); % initial value
     
%****************  begin bulk loop ********************************************
     
     %***********  first guess ************************************************
     du = u-us;
     dt = ts-t-.0098*zt;
     dq = Qs-Q;
     ta = t+tdk;
     ug = 0.5;
     dter  = 0.3; 
     ut    = sqrt(du.^2+ug.^2);
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
     Ribcu = -zu./zi/.004./Beta.^3;
     Ribu  = -grav.*zu./ta.*((dt-dter*jcool)+.61*ta.*dq)./ut.^2;
     zetu = CC.*Ribu.*(1+27/9*Ribu./CC);
     k50=find(zetu>50); % stable with very thin M-O length relative to zu
     k=find(Ribu<0);if ~isempty(k); zetu(k)=CC(k).*Ribu(k)./(1+Ribu(k)./Ribcu(k)); clear k;end;
     L10 = zu./zetu;
     usr = ut.*von./(log(zu./zo10)-psiu_26(zu./L10));
     tsr = -(dt-dter*jcool).*von*fdg./(log(zt./zot10)-psit_26(zt./L10));
     qsr = -(dq-wetc.*dter*jcool)*von*fdg./(log(zq./zot10)-psit_26(zq./L10));
     tkt = 0.001*ones(N,1);
     charn = 0.011*ones(N,1);
     k=find(ut>10); charn(k)=0.011+(ut(k)-10)/(18-10)*(0.018-0.011); clear k;
     k=find(ut>18); charn(k)=.018; clear k;
     nits=6; % number of iterations 
     
  %**************  bulk loop **************************************************
     
  for i=1:nits 
        
     zet=von.*grav.*zu./ta.*(tsr +.61*ta.*qsr)./(usr.^2);
     zo=charn.*usr.^2./grav+0.11*visa./usr; % surface roughness
     rr=zo.*usr./visa;
     L=zu./zet;
     zoq=min(1.15e-4,5.5e-5./rr.^0.6); % mositure roughness 
     zot=zoq; % temperature roughness
     usr=ut*von./(log(zu./zo)-psiu_26(zu./L));
     tsr=-(dt-dter*jcool).*von.*fdg./(log(zt./zot)-psit_26(zt./L));
     qsr=-(dq-wetc.*dter*jcool)*von*fdg./(log(zq./zoq)-psit_26(zq./L));
     tvsr=tsr+0.61*ta.*qsr;
     tssr=tsr+0.51*ta.*qsr;
     Bf=-grav./ta.*usr.*tvsr;
     ug=0.2*ones(N,1);
     k=find(Bf>0);if ~isempty(k); ug(k)=Beta*(Bf(k).*zi(k)).^.333; clear k;end;
     ut=sqrt(du.^2+ug.^2);
     hsb=-rhoa*cpa.*usr.*tsr;
     hlb=-rhoa.*Le.*usr.*qsr;
     qout=Rnl+hsb+hlb;
     dels=Rns.*(0.065+11*tkt-6.6e-5./tkt.*(1-exp(-tkt/8.0e-4)));
     qcol=qout-dels;
     alq=Al.*qcol+be*hlb*cpw./Le; 
     xlamx=6.0*ones(N,1);
     tkt=min(0.01, xlamx.*visw./(sqrt(rhoa./rhow).*usr));
     k=find(alq>0); xlamx(k)=6./(1+(bigc(k).*alq(k)./usr(k).^4).^0.75).^0.333;  
     tkt(k)=xlamx(k).*visw./(sqrt(rhoa(k)./rhow).*usr(k)); clear k; 
     dter=qcol.*tkt./tcw; 
     dqer=wetc.*dter;
     Rnl=0.97*(5.67e-8*(ts-dter.*jcool+tdk).^4-Rl); % update dter   
     if i==1; % save first iteration solution for case of zetu>50;
         usr50=usr(k50);tsr50=tsr(k50);qsr50=qsr(k50);L50=L(k50);
         zet50=zet(k50);dter50=dter(k50);dqer50=dqer(k50);tkt50=tkt(k50);
     end
     
  end  
  
  % insert first iteration solution for case with zetu>50
     usr(k50)=usr50;tsr(k50)=tsr50;qsr(k50)=qsr50;L(k50)=L50;
     zet(k50)=zet50;dter(k50)=dter50;dqer(k50)=dqer50;tkt(k50)=tkt50;
  
%****************  compute fluxes  ********************************************
     tau=rhoa.*usr.*usr.*du./ut;  % wind stress
     hsb=rhoa.*cpa.*usr.*tsr;     % sensible heat flux
     hlb=rhoa.*Le.*usr.*qsr;      % latent heat flux
     hbb=rhoa.*cpa.*usr.*tvsr;    % buoyancy flux
     hsbb=rhoa.*cpa.*usr.*tssr;   % sonic heat flux

%*****  compute transfer coeffs relative to ut @ meas. ht  ********************
     Cd=tau./rhoa./ut./max(.1,du);
     Ch=-usr.*tsr./ut./(dt-dter*jcool);
     Ce=-usr.*qsr./(dq-dqer*jcool)./ut;
     
%***  compute 10-m neutral coeff relative to ut (output if needed) ************
%    Cdn_10=von.^2./log(10./zo).^2;
%    Chn_10=von.^2.*fdg./log(10./zo)./log(10./zot);
%    Cen_10=von.^2.*fdg./log(10./zo)./log(10./zoq);

Urf=u+usr./von.*(log(zref_u./zu)-psiu_26(zref_u./L)+psiu_26(zu./L));
Trf=t+tsr./von./fdg.*(log(zref_t./zt)-psit_26(zref_t./L)+psit_26(zt./L))+0.0098*(zt-zref_t);
Qrf=Q+qsr./von./fdg.*(log(zref_q./zq)-psit_26(zref_q./L)+psit_26(zq./L));
RHrf=RHcalc(Trf,P,Qrf);
Qrf=Qrf*1000;
U10n=u+usr./von.*(log(10./zu)+psiu_26(zu./L));

%****************  output  ****************************************************
     
A=[usr tau hsb hlb hbb hsbb tsr qsr zo zot zoq Cd Ch Ce L zet dter tkt Urf Trf Qrf RHrf];   
%   1    2  3   4   5   6    7   8   9   0  1   2  3  4  5 6   7   8    9   0   1  2



%******** rain heat flux (save to use if desired) *****************************
% dwat=2.11e-5*((t+tdk)/tdk)^1.94; %! water vapour diffusivity
% dtmp=(1.+3.309e-3*t-1.44e-6*t*t)*0.02411/(rhoa*cpa); %! heat diffusivity
% dqs_dt=Q.*Le/(Rgas*(t+tdk).^2); %! Clausius-Clapeyron
% alfac= 1/(1+0.622*(dqs_dt.*Le*.dwat)./(cpa.*dtmp)); %! wet bulb factor
% RF= rain.*alfac.*cpw.*((ts-t-dter*jcool)+(Qs-Q-dqer*jcool)*Le./cpa)/3600;
  
%------------------------------------------------------------------------------
function psi=psit_26(zet)
% computes temperature structure function
   dzet=min(50,0.35*zet); % stable
   psi=-((1+0.6667*zet).^1.5+0.6667*(zet-14.28).*exp(-dzet)+8.525); 
   k=find(zet<0); % unstable
   x=(1-15*zet(k)).^0.5;
   psik=2*log((1+x)./2);
   x=(1-34.15*zet(k)).^0.3333;
   psic=1.5*log((1+x+x.^2)./3)-sqrt(3)*atan((1+2*x)./sqrt(3))+4*atan(1)./sqrt(3);
   f=zet(k).^2./(1+zet(k).^2);
   psi(k)=(1-f).*psik+f.*psic; 

%------------------------------------------------------------------------------
function psi=psiu_26(zet)
% computes velocity structure function 
   dzet=min(50,0.35*zet); % stable
   psi=-((1+zet)+0.6667*(zet-14.28).*exp(-dzet)+8.525);
   k=find(zet<0); % unstable
   x=(1-15*zet(k)).^0.25;
   psik=2*log((1+x)/2)+log((1+x.*x)/2)-2*atan(x)+2*atan(1);
   x=(1-10.15*zet(k)).^0.3333;
   psic=1.5*log((1+x+x.^2)/3)-sqrt(3)*atan((1+2*x)./sqrt(3))+4*atan(1)./sqrt(3);
   f=zet(k).^2./(1+zet(k).^2);
   psi(k)=(1-f).*psik+f.*psic;
   
%------------------------------------------------------------------------------
function exx=bucksat(T,P)
% computes saturation vapor pressure [mb] 
% given T [degC] and P [mb]   
   exx=6.1121.*exp(17.502.*T./(T+240.97)).*(1.0007+3.46e-6.*P);

%------------------------------------------------------------------------------
function qs=qsat26sea(T,P)
% computes surface saturation specific humidity [g/kg] 
% given T [degC] and P [mb] 
   ex=bucksat(T,P);
   es=0.98*ex; % reduction at sea surface
   qs=622*es./(P-0.378*es);

%------------------------------------------------------------------------------
function qs=qsat26air(T,P,rh)
% computes specific humidity [g/kg] 
% given T [degC], RH [%] and P [mb] 
   es=bucksat(T,P);
   em=0.01*rh.*es;
   qs=622*em./(P-0.378*em);

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
   
%------------------------------------------------------------------------------
function RHrf=RHcalc(T,P,Q)
% computes relative humidity given T,P, & Q 

   es=6.1121.*exp(17.502.*T./(T+240.97)).*(1.0007+3.46e-6.*P);
   em=Q.*P./(0.378.*Q+0.622);
   RHrf=100*em./es;

