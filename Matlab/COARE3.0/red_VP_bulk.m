y=load('C:\Data\cwf\dropbox\matlabstf\cwf\whoi_bulk\VP_test_data_1.txt');
jdy=y(:,1)-y(1,1)+datenum(2000,1,8)-datenum(1999,12,31);%julian day
%COLUMN A: DATE & TIME (UTC)
P=y(:,2);%COLUMN B: Atmospheric Pressure mbs
ta=y(:,3);%COLUMN C: Air Temperature (2m)
ts=y(:,4);%COLUMN D: Sea temperature (deprth 3m)
u=y(:,5);%COLUMN E: Wind speed (m/sec, 2m)
qa=y(:,6);%COLUMN F: Specific Humidity (gr/Kgr, 10m)
cloudf=y(:,7);%COLUMN G: Total Cloud Coverage (8 times the same daily value)

%**************  set height of the input data  ********
zu=2;
zt=2;
zq=10;
%******************************   Bogus in downward solar, IR flux, and BL
%height
Rs=200*ones(size(jdy));
Rl=400*ones(size(jdy));
zi=600*ones(size(jdy));
lat=35;%latitude of the site
%**************  set reference height for output of mean variables (e.g., 10-m)  ********
zref_u=10;
zref_t=10;
zref_q=10;
%%************************

rh=relhum([ta qa P]);
A=coare30vn_ref(u,zu,ta,zt,rh,zq,P,ts,Rs,Rl,lat,zi,zref_u,zref_t,zref_q);

%outputs of A
%A=[usr tau hsb hlb hbb hsbb tsr qsr zo  zot zoq Cd Ch Ce L zet dter tkt Urf Trf Qrf RHrf U10n];   
%    1    2  3   4   5   6    7   8   9   0   1  2  3 4  5    6   7   8   9   0   1   2  3
usr=A(:,1);%%   usr = friction velocity (m/s)
%   tau = wind stress (N/m^2)
%   hsb = sensible heat flux into ocean (W/m^2)
%   hlb = latent heat flux into ocean (W/m^2)
%   hbb = bouyancy flux into ocean (W/m^2)
%   hsbb = sonic bouyancy flux into ocean (W/m^2)
%   tsr = t*
%   qsr = q*
%   zo= z_o for velcoity turbulencen (m)
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
U10n=A(:,23);%   U10n = 10-m neutral wind (m/s)

yy=load('C:\Data\cwf\dropbox\matlabstf\cwf\whoi_bulk\VP_test_results_1.txt');
%COLUMN A: DATE & TIME
Rsday=yy(:,2);%COLUMN B: DAILY INSOLATION (8 same values)
Rl=yy(:,3);%COLUMN C: LONG-WAVE RADIATION
Hs=yy(:,4);%COLUMN D: SENSIBLE HEAT
Hl=yy(:,5);%COLUMN E: LATENT HEAT
Tau=yy(:,6);%COLUMN F: WIND STRESS


figure;plot(jdy,A(:,2),jdy,Tau,'x');xlabel('Julian Day');ylabel('Stress (N/m^2)');
figure;plot(jdy,A(:,3),jdy,Hs,'x');xlabel('Julian Day');ylabel('Sensible Heat (W/m^2)');
figure;plot(jdy,A(:,4),jdy,Hl,'x');xlabel('Julian Day');ylabel('Latent heat (W/m^2)');

figure;plot(jdy,u,jdy,U10n,'-x');xlabel('Julian Day');ylabel('Wind Speed (m/s)');
legend('u input height','U10n');




