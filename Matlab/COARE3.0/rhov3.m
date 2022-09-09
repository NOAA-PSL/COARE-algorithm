
function rhov=rhov3(y)
p=y(:,1);%pressure in mb
t=y(:,2);%temperature in C
h=y(:,3);%relative humidity in %
es=6.112.*exp(17.502.*t./(t+241.0)).*(1.0007+3.46e-6*p).*h/100;
q=es*622./(p-.378*es);

tdk=273.16;
Rgas=287.1;
rhoa=p*100/(Rgas*(t+tdk)*(1+0.61*q/1000));
rhov=rhoa*q;% (g/m^3)
y=rhov;
