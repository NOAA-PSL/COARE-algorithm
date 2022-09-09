function [y tkt]=dcoolf(usrx,Ts,xlamy,QQ,Rns)
%usr, tangential u* on ocean side
%Ts, sst in C
%xlamy=stability effect
%visa kinematic viscosity seawater  nomninaly 9e-7
%Sc  schmidt number heat in seawater  nomninally 5.9
%QQ net cooling at interface = hs +hl+rlnet
%Rns  net solar heating at interface
%x=[t h p]
usrx=usrx(:);Ts=Ts(:);xlamy=xlamy(:);QQ=QQ(:);Rns=Rns(:);
y=-4:.1:-1;
z=10.^y;
lam=10;
rhow=1022;
cpw=4e3;
Le   = (2.501-.00237*Ts)*1e6;

visa=SW_Kviscosity2(Ts,35);
Sc=5e4./(Ts.^2+155*Ts+3700);
dw=visa./Sc;
usr=usrx./xlamy;
delu=lam.*visa./usr;
delu=min(0.05,delu);
aay=dw.*delu;
bby=dw;
ccy=.4*usr;
ddy=sqrt(4*dw.*delu*.4.*usr-dw.^2);
ggy=aay+bby.*z+ccy.*z.^2;size(ggy);
delof=zeros(length(usr),length(z));
[m n]=size(ggy);
delof=1/2./ccy.*log(ggy./aay)+2*(delu-bby/2./ccy)./ddy.*(atan((bby+2*ccy.*z)./ddy)-atan(bby./ddy));

%Paulson and Simpson solar absorption coeffs
LL=[34.8 2.27 3.15e-2 5.48e-3 8.32e-4 1.27e-4 3.13e-4 7.82e-5 1.44e-5];
FF=[.237 .36 .179 .087 .08 .0246 .025 .007 .0004];
FF(6:9)=0;%%   Wick et al adjustment to reduce warm skins
Gp=zeros(4,9);Gt=zeros(m,length(z),9);Gp0z=Gt;
aa=aay;
bb=bby;
cc=ccy;
dd=ddy;

for i=1:5% iteration not need for 6:9 since coeffs are set to zero
   L=LL(i);%8e-4;%L=2.2;
gm=dd;j=sqrt(-1);
A=exp((bb-j*gm)./(2*cc.*L))./(2*cc.*j.*gm);
z1=z;
B1=(bb+2*cc.*z1-j*gm)./(2.*cc.*L);
F1=(j*gm-(bb-2*cc.*delu)).*(-expint(B1)-j*pi*0);
B2=(bb+2*cc.*z1+j*gm)./(2*cc.*L);
F2=(j*gm+(bb-2*cc.*delu)).*(-expint(B2)-j*pi*0);
G=A.*(F1+exp(j*gm./(cc.*L)).*F2);

B10=(bb-j*gm)./(2*cc*L);
F10=(j*gm-(bb-2*cc.*delu)).*(-expint(B10)-j*pi*0);
B20=(bb+j*gm)./(2*cc.*L);
F20=(j*gm+(bb-2*cc.*delu)).*(-expint(B20)-j*pi*0);
G0=A.*(F10+exp(j*gm./(cc.*L)).*F20);

%Gp(:,i)=-(G(iz)-G0-delof(iz))./z(iz)*d;%f(del) at selected depths
Gt(:,:,i)=-FF(i)*(G-G0-delof);%Qo/rho/cp
%Gp0z(:,:,i)=FF(i)*(1-L./z.*(1-exp(-z./L)));%f(del) select for old coare function
end;
Gtp=permute(Gt,[1 3 2]);
y=[];mx=[];up=[];
for k=1:m
up(k,:)=real(sum(Gtp(k,:,:)));
%%%%%%%%%%%%%  set sublayer thickness as 5*delu
ii=find(z>5*delu(k));
if isempty(ii)
    mx(k)=n;
else
    mx(k)=ii(1);
end;
y(k)=QQ(k).*delof(k,mx(k))/rhow/cpw-Rns(k).*up(k,mx(k))/rhow/cpw;
end;
tkt=5*delu';
%QQ(1).*delof(1,mx(1))/rhow/cpw

%%%%%%%%%%%%%%%%%%%%  total coolskin from surface to 5*dewlu



