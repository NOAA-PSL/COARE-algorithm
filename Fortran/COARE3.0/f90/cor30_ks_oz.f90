SUBROUTINE cor30_ks_oz(x,y)
!version with shortened iteration
!
!x is a 14 element input array from the calling program and
!y is the 6 element output array that this subroutine calculates.
! the Bessel function used in this routine are from Numerical Recipes.
!
    
! Variables
implicit none
real, intent(in):: x(14)
real, intent(out)::y(6)
real u,ts,t,ta,Q,Rl,zi,P,zu,usr,hsb,hlb,Aoz,alph,scw,Beta,von,fdg,tdk
real grav,Rgas,Le,cpa,cpv,rhoa,visa,Al,be,cpw,rhow,visw,tcw,bigc,Rnl
real du,wt,wq,tsr,qsr,Bf,ug,ut,qout,dels,qcol,alq,xlamx,tkt,dter,wbar,Cd
real lam,A,sca,ha,hw,usw,b,d,rwo,zoo,rw,bes0,bes,ra,vtc,vtco

real, external :: bessk0_s,bessk1_s
u=x(1) !wind speed, m/s
ts=x(2) !water temperature, C
t=x(3)
ta=t !air temperature, C
Q=x(4)/1000 !specific humidity
Rl=x(5) !downward IR flux, W/m**2
zi=x(6) !atmospheric inversion height, m
P=x(7) !atmospheric pressure, mb
zu=x(8) !height of the wind speed data, m
usr=x(9) !friction velocity, m/s
hsb=x(10) !sensible heat flux, W/m**2
hlb=x(11) !latent heat flux, W/m**2
Aoz=x(12) !Ozone reaction rate time scale, s**-1
alph=x(13) !Ozone dimensionless solubility
scw=x(14) !Ozone schmidt number in water


     !***********   set constants *************
     Beta=1.25 
     von=0.4 
     fdg=1.00 
     tdk=273.16 
     grav=9.82 
     !*************  air constants ************
     Rgas=287.1 
     Le=(2.501-.00237*ts)*1e6 
     cpa=1004.67 
     cpv=cpa*(1+0.84*Q) 
     rhoa=P*100/(Rgas*(t+tdk)*(1+0.61*Q)) 
     visa=1.325e-5*(1+6.542e-3*t+8.301e-6*t*t-4.8e-9*t*t*t) 
     !************  cool skin constants  *******
     Al=2.1e-5*(ts+3.2)**0.79 
     be=0.026 
     cpw=4000 
     rhow=1022 
     visw=1e-6 
     tcw=0.6 
     bigc=16*grav*cpw*(rhow*visw)**3/(tcw*tcw*rhoa*rhoa) 
     
     !**************  compute aux stuff *******
    if (Rl>0) then
        Rnl=0.97*(5.56e-8*(ts+tdk)**4-Rl) 
    else
        Rnl=50 
    end if
     du=u 
     !***************   Begin bulk loop *******
     wt=hsb/rhoa/cpa 
     wq=hlb/rhoa/Le 
     
     tsr=-wt/usr 
     qsr=-wq/usr 
     Bf=-grav/ta*usr*(tsr+.61*ta*qsr) 
     if (Bf>0) then
     ug=Beta*(Bf*zi)**.333 
     else
     ug=.2 
     end if
     ut=sqrt(du*du+ug*ug) 
     qout=Rnl+hsb+hlb 
     dels=0 !ignore sw effect
     qcol=qout-dels 
     alq=Al*qcol+be*hlb*cpw/Le ! Eq. 7 Buoy flux water

     if (alq>0) then 
       xlamx=6/(1+(bigc*alq/usr**4)**.75)**.333 ! Eq 13 Saunders
     else
       xlamx=6 ! Eq 13 Saunders 
     end if
      tkt=xlamx*visw/(sqrt(rhoa/rhow)*usr) !Eq.11 Sub. thk
      dter=qcol*tkt/tcw !  Eq.12 Cool skin
   
     
     !****************   Webb et al. correection  ************
     wbar=1.61*hlb/rhoa/Le+(1+1.61*Q)*hsb/rhoa/cpa/ta 
     
     !**************   compute transfer coeffs relative to du @meas. ht **********
     Cd=(usr/du)**2 
     
lam=13.3 
A=1.85 
!CO2 variables
sca=1 
!Fairall et al. 1999 parameterization
ha=lam !neglect air side sublayer buoyancy effects
hw=lam/A/6*xlamx !includes water side buoyancy effect

usw=usr/sqrt(rhow/rhoa)*6/xlamx 
b=2/von/usw 
d=visw/scw 
rwo=1/sqrt(Aoz*d) 
zoo=b/rwo 
rw=rwo*bessk0_s(zoo)/bessk1_s(zoo) 
bes0 = bessk0_s(zoo)
bes = bessk1_s(zoo)

ra=(ha*sqrt(sca)+1./sqrt(Cd)-5+.5*log(sca)/von)/usr  !air side resistance, see fairall et al, 2000 BLM
vtc=1/((rw/alph)+ra) !non-bubble xfer velocity
vtco=1/((rwo/alph)+ra) !non-bubble xfer velocity

y=(/rwo, ra, rw, vtco, vtc, (6/xlamx)/) 
!   1  2    3   4  5   6  7   8  9   0   1    2    3   4   5      6  7  8   9   10      11      12     13     14
!
END !SUBROUTINE cor30_ks_oz
