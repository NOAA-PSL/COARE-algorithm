Program cor3_0af
!toga coare bulk flux model version 2.6
!***************************************
!uses following subroutines:
!	cor30a.F90
!	psiu_30.F90
!	psit_30.F90
!	qsee.F90
!	grv.F90
!***************************************

!*********** basic specifications  *****
!	zu=			height of wind measurement
!	zt=			height of air temperature measurement
!	zq=			height of air humidity measurement
!	ts_depth	depth of water temperature measurement
!	jwarm=		0=no warm layer calc, 1 =do warm layer calc
!	jcool=		0=no cool skin calc, 1=do cool skin calc
!   jwave=      0= Charnock, 1=Oost et al, 2=Taylor and Yelland

!***********   input data **************
!	YYYYMMHHMMSS=		date in toga coare format, Y2K version
!	u=			wind speed (m/s), height zu
!	us=			surface current (m/s)
!	ts=			bulk surface sea temp (cent)
!	t=			air temp (cent), height zt
!	qs=			sea surface sat specific humidity (g/kg)
!	q=			air specific humidity (g/kg), height zq
!	Rs=			downward solar flux (w/m^2)
!	Rl=			downward IR flux (w/m^2)
!	zi=			inversion height (m)
!	P=			air pressure (mb)
!	rain=		rain rate (mm/hr)
!	lon=		longitude (deg E=+)
!	lat=		latitude (deg N=+)


!********** output data  ***************
!	hsb=			sensible heat flux (w/m^2)
!	hlb=			latent heat flux (w/m^2)
!	RF=			rain heat flux(w/m^2)
!	wbar=	   	webb mean w (m/s)
!	tau=			stress (nt/m^2)
!	zo=			velocity roughness length (m)
!	zot			temperature roughness length (m)
!	zoq=			moisture roughness length (m)
!	L=			Monin_Obukhov stability length
!	usr=			turbulent friction velocity (m/s), including gustiness
!	tsr			temperature scaling parameter (K)
!	qsr			humidity scaling parameter (g/g)
!	dter=			cool skin temperature depression (K)
!	dqer=			cool skin humidity depression (g/g)
!	tkt=			cool skin thickness (m)
!	Cd=			velocity drag coefficient at zu, referenced to u
!	Ch=			heat transfer coefficient at zt
!	Ce=			moisture transfer coefficient at zq
!	Cdn_10=			10-m velocity drag coeeficient, including gustiness
!	Chn_10=			10-m heat transfer coeeficient, including gustiness
!	Cen_10=			10-m humidity transfer coeeficient, including gustiness
!
IMPLICIT NONE

!fclose('all') 
!clear 
!
integer xin, ibg, indx(30), jdx(30)
real arrout(116,13), qsx(30), tsx(30), locx(30), dt(30), hwave
real :: x(19), y(30), arnl(35), Hrain(30), hnet(30), hs(30), hl(30), tau(30),hl_webb(30)
real :: u,tsnk,ta,qa,rs,rl,org,lat,lon,msp
real :: jcool, jwave
real :: a,al,b,be,cd,cdn_10,ce,cen_10,ch,chktime,chn_10,cpa,cpv,cpw,ctd1,ctd2
real :: didread,dqer,dsea,dt_wrm,dter,dtime,fxp,grav,hl_old,hlb,hs_old,hsb
integer :: i,icount,iday,ihr,imin,isec,iyr,jamset,jump,jwarm,l,le,mon
real :: loc,locx,lonx,newtime,p,q,q_pwp,qcol_ac,qjoule,qr_out,qs,qsr,rain,rf,rf_old,rgas
real :: rhoa,rhow,rich,rnl,rns,t,tau_ac,tau_old,taub,tcw,tdk,time,intime,jtime,tk_pwp,tkt,ts,ts_depth,tsea
real :: tsr,twave,us,usr,visa,visw,wbar,wg,zi,zo,zoq,zot,zq,zt,zu
double precision :: jdy,st
real, external :: grv, qsee
!real :: jdy,st
xin = 20
didread=0 

 !  jdy=x(xin,1) !time in the form YYYYMMDDHHSS.SS
 !  U=x(xin,2)  !true wind speed, m/s  etl sonic anemometer
 !  tsnk=x(xin,3) !sea snake temperature, C (0.05 m depth)
 !  ta=x(xin,4) !air temperature, C (z=14.5 m)
 !  qa=x(xin,5) !air specific humidity, g/kg (z=14.5  m)
 !  rs=x(xin,6) !downward solar flux, W/m^2 (ETL units)
 !  rl=x(xin,7) !downward IR flux, W/m^2 (ETL units)
 !  org=x(xin,8) !rainrate, mm/hr (ETL STI optical rain gauge, uncorrected)
 !  lat=x(xin,9) !latitude, deg  (SCS pcode)
 !  lon=x(xin,10) !longitude, deg (SCS pcode)
 !  msp=x(xin,11) !6-m deotg T from MSP, C    
 
! ************* open the output file
   open(unit=12,file='test3_0_ah_out.dat')
7 format(i6,i9,5i2.2,3f9.2,2f10.5,6f9.2)
!
zu=15 !anemometer ht
zt=15 !air T height
zq=15 !humidity height
ts_depth=6. !bulk water temperature sensor depth, ETL sea snake&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

jcool=1 
jwarm=1 
jwave=0 
icount=1 
!*********************  housekeep variables  ********	
qcol_ac=0 
tau_ac=0 
jtime=0 
jamset=0 
tau_old=.06 
hs_old=10 
hl_old=100 
RF_old=0 
dsea=0 
dt_wrm=0 
tk_pwp=19 
fxp=.5 
q_pwp=0 
jump=1 
!*******************  set constants  ****************
    tdk=273.16 
    grav=grv(-2.) !9.72 
    Rgas=287.1 
    cpa=1004.67    
    be=0.026 
    cpw=4000 
    rhow=1022 
    visw=1e-6 
    tcw=0.6 
    dter=0.3 

!***********   set variables not in data base  ********
    P=1008                      !air pressure
    us=0                        !surface current
    zi=600                     !inversion ht
!******************  setup read data loop  **********
open(unit=3,file='test3_0.txt')
8 format(f17.2,f7.2,f9.2,2f8.2,f9.2,f8.2,f6.2,f9.2,f9.2,f7.2)
do ibg = 1,116 !major read loop
   read(3,8) jdy, u, tsnk, ta, qa, rs, rl, org, lat, lon, msp
!*******   decode date  ************         
   st=(jdy) 
   print '(f17.0)',st
    iyr=floor(st/1e10) 
    mon=floor(st/1e8)-iyr*100 
    iday=floor((st/1e6) -iyr*1e4 - mon*100) 
    ihr=floor((st/1e4)-iyr*1e6-mon*1e4-iday*100 )
    imin=floor((st/100)- iyr*1e8-mon*1e6-iday*1e4-ihr*100)
    isec=0 
    print *,iyr, mon, iday, ihr, imin, isec 
!********   decode bulk met data ****
    if(ibg .eq. 1) ts=msp 
    tsea=msp !bulk sea surface temp
    t=ta !air temp
    qs=qsee(tsea, P) !bulk sea surface humidity
    q=qa !air humidity
    Rs=rs !downward solar flux
    Rl=rl !doward IR flux
    rain=org !rain rate
    grav=grv(lat) !9.72 
    lonx=lon !longitude
   
!*****  variables for warm layer  ***
!    time=(float(ihr*3600)+float(imin*60))/24./3600. 
    time=((float(ihr*3600)+float(imin*60))/24.) /3600. 
    intime=time 
    loc=(lonx+7.5)/15 
    locx(ibg)=loc 
    Rnl=.97*(5.67e-8*(ts-dter*jcool+273.16)**4-Rl) !oceanic broadband emissivity=0.97
    arnl(ibg)=Rnl 
    Rns=.945*Rs !oceanic albedo=0.055 daily average
!*********   set condition dependent stuff ******
    Le=(2.501-.00237*tsea)*1e6 
    cpv=cpa*(1+0.84*q/1000) 
    rhoa=P*100/(Rgas*(t+tdk)*(1+0.61*q/1000)) 
    visa=1.326e-5*(1+6.542e-3*t+8.301e-6*t*t-4.84e-9*t*t*t) 
    Al=2.1e-5*(tsea+3.2)**0.79 

!**************   apply warm layer  *********** 
    if (jwarm .EQ. 1) then    !do warm layer
        chktime=loc+(intime*24.0) 
!        chktime=loc+intime*24 
        newtime=(chktime-24*floor(chktime/24))*3600 
      if (icount .GT. 1) then !not first time thru
            if ((newtime .GT. 21600) .AND. (jump .EQ. 1)) then
         else
               jump=0 
                if (newtime .LT. jtime) then !re-zero at midnight
                    jamset=0 
                    fxp=.5 
                    tk_pwp=19 
                    tau_ac=0 
                    qcol_ac=0 
                    dt_wrm=0 
                    jump=0                    !goto 16
                else
!***************   set warm layer constants  **************
                    rich=.65    !crit rich	
                    ctd1=sqrt(2*rich*cpw/(Al*grav*rhow)) 
                    ctd2=sqrt(2*Al*grav/(rich*rhow))/(cpw**1.5) 
!************************************************
                    dtime=newtime-jtime !delta time for integrals
                    qr_out=Rnl+hs_old+hl_old+RF_old !total cooling at surface
                    q_pwp=fxp*Rns-qr_out    !tot heat abs in warm layer
                    if (q_pwp .LT. 50 .AND. jamset .EQ. 0) then !check for threshold
                        !goto 16		
                    else
                        jamset=1    !indicates threshold crossed
                        tau_ac=tau_ac+max(.002,tau_old)*dtime   !momentum integral
                    if ((qcol_ac+q_pwp*dtime) .GT. 0) then  !check threshold for warm layer existence
                     do i=1,5  !loop 5 times for fxp
                        
                                fxp=1-(0.28*0.014*(1-exp(-tk_pwp/0.014))+0.27*0.357*(1-exp(-tk_pwp/0.357))+0.45*12.82*(1-exp(-tk_pwp/12.82)))/tk_pwp 
                        !fg=fpaul(tk_pwp) fxp=fg(1) 
                        qjoule=(fxp*Rns-qr_out)*dtime 
                                if (qcol_ac+qjoule .GT. 0) then
                                    tk_pwp=min(19.,ctd1*tau_ac/sqrt(qcol_ac+qjoule)) 
                                endif 
                            enddo !  end i loop
                        else    !warm layer wiped out
                            fxp=0.75 
                            tk_pwp=19 
                            qjoule=(fxp*Rns-qr_out)*dtime 
                        endif !   end sign check on qcol_ac
                            qcol_ac=qcol_ac+qjoule  !heat integral
!***********        compute dt_warm     **************
                        if (qcol_ac .GT. 0) then
                            dt_wrm=ctd2*(qcol_ac)**1.5/tau_ac 
                        else 
                            dt_wrm=0 
                        endif 
                    endif ! end threshold check
                endif ! end midnight reset
                if (tk_pwp .LT. ts_depth) then
                    dsea=dt_wrm 
                else
                    dsea=dt_wrm*ts_depth/tk_pwp 
                endif 
            endif ! end 6am start first time thru
        endif ! end first time thru check
        jtime=newtime 
    endif ! end jwarm,  end warm layer model appl check    		
    ts=tsea+dsea 
    qs=qsee(ts, P) 
    qsx(ibg)=qs 
    tsx(20)= 1. !ts 
    a=.018 
    b=.729 
    twave=b*u 
    hwave=a*u**2.*(1+.015*u)
    x=(/u, us, ts, t, qs, q, Rs, Rl, rain, zi,  P, zu, zt, zq, lat, jcool, jwave, twave, hwave/)!set data for basic flux alogithm
!********    call modified LKB routine *******
    call cor30a(x,y) 
!****************** output from routine  *****************************
        hsb=y(1)                    !sensible heat flux W/m/m
        hlb=y(2)                    !latent
        taub=y(3)                    !stress
        zo=y(4)                     !vel roughness
        zot=y(5)                    !temp "
        zoq=y(6)                    !hum  "
        L=y(7)                      !Ob Length
        usr=y(8)                    !ustar
        tsr=y(9)                    !tstar
        qsr=y(10)                   !qstar  [g/g]
        dter=y(11)                  !cool skin delta T
        dqer=y(12)                  !  "   "     "   q
        tkt=y(13)                   !thickness of cool skin
        RF=y(14)                    !rain heat flux
        wbar=y(15)                  !webb mean w     
        Cd=y(16)                    !drag @ zu
        Ch=y(17)                    !
        Ce=y(18)                    !Dalton
        Cdn_10=y(19)                !neutral drag @ 10 [includes gustiness]
        Chn_10=y(20)                !
        Cen_10=y(21)                !
        Wg=y(22) 
!        zax(1)=jd                   !julian day
!        zax(2:10)=x(1:9)            !
!        zax(4)=tsea                 !Tsea [no cool skin]
!        zax(11:32)=y(1:22)          !
!        zax(33:35)=(/dt_wrm, tk_pwp, ts/)   !warm layer deltaT, thickness, corrected Tsea
    !*******   previous values from cwf hp basic code *****
    
    Hrain(ibg)=RF 
   !**********  new values from this code
    hnet(ibg)=Rns-Rnl-hsb-hlb-Hrain(ibg) !total heat input to ocean
    hs(ibg)=hsb 
    hl(ibg)=hlb 
    tau(ibg)=taub 
    hl_webb=rhoa*Le*wbar*qa/1000 
   !********************  save various parts of data **********************************
   !*************  create Bradley type out file
    indx(ibg)=ibg 
    arrout(ibg,1)=ibg 
    arrout(ibg,2)=jdy !output old results
    arrout(ibg,3)=hsb 
    arrout(ibg,4)=hlb 
    arrout(ibg,5)=ts-dter*jcool 
    arrout(ibg,6)=taub 
    arrout(ibg,7)=wbar 
    arrout(ibg,8)=RF 
    arrout(ibg,9)=dter 
    arrout(ibg,10)=dt_wrm 
    arrout(ibg,11)=tk_pwp 
    arrout(ibg,12)=tkt*1e3 
    arrout(ibg,13)=Wg  
    hs_old=hsb 
    hl_old=hlb 
    RF_old=RF 
    tau_old=taub 
    icount=icount+1 
!    print *, arrout(ibg,:)
    write(12,7)ibg,iyr, mon, iday, ihr, imin, isec ,arrout(ibg,3:13) 
enddo  !  data line loop
!*****************   write output file  ******
close(unit=3) 
close(unit=15)
read *, i
end PROGRAM !cor3_0af.F90

