
function psiuo(zet)
    x=(1.-15.*zet)**.25 
    psik=2.*log((1.+x)/2.)+log((1.+x*x)/2.)-2.*atan(x)+2.*atan(1.) 
    x=(1.-10.15*zet)**.3333 
    psic=1.5*log((1.+x+x*x)/3.)-sqrt(3.)*atan((1.+2.*x)/sqrt(3.))+4.*atan(1.)/sqrt(3.) 
    f=zet*zet/(1+zet*zet) 
    psiuo=(1-f)*psik+f*psic                                                
    if(zet>0)then 
      c=min(50.,.35*zet) 
      psiuo=-((1+1.0*zet)**1.0+.667*(zet-14.28)/exp(c)+8.525)
    endif 
END FUNCTION psiuo 
