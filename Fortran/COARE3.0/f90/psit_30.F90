
function psit_30(zet)
    x=(1.-(15*zet))**.5 
    psik=2*log((1+x)/2) 
    x=(1.-(34.15*zet))**.3333 
    psic=1.5*log((1.+x+x*x)/3.)-sqrt(3.)*atan((1.+2.*x)/sqrt(3.))+4.*atan(1.)/sqrt(3.) 
    f=zet*zet/(1+zet*zet) 
    psit_30=(1-f)*psik+f*psic   
   
    if(zet>0)then 
      c=min(50.,.35*zet) 
      psit_30=-((1.+2./3.*zet)**1.5+.6667*(zet-14.28)/exp(c)+8.525)
   endif
end FUNCTION psit_30

