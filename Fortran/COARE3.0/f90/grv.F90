function grv(lat)
real lat
gamma=9.7803267715
c1=0.0052790414
c2=0.0000232718
c3=0.0000001262
c4=0.0000000007
pi=3.141593

phi=lat*pi/180
x=sin(phi)
grv=gamma*(1+(c1*x**2)+(c2*x**4)+(c3*x**6)+(c4*x**8))
!print *,'grav=',grv,lat
end function grv

