function qsat(y)
real :: y(2)
x=y(1) !temp
p=y(2) !pressure
es=6.112*exp(17.502*x/(x+241.0))*(1.0007+3.46e-6*p)
qsat=es*622./(p-.378*es)
end function qsat