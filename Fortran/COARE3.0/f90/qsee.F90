function qsee(ts,Pa)
real :: ts,Pa
x=ts
p=Pa
es=6.112*exp(17.502*x/(x+240.97))*.98*(1.0007+3.46e-6*p)
qsee=es*621.97/(p-.378*es)
end function