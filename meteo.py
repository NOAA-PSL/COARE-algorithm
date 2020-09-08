"""
Functions for meteorological variable conversions and functions.

Translated and vectorized from NOAA bulk model and flux processing MATLAB scripts.

Uncomment code at the end of this file and execute '%run flux.py' from the iPython command line to test functions.

Byron Blomquist, CU/CIRES, NOAA/ESRL/PSD3
v1: Oct 2014
"""

def qsat(t,p):
    """
    usage: es = qsat(t,p)
    Returns saturation vapor pressure es (mb) given t(C) and p(mb).

    After Buck, 1981: J.Appl.Meteor., 20, 1527-1532

    Returns ndarray float for any numeric object input.
    """
    from numpy import copy, asarray, exp

    t2 = copy(asarray(t, dtype=float))  # convert to ndarray float
    p2 = copy(asarray(p, dtype=float))
    es = 6.1121 * exp(17.502 * t2 / (240.97 + t2))
    es = es * (1.0007 + p2 * 3.46e-6)
    return es


def qsea(sst,p):
    """
    usage: qs = qsea(sst,p)
    Returns saturation specific humidity (g/kg) at sea surface
    given sst(C) and p(mb) input of any numeric type.

    Returns ndarray float for any numeric object input.
    """
    ex = qsat(sst,p) # returns ex as ndarray float
    es = 0.98 * ex
    qs = 622 * es /(p - 0.378 * es)
    return qs


def qair(t,p,rh):
    """
    usage: qa, em = qair(t,p,rh)
    Returns specific humidity (g/kg) and partial pressure (mb)
    given t(C), p(mb) and rh(%).

    Returns ndarray float for any numeric object input.
    """
    from numpy import copy, asarray

    rh2 = copy(asarray(rh,dtype=float))  # conversion to ndarray float
    rh2 /= 100.0                         # frational rh
    p2 = copy(asarray(p, dtype=float))
    t2 = copy(asarray(t, dtype=float))
    em = rh2 * qsat(t2,p2)
    qa = 621.97 * em / (p2 - 0.378 * em)
    return (qa, em)


def rhcalc(t,p,q):
    """
    usage: rh = rhcalc(t,p,q)
    Returns RH(%) for given t(C), p(mb) and specific humidity, q(kg/kg)

    Returns ndarray float for any numeric object input.
    """
    from numpy import copy, asarray

    q2 = copy(asarray(q, dtype=float))    # conversion to ndarray float
    p2 = copy(asarray(p, dtype=float))
    t2 = copy(asarray(t, dtype=float))
    es = qsat(t2,p2)
    em = p2 * q2 / (0.622 + 0.378 * q2)
    rh = 100.0 * em / es
    return rh


def rhoa(t,p,rh):
    """
    computes moist air density from temperature, pressure and RH

    usage: Ra = rhoa(t,p,rh)

    inputs: t (deg C), p (mb or hPa) and rh

    output: Ra = moist air density in kg/m3

    """
    Md = 0.028964                       # mol wt dry air, kg/mole
    Mv = 0.018016                       # mol wt water, kg/mole
    Tk = t + 273.15                     # deg Kelvin
    Pa = p*100.0                        # Pascals
    Rgas = 8.314                        # in m3 Pa/mol K
    Pv = (rh/100.0)*qsat(t,p)*100.0     # H2O vapor pressure in Pa
    Pd = Pa - Pv                        # pressure dry air
    Ra = (Pd*Md + Pv*Mv)/(Rgas*Tk)      # moist air density
    return Ra


def rhod(t,p):
    """
    computes dry air density from temperature and pressure

    usage: Rd = rhod(t,p)

    inputs: t (deg C), and p (mb or hPa)

    output: Rd = dry air density in kg/m3

    """
    Rd = 287.058        # gas const for dry air in J/kg K
    tk = t+273.15       # deg Kelvin
    Pa = p*100          # Pascals
    Rdry = Pa/(Rd*Tk)   # dry air density, kg/m3
    return Rd


def grv(latitude):
    """
    usage: g = grv(latitude)
    Computes gravity, g [m/sec^2] given latitude in deg.

    Ref??

    Returns ndarray float for any numeric object input.
    """
    from numpy import copy, pi, sin, asarray

    lat = copy(asarray(latitude, dtype=float))    # conversion to ndarray float
    gamma = 9.7803267715
    c1 = 0.0052790414
    c2 = 0.0000232718
    c3 = 0.0000001262
    c4 = 0.0000000007
    phi = lat * pi/180;
    x = sin(phi);
    g = gamma * (1 + c1*x**2 + c2*x**4 + c3*x**6 + c4*x**8)
    return g


def psit_26(z_L):
    """
    usage psi = psit_26(z_L)

    Computes the temperature structure function given z/L.
    """
    from numpy import exp, log, sqrt, arctan, asarray, copy
    from util import find

    zet = copy(asarray(z_L, dtype=float))    # conversion to ndarray float
    dzet = 0.35*zet
    dzet[dzet>50] = 50.           # stable
    psi = -((1 + 0.6667*zet)**1.5 + 0.6667*(zet - 14.28)*exp(-dzet) + 8.525)
    k = find(zet < 0)            # unstable
    x = (1 - 15*zet[k])**0.5
    psik = 2*log((1 + x)/2.)
    x = (1 - 34.15*zet[k])**0.3333
    psic = 1.5*log((1.+x+x**2)/3.) - sqrt(3)*arctan((1 + 2*x)/sqrt(3))
    psic += 4*arctan(1.)/sqrt(3.)
    f = zet[k]**2 / (1. + zet[k]**2.)
    psi[k] = (1-f)*psik + f*psic
    return psi


def psiu_26(z_L):
    """
    usage: psi = psiu_26(z_L)

    Computes velocity structure function given z/L
    """
    from numpy import exp, log, sqrt, arctan, min, asarray, copy
    from util import find

    zet = copy(asarray(z_L, dtype=float))    # conversion to ndarray float
    dzet = 0.35*zet
    dzet[dzet>50] = 50.           # stable
    a = 0.7
    b = 3./4.
    c = 5.
    d = 0.35
    psi = -(a*zet + b*(zet - c/d)*exp(-dzet) + b*c/d)
    k = find(zet < 0)         # unstable
    x = (1 - 15*zet[k])**0.25
    psik = 2.*log((1.+x)/2.) + log((1.+x*x)/2.) - 2.*arctan(x) + 2.*arctan(1.)
    x = (1 - 10.15*zet[k])**0.3333
    psic = 1.5*log((1.+x+x**2)/3.) - sqrt(3.)*arctan((1.+2.*x)/sqrt(3.))
    psic += 4*arctan(1.)/sqrt(3.)
    f = zet[k]**2 / (1.+zet[k]**2)
    psi[k] = (1-f)*psik + f*psic
    return psi


def psiu_40(z_L):
    """
    usage: psi = psiu_40(z_L)

    Computes velocity structure function given z/L
    """
    from numpy import exp, log, sqrt, arctan, min, asarray, copy
    from util import find

    zet = copy(asarray(z_L, dtype=float))    # conversion to ndarray float
    dzet = 0.35*zet
    dzet[dzet>50] = 50.           # stable
    a = 1.
    b = 3./4.
    c = 5.
    d = 0.35
    psi = -(a*zet + b*(zet - c/d)*exp(-dzet) + b*c/d)
    k = find(zet < 0)         # unstable
    x = (1. - 18.*zet[k])**0.25
    psik = 2.*log((1.+x)/2.) + log((1.+x*x)/2.) - 2.*arctan(x) + 2.*arctan(1.)
    x = (1. - 10.*zet[k])**0.3333
    psic = 1.5*log((1.+x+x**2)/3.) - sqrt(3.)*arctan((1.+2.*x)/sqrt(3.))
    psic += 4.*arctan(1.)/sqrt(3.)
    f = zet[k]**2 / (1.+zet[k]**2)
    psi[k] = (1-f)*psik + f*psic
    return psi



def Le_water(t, sal):
    """
    computes latent heat of vaporization for pure water and seawater
    reference:  M. H. Sharqawy, J. H. Lienhard V, and S. M. Zubair, Desalination
                and Water Treatment, 16, 354-380, 2010. (http://web.mit.edu/seawater/)
    validity: 0 < t < 200 C;   0 <sal <240 g/kg

    usage: Le_w, Le_sw = Le_water(t, sal)

    inputs: T in deg C
            sal in ppt

    output: Le_w, Le_sw in J/g (kJ/kg)

    """

    # polynomial constants
    a = [2.5008991412E+06, -2.3691806479E+03, 2.6776439436E-01,
        -8.1027544602E-03, -2.0799346624E-05]

    Le_w = a[0] + a[1]*t + a[2]*t**2 + a[3]*t**3 + a[4]*t**4
    Le_sw = Le_w*(1 - 0.001*sal)
    return (Le_w/1000.0, Le_sw/1000.0)


def uv2spd_dir(u,v):
    """
    converts u, v meteorological wind components to speed/direction
    where u is velocity from N and v is velocity from E (90 deg)

    usage spd, dir = uv2spd_dir(u, v)

    """
    import numpy as np

    spd = np.zeros_like(u)
    dir = np.zeros_like(u)

    spd = np.sqrt(u**2 + v**2)
    dir = np.arctan2(v, u)*180.0/np.pi

    return (spd, dir)


def spd_dir2uv(spd,dir):
    """
    converts wind speed / direction to u, v meteorological wind components
    where u is velocity from N and v is velocity from E (90 deg)

    usage u, v = uv2spd_dir(spd, dir)

    """
    import numpy as np

    dir2 = (np.copy(dir) + 180.0)*np.pi/180.0
    s = np.sin(dir2)
    c = np.cos(dir2)
    v = -spd*s
    u = -spd*c

    return (u, v)







# This code executes if 'run flux.py' is executed from iPython cmd line
# Uncomment lines below to test particular functions
if __name__ == '__main__':
    import numpy as np
    import matplotlib.pyplot as plt

    # TEST MET FUNCTIONS
#     print 'testing qsat(25,1015) = 31.8037'
#     print qsat(25,1015)

#     print 'testing qsea(18,1013) = 12.5615'
#     print qsea(18,1013)

#     print 'testing qair(10,1010,77) = 5.8665'
#     print qair(10,1010,77)

#     print 'testing rhcalc(14,1013,9.7121) = 98'
#     print rhcalc(14,1013,9.7121)

#     print 'testing grv(45) = 9.8062'
#     print grv(45)

    # TEST PSI FUNCTIONS
#     zet = np.arange(1,501,25)
#     psi = psit_26(zet)
#     print 'psit_26: MATLAB result *10^3'
#     print '-0.0044   -0.0870   -0.2156   -0.3799'
#     print '-0.5734   -0.7922   -1.0337   -1.2959'
#     print '-1.5772   -1.8765   -2.1927   -2.5250'
#     print '-2.8726   -3.2349   -3.6113   -4.0013'
#     print '-4.4044   -4.8202   -5.2484   -5.6886'
#     print psi
#     psi = psiu_26(zet)
#     print 'psiu_26: MATLAB result'
#     print '  -4.3926  -28.9153  -46.4143  -63.9143  -81.4143'
#     print ' -98.9143 -116.4143 -133.9143 -151.4143 -168.9143'
#     print '-186.4143 -203.9143 -221.4143 -238.9143 -256.4143'
#     print '-273.9143 -291.4143 -308.9143 -326.4143 -343.9143'
#     print psi
#     psi = psiu_40(zet)
#     print 'psiu_40: MATLAB result'
#     print '  -4.6926  -36.7153  -61.7143  -86.7143 -111.7143'
#     print '-136.7143 -161.7143 -186.7143 -211.7143 -236.7143'
#     print '-261.7143 -286.7143 -311.7143 -336.7143 -361.7143'
#     print '-386.7143 -411.7143 -436.7143 -461.7143 -486.7143'
#     print psi
