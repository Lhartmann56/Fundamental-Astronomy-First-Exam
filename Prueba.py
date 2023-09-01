import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
from numpy import sin, cos, tan, sqrt, arctan, pi
from astropy.time import Time
import astropy.units as u
import astropy.units.cgs as q
#from astropy.coordinates import solar_system_ephemeris, get_body

#Functions definition.

def kepler (E, l, ecc):
    return E-ecc*sin(E)-(l.to('rad')).value

#Initial conditions in time J 2000

t0=Time('2023-08-29T00:00:00', scale='tdb', format='isot')
t0=t0.utc

#Jupiter's ephemeris
iota  = 1.303557901344705E+00*u.deg       #inclination
Omega = 1.0051944937794557E+02*u.deg      #Longitude of ascending node
omega = 2.733820595603509E+02*u.deg       #argument of perihelion
a = (7.781329845940255E+08*u.km).to('au') #semi-major axis
ecc = 4.816820832087293E-02               #eccentricity
n = 9.616837525063120E-07*u.deg/u.s       #mean motion
l0 = 1.839850388299260E+01*u.deg          #mean anomaly

#unit convertion.

iota = iota.to('rad')
Omega = Omega.to('rad')
omega = omega.to('rad')
n = n.to('rad/s')
l0 = l0.to('rad')

#initial conditions

tf = Time('2035-6-30T00:00:00', scale='utc', format='isot')
#tf = t0 + 15
time_steps = 1000
time_spam = np.linspace(t0, tf, time_steps)
data=np.zeros([time_steps, 3])
for i in range(time_steps):
    l = l0 + n*(time_spam[i]-t0)
    sol = root_scalar(kepler, x0=0.0, x1=pi, args=(l,ecc))
    E = sol.root*u.rad
    f = 2*arctan(sqrt((1+ecc)/(1-ecc))*tan(E/2))
    phi = f + omega
    r = a*(1-ecc*cos(E))
    x = r*(cos(Omega)*cos(phi) - cos(iota)*sin(Omega)*sin(phi))
    y = r*(sin(Omega)*cos(phi) + cos(iota)*cos(Omega)*sin(phi))
    z = r*sin(iota)*sin(phi)
    data[i,0] = x.value
    data[i,1] = y.value
    data[i,2] = z.value

#Earth's ephemeris

iotaT  = 4.164075603038432E-03*u.deg       #inclination
OmegaT = 1.498625963929686E+02*u.deg      #Longitude of ascending node
omegaT = 3.146587763491455E+02*u.deg       #argument of perihelion
aT = (1.495582533630905E+08*u.km).to('au') #semi-major axis
eccT = 1.694863932474438E-02               #eccentricity
nT = 1.141204629731537E-05*u.deg/u.s       #mean motion
l0T = 1.817846947871890E+02*u.deg          #mean anomaly

#unit convertion.

iotaT = iotaT.to('rad')
OmegaT = OmegaT.to('rad')
omegaT = omegaT.to('rad')
nT = nT.to('rad/s')
l0T = l0T.to('rad')

tfT = Time('2035-6-30T00:00:00', scale='utc', format='isot')
#tf = t0 + 15
time_stepsB = 1000
time_spamB = np.linspace(t0, tfT, time_stepsB)
data2=np.zeros([time_stepsB, 3])
for j in range(time_stepsB):
    f = l0 + nT*(time_spamB[j]-t0)
    solT = root_scalar(kepler, x0=0.0, x1=pi, args=(f,eccT))
    R = solT.root*q.rad
    g = 2*arctan(sqrt((1+eccT)/(1-eccT))*tan(R/2))
    phiT = g + omegaT
    ñ = aT*(1-eccT*cos(R))
    t = ñ*(cos(OmegaT)*cos(phiT) - cos(iotaT)*sin(OmegaT)*sin(phiT))
    u = ñ*(sin(OmegaT)*cos(phiT) + cos(iotaT)*cos(OmegaT)*sin(phiT))
    v = ñ*sin(iotaT)*sin(phiT)
    data2[j,0] = t.value
    data2[j,1] = u.value
    data2[j,2] = v.value



fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(projection='3d')
ax.scatter(data[:,0], data[:,1], data[:,2], color='k', s=0.1)     #orbit planet
ax.scatter(data2[:,0], data2[:,1], data2[:,2], color='k', s=0.1)     #orbit earth
ax.scatter(0,0,0, color='orange', s=100)                          #sun
plt.show()