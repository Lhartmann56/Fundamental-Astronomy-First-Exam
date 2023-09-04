import numpy as np
import matplotlib.pyplot as plt
from scipy .optimize import root_scalar
from numpy import sin, cos, tan, sqrt, arctan, pi, arccos
from astropy.time import Time
import astropy.units as u
from astropy.visualization import time_support
#from astropy.coordinates import solar_system_ephemeris, get_body

#Functions definition.

def kepler (E, l, ecc):
    return E-ecc*sin(E)-(l.to('rad')).value

#Initial conditions in time J 2000

t0=Time('2023-07-09T00:00:00', scale='tdb', format='isot')
t0=t0.utc

#Eart's ephemeris
iota  = 4.164075603038432E-03*u.deg       #inclination
Omega = 1.498625963929686E+02*u.deg      #Longitude of ascending node
omega = 3.146587763491455E+02*u.deg       #argument of perihelion
a = (1.495582533630905E+08*u.km).to('au') #semi-major axis
ecc = 1.694863932474438E-02               #eccentricity
n = 1.141204629731537E-05*u.deg/u.s       #mean motion
l0 = 1.817846947871890E+01*u.deg          #mean anomaly

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

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(projection='3d')
ax.scatter(data[:,0], data[:,1], data[:,2], color='k', s=0.1)     #orbit planet
ax.scatter(0,0,0, color='orange', s=100)                          #sun

#times of observation
times = [t0,
         Time('2027-07-09T21:30:00', scale='utc', format='isot'),
         Time('2028-07-09T21:30:00', scale='utc', format='isot'),
         Time('2029-07-09T21:30:00', scale='utc', format='isot')]

for t in times:
    l = l0 + n*(t-t0)
    sol = root_scalar(kepler, x0=0.0, x1=pi, args=(l,ecc))
    E = sol.root*u.rad
    f = 2*arctan(sqrt((1+ecc)/(1-ecc))*tan(E/2))
    phi = f + omega
    r = a*(1-ecc*cos(E))
    x = r*(cos(Omega)*cos(phi) - cos(iota)*sin(Omega)*sin(phi))
    y = r*(sin(Omega)*cos(phi) + cos(iota)*cos(Omega)*sin(phi))
    z = r*sin(iota)*sin(phi)
    ax.scatter(x,y,z, color='red', s=30)

#coordinates of asteroid
t1 = ['2023-07-09T10:00:00', '2023-07-09T12:00:00', '2023-07-09T14:00:00', '2023-08-03T11:00:00', '2023-08-15T08:00:00', '2023-08-21T09:00:00']
t = Time(t1, format='isot', scale='utc')
x1 = [-136.5637561083746, -136.52944516728272, -136.49508281581126, -123.8254221225029, -115.92607005859432,  -111.38702385389293]
y1 = [-58.86174117642695, -58.8582744676107, -58.85478559561711, -56.75184521313076, -54.90754850330492, -53.724113572894545]
z1 = [ 0.006477038631585696, 0.006476657161271048, 0.006476273252160792, 0.006244869528365741, 0.0060419264825281455, 0.0059117034614445335]


ax.scatter(x1,y1,z1, color='green', s=5)     
ax.set_xlabel('x (AU)')
ax.set_ylabel('y (AU)')
ax.set_zlabel('z (AU)')
ax.set_zlim(-1,1)
plt.title('Earth\'s orbit around the sun (Ecliptic J2000.0)')

r = np.zeros(len(t1))
vx = np.zeros(len(t1))
vy = np.zeros(len(t1))
vz = np.zeros(len(t1))

#radius
for i in range(len(t1)):
    r[i] = sqrt((x1[i]**2)+(y1[i]**2)+(z1[i]**2))

#velocity ----------------------------------------------------------------------------------------------------------------------------------------------------------

h = t[len(t1)-1].jd-t[len(t1)-2].jd
vx[len(t1)-1] = (x1[len(t1)-1]-x1[len(t1)-2])/h
vy[len(t1)-1] = (y1[len(t1)-1]-y1[len(t1)-2])/h
vz[len(t1)-1] = (z1[len(t1)-1]-z1[len(t1)-2])/h
for i in range(len(t1)-1):
    h = t[i-1].jd-t[i].jd
    vx[i] = (x1[i+1]-x1[i])/h
    vy[i] = (y1[i+1]-y1[i])/h
    vz[i] = (z1[i+1]-z1[i])/h



#calculate the total angular momentum -------------------------------------------------------------------------------------------------------------------------------------------------


momentum = np.zeros(len(t1))
momentum_z = np.zeros(len(t1))
for i in range(len(t1)):
    momentum[i] = sqrt(((y1[i]*vz[i]-z1[i]*vy[i])**2)+((z1[i]*vx[i]-x1[i]*vz[i])**2)+((y1[i]*vx[i]-x1[i]*vy[i])**2))
    #z component of momentum------------------------------------------------------------------------------------------------------------------------
    momentum_z[i] = y1[i]*vx[i]-x1[i]*vy[i]
t1 = ['2023-07-09T10:00:00', '2023-07-09T12:00:00', '2023-08-03T11:00:00', '2023-08-15T08:00:00', '2023-08-21T09:00:00']
nt = Time(t1, format='isot', scale='utc')

momentum = np.delete(momentum,2)
momentum_z = np.delete(momentum_z,2)

promedio_momentum = 0
promedio_momentum_z = 0
N = len(momentum)
for i in range(len(momentum)):
    promedio_momentum = promedio_momentum + momentum[i]
    promedio_momentum_z = promedio_momentum_z + momentum_z[i]
promedio_momentum = promedio_momentum/N
promedio_momentum_z = promedio_momentum_z/N
#inclination-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
iota_a = arccos(promedio_momentum_z/promedio_momentum_z)
#graph momentum -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plt.figure(figsize=(10,10))
plt.axhline(promedio_momentum, linewidth=2, color="red", label='Average of the momentum' )
plt.scatter(nt[:].jd, momentum[:], label='Momentum in each time')
plt.title('Momentum of asteroid')
plt.ylabel('Momentum')
plt.xlabel('Time')
plt.rcParams['legend.fontsize']= 10
plt.legend(loc="upper right")
plt.grid(visible=True, which='major', color=(0,0,0), linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color=(0.2,0.2,0.2), linestyle='-', alpha=0.2)
#graph momentum in z-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plt.figure(figsize=(10,10))
plt.axhline(promedio_momentum_z, linewidth=2, color="red", label='Average of the momentum' )
plt.scatter(nt[:].jd, momentum_z[:], label='Momentum in each time')
plt.title('Momentum of asteroid in z')
plt.ylabel('Momentum')
plt.xlabel('Time')
plt.rcParams['legend.fontsize']= 10
plt.legend(loc="upper right")
plt.grid(visible=True, which='major', color=(0,0,0), linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color=(0.2,0.2,0.2), linestyle='-', alpha=0.2)


plt.show()
