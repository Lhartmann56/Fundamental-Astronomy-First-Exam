import numpy as np
import matplotlib.pyplot as plt
from scipy .optimize import root_scalar
from numpy import sin, cos, tan, sqrt, arctan, pi, arccos, arccosh, sinh, absolute
from astropy.time import Time
import astropy.units as u
from astropy.visualization import time_support
from astropy import constants as const
#from astropy.coordinates import solar_syst
#definiendo las funciones
def kepler (E, l, ecc):
    return E-ecc*sin(E)-(l.to('rad')).value


t1 = ['2023-07-09T10:00:00', '2023-07-09T12:00:00', '2023-07-09T14:00:00']
t = Time(t1, format='isot', scale='utc')
x1 = [0.7059032210909959, 0.7064970232872655, 0.7070906151926998]
y1 = [-1.769547135717722, -1.7682131326526134,-1.7668786032741868]
z1 = [0.0001947177390506007, 0.00019457094778662054, 0.0001944240986080829]
#radius---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
r = np.zeros(len(t1))
for i in range(len(t1)):
    r[i] = sqrt((x1[i]**2)+(y1[i]**2)+(z1[i]**2))
print('matriz de radio -------------------------------------------------------')
print('la matriz con la separación es', r)

#calculo de la velocidad a partir de la velocidad:
vx = np.zeros(len(t1))
vy = np.zeros(len(t1))
vz = np.zeros(len(t1))
V = np.zeros(len(t1))

#primera derivada
h = t[1].byear-t[0].byear
vx[0] = (x1[1]-x1[0])/h
vy[0] = (y1[1]-y1[0])/h
vz[0] = (z1[1]-z1[0])/h
#derivada central
h = t[2].byear-t[1].byear
vx[1] = (x1[2]-x1[1])/h
vy[1] = (y1[2]-y1[1])/h
vz[1] = (z1[2]-z1[1])/h
#ultima derivada
h = t[len(t1)-1].byear - t[len(t1)-2].byear
vx[len(t1)-1] = (x1[len(t1)-1]-x1[len(t1)-2])/h
vy[len(t1)-1] = (y1[len(t1)-1]-y1[len(t1)-2])/h
vz[len(t1)-1] = (z1[len(t1)-1]-z1[len(t1)-2])/h

print('La velocidad en el eje x es: ', vx)
print('La velocidad en el eje y es: ', vy)
print('La velocidad en el eje z es: ', vz)
for i in range(len(t1)):
    V[i]=sqrt(vx[i]**2+vy[i]**2+vz[i]**2)
print('la velocidad en magnitud por posición es:--------------------------------------------------------------')
print('La rapides para cada medición es ', V)

#Calculo de la velocidad radial -----------------------------------------------------------------------------------------------------------------------
Ra_v = np.zeros(len(t1))
for i in range(len(t1)):
    Ra_v[i] = (x1[i]*vx[i]+y1[i]*vy[i]+z1[i]*vz[i])/r[i]
print('La rapidez radial para cada toma de datos es -----------------------------------------------------------')
print('La rapidez radial para cada punto es:', Ra_v)
#Calculo del momento angular --------------------------------------------------------------------------------------------------------------------------

L = np.zeros(len(t1))
Lx = np.zeros(len(t1))
Ly = np.zeros(len(t1))
Lz = np.zeros(len(t1))
Lt = 0
Ltx = 0
Lty = 0
Ltz = 0
for i in range(len(t1)):
    #x component of momentum------------------------------------------------------------------------------------------------------------------------
    Lx[i] = y1[i]*vz[i]-z1[i]*vy[i]
    Ltx = Ltx + Lx[i]
    #y component of momentum------------------------------------------------------------------------------------------------------------------------
    Ly[i] = z1[i]*vx[i]-x1[i]*vz[i]
    Lty = Lty + Ly[i]
    #z component of momentum------------------------------------------------------------------------------------------------------------------------
    Lz[i] = x1[i]*vy[i]-y1[i]*vx[i]
    Ltz = Ltz + Lz[i]
    #total momentum-------------------------------------------------------------------------------------------------------------------------------
    L[i] = sqrt(((Lx[i])**2)+((Ly[i])**2)+((Lz[i])**2))
    Lt = Lt + L[i]
Lt = Lt/len(L)
Ltx = Ltx/len(L)
Lty = Lty/len(L)
Ltz = Ltz/len(L)
print('El momento rotacional para el asteroide por toma de dato es: ------------------------------------------')
#print('La matriz momenum es: ', L)
print('el promedio del momento es: ', Lt)
#---------------------------------------------------------------------------------------------------------------------------------------------------
#Calculo de los parametros orbitales:
#---------------------------------------------------------------------------------------------------------------------------------------------------
#iota angulo ---------------------------------------------------------------------------------------------------------------------------------------
iota = arccos(Ltz/Lt)
print('El ángulo de inclinación iota es: ', iota*180/pi)
#Omega Angulo --------------------------------------------------------------------------------------------------------------------------------------
N = sqrt((Lty**2)+(Ltx**2))
OMEGA = arccos(-Lty/N)
print('El angulo OMEGA es: ', OMEGA*180/pi)
#excentricidad -------------------------------------------------------------------------------------------------------------------------------------
# Constantes de movimiento masa solar M y Constante de Gravitaciónb G
mu = const.G * const.M_sun
#cambio del sistema de unidades.
mu = mu.to('au^3/yr^2')
#se quitan las unidades para operar mejor
mu = mu.value
#se usan dos métodos para calcular las excentricidades
#método no 1-----------------------------------------------------------------------------------------------------------------------------------------
ecc1 = np.zeros(len(t1))
ecca1t = 0
for i in range(len(t1)):
    ecc1[i] = sqrt((r[i]**2/mu**2)*(((V[i]**2)-(mu/r[i]))**2+(Ra_v[i]*V[i])**2))
    ecca1t = ecca1t + ecc1[i]
ecca1t = ecca1t / len(ecc1)
print('Excentricidad calculada con el método no 1 ------------------------------------------------')
#print('La eccentricidad por dato es: ', ecc1)
print('el promedio de la excentricidad es: ', ecca1t)
# metodo no 2 ---------------------------------------------------------------------------------------------------------------------------------------
ecc2 = np.zeros(len(t1))
ecca2t = 0
for i in range(len(t1)):
    ecc2[i] = (1/mu)*sqrt(((2*mu-r[i]*V[i]**2)*r[i]*Ra_v[i]**2)+(mu-r[i]*V[i]**2)**2)
    ecca2t = ecca2t + ecc2[i]
ecca2t = ecca2t / len(ecc2)
print('Excentricidad calculada con el método no 2 ------------------------------------------------')
#print('La eccentricidad por dato es: ', ecc2)
print('el promedio de la excentricidad es: ', ecca2t)
#Calculo del semieje mayor a------------------------------------------------------------------------------------------------------------------------
P = (Lt**2)/mu
a = P/(1-(ecca1t**2))
a = absolute(a)
print('el semieje mayor: ', a)
#calculo del movimiento medio -----------------------------------------------------------------------------------------------------------------------
n = sqrt(mu/a**3)
print('el movimiento medio del asteroide es: ', n)
#omega-----------------------------------------------------------------------------------------------------------------------------------------------
omega = np.zeros(len(t1))
omegat = 0
for i in range(len(t1)):
    C1 = ((V[i]**2/mu)-(1/r[i]))
    C2 = Ra_v[i]*r[i]/mu
    C3 = 1/(N*ecca1t)
    e_cross_N = C3*(-(C1*x1[i]-C2*vx[i])*Ly[i]+(C1*y1[i]-C2*vy[i])*Lx[i]+(C1*z1[i]-C2*vz[i]))
    omega[i] = arccos(e_cross_N )
    omegat = omegat + omega[i] 
omegat = omegat/len(t1)
print('el valor de el angulo omega o argumento del perielio es: --------------------------------')
#print('Los valores angulares para cada coordenada son: ', omega*180/pi)
print('el argumento del perielioi es: ', omegat*180/pi)
#anomalia media----------------------------------------------------------------------------------------------------------------------------------------
#primero se halla Eo
Eo = arccos((1-(r[0]/a))*(1/ecca1t))
Lo = Eo - ecca1t*sin(Eo)
print((1-(r[0]/a))*(1/ecca1t))
#Eo = arccosh((1-(r[0]/a))*(1/ecca1t))
#Lo = Eo - ecca1t*sinh(Eo)
print('La anomalía media resulta ser -------------------------------------------------------------')
print('Anomalia media: ', Lo)
#-------------------------------------------------------------------------------------------------------------------------------------------------------
#condiciones iniciales para la orbit del asteroide
#-------------------------------------------------------------------------------------------------------------------------------------------------------
t0=Time('2023-07-09T00:00:00', scale='tdb', format='isot')
t0=t0.utc
iota = iota*u.rad
OMEGA = OMEGA*u.rad
omegat = omegat*u.rad
a = a*u.au
e = ecca2t
n = n*u.rad/u.s
Lo = Lo*u.rad


#unit convertion.

iota = iota.to('rad')
Omega = OMEGA.to('rad')
omegat = omegat.to('rad')
n = n.to('rad/s')
Lo = Lo.to('rad')

#condicion inicial -------------------------------------------------------------------------------------------------------------------------------------
tf = Time('2035-07-30T00:00:00', scale='utc', format='isot')
#tf = t0 + 15
time_steps = 1000
time_spam = np.linspace(t0, tf, time_steps)

data=np.zeros([time_steps, 3])

for i in range(time_steps):
    l = Lo + n*(time_spam[i]-t0)
    sol = root_scalar(kepler, x0=0.0, x1=pi, args=(l,e))
    E = sol.root*u.rad
    f = 2*arctan(sqrt((1+e)/(1-e))*tan(E/2))
    phi = f + omegat
    r = a*(1-e*cos(E))
    x = r*(cos(OMEGA)*cos(phi) - cos(iota)*sin(OMEGA)*sin(phi))
    y = r*(sin(OMEGA)*cos(phi) + cos(iota)*cos(OMEGA)*sin(phi))
    z = r*sin(iota)*sin(phi)
    data[i,0] = x.value
    data[i,1] = y.value
    data[i,2] = z.value

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(projection='3d')
ax.scatter(data[:,0], data[:,1], data[:,2], color='k', s=0.1)     #orbit planet
ax.scatter(0,0,0, color='orange', s=100)                          #sun

times = [t0,
         Time('2027-07-09T21:30:00', scale='utc', format='isot'),
         Time('2028-07-09T21:30:00', scale='utc', format='isot'),
         Time('2029-07-09T21:30:00', scale='utc', format='isot')]

for t in times:
    l = Lo + n*(t-t0)
    sol = root_scalar(kepler, x0=0.0, x1=pi, args=(l,e))
    E = sol.root*u.rad
    f = 2*arctan(sqrt((1+e)/(1-e))*tan(E/2))
    phi = f + omegat
    r = a*(1-e*cos(E))
    x = r*(cos(OMEGA)*cos(phi) - cos(iota)*sin(OMEGA)*sin(phi))
    y = r*(sin(OMEGA)*cos(phi) + cos(iota)*cos(OMEGA)*sin(phi))
    z = r*sin(iota)*sin(phi)
    ax.scatter(x,y,z, color='red', s=30)
plt.show()