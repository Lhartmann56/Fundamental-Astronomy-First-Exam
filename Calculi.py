import numpy as np
from numpy import sin, cos, tan, sqrt, arctan, pi
import astropy.units as u

x=90*u.deg
k= np.array([[cos(x), -sin(x), 0], [sin(x), cos(x), 0], [0, 0, 1]])
M= np.asmatrix(k)
N=np.linalg.inv(M)
P1=np.array([[-136.5637561083746*u.AU, -58.86174117642695*u.AU, 0.006477038631585696*u.AU]])
P2=np.array([[-136.52944516728272*u.AU, -58.8582744676107*u.AU, 0.0064766571612710486*u.AU]])
P3=np.array([[-136.49508281581126*u.AU, -58.85478559561711*u.AU, 0.006476273252160792*u.AU]])
P4=np.array([[-123.8254221225029*u.AU, -56.75184521313076*u.AU, 0.006244869528365741*u.AU]])
P5=np.array([[-115.92607005859432*u.AU, -54.90754850330492*u.AU, 0.0060419264825281455*u.AU]])
P6=np.array([[-111.38702385389293*u.AU, -53.724113572894545*u.AU, 0.0059117034614445335*u.AU]])
P6M=np.asmatrix(P6)


print(P6M)


