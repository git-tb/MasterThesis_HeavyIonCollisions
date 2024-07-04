import numpy as np
import scipy.integrate
from matplotlib import pyplot as plt
# from mayavi import mlab

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import PowerNorm

v0 = 2
mu2 = 1
lam = mu2/(v0**2)
eps = 0.1

def V(sigma,pi):
    return - mu2/2 * (sigma**2 + pi**2) + lam/4 * (sigma**2 + pi**2)**2 - eps * sigma

rmax = 1.5*v0
sigmagrid, pigrid = np.meshgrid(*np.transpose(np.linspace((-rmax,-rmax),(rmax,rmax),100)))
Vgrid = V(sigmagrid, pigrid)

def dy(tau, y):
    sigma, pi, Psigma, Ppi = y[0], y[1], y[2], y[3]

    dsigma = Psigma
    dpi = Ppi
    dPsigma = -1/tau * Psigma + (mu2 - lam* (pi**2 + sigma**2)) * sigma + eps
    dPpi = -1/tau * Ppi +  (mu2 - lam * (pi**2 + sigma**2)) * pi

    return (dsigma, dpi, dPsigma, dPpi)

def df(tau,y):
    sigma, pi, Psigma, Ppi = y[0], y[1], y[2], y[3]

    result = np.zeros((4,4))
    result[0,2] = 1
    result[1,3] = 1
    result[2,0] = - lam*(3*sigma**2 + pi**2)
    result[2,1] = -2 * lam * pi * sigma
    result[2,2] = -1/tau
    result[3,0] = -2 * lam * sigma * pi
    result[3,1] = -lam * (3 * pi**2 + sigma**2)
    result[3,3] = -1/tau

    return result

y0 = (-1.5,0.01,0,0)
tstart = 0.01
tend = 1000

result = scipy.integrate.solve_ivp(dy, (tstart,tend), y0, jac=df,method="Radau")
# result = scipy.integrate.solve_ivp(dy, (0.01,100), y0)

ts = result["t"]
ys = result["y"]
sigmas, pis, Psigmas, Ppis = ys[0],ys[1],ys[2],ys[3]
plt.plot(ts,sigmas)
plt.plot(ts,pis)
plt.show()

plt.plot(sigmas,pis)
plt.pcolor(sigmagrid,pigrid,Vgrid)
plt.show()