import numpy as np
import scipy.integrate
from matplotlib import pyplot as plt
from mayavi import mlab

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import PowerNorm

v0 = 2
mu2 = 1
lam = mu2/(v0**2)
eps = 0.2

def V(sigma,pi):
    return - mu2/2 * (sigma**2 + pi**2) + lam/4 * (sigma**2 + pi**2)**2 - eps * sigma
rmax = 1.5*v0
sigmagrid, pigrid = np.meshgrid(*np.transpose(np.linspace((-rmax,-rmax),(rmax,rmax),100)))
Vgrid = V(sigmagrid, pigrid)

def dy(tau, y):
    sigma, pi, Psigma, Ppi = y[0], y[1], y[2], y[3]

    dsigma = Psigma
    dpi = Ppi
    dPsigma = -1/tau * Psigma + (mu2 - (lam / 2)* (pi**2 + sigma**2)) * sigma + eps/2
    dPpi = -1/tau * Ppi +  (mu2 - lam * (pi**2 + sigma**2)) * pi

    return (dsigma, dpi, dPsigma, dPpi)

y0 = (-0.2,0.1,0,0)

result = scipy.integrate.solve_ivp(dy, (0.01,100), y0)

ts = result["t"]
ys = result["y"]
sigmas, pis, Psigmas, Ppis = ys[0],ys[1],ys[2],ys[3]
# plt.plot(ts,sigmas)
# plt.plot(ts,pis)
# plt.plot(ts,Psigmas)
# plt.plot(ts,Ppis)
# plt.show()
Vtaus = V(sigmas,pis)
# plt.plot(ts,Vtaus)
# plt.show()
# plt.pcolor(sigmagrid,pigrid,Vs,cmap="viridis",norm=PowerNorm(gamma=1/4))
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# surface = ax.plot_surface(sigmagrid, pigrid, Vgrid, cmap='viridis', norm=PowerNorm(gamma=1/4))
# ax.plot(sigmas,pis,Vtaus,c="r")
# plt.show()

fig = mlab.figure(bgcolor=(1,1,1))
su = mlab.surf(sigmagrid.T, pigrid.T, Vgrid.T)
sc = mlab.points3d(sigmas,pis,Vtaus, scale_factor=0.1, scale_mode='none',
                   opacity=1.0, resolution=20, color=(1,0,0))