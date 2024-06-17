import glob
import pandas as pd
from matplotlib import pyplot as plt
import scipy.interpolate
import numpy as np
import scipy.integrate
import scipy.special
import functools

print("[ DEBUG ] Loading data...")
# print(glob.glob("./../*"))
df = pd.read_csv("./../Mathematica/data/ExampleFreezeOut.csv")
# print(df)

####################################
# READ IN AND PROCESS THE INPUT DATA
####################################
alpha_arr = df["alpha"].to_numpy()
tau_arr = df["tau"].to_numpy()
r_arr = df["r"].to_numpy()
Dtau_arr = df["Dtau"].to_numpy()
Dr_arr = df["Dr"].to_numpy()
ur_arr = df["ur"].to_numpy()
utau_arr = df["utau"].to_numpy()

####################################
# interpolate to get smooth functions
tau = scipy.interpolate.CubicSpline(alpha_arr, tau_arr,extrapolate=True)
r = scipy.interpolate.CubicSpline(alpha_arr, r_arr,extrapolate=True)
Dtau = scipy.interpolate.CubicSpline(alpha_arr, Dtau_arr,extrapolate=True)
Dr = scipy.interpolate.CubicSpline(alpha_arr, Dr_arr,extrapolate=True)
ur = scipy.interpolate.CubicSpline(alpha_arr, ur_arr,extrapolate=True)
utau = scipy.interpolate.CubicSpline(alpha_arr, utau_arr,extrapolate=True)

####################################
# from interpolation function, define smoother sampled data
alphas = np.linspace(0,np.pi/2,1500)
dalphas = alphas[1] - alphas[0]
taus = tau(alphas)
rs = r(alphas)
Dtaus = Dtau(alphas)
Drs = Dr(alphas)
urs = ur(alphas)
utaus = utau(alphas)

print("[ DEBUG ] Data loaded.")
####################################
# COMPUTE FIELDS ON FREEZOUT SURFACE
####################################
print("[ DEBUG ] Define field functions on FO-surface...")

GeVtoIfm = 5.0677# 1GeV=5.0677(fm^-1)\[hbar]c
IfmtoGeV = 1/GeVtoIfm
fmtoIGeV = GeVtoIfm
IGeVtofm = 1/fmtoIGeV
msigma = 0.6;# mass of the sigma field ~ 400-800 MeV in GeV
epsilon = 0.01*0.160054 # energy density at freeze out in GeV/(fm^3)
vev = 0.15 # ??? vev of sigma field ??? pion decay constant ???, in \GeV
chi2 = (-msigma**2 + 2* msigma * np.sqrt(3*epsilon/vev**2 * IfmtoGeV**3 + msigma**2))/6 # in GeV^2
chi = np.sqrt(chi2) # in GeV

rho = np.sqrt(vev*(2*chi2+msigma**2)/(msigma**2))
thetas = chi * np.cumsum((Dtaus*utaus - Drs*urs) * fmtoIGeV * dalphas)

# We need function that return the pion and sigma field on the freezout surface. We could
#   1)  Define an interpolated function theta(alpha) from the integrated (np.cumsum) data 
#       above and define function pi0(alpha) and sigma0(alpha), referring back to theta(alpha).
#   2)  Compute arrays pi0s and sigma0s from the integrated (np.cumsum) data and interpolate
#       these to find functions pi0(alpha) and sigma0(alpha).

# following option 1):
theta = scipy.interpolate.CubicSpline(alphas, thetas)
pi0 = lambda alpha: np.sqrt(2) * rho * np.cos(theta(alpha))
sigma0 = lambda alpha: np.sqrt(2) * rho * np.sin(theta(alpha))

# following option 2)
# phi0s = rho * np.exp(1j*thetas)
# pi0s = np.sqrt(2) * np.imag(phi0s)
# sigma0s = np.sqrt(2) * np.real(phi0s)

# pi0 = scipy.interpolate.CubicSpline(alphas,pi0s)
# sigma0 = scipy.interpolate.CubicSpline(alphas,sigma0s)

alphatest = np.linspace(0,np.pi/2,10000)
plt.plot(alphatest,pi0(alphatest),label=r"$\pi_0$")
plt.plot(alphatest,sigma0(alphatest),label=r"$\sigma_0$")
plt.legend()
plt.show()

print("[ DEBUG ] Field functions done.")
####################################
# COMPUTE FOURIER SPECTRUM OF SOURCE
####################################
print("[ DEBUG ] Define Bessel and helper functions...")

mpi = 0.14 # in  GeV, m\[Pi]0 = 134.9768 MeV, m\[Pi] + -= 139.57039  MeV 
def omega(p):
    return np.sqrt(p**2 + mpi**2)

####################################
# Bessel functions
def J0rp(alpha, p):
    return scipy.special.j0(r(alpha)*p*GeVtoIfm)
def J1rp(alpha, p):
    return scipy.special.j1(r(alpha)*p*GeVtoIfm)
def Y0tw(alpha, p):
    return scipy.special.y0(tau(alpha)*omega(p)*GeVtoIfm)
def Y1tw(alpha, p):
    return scipy.special.y1(tau(alpha)*omega(p)*GeVtoIfm)
def J0tw(alpha, p):
    return scipy.special.j0(tau(alpha)*omega(p)*GeVtoIfm)
def J1tw(alpha, p):
    return scipy.special.j1(tau(alpha)*omega(p)*GeVtoIfm)

####################################
# Helper functions
def H1(alpha, p):
    return J0rp(alpha, p) * (-Y0tw(alpha, p) + 1j*J0tw(alpha,p))
def H2(alpha, p):
    return Dtau(alpha) * fmtoIGeV * p * J1rp(alpha,p) * (-Y0tw(alpha, p) + 1j*J0tw(alpha,p)) +\
        Dr(alpha) * fmtoIGeV * J0rp(alpha,p) * omega(p) * (-Y1tw(alpha,p) + 1j*J1tw(alpha,p))

print("[ DEBUG ] Bessel and helper functions done.")

####################################
# full fourier trafo of soruce
def ComputeSourceSpectr(p, pi, Dpi):
    # f = lambda alpha: tau(alpha) * r(alpha) * (fmtoIGeV**2) * (Dpi(alpha) * H1(alpha,p) + pi(alpha) * H2(alpha,p))
    # plt.plot(np.linspace(0,np.pi/2,100), f(np.linspace(0,np.pi/2,100)))
    # plt.show()
    integrand = lambda alpha: tau(alpha) * r(alpha) * (fmtoIGeV**2) * (Dpi(alpha) * H1(alpha,p) + pi(alpha) * H2(alpha,p))
    val, err = scipy.integrate.quad(
        integrand,
        0, np.pi/2, limit=200
        )
    return 2 * (np.pi**2) * val

ps = np.linspace(0,2,50)
myspectr = np.array([
    ComputeSourceSpectr(p, pi0, lambda alpha: chi * sigma0(alpha) * (-Dr(alpha)*utau(alpha) + Dtau(alpha)*ur(alpha)))
    for p in ps
])

plt.scatter(ps,np.abs(myspectr)**2)
plt.yscale("log")
plt.show()

# licserv5.rz.tu-ilmenau.de
