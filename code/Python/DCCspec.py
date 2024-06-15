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

alpha_arr = df["alpha"].to_numpy()
tau_arr = df["tau"].to_numpy()
r_arr = df["r"].to_numpy()
Dtau_arr = df["Dtau"].to_numpy()
Dr_arr = df["Dr"].to_numpy()
ur_arr = df["ur"].to_numpy()
utau_arr = df["utau"].to_numpy()

tau = scipy.interpolate.CubicSpline(alpha_arr, tau_arr,extrapolate=True)
r = scipy.interpolate.CubicSpline(alpha_arr, r_arr,extrapolate=True)
Dtau = scipy.interpolate.CubicSpline(alpha_arr, Dtau_arr,extrapolate=True)
Dr = scipy.interpolate.CubicSpline(alpha_arr, Dr_arr,extrapolate=True)
ur = scipy.interpolate.CubicSpline(alpha_arr, ur_arr,extrapolate=True)
utau = scipy.interpolate.CubicSpline(alpha_arr, utau_arr,extrapolate=True)

alphatest = np.linspace(0,np.pi/2,50)
# plt.plot(alphatest,r(alphatest),c="b")
# plt.scatter(alpha_arr[::],r_arr[::],c="r",marker="+")
# plt.show()

print("[ DEBUG ] Data loaded.")
####################################
# COMPUTE FIELD ON FREEZOUT SURFACE
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
def theta(alpha):
    return chi * scipy.integrate.quad(
        lambda s: (Dtau(s)*utau(s) - Dr(s)*ur(s)) * fmtoIGeV, 
        0, alpha
    )[0]
def phifreeze(alpha):
    return rho*np.exp(1j*theta(alpha))
def pi0(alpha):
    return np.sqrt(2)*np.imag(phifreeze(alpha))
def sigma0(alpha):
    return np.sqrt(2)*np.real(phifreeze(alpha))

# plt.scatter(alphatest,[theta(a) for a in alphatest])
# plt.show()

print("[ DEBUG ] Field functions done.")
####################################
# COMPUTE FOURIER SPECTRUM OF SOURCE
####################################
print("[ DEBUG ] Define Bessel and helper functions...")

mpi = 0.14 # in  GeV, m\[Pi]0 = 134.9768 MeV, m\[Pi] + -= 139.57039  MeV 
def omega(p):
    return np.sqrt(p**2 + mpi**2)

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

# Helper functions
def H1(alpha, p):
    return J0rp(alpha, p) * (-Y0tw(alpha, p) + 1j*J0tw(alpha,p))
def H2(alpha, p):
    return Dtau(alpha) * fmtoIGeV * p * J1rp(alpha,p) * (-Y0tw(alpha, p) + 1j*J0tw(alpha,p)) +\
        Dr(alpha) * fmtoIGeV * J0rp(alpha,p) * omega(alpha,p) * (-Y1tw(alpha,p) + 1j*J1tw(alpha,p))

print("[ DEBUG ] Bessel and helper functions done.")

#### introduce momentum grid
print("[ DEBUG ] Compute contributions to spectrum on fixed momentum grid...")
ps = np.linspace(0.01,2,2)

# H1ps = functools.partial(H1,p=ps)
# H2ps = functools.partial(H2,p=ps)

# "C(1/2)(s/p)ps = contribution of (1st/2nd) helper function with (sigma/pion) field integrated over FO surface, evaluated on ps"
C1sps = np.array([
    2*np.pi**2 * scipy.integrate.quad( 
    lambda alpha: 
        tau(alpha)*r(alpha) * fmtoIGeV**2 * (-1) * pi0(alpha) * chi * (-Dr(alpha)*utau(alpha) + Dtau(alpha)*ur(alpha)) * H1(alpha,p),
    0, np.pi/2
) 
for p in ps])
print("[ DEBUG ] C1sps done.")
C1pps = np.array([
    2*np.pi**2 * scipy.quad( 
    lambda alpha: 
        tau(alpha)*r(alpha) * fmtoIGeV**2 * sigma0(alpha) * chi * (-Dr(alpha)*utau(alpha) + Dtau(alpha)*ur(alpha)) * H1(alpha,p),
    0, np.pi/2
)
for p in ps])
print("[ DEBUG ] C1pps done.") 
C2sps = np.array([
    2*np.pi**2 * scipy.integrate.quad( 
    lambda alpha: 
        tau(alpha)*r(alpha) * fmtoIGeV**2 * sigma0(alpha) *H2ps(alpha),
    0, np.pi/2
) 
for p in ps])
print("[ DEBUG ] C2sps done.")
C2pps = np.array([
    2*np.pi**2 * scipy.integrate.quad( 
    lambda alpha: 
        tau(alpha)*r(alpha) * fmtoIGeV**2 * pi0(alpha) *H2ps(alpha),
    0, np.pi/2
)
for p in ps])
print("[ DEBUG ] C2pps done.")
            

