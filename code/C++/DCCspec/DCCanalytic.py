#%%

import numpy as np
from matplotlib import pyplot as plt
import scipy
import scipy.special as spc
import pandas as pd
from matplotlib.widgets import Button, Slider
import glob

get_ipython().run_line_magic("matplotlib","qt")
plt.style.use("mplstyles/myclassic_white.mplstyle")

#%%

GeVtoIfm = 5.0677
fmtoIGeV = GeVtoIfm
m = 0.14
tau0 = 0.4

def w(m,p):
    return np.sqrt(m**2 + p**2)

def JT(p,Tau,R, A,B):
    return 2 * (np.pi**2) * Tau * R * (fmtoIGeV)**2 / p * (
        B * (-spc.y0(Tau * w(m,p) * GeVtoIfm) + 1j * spc.j0(Tau * w(m,p) * GeVtoIfm)) + 
        A * (-spc.y1(Tau * w(m,p) * GeVtoIfm) + 1j * spc.j1(Tau * w(m,p) * GeVtoIfm)) * w(m,p)
    ) * spc.j1(R * p * GeVtoIfm)

def JR(p,Tau,R, A,B):
    return -2 * (np.pi**2) * R * (fmtoIGeV)**2 / w(m,p) * (
        B * spc.j0(R * p * GeVtoIfm) + 
        A * spc.j1(R * p * GeVtoIfm) * p
    ) * (
        Tau * (-spc.y1(Tau * w(m,p) * GeVtoIfm) + 1j * spc.j1(Tau * w(m,p) * GeVtoIfm)) -
        tau0 * (-spc.y1(tau0 * w(m,p) * GeVtoIfm) + 1j * spc.j1(tau0 * w(m,p) * GeVtoIfm))
    )

def specan(p, Tau, R, AT, BT, AR, BR):
    return 1/(2*np.pi)**3 * np.abs(JT(p, Tau, R, AT, BT) + JR(p, Tau, R, AR, BR))**2

#%%

# spectra = glob.glob("data/spectra_real_m140_constfield_20241029_110519/*/")
spectra = glob.glob("data/spectra_real_m140_consteps_20241029_110232/*/")
# spectra = glob.glob("data/spectra_real_m140_taudep_20241029_132506/*/")
spec = spectra[0]

df_init1 = pd.read_csv(spec+"/field0.txt",comment="#")
df_init2 = pd.read_csv(spec+"/field0_deriv.txt",comment="#")

fig, ax = plt.subplots()
ax.plot(df_init1["alpha"].to_numpy(),df_init1["field0Re"].to_numpy(),marker="",label="field")
ax.plot(df_init2["alpha"].to_numpy(),df_init2["Dfield0Re"].to_numpy(),marker="",label="deriv")
ax.legend()
plt.show()

df_spec = pd.read_csv(spec+"/spectr.txt", comment="#")

T0, R0 = 12.2, 8.5

Amax = 0.2
Bmax = 4/10 # value of (n^\mu d_\mu phi)/(typical value of r', tau')

ps = np.linspace(0,2,500)
myspec = specan(ps,T0, R0,0.1,0.1,0.1,0.1)

fig, ax = plt.subplots()
figsl = plt.figure()
ax.plot(df_spec["pT"].to_numpy(),df_spec["abs2Re"].to_numpy())

line, = ax.plot(ps,myspec)
ax.set_yscale("log")

ax1 = figsl.add_axes([0.1, 0.1, 0.03, 0.8])
sl1 = Slider(
    ax=ax1,
    valmin=0,
    valmax=20,
    label="Tau",
    valinit=T0,
    orientation="vertical"
)

ax2 = figsl.add_axes([0.2, 0.1, 0.03, 0.8])
sl2 = Slider(
    ax=ax2,
    valmin=0,
    valmax=10,
    label="R",
    valinit=R0,
    orientation="vertical"
)

ax3 = figsl.add_axes([0.3, 0.1, 0.03, 0.8])
sl3 = Slider(
    ax=ax3,
    valmin=-5,
    valmax=5,
    label="scale",
    valinit=0,
    orientation="vertical"
)

ax4 = figsl.add_axes([0.4, 0.1, 0.03, 0.8])
sl4 = Slider(
    ax=ax4,
    valmin=-Amax,
    valmax=Amax,
    label="A1",
    valinit=1,
    orientation="vertical"
)

ax5 = figsl.add_axes([0.5, 0.1, 0.03, 0.8])
sl5 = Slider(
    ax=ax5,
    valmin=-Bmax,
    valmax=Bmax,
    label="B1",
    valinit=1,
    orientation="vertical"
)

ax6 = figsl.add_axes([0.6, 0.1, 0.03, 0.8])
sl6 = Slider(
    ax=ax6,
    valmin=-Amax,
    valmax=Amax,
    label="A2",
    valinit=1,
    orientation="vertical"
)

ax7 = figsl.add_axes([0.7, 0.1, 0.03, 0.8])
sl7 = Slider(
    ax=ax7,
    valmin=-Bmax,
    valmax=Bmax,
    label="B2",
    valinit=1,
    orientation="vertical"
)

# The function to be called anytime a slider's value changes
def update(val):

    myspec = specan(ps,sl1.val, sl2.val,sl4.val, sl5.val, sl6.val, sl7.val)
    line.set_ydata(np.exp(sl3.val)*myspec)

    fig.canvas.draw_idle()

sl1.on_changed(update)
sl2.on_changed(update)
sl3.on_changed(update)
sl4.on_changed(update)
sl5.on_changed(update)
sl6.on_changed(update)
sl7.on_changed(update)

plt.show()

#%%

###
### COMPUTE A RECTANGULAR APPROXIMATION OF THE FREEZEOUT GEOMETRY
###

from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit
from scipy.integrate import quad
from scipy.optimize import minimize

df = pd.read_csv("./../../Mathematica/data/ExampleFreezeOutCorrected.csv")

alphas = df["alpha"].to_numpy()
taus = df["tau"].to_numpy()
rs = df["r"].to_numpy()

tau = CubicSpline(alphas, taus)
r = CubicSpline(alphas, rs)

tau0 = tau(np.pi/2)

def d(alpha, R, T):
    if(alpha <= np.arctan2(R,T)):
        return T/np.cos(alpha)
    if(np.arctan2(R,T) < alpha):
        return R/np.sin(alpha)
    
def Dd(alpha, R, T):
    if(alpha <= np.arctan2(R,T)):
        return np.tan(alpha) * T/np.cos(alpha)
    if(np.arctan2(R,T) < alpha):
        return -R/(np.sin(alpha)*np.tan(alpha))

def myfit_tau(alpha, R, T):
    return d(alpha, R, T) * np.cos(alpha)

def myfit_Dtau(alpha, R, T):
    # return Dd(alpha, R, T) * np.cos(alpha) - d(alpha, R, T) * np.sin(alpha)
    if(alpha <= np.arctan2(R,T)):
        return 0
    if(np.arctan2(R,T) < alpha):
        # return -R/(np.tan(alpha)**2) - R
        return -R/np.sin(alpha)**2

def myfit_r(alpha, R, T):
    x = d(alpha, R, T)
    return d(alpha, R, T) * np.sin(alpha)

def myfit_Dr(alpha, R, T):
    # return Dd(alpha, R, T) * np.sin(alpha) + d(alpha, R, T) * np.cos(alpha)
    if(alpha <= np.arctan2(R,T)):
        # return T * np.tan(alpha)**2 + T
        return T/np.cos(alpha)**2
    if(np.arctan2(R,T) < alpha):
        return 0

def local_cost(alpha, R, T):
    return (tau(alpha) - tau0 - myfit_tau(alpha, R, T))**2 + (r(alpha) - myfit_r(alpha, R, T))**2

def mycost(x):
    R, T = x
    return quad(local_cost, 0, np.pi/2,args=(R,T))[0]

# FIND OPTIMAL PARAMETERS FOR RECTANGULAR REGION
result = minimize(mycost, (9, 13))
Ropt, Topt = result.x

mytaus = np.array([myfit_tau(a,Ropt, Topt) + tau0 for a in alphas])
myrs = np.array([myfit_r(a,Ropt, Topt) for a in alphas])

myDtaus = np.array([myfit_Dtau(a,Ropt, Topt) for a in alphas])
myDrs = np.array([myfit_Dr(a,Ropt, Topt) for a in alphas])

numDtaus = (mytaus[1:] - mytaus[:-1])/(alphas[1]-alphas[0])
numDrs = (myrs[1:] - myrs[:-1])/(alphas[1]-alphas[0])

dataDtaus = df["Dtau"].to_numpy()
dataDrs = df["Dr"].to_numpy()

datanumDtaus = (taus[1:] - taus[:-1])/(alphas[1]-alphas[0])
datanumDrs = (rs[1:] - rs[:-1])/(alphas[1]-alphas[0])

fig, ax = plt.subplots()
ax.plot(alphas, taus,c="b",marker="",label=r"$\tau$")
ax.plot(alphas, mytaus,c="b",ls="--",marker="")
ax.plot(alphas, rs,c="r",marker="",label=r"$r$")
ax.plot(alphas, myrs,c="r",ls="--",marker="")
ax.legend()
plt.show()

fig, ax = plt.subplots()
ax.plot(alphas[:-1], numDtaus,c="b",marker="",label=r"$D\tau$, finite diff")
ax.plot(alphas, myDtaus,c="b",ls="--",marker="",label=r"$D\tau$, analytic")
ax.plot(alphas[:-1], datanumDtaus,c="b",ls="--",marker="x",label=r"$D\tau$, data finite diff")
ax.plot(alphas, dataDtaus,c="b",ls="-.",marker="",label=r"$D\tau$, data")

ax.plot(alphas[:-1], numDrs,c="r",marker="",label=r"$Dr$, finite diff")
ax.plot(alphas, myDrs,c="r",ls="--",marker="",label=r"$Dr$, analytic")
ax.plot(alphas[:-1], datanumDrs,c="r",ls="--",marker="x",label=r"$Dr$, data finite diff")
ax.plot(alphas, dataDrs,c="r",ls="-.",marker="",label=r"$Dr$, data")
ax.legend()
plt.show()

fig, ax = plt.subplots()
ax.plot(rs,taus,marker="",label="exact FO geometry",c="b")
ax.plot(myrs, mytaus,marker="",label="approximated FO geometry",c="r")
ax.legend()
plt.show()

# AND SAVE

# newdf = pd.DataFrame()
# keys = df.keys()

# newdf["alpha"] = alphas
# newdf["tau"] = mytaus
# newdf["r"] = myrs
# newdf["Dtau"] = myDtaus
# newdf["Dr"] = myDrs
# newdf["ur"] = df["ur"]
# newdf["utau"] = df["utau"]
# newdf.to_csv("./../../../code/Mathematica/data/ExampleFreezeOutRectangular.csv",index=False)

#%%

# CALCULATE CORRESPONDING FIELDS WHERE phi = const AND D_phi = const RESPECTIVELY

A_T = 1 # value of phi where tau = T = const
B_T = -1 # value of d_tau phi where tau = T = const

A_R = 1 # value of phi where tau = T = const
B_R = 1 # value of d_r phi where tau = T = const

def maskT(alpha,R,T):
    # return 1 when we are on the tau = T = const slice, else 0
    if(alpha <= np.arctan2(R,T)):
        return 1
    if(np.arctan2(R,T) < alpha):
        return 0  
    
def maskR(alpha,R,T):
    # return 1 when we are on the r = R = const slice, else 0
    return 1 - maskT(alpha, R, T)

field = np.array([maskT(a, Ropt, Topt) * A_T + maskR(a, Ropt, Topt) * A_R for a in alphas])
fieldderiv = np.array(myDtaus * B_T + myDrs * B_R)

plt.plot(fieldderiv)
plt.show()





