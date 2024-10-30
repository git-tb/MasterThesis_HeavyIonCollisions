#%%

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import scipy
import scipy.integrate

from IPython import get_ipython
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap

import datetime
import os


get_ipython().run_line_magic("matplotlib","qt")
plt.style.use("mplstyles/myclassic_white.mplstyle")

#%%

df_freezeout = pd.read_csv("./../../Mathematica/data/ExampleFreezeOutCorrected.csv")

alphas = df_freezeout["alpha"].to_numpy()
taus = df_freezeout["tau"].to_numpy()
rs = df_freezeout["r"].to_numpy()
Dtaus = df_freezeout["Dtau"].to_numpy()
Drs = df_freezeout["Dr"].to_numpy()
utaus = df_freezeout["utau"].to_numpy()
urs = df_freezeout["ur"].to_numpy()

# plt.plot(alphas,taus)
# plt.plot(alphas,rs)
# plt.plot(alphas,Dtau)
# plt.plot(alphas,Dr)
# plt.plot(alphas,utau)
# plt.plot(alphas,ur)
# plt.show()

tau = scipy.interpolate.CubicSpline(alphas,taus)
r = scipy.interpolate.CubicSpline(alphas,rs)
Dtau = scipy.interpolate.CubicSpline(alphas,Dtaus)
Dr = scipy.interpolate.CubicSpline(alphas,Drs)
utau = scipy.interpolate.CubicSpline(alphas,utaus)
ur = scipy.interpolate.CubicSpline(alphas,urs)

#%%
# ==========================================
# ============= CONSTANT FIELD =============
# ==========================================

# m = 0.14
Rs = np.linspace(0,1,11,endpoint=True)
epsmax = 0.160054 * 1e-3

sign = -1
ms = [0.28]
for m in ms:
    phimax = np.sqrt(2*epsmax)/m

    timestamp = "{:%Y%m%d_%H%M%S}".format(datetime.datetime.now())
    newpath = "data/init_real_m%d_s%d_constfield_"%(1000*m,sign)+timestamp
    os.makedirs(newpath)

    for R in Rs:
        phi = R * phimax
        dphi = sign*(1-R) * phimax

        ## SAVE
        outfilename = "R%d"%(10*R)
        with open(newpath+"/"+outfilename+".csv","w") as outfile:
            outfile.write("# phimax: "+str(phimax)+"\n")
            outfile.write("# mass: "+str(m)+"\n")
            outfile.write("alpha,RePhi,ImPhi,ReDPhi,ImDPhi"+"\n")
            outfile.write("0,"+str(phi)+",0,"+str(dphi)+",0"+"\n")
            outfile.write("1.570797,"+str(phi)+",0,"+str(dphi)+",0"+"\n")

#%%
# ==========================================
# ======== CONSTANT ENERGY DENSITY =========
# ==========================================

epsilon = 0.160054 * 1e-3

def uds(t):
    return -utau(t) * Dtau(t) + ur(t) * Dr(t) # utau === u^tau = - u_tau

sign = -1
ms = [0.28]
for m in ms:
    timestamp = "{:%Y%m%d_%H%M%S}".format(datetime.datetime.now())
    newpath = "data/init_real_m%d_s%d_consteps_"%(1000*m, sign)+timestamp
    os.makedirs(newpath)

    def f(y, t, *args):
        phi, chi = y[0], y[1]

        dphi = chi * uds(t)
        dchi = -m**2 * phi * uds(t)

        return (dphi, dchi)

    Rs = np.linspace(0,1,11,endpoint=True)
    phis = []
    dphis = []
    lines_init1 = []
    lines_init2 = []

    ts = np.linspace(0,np.pi/2,100,endpoint=True)
    for R in Rs:
        epot0 = R * epsilon
        ekin0 = epsilon - epot0

        phi0 = np.sqrt(2*epot0)/m
        chi0 = sign * np.sqrt(2*ekin0)

        y0 = (phi0, chi0)

        sol = scipy.integrate.odeint(f,y0,ts)

        phi = sol[:,0]
        chi = sol[:,1]
        dphi = chi*(-Dr(ts)*utau(ts) + Dtau(ts)*ur(ts))

        phis.append(phi)
        dphis.append(dphi)

        line_init1 = np.column_stack((ts, m*phi))
        line_init2 = np.column_stack((ts, dphi))

        lines_init1.append(line_init1)
        lines_init2.append(line_init2)

        ### SAVE

        timestamp = "{:%Y%m%d_%H%M%S}".format(datetime.datetime.now())
        outfilename = "R%d"%(10*R)
        with open(newpath+"/"+outfilename+".csv","w") as outfile:
            outfile.write("# epsilon: "+str(epsilon)+"\n")
            outfile.write("# mass: "+str(m)+"\n")
            outfile.write("alpha,RePhi,ImPhi,ReDPhi,ImDPhi"+"\n")

            for (n,t) in enumerate(ts):
                outfile.write(str(t)+","+str(phi[n])+",0,"+str(dphi[n])+",0"+"\n")

    fig, (ax1, ax2) = plt.subplots(figsize=(14,7),ncols=2)

    CMAP = LinearSegmentedColormap.from_list("custom", ["blue","red"])
    CMAP_LBWH = [0.025, 0.025, 0.05, 0.45]

    lc_init1 = LineCollection(lines_init1,array=Rs,cmap=CMAP)
    lc_init2 = LineCollection(lines_init2,array=Rs,cmap=CMAP)

    ax1.add_collection(lc_init1)
    ax2.add_collection(lc_init2)

    ax1.autoscale_view()
    ax2.autoscale_view()

    cax = ax1.inset_axes(CMAP_LBWH)
    cbar = fig.colorbar(lc_init1, cax=cax)

    plt.show()
    plt.close()

#%%
# ==========================================
# ===== TAU-INCREASING FIELD AMPLITUDE =====
# ==========================================

taumax = np.max(taus)
taumin = np.min(taus)
taumax = np.max(taus)

epsmax = 0.160054 * 1e-3

exponents = np.concatenate((1/np.arange(3,1,-1),np.arange(1,4,1)))
alphas = np.linspace(0,np.pi/2,100,endpoint=True)

phis = []
dphis = []
lines_init1 = []
lines_init2 = []
ms = [0.28]
for m in ms:
    timestamp = "{:%Y%m%d_%H%M%S}".format(datetime.datetime.now())
    newpath = "data/init_real_m%d_taudep_"%(1000*m)+timestamp
    os.makedirs(newpath)

    phimax = np.sqrt(2*epsmax)/m

    for (nvar,a) in enumerate(exponents):
        def phi_restr(alpha):
            return phimax * ((tau(alpha)-taumin)/(taumax-taumin))**a

        def dphi_restr(alpha):
            return Dr(alpha) * phimax * a * ((tau(alpha)-taumin)/(taumax-taumin))**(a-1) * 1/(taumax-taumin)
        
        phi = phi_restr(alphas)
        dphi = dphi_restr(alphas)

        line_init1 = np.column_stack((alphas, phi))
        line_init2 = np.column_stack((alphas, dphi))

        lines_init1.append(line_init1)
        lines_init2.append(line_init2)

        ### SAVE

        outfilename = "v%d"%(nvar)
        with open(newpath+"/"+outfilename+".csv","w") as outfile:
            outfile.write("# epsilon: "+str(epsilon)+"\n")
            outfile.write("# mass: "+str(m)+"\n")
            outfile.write("# exponent: "+str(a)+"\n")
            outfile.write("alpha,RePhi,ImPhi,ReDPhi,ImDPhi"+"\n")

            for (n,alpha) in enumerate(alphas):
                outfile.write(str(alpha)+","+str(phi[n])+",0,"+str(dphi[n])+",0"+"\n")

    fig, (ax1, ax2) = plt.subplots(figsize=(14,7),ncols=2)

    CMAP = LinearSegmentedColormap.from_list("custom", ["blue","red"])
    CMAP_LBWH = [0.025, 0.025, 0.05, 0.45]

    lc_init1 = LineCollection(lines_init1,array=np.arange(len(exponents)))
    lc_init2 = LineCollection(lines_init2,array=np.arange(len(exponents)))

    ax1.add_collection(lc_init1)
    ax2.add_collection(lc_init2)

    ax1.autoscale_view()
    ax2.autoscale_view()

    plt.show()
    # plt.close()