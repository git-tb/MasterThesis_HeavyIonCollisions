#%%
# ===================== LOAD LIBRARIES ======================
# ===========================================================

import numpy as np
import grids as gr
import field_conditions as fc
import integrators as it
from matplotlib import pyplot as plt
import operators as op
import functools as ft
import AnimPlayer
import matplotlib.animation as animation
from matplotlib.gridspec import GridSpec

plt.style.use("mplstyles/myclassic_white.mplstyle")

#%%
# =============== DEFINE PHYSICAL PARAMETERS ================
# ===========================================================

%matplotlib inline

fpi = 94.5 # MeV
mpi = 140 # MeV
msigma = 600 # MeV

# fpi = 2 # MeV
# mpi = 0.05*2 # MeV
# msigma = 1 # MeV

lam = (msigma**2/2 + mpi**2)/(fpi**2)

def V(sigma, pi1, pi2, pi3,T):
    rho2 = (sigma**2 + pi1**2 + pi2**2 + pi3**3)
    return lam/4 * (rho2**2) -lam/2 * rho2 * (fpi**2 -mpi**2/lam - T**2/2) - fpi*mpi**2 * sigma


### TEMPERATURE DEPENDENCE OF POTENTIAL
fig, (axsig, axpi) = plt.subplots(ncols=2,figsize=(20,10))
Ts = np.linspace(0,200,5)
xs = np.linspace(0,200,1000)

scale = 1e-9

for thisT in Ts:
    axsig.plot(xs,scale*V(xs,0,0,0,thisT),
        label=r"$T={%.0f}\ \mathrm{MeV}$"%(thisT),
        marker="")
    axpi.plot(xs,scale*V(0,xs,0,0,thisT),
        label=r"$T={%.0f}\ \mathrm{MeV}$"%(thisT),
        marker="")

axsig.set_ylim(-0.7,1.5)
axpi.set_ylim(-0.7,1.5)

axsig.set_ylabel(r"$V_\mathrm{pot}\ [\mathrm{GeV}]$")
axpi.set_ylabel(r"$V_\mathrm{pot}\ [\mathrm{GeV}]$")

axsig.legend()
axpi.legend()
plt.show()

### TEMPERATURE AS A FUNCTION OF TAU
T0 = 200

def T(tau):
    return 0

taugrid = np.linspace(0,10,100)
Tgrid = np.array([T(tau) for tau in taugrid])

fig, ax = plt.subplots()
plt.plot(taugrid, Tgrid)

# %%
# ================== DEFINE ODE SYSTEM ======================
# ===========================================================

def driver(tau, u):
    # unpack state vector
    sigma, pi0, pi1, pi2, Psigma, Ppi0, Ppi1, Ppi2 = u

    fields = (sigma, pi0, pi1, pi2)
    Pfields = (Psigma, Ppi0, Ppi1, Ppi2)

    rho2 = np.sum(np.array([f**2 for f in fields]), axis=0)

    # calculate all necessary derivatives
    Dfields, DPfields = [], []

    for (idx, (field, Pfield)) in enumerate(list(zip(fields, Pfields))):
        # time step
        Dfield = Pfield
        DPfield = -1/tau * Pfield + lam * (fpi**2 -mpi**2/lam - T(tau)**2/2 - rho2) * field + (idx == 0) * fpi*mpi**2

        Dfields.append(Dfield)
        DPfields.append(DPfield)

    # return update step
    res = (
        *Dfields,
        *DPfields
    )
    return np.vstack(res)


# %%
# ==================== INITIALIZATION =======================
# ===========================================================

tau0 = 0.1
tauf = 2

dtau = 0.0005
N_t = (tauf-tau0)/dtau
N_save = np.clip(100000,0,N_t)
skip_save = N_t//N_save
print("skip %.d before save"%(skip_save-1))
dtau = (tauf - tau0) / N_t

# =======================================
# initial conditions
sigma = -10
pi0 = 10
pi1 = 0
pi2 = 0
Psigma = 0
Ppi0 = 0
Ppi1 = 0
Ppi2 = 0

fields = (sigma, pi0, pi1, pi2)
Pfields = (Psigma, Ppi0, Ppi1, Ppi2)

# =======================================
# set up state vector
u_i = np.vstack((*fields, *Pfields))


# %%
# ================== PERFORM EVOLUTION ======================
# ===========================================================

driver_func = driver
# =======================================
# evolve
u_c = u_i
hist = [u_i]

tau = tau0
taus = [tau]
for t_idx in range(int(N_t)):
    print("tau=%.3f (%.d/%.d)"%(tau,t_idx,N_t))
    
    u_c = it.ev_step_RK4(tau, u_c, dtau, driver_func)

    tau = tau + dtau
    if(t_idx % skip_save == 0):
        taus.append(tau)
        hist.append(u_c)

hist = np.array(hist)

# %%
# ===================== PLOT THE DATA =======================
# ===========================================================
%matplotlib qt

sigmahist, pi0hist, pi1hist, pi2hist, Psigmahist, Ppi0hist, Ppi1hist, Ppi2hist = np.transpose(hist,(2,1,0))[0]

fieldhists = (sigmahist, pi0hist, pi1hist, pi2hist)
Pfieldhists = (Psigmahist, Ppi0hist, Ppi1hist, Ppi2hist)

fig = plt.figure(figsize=(19,10))
gs = GridSpec(2,4,figure=fig)

labels = [r"$\sigma$",r"$\pi^0$",r"$\pi^1$",r"$\pi^2$"]

for (idx, _) in enumerate(fieldhists):
    axfield = fig.add_subplot(gs[0,idx])
    axPfield = fig.add_subplot(gs[1,idx])

    axfield.plot(taus, fieldhists[idx])
    axPfield.plot(taus, Pfieldhists[idx])

    axfield.set_title(labels[idx])
    axPfield.set_title(r"$\Pi$"+labels[idx])

plt.show()

#%%
# =================== IN SIGMA-PI SPACE =====================
# ===========================================================

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

x = sigmahist
# y = np.sign(pi0hist*pi1hist*pi2hist)*np.sqrt(pi0hist**2+pi1hist**2+pi2hist**2)
y = pi0hist

points = np.array([x, y]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

fig, ax = plt.subplots()
lc = LineCollection(segments, cmap='jet')
lc.set_array(taus)
lc.set_linewidth(2)
line = ax.add_collection(lc)
fig.colorbar(line,ax=ax)

ax.set_xlim(np.min(x),np.max(x))
ax.set_ylim(np.min(y),np.max(y))
plt.show()


# %%
# =================== FOURIER ANALYSIS ======================
# ===========================================================
%matplotlib qt

fig = plt.figure(figsize=(20,7))
gs = GridSpec(1,4,figure=fig)

labels = [r"$\sigma$",r"$\pi^0$",r"$\pi^1$",r"$\pi^2$"]

for (idx, _) in enumerate(fieldhists):
    axfield = fig.add_subplot(gs[0,idx])

    spec = np.fft.fft(fieldhists[idx])
    freqs = 2*np.pi * np.fft.fftfreq(len(taus),(tauf-tau0)/N_t)

    spec = np.fft.fftshift(spec)
    freqs = np.fft.fftshift(freqs)

    axfield.plot(freqs, spec)

    axfield.set_title(labels[idx])

plt.show()