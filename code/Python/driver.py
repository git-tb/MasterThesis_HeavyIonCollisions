# %%

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

v0 = 1.5
mu2 = 2
lam = mu2/(v0**2)
eps = 0.05

# %%
# =======================================

def driver(tau, u, r, nghost=1):
    # unpack state vector
    sigma, pi, Psigma, Ppi = u

    # calculate all necessary derivatives
    dr_sigma = gr.ghost_pad(op.D1(sigma, r, nghost=nghost), nghost=nghost)
    ddr_sigma = gr.ghost_pad(op.D2(sigma, r, nghost=nghost), nghost=nghost)
    dr_pi = gr.ghost_pad(op.D1(pi, r, nghost=nghost), nghost=nghost)
    ddr_pi = gr.ghost_pad(op.D2(pi, r, nghost=nghost), nghost=nghost)

    fc.BC_parity(dr_sigma, is_even=False, left_side=True, nghost=nghost)
    fc.BC_parity(ddr_sigma, is_even=True, left_side=True, nghost=nghost)
    fc.BC_parity(dr_pi, is_even=False, left_side=True, nghost=nghost)
    fc.BC_parity(ddr_pi, is_even=True, left_side=True, nghost=nghost)

    # impose extrapolation (outflow) conditions
    fc.BC_outflow(dr_sigma, nghost=nghost, order=1)
    fc.BC_outflow(ddr_sigma, nghost=nghost, order=1)
    fc.BC_outflow(dr_pi, nghost=nghost, order=1)
    fc.BC_outflow(ddr_pi, nghost=nghost, order=1)

    # time evolution
    dsigma = Psigma
    dpi = Ppi
    dPsigma = -1/tau * Psigma + ddr_sigma + 1/r * dr_sigma  + (mu2 - lam* (pi**2 + sigma**2)) * sigma + eps
    dPpi = -1/tau * Ppi + ddr_pi + 1/r * dr_pi  + (mu2 - lam * (pi**2 + sigma**2)) * pi

    # impose parity conditions
    fc.BC_parity(dsigma, is_even=True, left_side=True, nghost=nghost)
    fc.BC_parity(dpi, is_even=True, left_side=True, nghost=nghost)
    fc.BC_parity(dPsigma, is_even=True, left_side=True, nghost=nghost)
    fc.BC_parity(dPpi, is_even=True, left_side=True, nghost=nghost)

    # impose extrapolation (outflow) conditions
    fc.BC_outflow(dsigma, nghost=nghost, order=1)
    fc.BC_outflow(dpi, nghost=nghost, order=1)
    fc.BC_outflow(dPsigma, nghost=nghost, order=1)
    fc.BC_outflow(dPpi, nghost=nghost, order=1)

    # return update step
    res = (
        dsigma,
        dpi,
        dPsigma,
        dPpi
    )

    return np.vstack(res)


# %%
# =======================================
# grid parameters
nghost = 2
Nr = 500
r_a, r_b = 0, 30

r = gr.gr_CC(Nr, r_a, r_b, nghost=nghost)
dr = r[1] - r[0]

tau_i, tau_f = 0.01, 50
CFL = 0.1
N_t = it.N_t_from_CFL(CFL, r[1] - r[0], tau_i, tau_f)
print(N_t)
dtau = (tau_f - tau_i) / N_t

# =======================================
# initial conditions

sigma0 = (v0+eps/(2*mu2)) * np.ones_like(r)
pi0 = 0.001*(r**2)*np.exp(-(r-10)**2)
Psigma0 = np.zeros_like(r)
Ppi0 = np.zeros_like(r)

# impose parity conditions
fc.BC_parity(sigma0, is_even=True, left_side=True, nghost=nghost)
fc.BC_parity(pi0, is_even=True, left_side=True, nghost=nghost)
fc.BC_parity(Psigma0, is_even=True, left_side=True, nghost=nghost)
fc.BC_parity(Ppi0, is_even=True, left_side=True, nghost=nghost)

# impose extrapolation (outflow) conditions
fc.BC_outflow(sigma0, nghost=nghost, order=1)
fc.BC_outflow(pi0, nghost=nghost, order=1)
fc.BC_outflow(Psigma0, nghost=nghost, order=1)
fc.BC_outflow(Ppi0, nghost=nghost, order=1)

# =======================================
# set up state vector
u_i = sigma0, pi0, Psigma0, Ppi0


# %%

driver_func = ft.partial(driver, r=r, nghost=nghost)
# =======================================
# evolve


u_c = u_i
hist = [u_i]

tau = tau_i
taus = [tau]
for t_idx in range(int(N_t)):
    print("tau=", tau)
    
    u_c = it.ev_step_RK4(tau, u_c, dtau, driver_func)

    tau = tau + dtau
    taus.append(tau)

    # sigma, pi, Psigma, Ppi = u_c
    hist.append(u_c)
    # print(u_c)

hist = np.array(hist)
sigmahist, pihist, Psigmahist, Ppihist = np.transpose(hist,(1,0,2))


# %%

skip = 5
taus = taus[::skip]
sigmahist, pihist, Psigmahist, Ppihist = sigmahist[::skip], pihist[::skip], Psigmahist[::skip], Ppihist[::skip]


%matplotlib qt
plt.style.use("mplstyles/myclassic_white.mplstyle")

fig = plt.figure(figsize=(19,10))
gs = GridSpec(2,2,figure=fig)

axsigma = fig.add_subplot(gs[0,0])
axpi = fig.add_subplot(gs[1,0])
axPsigma = fig.add_subplot(gs[0,1])
axPpi = fig.add_subplot(gs[1,1])

plotsigma = axsigma.plot(r,sigmahist[0], marker="o",markersize=2)[0]
plotpi = axpi.plot(r,pihist[0], marker="o",markersize=2)[0]
plotPsigma = axPsigma.plot(r,Psigmahist[0], marker="o",markersize=2)[0]
plotPpi = axPpi.plot(r,Ppihist[0], marker="o",markersize=2)[0]

axsigma.set_title(r"$\sigma$")
axpi.set_title(r"$\pi$")
axPsigma.set_title(r"$\Pi\sigma$")
axPpi.set_title(r"$\Pi\pi$")

def getlims(a):
    return np.min(a)-0.05*np.ptp(a), np.max(a)+0.05*np.ptp(a)

axsigma.set_ylim(getlims(sigmahist))
axpi.set_ylim(getlims(pihist))
axPsigma.set_ylim(getlims(Psigmahist))
axPpi.set_ylim(getlims(Ppihist))

def update(frame):
    fig.suptitle(r"$tau=$"+str(taus[frame]))

    plotsigma.set_ydata(sigmahist[frame])
    plotpi.set_ydata(pihist[frame])
    plotPsigma.set_ydata(Psigmahist[frame])
    plotPpi.set_ydata(Ppihist[frame])

    # axsigma.relim()
    # axpi.relim()
    # axPsigma.relim()
    # axPpi.relim()

    # axsigma.autoscale_view()
    # axpi.autoscale_view()
    # axPsigma.autoscale_view()
    # axPpi.autoscale_view()

    # axsigma.set_ylim(-1,1)


ani = AnimPlayer.Player(fig=fig, func=update,maxi=len(taus)-1)
plt.show()


# %%

sigmaavg, piavg = np.mean(sigmahist, axis=1), np.mean(pihist, axis=1)
plt.plot(taus,sigmaavg)
plt.plot(taus,piavg)
plt.hlines(v0,np.min(taus),np.max(taus))
# plt.plot(sigmaavg,piavg)
plt.show()
