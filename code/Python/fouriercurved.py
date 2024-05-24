import numpy as np
import matplotlib.pyplot as plt

Npoints = 256
L = 10
zs = np.linspace(0,L,Npoints) # somehow the left bound has to be 0???
dz = zs[1]-zs[0]

signal = np.exp(1j*3*2*np.pi/L*zs)
# plt.plot(zs,signal)
# plt.show()

spect = np.fft.fftshift(np.fft.fft(signal))
ps = 2*np.pi*np.fft.fftshift(np.fft.fftfreq(Npoints,dz)) # fftfreq gives ps as in np.exp(1j*2pi*p*x)
# plt.plot(ps,spect.real,c="b")
# plt.plot(ps,spect.imag,c="g",ls="-.")
# plt.show()


m = 1
def omega(p):
    return np.sqrt(m**2 + p**2)

# def phi(z,t):    
#     return np.sum(spect/Npoints * np.exp(1j*(ps*z)))

for item in list(zip(ps,omega(ps))):
    print(item)

def phi(t,zs,ps,coeffs):
    tw = np.outer(t*np.ones_like(zs),omega(ps))
    pz = np.outer(zs,ps)
    expmat = np.exp(1j*(pz-tw))
    return np.matmul(expmat,coeffs/Npoints)

# plt.plot(zs,phis(0,zs,ps,spect).real,c="g",ls="--")
# plt.plot(zs,np.array([phi(z,0) for z in zs]).real,c="r",ls="-.")
# plt.plot(zs,signal.real,c="b",ls=":")
# plt.show()

from matplotlib import animation
from IPython.display import HTML

tsteps = 100
tmax = 5
ts = np.linspace(0,tmax,tsteps)

phis = []
for t in ts:
    phis.append(phi(t,zs,ps,spect))

fig, ax = plt.subplots()
line_re, = ax.plot(zs,phis[0].real,c="b")
line0, = ax.plot(zs,signal.real,c="cyan",ls="-.")

# def animate(nt):
#     line_re.set_data(zs, np.real(phis[nt]))
#     # ax.relim()
#     # ax.autoscale_view()
#     return (line_re,)

# anim = animation.FuncAnimation(fig, animate,
#                                frames=np.arange(tsteps), # t-values
#                                interval=50, # wait time before displaying new frame in ms
#                                blit=True)

# HTML(anim.to_jshtml())
# plt.show()

import AnimPlayer

def update(frame):
    fig.suptitle(r"$t=$"+str(ts[frame]))
    line_re.set_ydata(phis[frame].real)
    line_re.relim()
    line_re.autoscale_view()

ani = AnimPlayer.Player(fig=fig, func=update,maxi=ts.size-1)
plt.show()