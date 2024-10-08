#%%
import numpy as np
from matplotlib import pyplot as plt
from IPython import get_ipython

get_ipython().run_line_magic("matplotlib","qt")
plt.style.use("mplstyles/myclassic_white.mplstyle")
#%%
###
### PLOT INTEGRATION DOMAIN FOR DECAY COMPUTATION (RESOLVING DIRAC DELTA)
###

###
# get_ipython().run_line_magic("matplotlib","inline")
get_ipython().run_line_magic("matplotlib","qt")

### TOGGLE WHETHER FIGURES SHOULD BE DISPLAYED IN INTERACTIVE MODE
plt.ioff()
# plt.ion()

ma, mb = 0.6, 0.14
mc = mb
Q = 1
pabc = (1/(2*ma))*np.sqrt(((ma+mb)**2 - mc**2)*((ma-mb)**2-mc**2))
Eabc = np.sqrt(mb**2 + pabc**2)

def w(m,p):
    return np.sqrt(m**2 + p**2)

def v1(t, p):
    return (w(mb, p)*w(ma, t*Q) - ma*Eabc)/(t*Q*p)

def vmin(p):
    if(p > pabc):
        return np.sqrt(p**2-pabc**2)/p
    return -1

def tmin(p):
    return np.clip(ma*(p*Eabc - w(mb,p)*np.sqrt(Eabc**2 - mb**2))/(Q*mb**2),0,np.inf)

def tmax(p):
    return ma*(p*Eabc + w(mb,p)*np.sqrt(Eabc**2 - mb**2))/(Q*mb**2)

def vlow(t,p):
    return np.clip(v1(t,p),-1, 1)

fig, ax = plt.subplots(figsize=(6*1.6,6))

# Select length of axes and the space between tick labels
xmin, xmax, ymin, ymax = -5, 5, -5, 5
ticks_frequency = 1

# Set bottom and left spines as x and y axes of coordinate system
ax.spines['bottom'].set_position('zero')
ax.spines['left'].set_position('zero')

# Remove top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.set_xticks([], minor=False)
ax.set_yticks([-1,0,1], minor=False)

# Draw major and minor grid lines
ax.grid(which='major', color='black', linewidth=1.5, linestyle='--', alpha=1)

# Draw arrows
arrow_fmt = dict(markersize=10, color='black', clip_on=False)
ax.plot((1), (0), marker='>', transform=ax.get_yaxis_transform(), **arrow_fmt)
ax.plot((0), (1), marker='^', transform=ax.get_xaxis_transform(), **arrow_fmt)

ts = np.linspace(0,8,1000)

# for i in range(3):
#     t = i/(3-1)
#     pmin = 0.3*pabc
#     pmax = 0.8*pabc
#     p = pmin + t*(pmax-pmin)

#     col = (0,0,0.5+t*(1-0.5))

#     ax.plot(ts,v1(ts,p),marker="",lw=3,c=col)
#     ax.fill_between(ts,vlow(ts,p),1,facecolor=(*col,0.2))

c1 = (0,0.4,1)
p1 = 0.5*pabc
ax.plot(ts,v1(ts,p1),marker="",lw=3,c=c1)
ax.fill_between(ts,vlow(ts,p1),1,facecolor=(*c1,0.2))

c2 = (1,0,0)
p2 = 1.4*pabc
ax.plot(ts,v1(ts,p2),marker="",lw=3,c=c2)
ax.fill_between(ts,vlow(ts,p2),1,facecolor=(*c2,0.2))

ax.vlines(tmin(p1),vmin(p1),1,lw=2.5,ls="--",color=c1)
ax.vlines(tmax(p1),vmin(p1),1,lw=2.5,ls="--",color=c1)
ax.hlines(vmin(p1),tmin(p1),tmax(p1),lw=2.5,ls="--",color=c1)
ax.hlines(1,tmin(p1),tmax(p1),lw=2.5,ls="--",color=c1)

ax.vlines(tmin(p2),vmin(p2),1,lw=2.5,ls="--",color=c2)
ax.vlines(tmax(p2),vmin(p2),1,lw=2.5,ls="--",color=c2)
ax.hlines(vmin(p2),tmin(p2),tmax(p2),lw=2.5,ls="--",color=c2)
ax.hlines(1,tmin(p2),tmax(p2),lw=2.5,ls="--",color=c2)

ax.text(4,0.8,r"$\omega_{p,T}>E^a_{b\vert c}$",fontsize=25,c=c2)
ax.text(0.7,-0.3,r"$\omega_{p,T}<E^a_{b\vert c}$",fontsize=25,c=c1)

ax.set_ylim((-1.5,1.5))
ax.set_xlim((0,np.max(ts)))

ax.set_ylabel(r"$v$",rotation=0)
ax.set_xlabel(r"$t$")

ax.xaxis.set_label_coords(1,0.45)
ax.yaxis.set_label_coords(-0.05,0.97)

fig.tight_layout()
fig.savefig("otherfiles/DecayCalc_IntegrDom.png")
plt.show()


