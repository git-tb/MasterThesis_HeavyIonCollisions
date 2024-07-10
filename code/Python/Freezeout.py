import glob
import pandas as pd
from matplotlib import pyplot as plt
import scipy.interpolate
import numpy as np
import scipy.integrate
import scipy.optimize
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
Nalphas = 1500
alphas = np.linspace(0,np.pi/2,1500)
dalphas = (np.pi/2)/Nalphas
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

pi0s = np.cumsum((-utaus*Dtaus+urs*Drs)*dalphas)

xsamples = np.linspace(0,np.pi/2,10)


def func(x, ysamples):
    intermediatefunc = scipy.interpolate.CubicSpline(xsamples,ysamples)  
    return intermediatefunc(x)  

# def target(ysamples):
#     chis = func(alphas, ysamples)
#     pi0s = np.cumsum(chis * (-utaus*Dtaus+urs*Drs)*dalphas)
#     return np.sum(np.abs((0.14*pi0s)**2 + chis**2)*dalphas)

# chis = func(alphas,np.cos(10*xsamples))
# plt.plot(alphas,chis)
# plt.plot(alphas,np.cos(10*alphas))
# plt.show()

# pi0s = np.cumsum(chis * (-utaus*Dtaus+urs*Drs)*dalphas)
# plt.plot(alphas, chis**2 + 0.14**2 * pi0s**2)
# plt.show()

from matplotlib.widgets import Button, Slider

# The parametrized function to be plotted
def f(t, a, b, c, d, e, f, g, h, i, j):
    amps = (a,b,c,d,e,f,g,h,i,j)
    return func(t,amps)

t = np.linspace(0, 1, 1000)

# Define initial parameters
init_amplitude = 5
init_frequency = 3

# Create the figure and the line that we will manipulate
fig, ax = plt.subplots()
line, = ax.plot(t, f(t, init_amplitude, init_frequency), lw=2)
ax.set_xlabel('Time [s]')

# adjust the main plot to make room for the sliders
fig.subplots_adjust(left=0.25, bottom=0.25)

# Make a horizontal slider to control the frequency.
axfreq = fig.add_axes([0.25, 0.1, 0.65, 0.03])
freq_slider = Slider(
    ax=axfreq,
    label='Frequency [Hz]',
    valmin=0.1,
    valmax=30,
    valinit=init_frequency,
)

# Make a vertically oriented slider to control the amplitude
axamp = fig.add_axes([0.1, 0.25, 0.0225, 0.63])
amp_slider = Slider(
    ax=axamp,
    label="Amplitude",
    valmin=0,
    valmax=10,
    valinit=init_amplitude,
    orientation="vertical"
)


# The function to be called anytime a slider's value changes
def update(val):
    line.set_ydata(f(t, amp_slider.val, freq_slider.val))
    fig.canvas.draw_idle()


# register the update function with each slider
freq_slider.on_changed(update)
amp_slider.on_changed(update)

# Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
resetax = fig.add_axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', hovercolor='0.975')


def reset(event):
    freq_slider.reset()
    amp_slider.reset()
button.on_clicked(reset)

plt.show()
