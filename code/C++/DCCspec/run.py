#%%
import subprocess
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from IPython import get_ipython
import glob

get_ipython().run_line_magic("matplotlib","qt")
plt.style.use("mplstyles/myclassic_white.mplstyle")

#%%
# COMPUTE A SPECTRUM FROM INITIALDATA

pTmax = 2
NpT = 200

initpaths = [
    "data/init_real/init0.csv",
    "data/init_real/init1.csv",
    "data/init_real/init2.csv",
    "data/init_real/init3.csv",
    "data/init_real/init4.csv",
    "data/init_real/init5.csv",
    "data/init_real/init6.csv",
    "data/init_real/init7.csv",
    "data/init_real/init8.csv",
    "data/init_real/init9.csv",
    "data/init_real/init10.csv",
    ""][:-1]

for initpath in initpaths:
    result = subprocess.run(args=[
        "./bin/spec",
        "--pTmax=%f"%(pTmax),
        "--NpT=%d"%(NpT),
        "--epsabs=0",
        "--epsrel=1e-5",
        "--iter=10000",
        "--initpath=%s"%(initpath)
    ])
    print(result)

#%%
# PLOT A SINGLE SPECTRUM

path = "data/sigmadecay/spec_20240807_223737"
# path = "data/realfield_inittest/spec_20240807_212830"

df_field0 = pd.read_csv(path+"/field0.txt",comment="#")
df_Dfield0 = pd.read_csv(path+"/field0_deriv.txt",comment="#")
df_spec = pd.read_csv(path+"/spectr.txt",comment="#")

fig_init, (ax_init1, ax_init2) = plt.subplots(nrows=1,ncols=2,figsize=(15,7))
fig_spec, ax_spec = plt.subplots(figsize=(7,7))

ax_init1.plot(
    df_field0["alpha"].to_numpy(),
    df_field0["field0Re"].to_numpy())
ax_init2.plot(
    df_Dfield0["alpha"].to_numpy(),
    df_Dfield0["Dfield0Re"].to_numpy())
ax_spec.plot(
    df_spec["pT"].to_numpy(),
    df_spec["abs2Re"].to_numpy())

ax_init1.set_ylabel(r"$\phi$")
ax_init2.set_ylabel(r"$n^\mu\partial_\mu\phi$")
ax_spec.set_ylabel(r"$\frac{1}{2\pi p^\perp}\frac{\mathrm{d}N}{\mathrm{d}p^\perp\mathrm{d}\eta_p}$")

ax_init1.set_xlabel(r"$\alpha$")
ax_init2.set_xlabel(r"$\alpha$")
ax_spec.set_xlabel(r"$p^\perp$")

ax_spec.set_yscale("log")

plt.show()

#%%
# PLOT LIST OF SPECTRA

import scipy.interpolate

paths = glob.glob("data/realfield_inittest/*/")

fig_init, (ax_init1, ax_init2) = plt.subplots(nrows=1,ncols=2,figsize=(15,7))
fig_spec, ax_spec = plt.subplots(figsize=(7,7))
fig_fullspec, ax_fullspec = plt.subplots(figsize=(7,7))

ax_init1.set_ylabel(r"$\pi^0$")
ax_init2.set_ylabel(r"$n^\mu\partial_\mu\pi^0$")
ax_spec.set_ylabel(r"$\frac{1}{2\pi p^\perp}\frac{\mathrm{d}N}{\mathrm{d}p^\perp\mathrm{d}\eta_p}$")

ax_init1.set_xlabel(r"$\alpha$")
ax_init2.set_xlabel(r"$\alpha$")
ax_spec.set_xlabel(r"$p^\perp$")

ax_spec.set_yscale("log")
ax_fullspec.set_yscale("log")

decaypath = "data/sigmadecay/decayspec_20240807_230447"

df_decayspec = pd.read_csv(decaypath+"/decayspec.txt",comment="#")
df_decayprimespec = pd.read_csv(decaypath+"/primespec_interp.txt",comment="#")

decayfunc = scipy.interpolate.CubicSpline(
    df_decayspec["p"].to_numpy()[1:],
    0.6666667*df_decayspec["finalspecRe"].to_numpy()[1:])

fig_decay, (ax_decay1, ax_decay2) = plt.subplots(nrows=1,ncols=2,figsize=(15,7))
ax_decay1.plot(
    df_decayprimespec["q"].to_numpy(),
    df_decayprimespec["primespecRe"].to_numpy(),
)
ax_decay2.plot(
    df_decayspec["p"].to_numpy(),
    df_decayspec["finalspecRe"].to_numpy()
)

ax_decay1.set_yscale("log")
ax_decay2.set_yscale("log")

ax_decay1.set_xlabel(r"$p^\perp$")
ax_decay2.set_xlabel(r"$p^\perp$")

for (n,path) in enumerate(paths):
    df_spec = pd.read_csv(path+"spectr.txt",comment="#")
    df_field0 = pd.read_csv(path+"field0.txt",comment="#")
    df_Dfield0 = pd.read_csv(path+"field0_deriv.txt",comment="#")

    t = n / (len(paths)-1)
    col=(t,0,1-t)

    ax_init1.plot(
        df_field0["alpha"].to_numpy(),
        df_field0["field0Re"].to_numpy(),
        marker="",
        color=col)

    ax_init2.plot(
        df_Dfield0["alpha"].to_numpy(),
        df_Dfield0["Dfield0Re"].to_numpy(),
        marker="",
        color=col)

    ax_spec.plot(
        df_spec["pT"].to_numpy(),
        df_spec["abs2Re"].to_numpy(),
        marker="",
        color=col)
    
    ax_fullspec.plot(
        df_spec["pT"].to_numpy(),
        df_spec["abs2Re"].to_numpy()+decayfunc(df_spec["pT"].to_numpy()),
        marker="",
        color=col)

fig_init.tight_layout()
fig_spec.tight_layout()
plt.show()

# %%
# COMPUTE DECAY SPECTRUM FROM PRIMARY SPECTRUM

pTmax = 1
NpT = 100
# primespecpath = "data/spec_20240730_131734/spectr.txt"
primespecpath = "data/spec_20240806_151611/spectr.txt"

subprocess.run(args=[
    "./bin/decay",
    "--pTmax=%f"%(pTmax),
    "--NpT=%d"%(NpT),
    "--primespecpath=%s"%(primespecpath)
])
# %%
# PLOT SPECTRUM

plt.style.use("mplstyles/myclassic_white.mplstyle")

paths = [
"data/decaytest_1to10Gev/decayspec_20240806_161324/",
"data/decaytest_1to10Gev/decayspec_20240806_161542/",
"data/decaytest_1to10Gev/decayspec_20240806_162432/",
"data/decaytest_1to10Gev/decayspec_20240806_162634/",
"data/decaytest_1to10Gev/decayspec_20240806_163017/",
"data/decaytest_1to10Gev/decayspec_20240806_163439/",
"data/decaytest_1to10Gev/decayspec_20240806_163743/",
"data/decaytest_1to10Gev/decayspec_20240806_164059/",
"data/decaytest_1to10Gev/decayspec_20240806_164445/",
"data/decaytest_1to10Gev/decayspec_20240806_164831/",
""][:-1]

fig, ax = plt.subplots()

for path in paths:
    df_ps = pd.read_csv(path+"primespec_interp.txt",comment="#")
    qmax = max(df_ps["q"].to_numpy())

    df_ds = pd.read_csv(path+"decayspec.txt",comment="#")
    ax.plot(df_ds["p"].to_numpy(),df_ds["finalspecRe"].to_numpy(),
                label=r"$q_{\mathrm{max}}=%.1f\,\mathrm{GeV}$"%(qmax),ls="None",ms=10)

ax.set_yscale("log")

ax.set_ylabel(r"$\frac{1}{2\pi p^\perp}\frac{\mathrm{d}N}{\mathrm{d}p^\perp\mathrm{d}\eta_p}$")
ax.set_xlabel(r"$p^\perp$")

plt.legend()
plt.show()
# %%

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit

path = "./../../../code/Mathematica/data/ExampleFreezeOut.csv"
df = pd.read_csv(path)
print(df)

alpha = df["alpha"].to_numpy()

keys = ["tau","r","Dtau","Dr","ur","utau"]
tau = df["tau"].to_numpy()
r = df["r"].to_numpy()
Dtau = df["Dtau"].to_numpy()
Dr = df["Dr"].to_numpy()
ur = df["ur"].to_numpy()
utau = df["utau"].to_numpy()

newdf = pd.DataFrame()

newalpha = np.linspace(0,np.pi/2.0 + 0.001,500)
newdf["alpha"] = newalpha

for key in keys:
    popt = np.polyfit(alpha,df[key].to_numpy(),55)    
    fitfunc = np.poly1d(popt)

    offset = 0
    if(key == "Dtau" or key == "r" or key == "ur"):
        offset = fitfunc(0)
    
    newdf[key] = fitfunc(newalpha)-offset

    fig, ax = plt.subplots()
    ax.scatter(alpha, df[key].to_numpy())
    ax.plot(newalpha, fitfunc(newalpha)-offset)
    ax.set_title(key)

plt.show()
print(newdf)
newdf.to_csv("./../../../code/Mathematica/data/ExampleFreezeOutCorrected.csv",index=False)

# %%
# CREATE ESSENTIALLY RANDOM INITIAL DATA

def f0(alpha):
    return np.sin(11.4*alpha) * np.exp(-(alpha-1.4*1j)**2)

def Df0(alpha):
    return np.cos(12-7*alpha) * np.exp(-3*(alpha-0.7*1j)**2)

x = np.linspace(0,np.pi/2,100)
f0x = f0(x)
Df0x = Df0(x)
plt.plot(x,np.real(f0x))
plt.plot(x,np.imag(f0x))
plt.show()
plt.plot(x,np.real(Df0x))
plt.plot(x,np.imag(Df0x))
plt.show()

mydf = pd.DataFrame()
mydf["alpha"] = x
mydf["f0Re"] = np.real(f0x)
mydf["f0Im"] = np.imag(f0x)
mydf["Df0Re"] = np.real(Df0x)
mydf["Df0Im"] = np.imag(Df0x)
mydf.to_csv("initial1.csv",index=False,sep=",")