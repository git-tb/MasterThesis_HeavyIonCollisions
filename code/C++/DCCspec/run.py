#%%
import subprocess
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from IPython import get_ipython

get_ipython().run_line_magic("matplotlib","qt")
plt.style.use("mplstyles/myclassic_white.mplstyle")

#%%
# COMPUTE A SPECTRUM FROM INITIALDATA

pTmax = 10
NpT = 100
# initpath = "data/init_20240730_172100/initialfields_pi0.csv"
initpath = "initial1.csv"

result = subprocess.run(args=[
    "./bin/spec",
    "--pTmax=%f"%(pTmax),
    "--NpT=%d"%(NpT),
    "--epsabs=0",
    "--epsrel=1e-6",
    "--iter=10000",
    "--initpath=%s"%(initpath)
])
print(result)

#%%
# PLOT SPECTRUM

# path = "data/spec_20240731_142132/"
# path = "data/spec_20240731_143403/"
# path = "data/spec_20240731_150822/"
# path = "data/spec_20240731_151921/"
# path = "data/spec_20240731_130022/"
path = "data/spec_20240801_111153/"
prop = "ur"

df_spec = pd.read_csv(path+"spectr.txt",comment="#")
df_specanti = pd.read_csv(path+"spectr_anti.txt",comment="#")
df_init = pd.read_csv(path+prop+"_samp.txt",comment="#")
df_init2 = pd.read_csv(path+prop+"_interp.txt",comment="#")
df_field = pd.read_csv(path+"field0.txt",comment="#")

fig, ax = plt.subplots()
ax.plot(df_init2["alpha"].to_numpy(), df_init2[prop+"Re"].to_numpy(),c="r")
ax.scatter(df_init["alpha"].to_numpy(), df_init[prop+"Re"].to_numpy(),c="g")
ax.set_title(prop)
fig.suptitle(path)
plt.show()

fig, ax = plt.subplots()
ax.plot(df_spec["pT"].to_numpy(), df_spec["abs2Re"].to_numpy(),lw=3)
ax.plot(df_specanti["pT"].to_numpy(), df_specanti["abs2Re"].to_numpy(),lw=1)
ax.set_yscale("log")
fig.suptitle(path)
plt.show()

fig, ax = plt.subplots()
ax.plot(df_field["alpha"].to_numpy(), df_field["field0Re"].to_numpy())
ax.plot(df_field["alpha"].to_numpy(), df_field["field0Im"].to_numpy())
fig.suptitle(path)
plt.show()

# %%
# COMPUTE DECAY SPECTRUM FROM PRIMARY SPECTRUM

pTmax = 1
NpT = 100
# primespecpath = "data/spec_20240730_131734/spectr.txt"
primespecpath = "data/spec_20240730_131221/spectr.txt"

subprocess.run(args=[
    "./bin/decay",
    "--pTmax=%f"%(pTmax),
    "--NpT=%d"%(NpT),
    "--primespecpath=%s"%(primespecpath)
])
# %%
# PLOT SPECTRUM

path1="data/decayspec_20240730_131835/"
path2="data/decayspec_20240730_132116/"

df_ds1 = pd.read_csv(path1+"decayspec.txt",comment="#")
df_ps1 = pd.read_csv(path1+"primespec_interp.txt",comment="#")
df_ds2 = pd.read_csv(path2+"decayspec.txt",comment="#")
df_ps2 = pd.read_csv(path2+"primespec_interp.txt",comment="#")

fig, ax = plt.subplots(nrows=1,ncols=2)
ax[0].plot(df_ps1["q"].to_numpy(),df_ps1["primespecRe"].to_numpy())
ax[1].plot(df_ds1["p"].to_numpy(),df_ds1["finalspecRe"].to_numpy())
ax[0].plot(df_ps2["q"].to_numpy(),df_ps2["primespecRe"].to_numpy())
ax[1].plot(df_ds2["p"].to_numpy(),df_ds2["finalspecRe"].to_numpy())
ax[0].set_yscale("log")
ax[1].set_yscale("log")
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