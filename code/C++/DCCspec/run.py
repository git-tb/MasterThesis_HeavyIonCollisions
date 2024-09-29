#%%
import subprocess
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from IPython import get_ipython
import glob
import os
import scipy.interpolate
from matplotlib import gridspec

get_ipython().run_line_magic("matplotlib","qt")
plt.style.use("mplstyles/myclassic_white.mplstyle")

#%%
###############################################################
###############################################################
# COMPUTE SPECTRA FROM LIST OF INITIAL DATA
###############################################################
###############################################################

pTmax = 2
NpT = 200

# parentdir = "data/init_real_m140_consteps_20240826_143119/*"
# parentdir = "data/init_real_m226_consteps_20240826_143119/*"
# parentdir = "data/init_real_m312_consteps_20240826_143137/*"
# parentdir = "data/init_real_m398_consteps_20240826_143151/*"
# parentdir = "data/init_real_m484_consteps_20240826_143205/*"
# parentdir = "data/init_real_m570_consteps_20240826_143224/*"
# parentdir = "data/init_real_m656_consteps_20240826_143237/*"
# parentdir = "data/init_real_m742_consteps_20240826_143253/*"
# parentdir = "data/init_real_m828_consteps_20240826_143307/*"
# parentdir = "data/init_real_m914_consteps_20240826_143326/*"
# parentdir = "data/init_real_m1000_consteps_20240826_143341/*"

parentdir = "data/init_real_m140_s-1_consteps_20240926_124611/*"
# parentdir = "data/init_real_m226_s-1_consteps_20240926_124618/*"
# parentdir = "data/init_real_m312_s-1_consteps_20240926_124625/*"
# parentdir = "data/init_real_m398_s-1_consteps_20240926_124635/*"
# parentdir = "data/init_real_m484_s-1_consteps_20240926_124646/*"
# parentdir = "data/init_real_m570_s-1_consteps_20240926_124658/*"
# parentdir = "data/init_real_m656_s-1_consteps_20240926_124709/*"
# parentdir = "data/init_real_m742_s-1_consteps_20240926_124716/*"
# parentdir = "data/init_real_m828_s-1_consteps_20240926_124724/*"
# parentdir = "data/init_real_m914_s-1_consteps_20240926_124737/*"
# parentdir = "data/init_real_m1000_s-1_consteps_20240926_124748/*"

m = 0.140
# m = 0.226
# m = 0.312
# m = 0.398
# m = 0.484
# m = 0.570
# m = 0.656
# m = 0.742
# m = 0.828
# m = 0.914
m = 1.000

initpaths = glob.glob(parentdir)
initpaths.sort(key=lambda s:(len(s), s))

for initpath in initpaths:
    result = subprocess.run(args=[
        "./bin/spec",
        "--m=%f"%(m),
        "--pTmax=%f"%(pTmax),
        "--NpT=%d"%(NpT),
        "--epsabs=0",
        "--epsrel=1e-5",
        "--iter=10000",
        "--initpath=%s"%(initpath)
    ])
    print(result)

newdir = "data/newspectra"
idx = 0
while(os.path.isdir(newdir+str(idx))):
    idx += 1
subprocess.run(args=["mkdir",newdir+str(idx)])

lastspecs = glob.glob("data/spec_????????_??????")
for spec in lastspecs:
    subprocess.run(args=["mv",spec,newdir+str(idx)])

#%%
###############################################################
###############################################################
# COMPUTE SINGLE SPECTRUM
###############################################################
###############################################################

pTmax = 2
NpT = 200
m = 0.14

initpath = "data/init_real_pi_constfield_20240822_161734/init5.csv"
result = subprocess.run(args=[
        "./bin/spec",
        "--m=%f"%(m),
        "--pTmax=%f"%(pTmax),
        "--NpT=%d"%(NpT),
        "--epsabs=0",
        "--epsrel=1e-5",
        "--iter=10000",
        "--initpath=%s"%(initpath)
    ])
print(result)

#%%
###############################################################
###############################################################
# PLOT LIST OF REAL FIELD SPECTRA WITH SURFACE GEOMETRIES
###############################################################
###############################################################

###
# get_ipython().run_line_magic("matplotlib","inline")
get_ipython().run_line_magic("matplotlib","qt")

### TOGGLE WHETHER FIGURES SHOULD BE DISPLAYED IN INTERACTIVE MODE
plt.ioff()
# plt.ion()

parentdir = "data/fosurf_scaletest/*"

paths = glob.glob(parentdir)

fig_fo, ax_fo = plt.subplots(figsize=(7,7))
fig_spec, ax_spec = plt.subplots(figsize=(7,7))

ax_spec.set_yscale("log")

for (n,path) in enumerate(paths):
    df_spec = pd.read_csv(path+"/spectr.txt",comment="#")
    df_tau = pd.read_csv(path+"/tau_interp.txt",comment="#")
    df_r = pd.read_csv(path+"/r_interp.txt",comment="#")

    t = n / (len(paths)-1)
    col=(t,0,1-t)

    ax_fo.plot(df_r["rRe"].to_numpy(),
               df_tau["tauRe"].to_numpy(),
               marker="",
               color=col)

    ax_spec.plot(
        df_spec["pT"].to_numpy(),
        df_spec["abs2Re"].to_numpy(),
        marker="",
        color=col)

fig_fo.suptitle(parentdir)
fig_spec.suptitle(parentdir)

fig_fo.tight_layout()
fig_spec.tight_layout()

plt.show()
# plt.close()
#%%
###############################################################
###############################################################
# COMPUTE SPECTRA FOR DIFFERENT MASSES BUT SAME INITIAL DATA
###############################################################
###############################################################

pTmax = 2
NpT = 200
# initpath = "data/init_real_pi_constfield_20240822_161734/init5.csv"
# initpath = "data/init_real_pi_consteps_20240822_135440/init8.csv"
initpath = "data/init_real_pi_consteps_20240822_135426/init8.csv"
ms = np.linspace(0.14,0.8,10)

for m in ms:
    result = subprocess.run(args=[
        "./bin/spec",
        "--pTmax=%f"%(pTmax),
        "--NpT=%d"%(NpT),
        "--m=%.2f"%(m),
        "--epsabs=0",
        "--epsrel=1e-5",
        "--iter=10000",
        "--initpath=%s"%(initpath)
    ])
    print(result)

newdir = "data/newspectra"
idx = 0
while(os.path.isdir(newdir+str(idx))):
    idx += 1
subprocess.run(args=["mkdir",newdir+str(idx)])

lastspecs = glob.glob("data/spec_????????_??????")
for spec in lastspecs:
    subprocess.run(args=["mv",spec,newdir+str(idx)])

#%%
###############################################################
###############################################################
# PLOT A SINGLE SPECTRUM
###############################################################
###############################################################

path = "data/sigmadecay/spec_20240807_223737"
path = "data/complexfield_inittest/spec_20240809_110228"
path = "data/realfield_inittest/spec_20240807_212830"

df_field0 = pd.read_csv(path+"/field0.txt",comment="#")
df_Dfield0 = pd.read_csv(path+"/field0_deriv.txt",comment="#")
df_spec = pd.read_csv(path+"/spectr.txt",comment="#")
df_spec_anti = pd.read_csv(path+"/spectr_anti.txt",comment="#")

fig_init, (ax_init1, ax_init2) = plt.subplots(nrows=1,ncols=2,figsize=(14,6))
fig_spec, ax_spec = plt.subplots(figsize=(6,6))

fig_init.suptitle(path)
fig_spec.suptitle(path)

ax_init1.plot(
    df_field0["alpha"].to_numpy(),
    df_field0["field0Re"].to_numpy(),
    marker="",
    label=r"$\Re$")
ax_init1.plot(
    df_field0["alpha"].to_numpy(),
    df_field0["field0Im"].to_numpy(),
    c="r",
    marker="",
    label=r"$\Im$")
ax_init2.plot(
    df_Dfield0["alpha"].to_numpy(),
    df_Dfield0["Dfield0Re"].to_numpy(),
    marker="",
    label=r"$\Re$")
ax_init2.plot(
    df_Dfield0["alpha"].to_numpy(),
    df_Dfield0["Dfield0Im"].to_numpy(),
    c="r",
    marker="",
    label=r"$\Im$")
ax_spec.plot(
    df_spec["pT"].to_numpy(),
    df_spec["abs2Re"].to_numpy(),
    lw=2,
    c="b",
    marker="",
    label=r"$\mathrm{particle}$")
ax_spec.plot(
    df_spec_anti["pT"].to_numpy(),
    df_spec_anti["abs2Re"].to_numpy(),
    c="r",
    marker="",
    label=r"$\mathrm{anti particle}$")

ax_init1.set_ylabel(r"$\phi$")
ax_init2.set_ylabel(r"$n^\mu\partial_\mu\phi$")
ax_spec.set_ylabel(r"$\frac{1}{2\pi p^\perp}\frac{\mathrm{d}N}{\mathrm{d}p^\perp\mathrm{d}\eta_p}$")

ax_init1.set_xlabel(r"$\alpha$")
ax_init2.set_xlabel(r"$\alpha$")
ax_spec.set_xlabel(r"$p^\perp$")

ax_spec.set_yscale("log")
ax_spec.set_ylim(10**(-3),10**5)
ax_spec.set_xlim(0,2)

ax_init1.legend()
ax_init2.legend()
ax_spec.legend()

fig_init.tight_layout()
fig_spec.tight_layout()

plt.show()

#%%
###############################################################
###############################################################
# PLOT LIST OF REAL FIELD SPECTRA
###############################################################
###############################################################

###
# get_ipython().run_line_magic("matplotlib","inline")
get_ipython().run_line_magic("matplotlib","qt")

### TOGGLE WHETHER FIGURES SHOULD BE DISPLAYED IN INTERACTIVE MODE
plt.ioff()
# plt.ion()


# parentdir = "data/realfield_inittest/*/"
# parentdir = "data/specsreal_v3/*/"
# parentdir = "data/realfield_inittest_v2/*/"
# parentdir = "data/spectra_real_consteps_20240822_135426/*/"
# parentdir = "data/spectra_real_consteps_20240822_135440/*/"
# parentdir = "data/spectra_real_constfield_20240822_161734/*/"
# parentdir = "data/spectra_real_constfield_20240822_161738/*/"
# parentdir = "data/spectra_real_constfield_20240822_161734_masses/*/"
# parentdir = "data/spectra_real_consteps_20240822_135440_masses/*/"
# parentdir = "data/spectra_real_consteps_20240822_135426_masses/*/"
# parentdir = "data/spectra_real_consteps_20240826_143119_m226/*/"
# parentdir = "data/spectra_real_consteps_20240826_143137_m312/*/"
# parentdir = "data/spectra_real_consteps_20240826_143151_m398/*/"
# parentdir = "data/spectra_real_consteps_20240826_143205_m484/*/"
# parentdir = "data/spectra_real_consteps_20240826_143224_m570/*/"
# parentdir = "data/spectra_real_consteps_20240826_143237_m656/*/"
# parentdir = "data/spectra_real_consteps_20240826_143253_m742/*/"
# parentdir = "data/spectra_real_consteps_20240826_143307_m828/*/"
# parentdir = "data/spectra_real_consteps_20240826_143326_m914/*/"
# parentdir = "data/spectra_real_consteps_20240826_143341_m1000/*/"

# parentdir = "data/spectra_real_consteps_20240926_124618_m226_s-1/*/"
# parentdir = "data/spectra_real_consteps_20240926_124625_m312_s-1/*/"
# parentdir = "data/spectra_real_consteps_20240926_124635_m398_s-1/*/"
# parentdir = "data/spectra_real_consteps_20240926_124646_m484_s-1/*/"
# parentdir = "data/spectra_real_consteps_20240926_124658_m570_s-1/*/"
# parentdir = "data/spectra_real_consteps_20240926_124709_m656_s-1/*/"
# parentdir = "data/spectra_real_consteps_20240926_124716_m742_s-1/*/"
# parentdir = "data/spectra_real_consteps_20240926_124724_m828_s-1/*/"
# parentdir = "data/spectra_real_consteps_20240926_124737_m914_s-1/*/"
parentdir = "data/spectra_real_consteps_20240926_124748_m1000_s-1/*/"
paths = glob.glob(parentdir)

# fig_init, (ax_init1, ax_init2) = plt.subplots(nrows=1,ncols=2,figsize=(15,7))
fig_init = plt.figure(figsize=(7,7))
gs = gridspec.GridSpec(nrows=2, ncols=1, hspace=0)
ax_init1, ax_init2 = fig_init.add_subplot(gs[0]), fig_init.add_subplot(gs[1])

fig_spec, ax_spec = plt.subplots(figsize=(7,7))
fig_fullspec, ax_fullspec = plt.subplots(figsize=(7,7))

ax_init1.set_ylabel(r"$\pi^0\ [\mathrm{GeV}]$")
ax_init2.set_ylabel(r"$n^\mu\partial_\mu\pi^0\ [\mathrm{GeV}^2]$")
ax_spec.set_ylabel(r"$\frac{1}{2\pi p_T}\frac{\mathrm{d}N}{\mathrm{d}p_T\mathrm{d}\eta_p}\ [\mathrm{GeV}^{-2}]$")

ax_init1.set_xlabel(r"$\alpha$")
ax_init2.set_xlabel(r"$\alpha$")
ax_spec.set_xlabel(r"$p_T\ [\mathrm{GeV}]$")

ax_spec.set_yscale("log")
ax_fullspec.set_yscale("log")

###
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
    marker=""
)
ax_decay2.plot(
    df_decayspec["p"].to_numpy(),
    df_decayspec["finalspecRe"].to_numpy(),
    marker=""
)

ax_decay1.set_yscale("log")
ax_decay2.set_yscale("log")

ax_decay1.set_xlabel(r"$p^\perp$")
ax_decay2.set_xlabel(r"$p^\perp$")
ax_fullspec.set_xlabel(r"$p^\perp$")

ax_decay1.set_ylabel(r"$\frac{1}{2\pi p^\perp}\frac{\mathrm{d}N}{\mathrm{d}p^\perp\mathrm{d}\eta_p}$")
ax_decay2.set_ylabel(r"$\frac{1}{2\pi p^\perp}\frac{\mathrm{d}N}{\mathrm{d}p^\perp\mathrm{d}\eta_p}$")
ax_fullspec.set_ylabel(r"$\frac{1}{2\pi p^\perp}\frac{\mathrm{d}N}{\mathrm{d}p^\perp\mathrm{d}\eta_p}$")

fig_decay.tight_layout()
###

ax_fullspec.set_xlim(ax_decay2.get_xlim())

minval1, maxval1 = np.inf, -np.inf
minval2, maxval2 = np.inf, -np.inf
for (n,path) in enumerate(paths):
    df_spec = pd.read_csv(path+"/spectr.txt",comment="#")
    df_field0 = pd.read_csv(path+"/field0.txt",comment="#")
    df_Dfield0 = pd.read_csv(path+"/field0_deriv.txt",comment="#")

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
    
    mindata1 = np.min(df_field0["field0Re"].to_numpy())
    maxdata1 = np.max(df_field0["field0Re"].to_numpy())
    minval1 = np.min((minval1,mindata1))
    maxval1 = np.max((maxval1,maxdata1))

    mindata2 = np.min(df_Dfield0["Dfield0Re"].to_numpy())
    maxdata2 = np.max(df_Dfield0["Dfield0Re"].to_numpy())
    minval2 = np.min((minval2,mindata2))
    maxval2 = np.max((maxval2,maxdata2))

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

myrange1 = maxval1 - minval1
myrange2 = maxval2 - minval2
if(myrange1 > 0):
    ax_init1.set_ylim(minval1-0.05*myrange1,maxval1+0.05*myrange1)
if(myrange2 > 0):
    ax_init2.set_ylim(minval2-0.05*myrange2,maxval2+0.05*myrange2)

ax_init1.set_xticklabels(ax_init1.get_xticklabels(), visible=False)
labels = ax_init1.get_yticklabels()
labels[0] = labels[-1] = ""
ax_init1.set_yticklabels(labels)
labels = ax_init2.get_yticklabels()
labels[0] = labels[-1] = ""
ax_init2.set_yticklabels(labels)


xticks = [0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2]
xticklabels = [r"$0$",r"$\pi/8$",r"$\pi/4$",r"$3\pi/8$",r"$\pi/2$"]
ax_init1.set_xticks(xticks, xticklabels)
ax_init2.set_xticks(xticks, xticklabels)

ax_init1.set_xlim(0,np.pi/2)
ax_init2.set_xlim(0,np.pi/2)

fig_init.suptitle(parentdir)
fig_spec.suptitle(parentdir)
fig_decay.suptitle(decaypath)
fig_fullspec.suptitle(parentdir+"\n"+decaypath)

fig_init.tight_layout()
fig_spec.tight_layout()
fig_fullspec.tight_layout()

fig_init.savefig("data/images/"+parentdir.replace("data/","").replace("/*/","")+"_init.png",dpi=150)
fig_spec.savefig("data/images/"+parentdir.replace("data/","").replace("/*/","")+"_spec.png",dpi=150)

plt.close()

#%%
###############################################################
###############################################################
# PLOT LIST OF COMPLEX FIELD SPECTRA
###############################################################
###############################################################

###
# get_ipython().run_line_magic("matplotlib","inline")
get_ipython().run_line_magic("matplotlib","qt")

### TOGGLE WHETHER FIGURES SHOULD BE DISPLAYED IN INTERACTIVE MODE
plt.ioff()
# plt.ion()

# parentdir = "data/spectra_comp_composed_consteps_20240822_145641/*/"
# parentdir = "data/spectra_comp_composed_consteps_20240822_145647/*/"
parentdir = "data/spectra_comp_constfield_20240822_174903/*/"
# parentdir = "data/spectra_comp_constfield_20240822_162959/*/"
paths = sorted(glob.glob(parentdir))

fig_init, ((ax_init1, ax_init2),(ax_init1_im, ax_init2_im)) = plt.subplots(nrows=2,ncols=2,figsize=(12,12))
fig_spec, (ax_spec,ax_specanti) = plt.subplots(nrows=1,ncols=2,figsize=(12,6))
fig_fullspec, (ax_fullspec,ax_fullspecanti) = plt.subplots(nrows=1,ncols=2,figsize=(12,6))

ax_init1.set_ylabel(r"$\Re(\pi^+)$")
ax_init2.set_ylabel(r"$\Re(n^\mu\partial_\mu\pi^+$)")
ax_init1_im.set_ylabel(r"$\Im(\pi^+)$")
ax_init2_im.set_ylabel(r"$\Im(n^\mu\partial_\mu\pi^+$)")
ax_spec.set_ylabel(r"$\frac{1}{2\pi p^\perp}\frac{\mathrm{d}N}{\mathrm{d}p^\perp\mathrm{d}\eta_p}$")
ax_specanti.set_ylabel(r"$\frac{1}{2\pi p^\perp}\frac{\mathrm{d}\overline{N}}{\mathrm{d}p^\perp\mathrm{d}\eta_p}$")
ax_fullspec.set_ylabel(r"$\frac{1}{2\pi p^\perp}\frac{\mathrm{d}N}{\mathrm{d}p^\perp\mathrm{d}\eta_p}$")
ax_fullspecanti.set_ylabel(r"$\frac{1}{2\pi p^\perp}\frac{\mathrm{d}\overline{N}}{\mathrm{d}p^\perp\mathrm{d}\eta_p}$")

ax_init1.set_xlabel(r"$\alpha$")
ax_init2.set_xlabel(r"$\alpha$")
ax_init1_im.set_xlabel(r"$\alpha$")
ax_init2_im.set_xlabel(r"$\alpha$")
ax_spec.set_xlabel(r"$p^\perp$")
ax_specanti.set_xlabel(r"$p^\perp$")
ax_fullspec.set_xlabel(r"$p^\perp$")
ax_fullspecanti.set_xlabel(r"$p^\perp$")

ax_spec.set_yscale("log")
ax_specanti.set_yscale("log")
ax_fullspec.set_yscale("log")
ax_fullspecanti.set_yscale("log")

######
decaypath = "data/sigmadecay/decayspec_20240807_230447"

df_decayspec = pd.read_csv(decaypath+"/decayspec.txt",comment="#")
df_decayprimespec = pd.read_csv(decaypath+"/primespec_interp.txt",comment="#")

decayfunc = scipy.interpolate.CubicSpline(
    df_decayspec["p"].to_numpy()[1:],
    0.6666667*df_decayspec["finalspecRe"].to_numpy()[1:])

fig_decay, (ax_decay1, ax_decay2) = plt.subplots(nrows=1,ncols=2,figsize=(12,6))
fig_decay.suptitle(decaypath)

ax_decay1.plot(
    df_decayprimespec["q"].to_numpy(),
    df_decayprimespec["primespecRe"].to_numpy(),
    marker=""
)
ax_decay2.plot(
    df_decayspec["p"].to_numpy(),
    df_decayspec["finalspecRe"].to_numpy(),
    marker=""
)

ax_decay1.set_yscale("log")
ax_decay2.set_yscale("log")

ax_decay1.set_xlabel(r"$p^\perp$")
ax_decay2.set_xlabel(r"$p^\perp$")
ax_fullspec.set_xlabel(r"$p^\perp$")

ax_decay1.set_ylabel(r"$\frac{1}{2\pi p^\perp}\frac{\mathrm{d}N}{\mathrm{d}p^\perp\mathrm{d}\eta_p}$")
ax_decay2.set_ylabel(r"$\frac{1}{2\pi p^\perp}\frac{\mathrm{d}N}{\mathrm{d}p^\perp\mathrm{d}\eta_p}$")
ax_fullspec.set_ylabel(r"$\frac{1}{2\pi p^\perp}\frac{\mathrm{d}N}{\mathrm{d}p^\perp\mathrm{d}\eta_p}$")

fig_decay.tight_layout()
######

for (n,path) in enumerate(paths):
    df_spec = pd.read_csv(path+"spectr.txt",comment="#")
    df_specanti = pd.read_csv(path+"spectr_anti.txt",comment="#")
    df_field0 = pd.read_csv(path+"field0.txt",comment="#")
    df_Dfield0 = pd.read_csv(path+"field0_deriv.txt",comment="#")

    t = n / (len(paths)-1)
    col=(t,0,1-t)
    col2=(t,1,1-t)

    ax_init1.plot(
        df_field0["alpha"].to_numpy(),
        df_field0["field0Re"].to_numpy(),
        marker="",
        color=col)    
    ax_init1_im.plot(
        df_field0["alpha"].to_numpy(),
        df_field0["field0Im"].to_numpy(),
        marker="",
        color=col)

    ax_init2.plot(
        df_Dfield0["alpha"].to_numpy(),
        df_Dfield0["Dfield0Re"].to_numpy(),
        marker="",
        color=col)
    ax_init2_im.plot(
        df_Dfield0["alpha"].to_numpy(),
        df_Dfield0["Dfield0Im"].to_numpy(),
        marker="",
        color=col)

    ax_spec.plot(
        df_spec["pT"].to_numpy(),
        df_spec["abs2Re"].to_numpy(),
        marker="",
        color=col)
    ax_specanti.plot(
        df_specanti["pT"].to_numpy(),
        df_specanti["abs2Re"].to_numpy(),
        marker="",
        color=col)
    
    ax_fullspec.plot(
        df_spec["pT"].to_numpy(),
        df_spec["abs2Re"].to_numpy()+decayfunc(df_spec["pT"].to_numpy()),
        marker="",
        color=col)
    ax_fullspecanti.plot(
        df_specanti["pT"].to_numpy(),
        df_specanti["abs2Re"].to_numpy()+decayfunc(df_specanti["pT"].to_numpy()),
        marker="",
        color=col)

fig_init.suptitle(parentdir)
fig_spec.suptitle(parentdir)
fig_fullspec.suptitle(parentdir+"\n"+decaypath)

fig_init.tight_layout()
fig_spec.tight_layout()
fig_fullspec.tight_layout()

fig_init.savefig("data/images/"+parentdir.replace("data/","").replace("/*/","")+"_init.png",dpi=150)
fig_spec.savefig("data/images/"+parentdir.replace("data/","").replace("/*/","")+"_spec.png",dpi=150)

plt.close()

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
# PLOT DECAYTEST SPECTRA

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

fig_decayspec, (ax_primespec,ax_decayspec)= plt.subplots(ncols=2,figsize=(15,7))
fig_decayspec.suptitle("data/decaytest_1to10Gev")

df_ps = pd.read_csv(paths[-1]+"primespec_interp.txt",comment="#")
ax_primespec.plot(df_ps["q"],df_ps["primespecRe"],marker="")
ax_primespec.set_yscale("log")
ax_primespec.set_ylabel(r"$\frac{1}{2\pi q^\perp}\frac{\mathrm{d}N}{\mathrm{d}q^\perp\mathrm{d}\eta_q}$")
ax_primespec.set_xlabel(r"$q^\perp$")

istart, iend = 0, 200
axin = ax_primespec.inset_axes([0.35,0.45,0.6,0.5])
axin.plot(df_ps["q"][istart:iend],df_ps["primespecRe"][istart:iend],marker="")
axin.set_yscale("log")

for path in paths:
    df_ps = pd.read_csv(path+"primespec_interp.txt",comment="#")
    qmax = max(df_ps["q"].to_numpy())

    df_ds = pd.read_csv(path+"decayspec.txt",comment="#")
    ax_decayspec.plot(df_ds["p"].to_numpy(),df_ds["finalspecRe"].to_numpy(),
                label=r"$q_{\mathrm{max}}=%.1f\,\mathrm{GeV}$"%(qmax),ls="None",ms=10)

ax_decayspec.set_yscale("log")

ax_decayspec.set_ylabel(r"$\frac{1}{2\pi p^\perp}\frac{\mathrm{d}N}{\mathrm{d}p^\perp\mathrm{d}\eta_p}$")
ax_decayspec.set_xlabel(r"$p^\perp$")

fig_decayspec.tight_layout()

plt.legend()
plt.show()
# %%
#############################################################
########### COMPUTE SMOOTHER FREEZEOUT GEOMETRY #############
#############################################################

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit

path = "./../../../code/Mathematica/data/ExampleFreezeOut.csv"
df = pd.read_csv(path)
print(df)

factor = 10

df["tau"] *= factor
df["r"] *= factor
df["Dtau"] *= factor
df["Dr"] *= factor

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

# plt.show()
plt.close()
print(newdf)
newdf.to_csv("./../../../code/Mathematica/data/ExampleFreezeOutCorrected_"+str(factor)+"x.csv",index=False)

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