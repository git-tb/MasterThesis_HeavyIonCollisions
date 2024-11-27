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
from matplotlib.widgets import Button, Slider

GOLDEN = 1.61803

get_ipython().run_line_magic("matplotlib","qt")
plt.style.use("mplstyles/myclassic_white.mplstyle")

#%%

def data_to_bins(datax, datay, bins_lower, bins_upper):
    bins_x = ( bins_upper + bins_lower ) / 2

    bins_all = [[] for i in range(len(bins_x))]

    for (x,y) in list(zip(datax, datay)):
        binidx = np.where((x >= bins_lower).astype(int) * (x <= bins_upper).astype(int))[0][0]
        bins_all[binidx].append(y)

    bins_y = []
    for bin_y in bins_all:
        bins_y.append(np.mean(bin_y))

    return bins_x, bins_y

#%%
###############################################################
###############################################################
# COMPUTE SPECTRA FROM LIST OF INITIAL DATA
###############################################################
###############################################################

pTmax = 2
NpT = 200

# parentdir = "data/init_real_m140_s-1_consteps_20240926_124611/*"
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

# parentdir = "data/init_real_m140_consteps_20241029_110232/*"
# parentdir = "data/init_real_m140_constfield_20241029_110519/*"
# parentdir = "data/init_real_m140_taudep_20241029_132506/*"
# parentdir = "data/init_real_m280_consteps_20241029_140005/*"
# parentdir = "data/init_real_m280_constfield_20241029_140034/*"
# parentdir = "data/init_real_m140_constfield_20241029_163722/*"

# parentdir = "data/init_real_m140_s-1_consteps_20241030_104847/*"
# parentdir = "data/init_real_m280_s-1_constfield_20241030_111128/*"
# parentdir = "data/init_real_m280_s-1_consteps_20241030_111119/*"
parentdir = "data/init_real_m280_taudep_20241030_150721/*"
m = 0.28

# m = 0.140
# m = 0.226
# m = 0.312
# m = 0.398
# m = 0.484
# m = 0.570
# m = 0.656
# m = 0.742
# m = 0.828
# m = 0.914
# m = 1.000

initpaths = glob.glob(parentdir)
initpaths.sort(key=lambda s:(len(s), s))
initpaths.reverse()

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

initpath = "data/init_compl_taudep_20241023_171719/init.csv"
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

    scale = df_tau["tauRe"][0]/13.1

    ax_fo.plot(df_r["rRe"].to_numpy(),
               df_tau["tauRe"].to_numpy(),
               marker="",
               color=col)

    ax_spec.plot(
        df_spec["pT"].to_numpy(),
        1/(scale**4)*df_spec["abs2Re"].to_numpy(),
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

pTmax = 5
NpT = 500
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

pTmax = 5
NpT = 500
# files = glob.glob("data/init_real_sigmamasses_consteps_20241119_171318/*")
# files.sort() # sorts normally by alphabetical order
# files.sort(key=len) # sorts by descending length
# files = 40*["data/init_real_sigma_constfield/init0.csv"]
# masses = [312,329,346,364,381,399,416,434,451,468,486,503,521,538,555,573,590,608,625,643,660,677,695,712,730,747,764,782,799,817,834,851,869,886,904,921,939,956,973,991]
# masses = np.array(masses)/1000
masses = np.linspace(280,2,100,endpoint=True)
files = 100*["data/init_real_sigma_constfield/init0.csv"]

for (m,file) in list(zip(masses,files)):
    result = subprocess.run(args=[
        "./bin/spec",
        "--pTmax=%f"%(pTmax),
        "--NpT=%d"%(NpT),
        "--m=%.2f"%(m),
        "--epsabs=0",
        "--epsrel=1e-5",
        "--iter=10000",
        "--initpath=%s"%(file)
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

# path = "data/spec_20241023_153952"
# path = "data/spec_20241023_170656"
path = "data/spec_20241023_171739"

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
# parentdir = "data/spectra_real_consteps_20240926_124748_m1000_s-1/*/"

# paths = glob.glob(parentdir)

folders = [
 'data/spectra_real_consteps_20240822_135426',
 'data/spectra_real_consteps_20240826_143119_m226',
 'data/spectra_real_consteps_20240826_143137_m312',
 'data/spectra_real_consteps_20240826_143151_m398',
 'data/spectra_real_consteps_20240826_143205_m484',
 'data/spectra_real_consteps_20240826_143224_m570',
 'data/spectra_real_consteps_20240826_143237_m656',
 'data/spectra_real_consteps_20240826_143253_m742',
 'data/spectra_real_consteps_20240826_143307_m828',
 'data/spectra_real_consteps_20240826_143326_m914',
 'data/spectra_real_consteps_20240826_143341_m1000']
parentdir = "epsconst_masses"

ms = [
    140.0, 
    226.0, 
    312.0, 
    398.0, 
    484.0, 
    570.0, 
    656.0, 
    742.0, 
    828.0, 
    914.0, 
    1000.0]

folders.sort()
paths = []
for folder in folders:
    localpaths = glob.glob(folder+"/*")
    paths.append(localpaths[0])

parentdir = "epsconst_masses"

# fig_init, (ax_init1, ax_init2) = plt.subplots(nrows=1,ncols=2,figsize=(15,7))
fig_init = plt.figure(figsize=(7,7))
gs = gridspec.GridSpec(nrows=2, ncols=1, hspace=0)
ax_init1, ax_init2 = fig_init.add_subplot(gs[0]), fig_init.add_subplot(gs[1])

fig_spec, ax_spec = plt.subplots(figsize=(7,7))
fig_fullspec, ax_fullspec = plt.subplots(figsize=(7,7))

ax_init1.set_ylabel(r"$m\cdot\pi^0\ [\mathrm{GeV}^2]$")
ax_init2.set_ylabel(r"$n^\mu\partial_\mu\pi^0\ [\mathrm{GeV}^2]$")
ax_spec.set_ylabel(r"$\mathcal{N}\cdot\frac{1}{2\pi p_T}\frac{\mathrm{d}N}{\mathrm{d}p_T\mathrm{d}\eta_p}\ [\mathrm{GeV}^{-2}]$")

ax_init1.set_xlabel(r"$\alpha$")
ax_init2.set_xlabel(r"$\alpha$")
ax_spec.set_xlabel(r"$p_T\ [\mathrm{GeV}]$")

ax_spec.set_yscale("log")
ax_fullspec.set_yscale("log")

###
decaypath = "data/decay_inittest/decayspec_20241010_154911"

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

    m = ms[n]/1000

    ax_init1.plot(
        df_field0["alpha"].to_numpy(),
        m*df_field0["field0Re"].to_numpy(),
        marker="",
        color=col)

    ax_init2.plot(
        df_Dfield0["alpha"].to_numpy(),
        df_Dfield0["Dfield0Re"].to_numpy(),
        marker="",
        color=col)
    
    mindata1 = np.min(m*df_field0["field0Re"].to_numpy())
    maxdata1 = np.max(m*df_field0["field0Re"].to_numpy())
    minval1 = np.min((minval1,mindata1))
    maxval1 = np.max((maxval1,maxdata1))

    mindata2 = np.min(df_Dfield0["Dfield0Re"].to_numpy())
    maxdata2 = np.max(df_Dfield0["Dfield0Re"].to_numpy())
    minval2 = np.min((minval2,mindata2))
    maxval2 = np.max((maxval2,maxdata2))

    scale = 1000/df_spec["abs2Re"].to_numpy()[0]
    print("mass: ", m, " | scalefactor: ", scale, " | 0.63*m2: ", 0.6*m**2)

    ax_spec.plot(
        df_spec["pT"].to_numpy(),
        scale*df_spec["abs2Re"].to_numpy(),
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

# fig_init.savefig("data/images/"+parentdir.replace("data/","").replace("/*/","")+"_init.png",dpi=150)
# fig_spec.savefig("data/images/"+parentdir.replace("data/","").replace("/*/","")+"_spec.png",dpi=150)

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
###############################################################
###############################################################
# COMPUTE DECAY SPECTRA
###############################################################
###############################################################

ma, mb, mc = 0.656, 0.14, 0.14
pTmax = 2
NpT = 200
primespecpath = "data/decay_convtest/spec_20241009_154507_10GeV/spectr.txt"

subprocess.run(args=[
    "./bin/decay",
    "--ma=%f"%(ma),
    "--mb=%f"%(mb),
    "--mc=%f"%(mc),
    "--pTmax=%f"%(pTmax),
    "--NpT=%d"%(NpT),
    "--primespecpath=%s"%(primespecpath)
])

# %%
###############################################################
###############################################################
# COMPUTE LIST OF DECAY SPECTRA
###############################################################
###############################################################

mb, mc = 0.14, 0.14

pTmax = 2
NpT = 200
files = glob.glob("data/spectra_real_sigma_constfield/*/spectr.txt")
files.sort() # sorts normally by alphabetical order
masses = [312,329,346,364,381,399,416,434,451,468,486,503,521,538,555,573,590,608,625,643,660,677,695,712,730,747,764,782,799,817,834,851,869,886,904,921,939,956,973,991]
masses = np.array(masses)/1000

for (ma,file) in list(zip(masses,files)):
    result = subprocess.run(args=[
        "./bin/decay",
        "--ma=%f"%(ma),
        "--mb=%f"%(mb),
        "--mc=%f"%(mc),
        "--pTmax=%f"%(pTmax),
        "--NpT=%d"%(NpT),
        "--primespecpath=%s"%(file)
    ])
    print(result)

newdir = "data/newspectra"
idx = 0
while(os.path.isdir(newdir+str(idx))):
    idx += 1
subprocess.run(args=["mkdir",newdir+str(idx)])

lastspecs = glob.glob("data/decayspec_????????_??????")
for spec in lastspecs:
    subprocess.run(args=["mv",spec,newdir+str(idx)])

#%%
###############################################################
###############################################################
# PLOT SINGLE DECAY SPECTRUM
###############################################################
###############################################################

plt.style.use("mplstyles/myclassic_white.mplstyle")

dir = "decayspec_20241009_170755"
path = "data/decaydecay_convtest/"+dir+"/"

fig_decayspec, ax_decayspec = plt.subplots(figsize=(7,7))
fig_primespec, ax_primespec = plt.subplots(figsize=(7,7))
fig_decayspec.suptitle(path+"decayspec")
fig_primespec.suptitle(path+"primespec")

df_ps = pd.read_csv(path+"primespec_interp.txt",comment="#")
ax_primespec.plot(df_ps["q"],df_ps["primespecRe"],marker="")

ax_primespec.set_yscale("log")
ax_primespec.set_ylabel(r"$\frac{1}{2\pi q_T}\frac{\mathrm{d}N}{\mathrm{d}q_T\mathrm{d}\eta_q}\ [\mathrm{GeV}^{-2}]$")
ax_primespec.set_xlabel(r"$q_T\ [\mathrm{GeV}]$")

ylims = ax_primespec.get_ylim()
ax_primespec.set_ylim((ylims[0],ylims[1]*100))

istart, iend = 0, 200
pad = 0.05
x0, y0 = 0.4, 0.45
dx, dy = 1-x0-pad, 1-y0-pad
axin = ax_primespec.inset_axes([x0, y0, dx, dy])
axin.plot(df_ps["q"][istart:iend],df_ps["primespecRe"][istart:iend],marker="")
axin.set_yscale("log")

df_ds = pd.read_csv(path+"decayspec.txt",comment="#")
ax_decayspec.plot(df_ds["p"].to_numpy(),df_ds["finalspecRe"].to_numpy(),marker="")

ax_decayspec.set_yscale("log")

ax_decayspec.set_ylabel(r"$\frac{1}{2\pi p_T}\frac{\mathrm{d}N}{\mathrm{d}p_T\mathrm{d}\eta_p}\ [\mathrm{GeV}^{-2}]$")
ax_decayspec.set_xlabel(r"$p_T\ [\mathrm{GeV}]$")

ax_decayspec.legend(loc=3,fontsize=15)

fig_decayspec.tight_layout()
fig_primespec.tight_layout()

fig_primespec.savefig("data/images/"+dir+"_primespec.png")
fig_decayspec.savefig("data/images/"+dir+"_decayspec.png")

plt.show()

# %%
###############################################################
###############################################################
# PLOT DECAY SPECTRA FOR CONVTEST DATA
###############################################################
###############################################################

plt.style.use("mplstyles/myclassic_white.mplstyle")

paths = glob.glob("data/decay_convtest/decayspec_*/")
paths.sort()

fig_decayspec, ax_decayspec = plt.subplots(figsize=(7,7))
fig_primespec, ax_primespec = plt.subplots(figsize=(7,7))
fig_decayspec.suptitle("data/decaydecay_convtest/*/decayspec_*")
fig_primespec.suptitle("data/decaydecay_convtest/*/primespec_*")

df_ps = pd.read_csv(paths[-1]+"primespec_interp.txt",comment="#")
ax_primespec.plot(df_ps["q"],df_ps["primespecRe"],marker="")

ax_primespec.set_yscale("log")
ax_primespec.set_ylabel(r"$\frac{1}{2\pi q_T}\frac{\mathrm{d}N}{\mathrm{d}q_T\mathrm{d}\eta_q}\ [\mathrm{GeV}^{-2}]$")
ax_primespec.set_xlabel(r"$q_T\ [\mathrm{GeV}]$")

ylims = ax_primespec.get_ylim()
ax_primespec.set_ylim((ylims[0],ylims[1]*100))

istart, iend = 0, 200
pad = 0.05
x0, y0 = 0.4, 0.45
dx, dy = 1-x0-pad, 1-y0-pad
axin = ax_primespec.inset_axes([x0, y0, dx, dy])
axin.plot(df_ps["q"][istart:iend],df_ps["primespecRe"][istart:iend],marker="")
axin.set_yscale("log")

for (i,path) in enumerate(paths):
    df_ps = pd.read_csv(path+"primespec_interp.txt",comment="#")
    qmax = max(df_ps["q"].to_numpy())

    df_ds = pd.read_csv(path+"decayspec.txt",comment="#")
    ax_decayspec.plot(df_ds["p"].to_numpy(),df_ds["finalspecRe"].to_numpy(),
                label=r"$q_{\mathrm{max}}=%.1f\,\mathrm{GeV}$"%(qmax),
                marker="")
    
    # t = i/(len(paths)-1)
    # col=(1-t,0,t)
    # ax_decayspec.plot(df_ds["p"].to_numpy(),df_ds["finalspecRe"].to_numpy(),
    #             label=r"$q_{\mathrm{max}}=%.1f\,\mathrm{GeV}$"%(qmax),
    #             marker="",c=col)

ax_decayspec.set_yscale("log")

ax_decayspec.set_ylabel(r"$\frac{1}{2\pi p_T}\frac{\mathrm{d}N}{\mathrm{d}p_T\mathrm{d}\eta_p}\ [\mathrm{GeV}^{-2}]$")
ax_decayspec.set_xlabel(r"$p_T\ [\mathrm{GeV}]$")

ax_decayspec.legend(loc=3,fontsize=15)

fig_decayspec.tight_layout()
fig_primespec.tight_layout()

# fig_primespec.savefig("data/images/decaydecay_convtest_primespec.png")
# fig_decayspec.savefig("data/images/decaydecay_convtest_decayspecs.png")

plt.show()

# %%
###############################################################
###############################################################
# PLOT DECAY SPECTRA FOR INITTEST DATA
###############################################################
###############################################################

plt.style.use("mplstyles/myclassic_white.mplstyle")

paths = glob.glob("data/decay_inittest/decayspec_*/")
paths.sort()

fig_decayspec, ax_decayspec = plt.subplots(figsize=(7,7))
fig_primespec, ax_primespec = plt.subplots(figsize=(7,7))
fig_decayspec.suptitle("data/decaydecay_inittest/*/decayspec_*")
fig_primespec.suptitle("data/decaydecay_inittest/*/primespec_*")

for (i,path) in enumerate(paths):
    t = i/(len(paths)-1)
    col=(1-t,0,t)

    df_ps = pd.read_csv(path+"primespec_interp.txt",comment="#")
    df_ds = pd.read_csv(path+"decayspec.txt",comment="#")

    ax_primespec.plot(df_ps["q"].to_numpy(),df_ps["primespecRe"].to_numpy(), marker="",c=col)
    ax_decayspec.plot(df_ds["p"].to_numpy(),df_ds["finalspecRe"].to_numpy(), marker="",c=col)

ax_primespec.set_yscale("log")
ax_decayspec.set_yscale("log")

ax_primespec.set_ylabel(r"$\frac{1}{2\pi q_T}\frac{\mathrm{d}N}{\mathrm{d}q_T\mathrm{d}\eta_p}\ [\mathrm{GeV}^{-2}]$")
ax_primespec.set_xlabel(r"$q_T\ [\mathrm{GeV}]$")

ax_decayspec.set_ylabel(r"$\frac{1}{2\pi p_T}\frac{\mathrm{d}N}{\mathrm{d}p_T\mathrm{d}\eta_p}\ [\mathrm{GeV}^{-2}]$")
ax_decayspec.set_xlabel(r"$p_T\ [\mathrm{GeV}]$")

fig_decayspec.tight_layout()
fig_primespec.tight_layout()

fig_primespec.savefig("data/images/decaydecay_inittest_primespec.png")
fig_decayspec.savefig("data/images/decaydecay_inittest_decayspecs.png")

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
# newdf.to_csv("./../../../code/Mathematica/data/ExampleFreezeOutCorrected_"+str(factor)+"x.csv",index=False)

# %%
#############################################################
########### COMPARE FLUIDUM WITH ALICE DATASETS #############
#############################################################

df_alice = pd.read_csv("./../../Mathematica/data/HEPData-ins1222333-v1-Table_1.csv",comment="#")

pTs_alice = df_alice["PT [GEV]"].to_numpy()[:41].astype(float)
spec_alice = df_alice["(1/Nev)*(1/(2*PI*PT))*D2(N)/DPT/DYRAP [GEV**-2]"].to_numpy()[:41].astype(float)
bins_lower = df_alice["PT [GEV] LOW"].to_numpy()[:41].astype(float)
bins_upper = df_alice["PT [GEV] HIGH"].to_numpy()[:41].astype(float)

df = pd.read_csv("./../../Mathematica/data/pionListBig.txt")
pTs_fluidum= np.array(df.keys()).astype(float)
spec_fluidum_pi0 = df.iloc[0].to_numpy()
spec_fluidum_piplus = df.iloc[1].to_numpy()

fig_spec, ax_spec = plt.subplots(figsize=(7,7))
fig_diff, ax_diff = plt.subplots(figsize=(7,7))

ax_spec.set_ylabel(r"$(2\pi p_T)^{-1}dN/(dp_Td\eta_p)\ [{GeV}^{-2}]$")
ax_diff.set_ylabel(r"$(2\pi p_T)^{-1}dN/(dp_Td\eta_p)\ [{GeV}^{-2}]$")

ax_spec.set_xlabel(r"$p_T\ [GeV]$")
ax_diff.set_xlabel(r"$p_T\ [GeV]$")

ax_spec.plot(pTs_fluidum, spec_fluidum_piplus,label="FluiduM",marker="",lw=2,c="r")
ax_spec.plot(pTs_alice, spec_alice,label="ALICE",marker="",lw=2,c="b")
ax_spec.set_yscale("log")

spec_fluidum_piplus = np.exp(scipy.interpolate.interp1d(pTs_fluidum, np.log(spec_fluidum_piplus))(pTs_alice))

ax_diff.plot(pTs_alice, np.abs(spec_alice - spec_fluidum_piplus),label=r"ALICE$-$FluiduM",marker="",c="b",lw=2)
ax_diff.set_yscale("log")

ax_spec.legend()
ax_diff.legend()

fig_spec.tight_layout()
fig_diff.tight_layout()

# fig_spec.savefig("data/images/FluidumAliceCompare.png")
# fig_diff.savefig("data/images/FluidumAliceDiff.png")

plt.show()

#%%

# spec = "data/spectra_taudep/spec_20241023_171739"
# spec = "data/spectra_taudep/spec_20241023_170656"
# spec = "data/spectra_taudep/spec_20241023_153952"

# spectra = glob.glob("data/spectra_real_constfield_20240822_161734/*")
# spec = spectra[5]

# spectra = glob.glob("data/spectra_real_constfield_20240822_161734_masses/*")
# spec = spectra[3]

# spectra = glob.glob("data/spectra_real_consteps_20240822_135426_masses/*")
# spec = spectra[0]

# spectra = glob.glob("data/spectra_real_m140_taudep_20241029_132506/*")
# spec = spectra[0]

# spectra = glob.glob("data/spectra_real_m280_taudep_20241029_135948/*/")
# spectra = glob.glob("data/spectra_real_m280_consteps_20241029_140005/*/")
# spectra = glob.glob("data/spectra_real_m280_constfield_20241029_140034/*/")
# spectra = glob.glob("data/spectra_real_m280_s-1_constfield_20241030_111128/*")
# spec = spectra[0]

# spectra = glob.glob("data/spectra_real_m140_consteps_20241029_110232/*/")
# spectra = glob.glob("data/spectra_real_m140_taudep_20241029_132506/*/")
spectra = glob.glob("data/spectra_real_m140_constfield_20241029_110519/*/")
# spectra = glob.glob("data/spectra_real_m140_constfield_20241029_162545/*/")
# spectra = glob.glob("data/spectra_real_m140_constfield_20241029_163722/*/")
# spectra = glob.glob("data/spectra_real_m140_s-1_consteps_20241030_104847/*")
spec = spectra[0]

# spectra = glob.glob("data/spectra_real_consteps_20240822_135440/*")
# spectra = glob.glob("data/spectra_real_consteps_20240826_143119_m226/*")
# spectra = glob.glob("data/spectra_real_consteps_20240826_143137_m312/*")
# spectra = glob.glob("data/spectra_real_consteps_20240826_143151_m398/*")
# spectra = glob.glob("data/spectra_real_consteps_20240826_143205_m484/*")
# spectra = glob.glob("data/spectra_real_consteps_20240826_143224_m570/*")
# spectra = glob.glob("data/spectra_real_consteps_20240826_143237_m656/*")
# spectra = glob.glob("data/spectra_real_consteps_20240826_143253_m742/*")
# spectra = glob.glob("data/spectra_real_consteps_20240826_143307_m828/*")
# spectra = glob.glob("data/spectra_real_consteps_20240826_143326_m914/*")
# spectra = glob.glob("data/spectra_real_consteps_20240826_143341_m1000/*")
# spectra = glob.glob("data/spectra_real_consteps_20240926_124618_m226_s-1/*")
# spectra = glob.glob("data/spectra_real_consteps_20240926_124625_m312_s-1/*")
# spectra = glob.glob("data/spectra_real_consteps_20240926_124635_m398_s-1/*")
# spectra = glob.glob("data/spectra_real_consteps_20240926_124646_m484_s-1/*")
# spectra = glob.glob("data/spectra_real_consteps_20240926_124658_m570_s-1/*")
# spectra = glob.glob("data/spectra_real_consteps_20240926_124709_m656_s-1/*")
# spectra = glob.glob("data/spectra_real_consteps_20240926_124716_m742_s-1/*")
# spectra = glob.glob("data/spectra_real_consteps_20240926_124724_m828_s-1/*")
# spectra = glob.glob("data/spectra_real_consteps_20240926_124737_m914_s-1/*")
# spectra = glob.glob("data/spectra_real_consteps_20240926_124748_m1000_s-1/*")
# spec = spectra[5]

FIGSIZE = (7,7)
AXISLABELSIZE = 20
TICKLABELSIZE = 15

FIELD_XLABEL = r"$\alpha$"
DFIELD_XLABEL = r"$\alpha$"
SPEC_XLABEL = r"$p_T\ [GeV]$"

FIELD_YLABEL = r"$\pi^0\ [GeV]$"
DFIELD_YLABEL = r"$n^\mu\partial_\mu\pi^0\ [GeV^2]$"
SPEC_YLABEL = r"$(2\pi p_T)^{-1}dN_{coherent}/(dp_Td\eta_p)\ [GeV^{-2}]$"

with open(spec+"/spectr.txt", "r") as file:
    for l in file.readlines():
        if(l[0]!="#"):
            break
        print(l.replace("\n",""))

df_field0 = pd.read_csv(spec+"/field0.txt",comment="#")
df_Dfield0 = pd.read_csv(spec+"/field0_deriv.txt",comment="#")
df_spec = pd.read_csv(spec+"/spectr.txt",comment="#")
df_specanti = pd.read_csv(spec+"/spectr_anti.txt",comment="#")

fig_init = plt.figure(figsize=FIGSIZE)
gs = gridspec.GridSpec(nrows=2, ncols=1, hspace=0)
ax_init1, ax_init2 = fig_init.add_subplot(gs[0]), fig_init.add_subplot(gs[1])

fig_spec, ax_spec = plt.subplots(figsize=(7,7))
fig_spec.subplots_adjust(left=0.25, bottom=0.25)

x_spec = df_spec["pT"].to_numpy()
x_specanti = df_specanti["pT"].to_numpy()

y_spec = df_spec["abs2Re"].to_numpy()
y_specanti = df_specanti["abs2Re"].to_numpy()

ax_init1.plot(
        df_field0["alpha"].to_numpy(),
        df_field0["field0Re"].to_numpy(),
        marker="",
        lw=2)

ax_init2.plot(
        df_Dfield0["alpha"].to_numpy(),
        df_Dfield0["Dfield0Re"].to_numpy(),
        marker="",
        lw=2)

ax_init1.plot(
        df_field0["alpha"].to_numpy(),
        df_field0["field0Im"].to_numpy(),
        ls="-.",
        lw=2,
        marker="")

ax_init2.plot(
        df_Dfield0["alpha"].to_numpy(),
        df_Dfield0["Dfield0Im"].to_numpy(),
        ls="-.",
        lw=2,
        marker="")

ax_spec.plot(ps, np.abs(Jps)**2)

line_spec, = ax_spec.plot(x_spec, y_spec,
    label=r"model (particle)",
    marker="",
    lw=2)

line_specanti, = ax_spec.plot(x_specanti, y_specanti,
    label=r"model (antiparticle)",
    marker="",
    lw=2,
    ls="-.")

ax_spec.plot(pTs_alice, spec_alice - spec_fluidum_piplus,
    label=r"ALICE$-$FluiduM",
    marker="",
    lw=2)

ax_init1.set_ylabel(FIELD_YLABEL, fontsize=AXISLABELSIZE)
ax_init2.set_ylabel(DFIELD_YLABEL, fontsize=AXISLABELSIZE)
ax_spec.set_ylabel(SPEC_YLABEL, fontsize=AXISLABELSIZE)

ax_init1.set_xlabel(FIELD_XLABEL, fontsize=AXISLABELSIZE)
ax_init2.set_xlabel(DFIELD_XLABEL, fontsize=AXISLABELSIZE)
ax_spec.set_xlabel(SPEC_XLABEL, fontsize=AXISLABELSIZE)

# ax_init1.set_xticklabels(ax_init1.get_xticklabels(), visible=False)
# labels = ax_init1.get_yticklabels()
# labels[0] = labels[-1] = ""
# ax_init1.set_yticklabels(labels)
# labels = ax_init2.get_yticklabels()
# labels[0] = labels[-1] = ""
# ax_init2.set_yticklabels(labels)

xticks = [0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2]
xticklabels = [r"$0$",r"$\pi/8$",r"$\pi/4$",r"$3\pi/8$",r"$\pi/2$"]
ax_init1.set_xticks(xticks, xticklabels)
ax_init2.set_xticks(xticks, xticklabels)

ax_init1.tick_params(axis="both",labelsize=TICKLABELSIZE)
ax_init2.tick_params(axis="both",labelsize=TICKLABELSIZE)
ax_spec.tick_params(axis="both",labelsize=TICKLABELSIZE)

ax_init1.set_xlim(0,np.pi/2)
ax_init2.set_xlim(0,np.pi/2)

ax_init1.grid(False,which="both")
ax_init2.grid(False,which="both")
ax_spec.grid(False,which="both")

ax_spec.set_yscale("log")
ax_spec.legend()

fig_init.tight_layout()
fig_spec.tight_layout()

# Make a horizontal slider
ax_sx = fig_spec.add_axes([0.25, 0.1, 0.65, 0.03])
sx_slider = Slider(
    ax=ax_sx,
    label='scale x',
    valmin=0.5,
    valmax=2,
    valinit=1,
)

# Make a vertically oriented slider
ax_sy = fig_spec.add_axes([0.1, 0.25, 0.0225, 0.63])
sy_slider = Slider(
    ax=ax_sy,
    label="ln(scale y)",
    valmin=-5,
    valmax=5,
    valinit=0,
    orientation="vertical"
)

# The function to be called anytime a slider's value changes
def update(val):
    line_spec.set_xdata(sx_slider.val * x_spec)
    line_spec.set_ydata(np.exp(sy_slider.val) * y_spec)
    line_specanti.set_xdata(sx_slider.val * x_specanti)
    line_specanti.set_ydata(np.exp(sy_slider.val) * y_specanti)
    fig_spec.canvas.draw_idle()

sx_slider.on_changed(update)
sy_slider.on_changed(update)

plt.show()