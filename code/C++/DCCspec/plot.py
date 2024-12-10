#%%
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from IPython import get_ipython
import glob
import scipy.interpolate
from matplotlib import gridspec
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
import scipy.optimize

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

pT_compare = pTs_alice
spec_compare = np.abs(spec_alice - spec_fluidum_piplus)

ax_diff.plot(pT_compare, spec_compare,label=r"ALICE$-$FluiduM",marker="",c="b",lw=2)
ax_diff.set_yscale("log")

ax_spec.legend(fontsize=LEGENDSIZE,fancybox=True, framealpha=0.85,shadow=False)
ax_diff.legend(fontsize=LEGENDSIZE,fancybox=True, framealpha=0.85,shadow=False)

ax_spec.grid(False, which="both")
ax_diff.grid(False, which="both")

fig_spec.tight_layout()
fig_diff.tight_layout()

# fig_spec.savefig(dpi=150,"images/FluidumAliceCompare.png")
# fig_diff.savefig(dpi=150,"images/FluidumAliceDiff.png")

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
# plt.ioff()
plt.ion()

# parentdir = "data/realfield_inittest/*/"
# parentdir = "data/specsreal_v3/*/"
# parentdir = "data/realfield_inittest_v2/*/"
# parentdir = "data_old/spectra_real_consteps_20240822_135426/*/"
# parentdir = "data/spectra_real_consteps_20240822_135440/*/"
# parentdir = "data/spectra_real_constfield_20240822_161734/*/"
# parentdir = "data/spectra_real_constfield_20240822_161738/*/"
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
# ms = 140*np.ones(len(paths))
# ms = 226*np.ones(len(paths))
# ms = 312*np.ones(len(paths))
# ms = 398*np.ones(len(paths))
# ms = 484*np.ones(len(paths))
# ms = 570*np.ones(len(paths))
# ms = 656*np.ones(len(paths))
# ms = 742*np.ones(len(paths))
# ms = 828*np.ones(len(paths))
# ms = 914*np.ones(len(paths))
# ms = 1000*np.ones(len(paths))

# parentdir = "data_old/spectra_real_constfield_20240822_161734_masses/*/"
# paths = glob.glob(parentdir)
# ms = np.linspace(140,800,10,endpoint=True)

# spectra = glob.glob("data/spectra_real_m140_consteps_20241029_110232/*/")
# spectra = glob.glob("data/spectra_real_m140_taudep_20241029_132506/*/")
# spectra = glob.glob("data/spectra_real_m140_constfield_20241029_110519/*/")
# spectra = glob.glob("data/spectra_real_m140_constfield_20241029_162545/*/")
# spectra = glob.glob("data/spectra_real_m140_constfield_20241029_163722/*/")
# spectra = glob.glob("data/spectra_real_m140_s-1_constfield_20241030_102134/*/")
# spec = spectra[0]

# parentdir = "data/spectra_real_constfield_20240822_161734/*/"
# paths = glob.glob(parentdir)
# ms = 140*np.ones(len(paths))

# parentdir = "data/spectra_real_m140_consteps_20241029_110232/*/"
parentdir = "data/spectra_real_m140_taudep_20241029_132506/*/"
# parentdir = "data/spectra_real_m140_constfield_20241029_110519/*/"
# parentdir = "data/spectra_real_m140_constfield_20241029_162545/*/"
# parentdir = "data/spectra_real_m140_constfield_20241029_163722/*/"
# parentdir = "data/spectra_real_m140_s-1_constfield_20241030_102134/*/"
# parentdir = "data/spectra_real_m140_s-1_consteps_20241030_104847/*/"
paths = glob.glob(parentdir)
paths.sort()
paths.reverse()
ms = 140*np.ones(len(paths))

# parentdir = "data/spectra_real_m280_taudep_20241029_135948/*/"
# parentdir = "data/spectra_real_m280_consteps_20241029_140005/*/"
# parentdir = "data/spectra_real_m280_constfield_20241029_140034/*/"
# parentdir = "data/spectra_real_m280_s-1_constfield_20241030_111128/*/"
# paths = glob.glob(parentdir)
# ms = 280*np.ones(len(paths))

## epsconstmasses
# folders = [
#  'data_old/spectra_real_consteps_20240822_135426',
#  'data_old/spectra_real_consteps_20240826_143119_m226',
#  'data_old/spectra_real_consteps_20240826_143137_m312',
#  'data_old/spectra_real_consteps_20240826_143151_m398',
#  'data_old/spectra_real_consteps_20240826_143205_m484',
#  'data_old/spectra_real_consteps_20240826_143224_m570',
#  'data_old/spectra_real_consteps_20240826_143237_m656',
#  'data_old/spectra_real_consteps_20240826_143253_m742',
#  'data_old/spectra_real_consteps_20240826_143307_m828',
#  'data_old/spectra_real_consteps_20240826_143326_m914',
#  'data_old/spectra_real_consteps_20240826_143341_m1000']

# ms = [
#     140.0, 
#     226.0, 
#     312.0, 
#     398.0, 
#     484.0, 
#     570.0, 
#     656.0, 
#     742.0, 
#     828.0, 
#     914.0, 
#     1000.0]

# folders.sort()
# paths = []
# for folder in folders:
#     localpaths = glob.glob(folder+"/*")
#     paths.append(localpaths[0])

# parentdir = "epsconst_masses"
### \epsconstmasses

# paths = [
#     glob.glob("data/spectra_real_m140_constfield_20241029_110519/*")[0],
#     glob.glob("data/spectra_real_m140_consteps_20241029_110232/*")[0],
#     # glob.glob("data/spectra_real_m140_s-1_consteps_20241030_104847/*")[5],
#     glob.glob("data/spectra_real_m140_taudep_20241029_132506/*")[0]
# ]
# parentdir = "m140_compareICs"

# paths = [
#     glob.glob("data/spectra_real_m280_constfield_20241029_140034/*")[5],
#     glob.glob("data/spectra_real_m280_consteps_20241029_140005/*")[0],
#     # glob.glob("data/spectra_real_m280_s-1_consteps_20241030_111119/*")[5],
#     glob.glob("data/spectra_real_m280_taudep_20241030_150721/*")[2]
# ]
# parentdir = "m280_compareICs"

# SAVE = True
SAVE = False
SAVETITLE = parentdir.replace("data/","").replace("/*/","").replace("data_old/","").replace("/*/","")

# COMPARE = True
COMPARE = False

FIELD_XLABEL = r"$\alpha$"
DFIELD_XLABEL = r"$\alpha$"
SPEC_XLABEL = r"$p_T\ [GeV]$"

FIELD_YLABEL = r"$\pi^0\ [GeV]$"
DFIELD_YLABEL = r"$n^\mu\partial_\mu\pi^0\ [GeV]$"
# SPEC_YLABEL = r"$\mathcal{N}\cdot\frac{1}{2\pi p_T}\frac{dN_{coherent}}{dp_Td\eta_p}\ [GeV^{-2}]$"
SPEC_YLABEL = r"$(2\pi p_T)^{-1}dN_{coherent}/(dp_Td\eta_p)\ [GeV^{-2}]$"

# SCALEBYMASS = True
SCALEBYMASS = False
if(SCALEBYMASS):
    FIELD_YLABEL = r"$m\cdot\pi^0\ [GeV^2]$"
# NORMALIZE_SPECS = True
NORMALIZE_SPECS = False
if(NORMALIZE_SPECS):
    SPEC_YLABEL = r"$\mathcal{N}\cdot(2\pi p_T)^{-1}dN_{coherent}/(dp_Td\eta_p)\ [GeV^{-2}]$"

# TOBINS = True
TOBINS = False

FIGSIZE = (7,7)
ADJUSTLABELS = False
AXISLABELSIZE = 20
TICKLABELSIZE = 15
LINEWIDTH = 2
LEGENDSIZE = 15

# LINECOLLECTION = True
LINECOLLECTION = False
WITHCMAP = True
# WITHCMAP = False
CMAP = LinearSegmentedColormap.from_list("custom", ["blue","red"])
CMAP_LBWH = [0.025, 0.025, 0.05, 0.45]
CMAP_LABELSIZE = 15
CMAP_TICKSIZE = 15
LC_LABEL = r"$m\ [GeV]$"
LC_ARRAY = np.array(ms)/1000
# LC_ARRAY = np.arange(len(paths))
lines_init1 = []
lines_init2 = []
lines_spec = []

COMPLEX_FIELD = False

LEGEND = True
# LEGEND = False
# labels = [r"$m=%i\ MeV$"%(m) for m in ms]
labels = [r"$a=1/3$",r"$a=1/2$",r"$a=1$",r"$a=2$",r"$a=3$"]
# labels.reverse()
# labels = ["1","2","3","4"]
# labels = [r"Type A, $m=140\ MeV$",r"Type B, $m=140\ MeV$",r"Type C, $m=140\ MeV$"]
# labels = [r"Type A, $m=280\ MeV$",r"Type B, $m=280\ MeV$",r"Type C, $m=280\ MeV$"]
if(not LEGEND):
    labels = ["" for n in range(len(paths))]

# cols = [(t, 0, 1-t) for t in np.linspace(0,1,len(paths),endpoint=True)]
cols = ["r", "g", "b","purple", "orange"]

fig_init = plt.figure(figsize=FIGSIZE)
gs = gridspec.GridSpec(nrows=2, ncols=1, hspace=0)
ax_init1, ax_init2 = fig_init.add_subplot(gs[0]), fig_init.add_subplot(gs[1])

fig_spec, ax_spec = plt.subplots(figsize=FIGSIZE)

minval1, maxval1 = np.inf, -np.inf
minval2, maxval2 = np.inf, -np.inf

for (n,path) in enumerate(paths):
    df_spec = pd.read_csv(path+"/spectr.txt",comment="#")
    df_specanti = pd.read_csv(path+"/spectr_anti.txt",comment="#")
    df_field0 = pd.read_csv(path+"/field0.txt",comment="#")
    df_Dfield0 = pd.read_csv(path+"/field0_deriv.txt",comment="#")

    col = cols[n]

    scale_specs = 1
    if(NORMALIZE_SPECS):
        scale_specs = 5000/df_spec["abs2Re"].to_numpy()[0]
    m = ms[n]/1000

    x_init1 = df_field0["alpha"].to_numpy()
    x_init2 = df_Dfield0["alpha"].to_numpy()
    x_spec = df_spec["pT"].to_numpy()

    y_init1 = df_field0["field0Re"].to_numpy()
    if(SCALEBYMASS):
        y_init1 = m * df_field0["field0Re"].to_numpy()
    y_init2 = df_Dfield0["Dfield0Re"].to_numpy()
    y_spec = scale_specs * df_spec["abs2Re"].to_numpy()

    mindata1 = np.min(y_init1)
    maxdata1 = np.max(y_init1)
    minval1 = np.min((minval1,mindata1))
    maxval1 = np.max((maxval1,maxdata1))

    mindata2 = np.min(y_init2)
    maxdata2 = np.max(y_init2)
    minval2 = np.min((minval2,mindata2))
    maxval2 = np.max((maxval2,maxdata2))

    ax_init1.plot(
        x_init1,
        y_init1,
        marker="",
        color=col,
        lw=LINEWIDTH)

    ax_init2.plot(
        x_init2,
        y_init2,
        marker="",
        color=col,
        lw=LINEWIDTH)

    if(not TOBINS):
        ax_spec.plot(
            x_spec,
            y_spec,
            marker="",
            color=col,
            label=labels[n],
            lw=LINEWIDTH)
    
    if(TOBINS):
        idx = np.where(x_spec >= bins_lower[0])[0][0]

        ax_spec.plot(
            x_spec[:idx],
            y_spec[:idx],
            marker="",
            color=col,
            label=labels[n],
            lw=LINEWIDTH)

        x_spec,y_spec = data_to_bins(x_spec[idx:], y_spec[idx:], bins_lower, bins_upper) 

        ax_spec.plot(
            x_spec,
            y_spec,
            color=col,
            marker="",
            ls = "--",
            # label=labels[n],
            lw=LINEWIDTH)


    if LINECOLLECTION:
        line_init1 = np.column_stack((x_init1, y_init1))
        line_init2 = np.column_stack((x_init2, y_init2))
        line_spec = np.column_stack((x_spec, y_spec))

        lines_spec.append(line_spec)
        lines_init1.append(line_init1)
        lines_init2.append(line_init2)    
    
    if(COMPLEX_FIELD):
        ax_init1.plot(
            df_field0["alpha"].to_numpy(),
            m*df_field0["field0Im"].to_numpy(),
            marker="",
            ls="-.",
            color=col)

        ax_init2.plot(
            df_Dfield0["alpha"].to_numpy(),
            df_Dfield0["Dfield0Im"].to_numpy(),
            marker="",
            ls="-.",
            color=col)

        ax_spec.plot(
            df_specanti["pT"].to_numpy(),
            scale*df_specanti["abs2Re"].to_numpy(),
            marker="",
            ls="-.",
            color=col)

if(LINECOLLECTION):
        # print(float(file.readlines()[3].replace("\n","").split("\t")[1]))
    ax_init1.clear()
    ax_init2.clear()
    ax_spec.clear()

    lc_init1 = LineCollection(lines_init1,array=LC_ARRAY,cmap=CMAP)
    lc_init2 = LineCollection(lines_init2,array=LC_ARRAY,cmap=CMAP)
    lc_spec = LineCollection(lines_spec,array=LC_ARRAY,cmap=CMAP)

    ax_init1.add_collection(lc_init1)
    ax_init2.add_collection(lc_init2)
    ax_spec.add_collection(lc_spec)

    ax_init1.autoscale_view()
    ax_init2.autoscale_view()
    ax_spec.autoscale_view()

    if(WITHCMAP):
        cax = ax_spec.inset_axes(CMAP_LBWH)
        cbar = fig_spec.colorbar(lc_spec, cax=cax)
        cbar.set_label(LC_LABEL, fontsize=CMAP_LABELSIZE)
        cbar.ax.tick_params(labelsize=CMAP_TICKSIZE)

        # fig_spec.colorbar(lc_spec,label=LC_LABEL)

ax_spec.set_yscale("log")

ax_init1.set_ylabel(FIELD_YLABEL, fontsize=AXISLABELSIZE)
ax_init2.set_ylabel(DFIELD_YLABEL, fontsize=AXISLABELSIZE)
ax_spec.set_ylabel(SPEC_YLABEL, fontsize=AXISLABELSIZE)

ax_init1.set_xlabel(FIELD_XLABEL, fontsize=AXISLABELSIZE)
ax_init2.set_xlabel(DFIELD_XLABEL, fontsize=AXISLABELSIZE)
ax_spec.set_xlabel(SPEC_XLABEL, fontsize=AXISLABELSIZE)

myrange1 = maxval1 - minval1
myrange2 = maxval2 - minval2
if(myrange1 > 0):
    ax_init1.set_ylim(minval1-0.05*myrange1,maxval1+0.05*myrange1)
if(myrange2 > 0):
    ax_init2.set_ylim(minval2-0.05*myrange2,maxval2+0.05*myrange2)

ax_init1.set_xticklabels(ax_init1.get_xticklabels(), visible=False)

if(ADJUSTLABELS):
    labels = ax_init1.get_yticklabels()
    labels[0] = labels[-1] = ""
    ax_init1.set_yticklabels(labels)
    labels = ax_init2.get_yticklabels()
    labels[0] = labels[-1] = ""
    ax_init2.set_yticklabels(labels)

xticks = [0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2]
xticklabels = [r"$0$",r"$\pi/8$",r"$\pi/4$",r"$3\pi/8$",r"$\pi/2$"]
ax_init1.set_xticks([],[])
ax_init2.set_xticks(xticks, xticklabels)

ax_init1.tick_params(axis="both",labelsize=TICKLABELSIZE)
ax_init2.tick_params(axis="both",labelsize=TICKLABELSIZE)
ax_spec.tick_params(axis="both",labelsize=TICKLABELSIZE)

ax_init1.xaxis.set_ticks_position("bottom")
ax_init1.yaxis.set_ticks_position("left")
ax_init2.xaxis.set_ticks_position("bottom")
ax_init2.yaxis.set_ticks_position("left")
ax_spec.xaxis.set_ticks_position("bottom")
ax_spec.yaxis.set_ticks_position("left")

ax_init1.set_xlim(0,np.pi/2)
ax_init2.set_xlim(0,np.pi/2)
ax_spec.set_xlim(0,2)

ax_init1.grid(False,which="both")
ax_init2.grid(False,which="both")
ax_spec.grid(False,which="both")

if(COMPARE):
    ax_spec.plot(pT_compare, spec_compare,lw=2,label=r"ALICE$-$FluiduM",c="k",ls="--",marker="")
    ax_spec.legend(fontsize=LEGENDSIZE,fancybox=True, framealpha=0.85,shadow=False)

if(LEGEND):
    ax_spec.legend(fontsize=LEGENDSIZE,fancybox=True, framealpha=0.85,shadow=False)

# fig_init.suptitle(parentdir)
# fig_spec.suptitle(parentdir)

fig_init.tight_layout()
fig_spec.tight_layout()

if(SAVE):
    fig_init.savefig("images/"+SAVETITLE+"_init.png",dpi=150)
    fig_spec.savefig("images/"+SAVETITLE+"_spec.png",dpi=150)

# fig_init.show()
# fig_spec.show()

plt.show()

# %%
###############################################################
###############################################################
# PLOT DECAY SPECTRA FOR INITTEST DATA
###############################################################
###############################################################

paths = glob.glob("data_old/decay_inittest/decayspec_*/")
paths.sort()

fig_decayspec, ax_decayspec = plt.subplots(figsize=(7,7))
fig_primespec, ax_primespec = plt.subplots(figsize=(7,7))
# fig_decayspec.suptitle("data/decaydecay_inittest/*/decayspec_*")
# fig_primespec.suptitle("data/decaydecay_inittest/*/primespec_*")

lines_ps = []
lines_ds = []

SPEC_XLABEL_PS = r"$q_T\ [GeV]$"
SPEC_XLABEL_DS = r"$p_T\ [GeV]$"
SPEC_YLABEL_PS = r"$(2\pi q_T)^{-1}dN_{primary}/(dq_Td\eta_q)\ [GeV^{-2}]$"
SPEC_YLABEL_DS = r"$(2\pi p_T)^{-1}dN_{decay}/(dp_Td\eta_p)\ [GeV^{-2}]$"

LC_ARRAY = np.arange(len(paths))

cols = [(t, 0, 1-t) for t in np.linspace(0,1,len(paths),endpoint=True)]
for (n,path) in enumerate(paths):
    df_ps = pd.read_csv(path+"primespec_interp.txt",comment="#")
    df_ds = pd.read_csv(path+"decayspec.txt",comment="#")

    col = cols[n]

    x_ps = df_ps["q"].to_numpy()
    y_ps = df_ps["primespecRe"].to_numpy()
    x_ds = df_ds["p"].to_numpy()
    y_ds = df_ds["finalspecRe"].to_numpy()

    # ax_primespec.plot(x_ps, y_ps, marker="",c=col)
    # ax_decayspec.plot(x_ds, y_ds, marker="",c=col)
        
    line_ps = np.column_stack((x_ps, y_ps))
    line_ds = np.column_stack((x_ds, y_ds))

    lines_ps.append(line_ps)
    lines_ds.append(line_ds)

lc_ps = LineCollection(lines_ps,array=LC_ARRAY,cmap=CMAP)
lc_ds = LineCollection(lines_ds,array=LC_ARRAY,cmap=CMAP)

ax_primespec.add_collection(lc_ps)
ax_decayspec.add_collection(lc_ds)

ax_primespec.autoscale_view()
ax_decayspec.autoscale_view()

ax_primespec.set_yscale("log")
ax_decayspec.set_yscale("log")

ax_primespec.grid(False, which="both")
ax_decayspec.grid(False, which="both")

ax_primespec.tick_params(axis="both",labelsize=TICKLABELSIZE)
ax_decayspec.tick_params(axis="both",labelsize=TICKLABELSIZE)

ax_primespec.xaxis.set_ticks_position("bottom")
ax_primespec.yaxis.set_ticks_position("left")
ax_decayspec.xaxis.set_ticks_position("bottom")
ax_decayspec.yaxis.set_ticks_position("left")

ax_primespec.set_ylabel(SPEC_YLABEL_PS)
ax_primespec.set_xlabel(SPEC_XLABEL_PS)

ax_decayspec.set_ylabel(SPEC_YLABEL_DS)
ax_decayspec.set_xlabel(SPEC_XLABEL_DS)

fig_decayspec.tight_layout()
fig_primespec.tight_layout()

fig_primespec.savefig("images/decaydecay_inittest_primespec.png")
fig_decayspec.savefig("images/decaydecay_inittest_decayspecs.png")

plt.show()

# %%
###############################################################
###############################################################
# PLOT FREEZEOUT GEOMETRY
###############################################################
###############################################################

ARROW_WIDTH = 0.3
ARROW_LENGTH = 0.2
ARROW_OVERHANG = 0.1
rmax = 11
taumax = 15

path = "./../../../code/Mathematica/data/ExampleFreezeOutCorrected.csv"
df = pd.read_csv(path)

fig, ax = plt.subplots(figsize=(7,7))

ax.plot(df["r"].to_numpy(), df["tau"].to_numpy(),
        marker="",
        label="exact freezeout geometry")
ax.hlines(0.4,0,10,color="k",ls="--",lw=2)
ax.hlines(12.5,0,8.5,color="b",lw=2,ls="--")
ax.vlines(8.5,0.4,12.5,color="b",lw=2,ls="--")

ax.set_xlim(-ARROW_WIDTH,rmax)
ax.set_ylim(-ARROW_WIDTH,taumax)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)

# Set bottom and left spines as x and y axes of coordinate system
ax.spines['bottom'].set_position('zero')
ax.spines['left'].set_position('zero')

ax.grid(False)
ax.set_xticks([8.5])
ax.set_yticks([0.4,12.5])
ax.set_xticklabels([r"$\mathcal{R}$"])
ax.set_yticklabels([r"$\tau_0$",r"$\mathcal{T}$"])
ax.xaxis.set_tick_params(width=2,length=7,direction='out')
ax.yaxis.set_tick_params(width=2,length=7,direction='out')
ax.xaxis.set_ticks_position("bottom")
ax.yaxis.set_ticks_position("left")

ax.set_xlabel(r"$r$",rotation=0,)
ax.set_ylabel(r"$\tau$",rotation=0)

ax.xaxis.set_label_coords(0.95,0)
ax.yaxis.set_label_coords(-0.03,0.95)

# ax.plot([0,8.5/2],[0.4,12.5],c="grey",lw=1,marker="")
# ax.plot([0,8.5],[0.4,12.5],c="grey",lw=2,marker="")

# Draw arrows
ax.axes.arrow(0,0,rmax,0,
              color="k",
              head_width=ARROW_WIDTH,
              head_length=ARROW_LENGTH,
              overhang = ARROW_OVERHANG,
              length_includes_head=True)
ax.axes.arrow(0,0,0,taumax,
              color="k",
              head_width=ARROW_WIDTH*rmax/taumax,
              head_length=ARROW_LENGTH*taumax/rmax,
              overhang = ARROW_OVERHANG,
              length_includes_head=True)

# ax.legend()

fig.tight_layout()

fig.savefig("images/freezeoutsurface_simplified.png")

plt.show()

#%%
###############################################################
###############################################################
# PLOT DECAY SPECTRA FOR SIGMA RESONANCE
###############################################################
###############################################################

# SAVE = True
SAVE = False
SAVETITLE = "decay_sigmaresonance_pi_m140"

COMPARE = True
# COMPARE = False

PSPEC_XLABEL = r"$q_T\ [GeV]$"
DSPEC_XLABEL = r"$p_T\ [GeV]$"

PSPEC_YLABEL = r"$(2\pi q_T)^{-1}dN_{\sigma}/(dq_Td\eta_q)\ [GeV^{-2}]$"
DSPEC_YLABEL = r"$(2\pi p_T)^{-1}dN_{\pi}/(dp_Td\eta_p)\ [GeV^{-2}]$"

FIGSIZE = (7,7)
AXISLABELSIZE = 20
TICKLABELSIZE = 15
LINEWIDTH = 2
LEGENDSIZE = 10

CMAP = LinearSegmentedColormap.from_list("custom", ["blue","red"])
CMAP_LBWH = [0.025, 0.025, 0.05, 0.45]
CMAP_LABELSIZE = 15
CMAP_TICKSIZE = 15
LC_LABEL = r"$m\ [GeV]$"

# folders = glob.glob("data/decayspectra_sigma_constfield_pi_m280/decayspec*")
folders = glob.glob("data/newspectra1/decayspec*")
# folders = glob.glob("data/decayspectra_sigma_20241119_171318_pi_m140/decayspec*")
# folders = glob.glob("data/decayspectra_sigma_20241119_171318_pi_m280/decayspec*")
# folders = glob.glob("data/decayspectra_sigma_p5GeV_20241119_171318_pi_m280/decayspec*")
folders.sort()

masses = []
for f in folders:
    with open(f+"/decayspec.txt") as file:
        masses.append(float(file.readlines()[3].replace("\n","").split("\t")[1]))
masses = 1000*np.array(masses)

mpi = 140

Mpole = 479
Gpole = 2* 279

msigma = np.sqrt(1/4 * (16 * mpi**2 + 
                        np.sqrt(16 * Gpole**2 * Mpole**2 + 
                                (-16 * mpi**2 - Gpole**2 + 4*Mpole**2)**2)))
Gam = np.sqrt(1/2 * (16 * mpi**2 + Gpole**2 - 
     4*  Mpole**2 + np.sqrt(16*Gpole**2 * Mpole**2 + (-16 * mpi**2 - Gpole**2 + 4*Mpole**2)**2)))

def Delta(s):
    return 1/(s-msigma**2+1j*Gam*np.sqrt(s-(2*mpi)**2))

def S(k):
    return -1/np.pi*np.imag(Delta(k**2))*(masses<=1000)

weights = 2*masses*S(masses)*np.ptp(masses)/len(masses)
weights /= np.sum(weights)

fig_weights, ax_weights = plt.subplots(figsize=(7,7))
ax_weights.plot(masses,weights)

ax_weights.legend(fontsize=LEGENDSIZE,fancybox=True, framealpha=0.85,shadow=False)
ax_weights.tick_params(axis="both",labelsize=TICKLABELSIZE)
ax_weights.xaxis.set_ticks_position("bottom")
ax_weights.yaxis.set_ticks_position("left")
ax_weights.set_xlabel(r"$\mu$")
ax_weights.set_ylabel(r"$2\mu\text{d}\mu\rho(\mu^2)$")
ax_weights.grid(False,which="both")
fig_weights.tight_layout()


lines_ps, lines_ds = [], []

fig_ps, ax_ps = plt.subplots(figsize=(7,7))
fig_ds, ax_ds = plt.subplots(figsize=(7,7))
fig_fs, ax_fs = plt.subplots(figsize=(7,7))

x_fs, y_fs = np.zeros(shape=(2,1))

for (idx,(f,m)) in enumerate(list(zip(folders,masses))):
    df_ps = pd.read_csv(f+"/primespec_interp.txt",comment="#")
    df_ds = pd.read_csv(f+"/decayspec.txt",comment="#")

    x_ps, y_ps = df_ps["q"].to_numpy(),df_ps["primespecRe"].to_numpy()
    x_ds, y_ds = df_ds["p"].to_numpy(),df_ds["finalspecRe"].to_numpy()

    line_ps = np.column_stack((x_ps, y_ps))
    line_ds = np.column_stack((
        x_ds[(~np.isnan(y_ds))*(~np.isinf(y_ds))],
        y_ds[(~np.isnan(y_ds))*(~np.isinf(y_ds))]))

    y_ds[np.isnan(y_ds)+np.isinf(y_ds)] = 0
    y_ds[np.isnan(y_ds)+np.isinf(y_ds)] = 0

    x_fs = x_ds
    y_fs = y_fs + y_ds * weights[idx]

    lines_ps.append(line_ps)
    lines_ds.append(line_ds)

x_fs = x_fs[(~np.isnan(y_fs))*(~np.isinf(y_fs))]
y_fs = y_fs[(~np.isnan(y_fs))*(~np.isinf(y_fs))]
ax_fs.plot(x_fs, y_fs,lw=2,c="b",marker="",label=r"DCC decay spectrum $\sigma\to\pi\pi$")
ax_fs.plot(pT_compare, spec_compare,lw=2,c="k",ls="--",marker="",label=r"ALICE$-$FluiduM")

lc_ps = LineCollection(lines_ps,array=masses,cmap=CMAP)
lc_ds = LineCollection(lines_ds,array=masses,cmap=CMAP)

ax_ps.add_collection(lc_ps)
ax_ds.add_collection(lc_ds)
ax_ps.autoscale_view()
ax_ds.autoscale_view()

# STYLE

ax_fs.legend(fontsize=LEGENDSIZE,fancybox=True, framealpha=0.85,shadow=False)

cax = ax_ds.inset_axes(CMAP_LBWH)
cbar = fig_spec.colorbar(lc_ds, cax=cax)
cbar.set_label(LC_LABEL, fontsize=CMAP_LABELSIZE)
cbar.ax.tick_params(labelsize=CMAP_TICKSIZE)

ax_ps.tick_params(axis="both",labelsize=TICKLABELSIZE)
ax_ds.tick_params(axis="both",labelsize=TICKLABELSIZE)
ax_fs.tick_params(axis="both",labelsize=TICKLABELSIZE)

ax_ps.xaxis.set_ticks_position("bottom")
ax_ds.xaxis.set_ticks_position("bottom")
ax_fs.xaxis.set_ticks_position("bottom")
ax_ps.yaxis.set_ticks_position("left")
ax_ds.yaxis.set_ticks_position("left")
ax_fs.yaxis.set_ticks_position("left")

ax_ps.set_yscale("log")
ax_ds.set_yscale("log")
ax_fs.set_yscale("log")

ax_ps.set_xlabel(PSPEC_XLABEL)
ax_ps.set_ylabel(PSPEC_YLABEL)
ax_ds.set_xlabel(DSPEC_XLABEL)
ax_ds.set_ylabel(DSPEC_YLABEL)
ax_fs.set_xlabel(DSPEC_XLABEL)
ax_fs.set_ylabel(DSPEC_YLABEL)

ax_ds.grid(False,which="both")
ax_ps.grid(False,which="both")
ax_fs.grid(False,which="both")

fig_ds.tight_layout()
fig_ps.tight_layout()
fig_fs.tight_layout()

if(SAVE):
    print("SAVED!")
    fig_ps.savefig("images/"+SAVETITLE+"_primespec.png")
    fig_ds.savefig("images/"+SAVETITLE+"_decayspecs.png")
    fig_fs.savefig("images/"+SAVETITLE+"_decayspec_weighted.png")
else:
    print("NOT SAVED!")


#%%
###############################################################
###############################################################
# DISCRETIZE SIGMA RESONANCE
###############################################################
###############################################################

# df_sig600 = pd.read_csv("sigmaresonance600.txt",header=None,sep=" ")

# m_sig = df_sig600.iloc[:,0].to_numpy()

# imax = np.where(m_sig > 1000)[0][0]

# m_sig = m_sig[:imax]
# w_sig = df_sig600.iloc[:,1].to_numpy()[:imax]

# m_min, m_max = np.min(m_sig), np.max(m_sig)

# def myfitspec(k,g,msig,mpi):
#     c = np.abs(1 - 4*mpi**2/k**2)
#     I0 = 1/(2*np.sqrt(c))*np.log(np.abs(
#         (np.sqrt(c) - 1)/(np.sqrt(c) + 1)
#     ))
#     ReSigK =-3*g**2*k**2/(64*np.pi**2) * (
#         2/3 - c+(c-c**2)*I0
#     )
#     ImSigK = -3*g**2/(32*np.pi) * mpi**2 * (1-4*mpi**2/k**2)**(1/2)
#     return -2*np.pi*ImSigK / ( (k*2 - msig**2 -ReSigK)**2 + ImSigK**2)

# def myfit(xs,*allcoeffs):
#     allcoeffs = np.array(allcoeffs)
#     a, x0 = allcoeffs[0:2]
#     coeffs = allcoeffs[2:]
#     xs = np.array(xs)
#     exps = np.arange(len(coeffs))
#     return np.array([np.sum(coeffs*(a*(x-x0))**exps)*np.exp(-a**2*(x-x0)**2) for x in xs])

# popt, pcov = scipy.optimize.curve_fit(myfit, m_sig, w_sig,p0=(0.01,600,*np.ones(9)))
# def fitfunc(k):
#     return myfit(k,*popt)

mpi = 140
msigma = 479
Gam = 2*279

def Delta(s):
    return 1/(s-msigma**2+1j*Gam*np.sqrt(s-(2*mpi)**2))

def S(k):
    return -1/np.pi*np.imag(Delta(k**2))

# def spectralfunc(k):
#     return 1/(k**2-msigma**2+1j*np.sqrt(k**2-2*mpi**2))

m_min, m_max, m_mid = 280, 1500, 600
xs = np.linspace(m_min,m_max,10000)
ys = 2*xs*S(xs)


fig, ax = plt.subplots(figsize=(7,7))

# ax.plot(m_sig, w_sig)
# ax.plot(xs, fitfunc(xs))

### spectral bins
# bins_edges = np.concatenate((np.linspace(m_min, m_mid, 50, endpoint=False),np.linspace(m_mid, m_max, 20, endpoint=True)))
bins_edges = np.linspace(m_min,m_max,100,endpoint=True)
bins_lower = bins_edges[:-1]
bins_upper = bins_edges[1:]
bins_x, bins_y = data_to_bins(xs, ys,bins_lower, bins_upper)

# bins_y /= np.sum(bins_y)

ax.plot(xs,ys,marker="")
ax.plot(bins_x, bins_y)
print(bins_x.astype(int))
print(bins_y)
print("".join(["%d,"%(b) for b in bins_x.astype(int)]))
print("".join(["%.5f,"%(b) for b in bins_y]))

ax.set_xlabel(r"$k\ [MeV]$")
ax.set_ylabel(r"$S(k^2)\ [Mev^{-2}]$")

ax.tick_params(axis="both",labelsize=TICKLABELSIZE)

ax.xaxis.set_ticks_position("bottom")
ax.yaxis.set_ticks_position("left")

ax.grid(False,which="both")

fig.tight_layout()

# fig.savefig("images/sigmaresonance_spectral.png")

"""
[312 329 346 364 381 399 416 434 451 468 486 503 521 538 555 573 590 608
 625 643 660 677 695 712 730 747 764 782 799 817 834 851 869 886 904 921
 939 956 973 991]

[0.0113022  0.01342601 0.01521784 0.01710376 0.01888684 0.02073581
 0.02292558 0.02557391 0.02864946 0.03212066 0.03605239 0.04055638
 0.04562848 0.05098363 0.056002   0.05983979 0.0616738  0.06098072
 0.05773101 0.05241178 0.04586667 0.03902149 0.03261387 0.02703996
 0.02236943 0.01849156 0.01528533 0.0127     0.01070284 0.00916261
 0.00781289 0.00639653 0.00491597 0.00372084 0.00318787 0.00316915
 0.00296755 0.00237958 0.00227932 0.00211452]

[312,329,346,364,381,399,416,434,451,468,486,503,521,538,555,573,590,608,625,643,660,677,695,712,730,747,764,782,799,817,834,851,869,886,904,921,939,956,973,991]

0.01130,0.01343,0.01522,0.01710,0.01889,0.02074,0.02293,0.02557,0.02865,0.03212,0.03605,0.04056,0.04563,0.05098,0.05600,0.05984,0.06167,0.06098,0.05773,0.05241,0.04587,0.03902,0.03261,0.02704,0.02237,0.01849,0.01529,0.01270,0.01070,0.00916,0.00781,0.00640,0.00492,0.00372,0.00319,0.00317,0.00297,0.00238,0.00228,0.00211,
"""