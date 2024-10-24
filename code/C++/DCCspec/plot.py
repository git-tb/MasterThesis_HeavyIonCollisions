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

get_ipython().run_line_magic("matplotlib","qt")
plt.style.use("mplstyles/myclassic_white.mplstyle")

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
parentdir = "data/spectra_real_consteps_20240822_135426/*/"
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

# parentdir = "data/spectra_real_constfield_20240822_161734_masses/*/"
# paths = glob.glob(parentdir)
# ms = np.linspace(140,800,10,endpoint=True)

# parentdir = "data/spectra_real_constfield_20240822_161734/*/"
# paths = glob.glob(parentdir)
# ms = 140*np.ones(len(paths))

parentdir = "data/spectra_taudep/*/"
paths = glob.glob(parentdir)
ms = 140*np.ones(len(paths))

### epsconstmasses
# folders = [
#  'data/spectra_real_consteps_20240822_135426',
#  'data/spectra_real_consteps_20240826_143119_m226',
#  'data/spectra_real_consteps_20240826_143137_m312',
#  'data/spectra_real_consteps_20240826_143151_m398',
#  'data/spectra_real_consteps_20240826_143205_m484',
#  'data/spectra_real_consteps_20240826_143224_m570',
#  'data/spectra_real_consteps_20240826_143237_m656',
#  'data/spectra_real_consteps_20240826_143253_m742',
#  'data/spectra_real_consteps_20240826_143307_m828',
#  'data/spectra_real_consteps_20240826_143326_m914',
#  'data/spectra_real_consteps_20240826_143341_m1000']
# parentdir = "epsconst_masses"

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
#     paths.append(localpaths[-1])

# parentdir = "epsconst_masses"
### \epsconstmasses

SAVETITLE = parentdir.replace("data/","").replace("/*/","")

FIELD_XLABEL = r"$\alpha$"
DFIELD_XLABEL = r"$\alpha$"
SPEC_XLABEL = r"$p_T\ [GeV]$"

FIELD_YLABEL = r"$\pi^0\ [GeV]$"
DFIELD_YLABEL = r"$n^\mu\partial_\mu\pi^0\ [GeV^2]$"
# SPEC_YLABEL = r"$\mathcal{N}\cdot\frac{1}{2\pi p_T}\frac{dN_{coherent}}{dp_Td\eta_p}\ [GeV^{-2}]$"
SPEC_YLABEL = r"$(2\pi p_T)^{-1}dN_{coherent}/(dp_Td\eta_p)\ [GeV^{-2}]$"

SCALEBYMASS = False
if(SCALEBYMASS):
    FIELD_YLABEL = r"$m\cdot\pi^0\ [GeV^2]$"
NORMALIZE_SPECS = False
if(NORMALIZE_SPECS):
    SPEC_YLABEL = r"$\mathcal{N}\cdot(2\pi p_T)^{-1}dN_{coherent}/(dp_Td\eta_p)\ [GeV^{-2}]$"


FIGSIZE = (7,7)
AXISLABELSIZE = 20
TICKLABELSIZE = 15

LINECOLLECTION = False
WITHCMAP = False
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

LEGEND = False
LEGENDSIZE = 10
labels = [r"$m=%i\ MeV$"%(m) for m in ms]
if(not LEGEND):
    labels = ["" for n in range(len(paths))]

fig_init = plt.figure(figsize=FIGSIZE)
gs = gridspec.GridSpec(nrows=2, ncols=1, hspace=0)
ax_init1, ax_init2 = fig_init.add_subplot(gs[0]), fig_init.add_subplot(gs[1])

fig_spec, ax_spec = plt.subplots(figsize=FIGSIZE)

minval1, maxval1 = np.inf, -np.inf
minval2, maxval2 = np.inf, -np.inf

cols = [(t, 0, 1-t) for t in np.linspace(0,1,len(paths),endpoint=True)]
for (n,path) in enumerate(paths):
    df_spec = pd.read_csv(path+"/spectr.txt",comment="#")
    df_specanti = pd.read_csv(path+"/spectr_anti.txt",comment="#")
    df_field0 = pd.read_csv(path+"/field0.txt",comment="#")
    df_Dfield0 = pd.read_csv(path+"/field0_deriv.txt",comment="#")

    col = cols[n]

    scale_specs = 1
    if(NORMALIZE_SPECS):
        scale_specs = 1000/df_spec["abs2Re"].to_numpy()[0]
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
        color=col)

    ax_init2.plot(
        x_init2,
        y_init2,
        marker="",
        color=col)

    ax_spec.plot(
        x_spec,
        y_spec,
        marker="",
        color=col,
        label=labels[n])

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

ax_init1.tick_params(axis="both",labelsize=TICKLABELSIZE)
ax_init2.tick_params(axis="both",labelsize=TICKLABELSIZE)
ax_spec.tick_params(axis="both",labelsize=TICKLABELSIZE)

ax_init1.set_xlim(0,np.pi/2)
ax_init2.set_xlim(0,np.pi/2)

ax_init1.grid(False,which="both")
ax_init2.grid(False,which="both")
ax_spec.grid(False,which="both")

if(LEGEND):
    ax_spec.legend(fontsize=LEGENDSIZE,fancybox=True, framealpha=0.85,shadow=False)

# fig_init.suptitle(parentdir)
# fig_spec.suptitle(parentdir)

fig_init.tight_layout()
fig_spec.tight_layout()

fig_init.savefig("data/images/"+SAVETITLE+"_init.png",dpi=150)
fig_spec.savefig("data/images/"+SAVETITLE+"_spec.png",dpi=150)

# fig_init.show()
# fig_spec.show()

plt.show()

# %%
###############################################################
###############################################################
# PLOT DECAY SPECTRA FOR INITTEST DATA
###############################################################
###############################################################

paths = glob.glob("data/decay_inittest/decayspec_*/")
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

ax_primespec.set_ylabel(SPEC_YLABEL_PS)
ax_primespec.set_xlabel(SPEC_XLABEL_PS)

ax_decayspec.set_ylabel(SPEC_YLABEL_DS)
ax_decayspec.set_xlabel(SPEC_XLABEL_DS)

fig_decayspec.tight_layout()
fig_primespec.tight_layout()

fig_primespec.savefig("data/images/decaydecay_inittest_primespec.png")
fig_decayspec.savefig("data/images/decaydecay_inittest_decayspecs.png")

plt.show()