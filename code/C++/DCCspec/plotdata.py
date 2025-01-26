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

# %%
#############################################################
########### COMPARE FLUIDUM WITH ALICE DATASETS #############
#############################################################

SAVE = False
# SAVE = True

TITLESPEC = "Images/AliceFluidumCompare_specs"
TITLEDIFF = "Images/AliceFluidumCompare_diff"

FIGSIZE = (7,7)
LINEWIDTH = 2
MARKERSIZE = 50
LEGENDSIZE = 20
TICKLABELSIZE = 20
AXISLABELSIZE = 30
TICKSIZE = 10
TICKWIDTH = 2

df_alice = pd.read_csv("./../../Mathematica/data/HEPData-ins1222333-v1-Table_1.csv",comment="#")

pTs_alice = df_alice["PT [GEV]"].to_numpy()[:41].astype(float)
spec_alice = df_alice["(1/Nev)*(1/(2*PI*PT))*D2(N)/DPT/DYRAP [GEV**-2]"].to_numpy()[:41].astype(float)
spec_alice_err = df_alice["sys +"].to_numpy()[:41].astype(float)+df_alice["stat +"].to_numpy()[:41].astype(float)
bins_lower = df_alice["PT [GEV] LOW"].to_numpy()[:41].astype(float)
bins_upper = df_alice["PT [GEV] HIGH"].to_numpy()[:41].astype(float)

df = pd.read_csv("./../../Mathematica/data/pionListBig.txt")
pTs_fluidum= np.array(df.keys()).astype(float)
spec_fluidum_pi0 = df.iloc[0].to_numpy()
spec_fluidum_piplus = df.iloc[1].to_numpy()
spec_fluidum_piplus_compare = np.exp(scipy.interpolate.interp1d(pTs_fluidum, np.log(spec_fluidum_piplus))(pTs_alice))
pT_compare = pTs_alice
spec_compare = np.abs(spec_alice - spec_fluidum_piplus_compare)


fig_spec, ax_spec = plt.subplots(figsize=FIGSIZE)
fig_diff, ax_diff = plt.subplots(figsize=FIGSIZE)

ax_spec.plot(pTs_fluidum, spec_fluidum_piplus,label="FluiduM",marker="",lw=LINEWIDTH,c="r")
# ax_spec.plot(pTs_alice, spec_alice,label="ALICE",marker="",lw=LINEWIDTH,c="b")
# ax_spec.scatter(pTs_alice, spec_alice,label="ALICE",marker="x",c="b")
ax_spec.errorbar(pTs_alice, spec_alice, spec_alice_err,label="ALICE",c="b",fmt="o",markersize=MARKERSIZE/12,lw=LINEWIDTH)
# ax_diff.scatter(pT_compare, spec_compare,label=r"ALICE$-$FluiduM",marker="o",c="b",s=MARKERSIZE)
ax_diff.plot(pT_compare, spec_compare,label=r"ALICE$-$FluiduM",marker="",c="b",lw=LINEWIDTH)

ax_spec.set_yscale("log")
ax_diff.set_yscale("log")

ax_spec.grid(False, which="both")
ax_diff.grid(False, which="both")

ax_spec.set_ylabel(r"$(2\pi p_T)^{-1}dN/(dp_Td\eta_p)\ [{GeV}^{-2}]$",fontsize=AXISLABELSIZE)
ax_diff.set_ylabel(r"$(2\pi p_T)^{-1}dN/(dp_Td\eta_p)\ [{GeV}^{-2}]$",fontsize=AXISLABELSIZE)

ax_spec.set_xlabel(r"$p_T\ [GeV]$",fontsize=AXISLABELSIZE)
ax_diff.set_xlabel(r"$p_T\ [GeV]$",fontsize=AXISLABELSIZE)

ax_spec.legend(fontsize=LEGENDSIZE,fancybox=True, framealpha=0.85,shadow=False)
ax_diff.legend(fontsize=LEGENDSIZE,fancybox=True, framealpha=0.85,shadow=False)

ax_diff.tick_params(size=TICKSIZE,width=TICKWIDTH,labelsize=TICKLABELSIZE)
ax_spec.tick_params(size=TICKSIZE,width=TICKWIDTH,labelsize=TICKLABELSIZE)

ax_diff.tick_params(bottom=True,left=True,top=False,right=False)
ax_spec.tick_params(bottom=True,left=True,top=False,right=False)

fig_spec.tight_layout()
fig_diff.tight_layout()

if(SAVE):
    fig_spec.savefig(TITLESPEC+".svg")
    fig_diff.savefig(TITLEDIFF+".svg")
    print("SAVED!")
else:
    print("NOT SAVED!")

fig_spec.show()
fig_diff.show()

# %%
#############################################################
######### PLOT 1-PARAMETER FAMILY OF REAL FIELDS ############
#############################################################

parentdir = "Data/spec_real_constfield_m140_20250116_115949"
# parentdir = "Data/spec_real_constfield_varm"
# parentdir = "Data/spec_real_consteps_m140_20250117_113521"
# parentdir = "Data/spec_real_consteps_varm_20250116_114850"
# parentdir = "Data/spec_real_taudep_m140_20250117_162422"
folders = sorted(glob.glob(parentdir+"/*/"))

TITLEINIT = parentdir.replace("Data/","Images/") + "_init"
TITLESPEC = parentdir.replace("Data/","Images/") + "_spec"

SAVE = False
# SAVE = True

FIGSIZE = (7,7)
LINEWIDTH = 1.5
TICKLABELSIZE = 20
AXISLABELSIZE = 25
TICKSIZE = 10
TICKWIDTH = 2

DFIELD_XLABEL = r"$\alpha$"
SPEC_XLABEL = r"$p_T\ [GeV]$"

FIELD_YLABEL = r"$\phi_{DCC}\ [GeV]$"
DFIELD_YLABEL = r"$n^\mu\partial_\mu\phi_{DCC}\ [GeV]$"
SPEC_YLABEL = r"$(2\pi p_T)^{-1}dN_{DCC}/(dp_Td\eta_p)\ [GeV^{-2}]$"

CMAP = LinearSegmentedColormap.from_list("custom", [(0,0,1,0.6),(1,0,0,0.6)])
# CMAP = LinearSegmentedColormap.from_list("custom", [(0,0,1),(1,0,0)])
CMAP_LBWH = [0.025, 0.025, 0.05, 0.45]
CMAP_LABELSIZE = 20
CMAP_TICKSIZE = 15

LABELBYMASS = False
# LABELBYMASS = True

COMPARE = False
# COMPARE = True
LS_COMPARE = "-"



LC_ARRAY = np.arange(len(folders))
masses = np.zeros(len(folders))
if(LABELBYMASS):
    LC_ARRAY = masses
    LC_LABEL = r"$m\ [GeV]$"

fig_init = plt.figure(figsize=FIGSIZE)
gs = gridspec.GridSpec(nrows=2, ncols=1, hspace=0.0)
ax_init_f, ax_init_df = fig_init.add_subplot(gs[0]), fig_init.add_subplot(gs[1])
ax_init_f.sharex(ax_init_df)
fig_spec, ax_spec = plt.subplots(figsize=FIGSIZE)

lines_init_f= []
lines_init_df = []
lines_spec = []

cols = [(t, 0, 1-t) for t in np.linspace(0,1,len(folders),endpoint=True)]
for (n,folder) in enumerate(folders):
    with open(folder+"spectr.txt") as f:
        lines = f.readlines()
        j = [k for (k,x) in enumerate(["mass" in l for l in lines]) if x][0]
        masses[n]=float(lines[j].replace("# particle mass:\t","").replace("\n",""))

    df_spec = pd.read_csv(folder+"spectr.txt",comment="#") 
    df_init_f = pd.read_csv(folder+"field0.txt",comment="#")
    df_init_df = pd.read_csv(folder+"field0_deriv.txt",comment="#") 

    x_init_f = df_init_f["alpha"].to_numpy()
    x_init_df = df_init_df["alpha"].to_numpy()
    x_spec = df_spec["pT"].to_numpy() 

    y_init_f = df_init_f["field0Re"].to_numpy()
    y_init_df = df_init_df["Dfield0Re"].to_numpy()
    y_spec = df_spec["abs2Re"].to_numpy()

    col = cols[n]

    line_init_f = np.column_stack((x_init_f, y_init_f))
    line_init_df = np.column_stack((x_init_df, y_init_df))
    line_spec = np.column_stack((x_spec, y_spec))

    lines_spec.append(line_spec)
    lines_init_f.append(line_init_f)
    lines_init_df.append(line_init_df)  
    
lc_init_f = LineCollection(lines_init_f,array=LC_ARRAY,cmap=CMAP,linewidths=LINEWIDTH)
lc_init_df = LineCollection(lines_init_df,array=LC_ARRAY,cmap=CMAP,linewidths=LINEWIDTH)
lc_spec = LineCollection(lines_spec,array=LC_ARRAY,cmap=CMAP,linewidths=LINEWIDTH)

ax_init_f.add_collection(lc_init_f)
ax_init_df.add_collection(lc_init_df)
ax_spec.add_collection(lc_spec)

ax_init_f.autoscale_view()
ax_init_df.autoscale_view()
ax_spec.autoscale_view()

if(COMPARE):
    ax_spec.plot(pT_compare, spec_compare,lw=2*LINEWIDTH,label=r"ALICE$-$FluiduM",c="k",ls=LS_COMPARE,marker="")
    ax_spec.legend(fontsize=LEGENDSIZE,fancybox=True, framealpha=0.85,shadow=False)

if(LABELBYMASS):
    cax = ax_spec.inset_axes(CMAP_LBWH)
    cbar = fig_spec.colorbar(lc_spec, cax=cax)
    cbar.set_label(LC_LABEL, fontsize=CMAP_LABELSIZE)
    cbar.ax.tick_params(labelsize=CMAP_TICKSIZE)

ax_spec.set_yscale("log")

ax_init_f.set_ylabel(FIELD_YLABEL, fontsize=AXISLABELSIZE)
ax_init_df.set_ylabel(DFIELD_YLABEL, fontsize=AXISLABELSIZE)
ax_spec.set_ylabel(SPEC_YLABEL, fontsize=AXISLABELSIZE)

ax_init_df.set_xlabel(DFIELD_XLABEL, fontsize=AXISLABELSIZE)
ax_spec.set_xlabel(SPEC_XLABEL, fontsize=AXISLABELSIZE)

xticks = [0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2]
xticklabels = [r"$0$",r"$\pi/8$",r"$\pi/4$",r"$3\pi/8$",r"$\pi/2$"]
ax_init_f.set_xticks(xticks, xticklabels)
# ax_init_df.set_xticks(xticks, xticklabels)
ax_init_f.set_xticklabels(xticklabels,visible=False)

ax_spec.set_xlim(0,2)

yticks = ax_init_df.get_yticks()
yticklabels = ax_init_df.get_yticklabels()
yticklabels[-1] = ""
ax_init_df.set_yticks(yticks,yticklabels)

ax_init_f.tick_params(size=TICKSIZE,width=TICKWIDTH,labelsize=TICKLABELSIZE)
ax_init_df.tick_params(size=TICKSIZE,width=TICKWIDTH,labelsize=TICKLABELSIZE)
ax_spec.tick_params(size=TICKSIZE,width=TICKWIDTH,labelsize=TICKLABELSIZE)

ax_init_f.tick_params(bottom=True,left=True,top=True,right=True)
ax_init_df.tick_params(bottom=True,left=True,top=True,right=True)
ax_spec.tick_params(bottom=True,left=True,top=True,right=True)

ax_init_f.set_xlim(0,np.pi/2)
# ax_init_df.set_xlim(0,np.pi/2)

ax_init_f.grid(False,which="both")
ax_init_df.grid(False,which="both")
ax_spec.grid(False,which="both")

fig_init.tight_layout()
fig_spec.tight_layout()

if(SAVE):
    fig_init.savefig(TITLEINIT+".svg")
    fig_spec.savefig(TITLESPEC+".svg")
    print("SAVED!")
else:
    print("NOT SAVED!")

fig_init.show()
fig_spec.show()

# %%
#############################################################
############# COMPARE 1-PARAMETER FAMILIES #################
#############################################################

folders = [
    sorted(glob.glob("Data/spec_real_constfield_m140_20250116_115949/*/"))[0],
    sorted(glob.glob("Data/spec_real_consteps_m140_20250117_113521/*/"))[-1],
    sorted(glob.glob("Data/spec_real_taudep_m140_20250117_162422/*/"))[-1]
]

# folders = [    
#     sorted(glob.glob("Data/spec_real_constfield_m280_20250122_110401/*/"))[0],
#     sorted(glob.glob("Data/spec_real_consteps_m280_20250122_110422/*/"))[-1],
#     sorted(glob.glob("Data/spec_real_taudep_m280_20250122_110439/*/"))[-1]
# ]

# folders = [
#     sorted(glob.glob("Data/spec_real_constfield_m420_20250122_111331/*/"))[0],
#     sorted(glob.glob("Data/spec_real_consteps_m420_20250122_111354/*/"))[-1],
#     sorted(glob.glob("Data/spec_real_taudep_m420_20250122_111415/*/"))[-1]
# ]

TITLEINIT = "Images/spec_real_compare_m140_init"
TITLESPEC = "Images/spec_real_compare_m140_spec"

LABELS = [
    "Type A", 
    "Type B",
    "Type C"
]

# SAVE = False
SAVE = True

FIGSIZE = (7,7)
LINEWIDTH = 2
TICKLABELSIZE = 20
AXISLABELSIZE = 25
TICKSIZE = 10
TICKWIDTH = 2

DFIELD_XLABEL = r"$\alpha$"
SPEC_XLABEL = r"$p_T\ [GeV]$"

FIELD_YLABEL = r"$\phi_{DCC}\ [GeV]$"
DFIELD_YLABEL = r"$n^\mu\partial_\mu\phi_{DCC}\ [GeV]$"
SPEC_YLABEL = r"$(2\pi p_T)^{-1}dN_{DCC}/(dp_Td\eta_p)\ [GeV^{-2}]$"

# COMPARE = False
COMPARE = True
LS_COMPARE = "-"

COLS = ['r','b','g']

fig_init = plt.figure(figsize=FIGSIZE)
gs = gridspec.GridSpec(nrows=2, ncols=1, hspace=0.0)
ax_init_f, ax_init_df = fig_init.add_subplot(gs[0]), fig_init.add_subplot(gs[1])
ax_init_f.sharex(ax_init_df)
fig_spec, ax_spec = plt.subplots(figsize=FIGSIZE)

for (n,folder) in enumerate(folders):

    df_spec = pd.read_csv(folder+"spectr.txt",comment="#") 
    df_init_f = pd.read_csv(folder+"field0.txt",comment="#")
    df_init_df = pd.read_csv(folder+"field0_deriv.txt",comment="#") 

    x_init_f = df_init_f["alpha"].to_numpy()
    x_init_df = df_init_df["alpha"].to_numpy()
    x_spec = df_spec["pT"].to_numpy() 

    y_init_f = df_init_f["field0Re"].to_numpy()
    y_init_df = df_init_df["Dfield0Re"].to_numpy()
    y_spec = df_spec["abs2Re"].to_numpy()

    ax_init_f.plot(x_init_f, y_init_f,lw=LINEWIDTH,c=COLS[n],marker="")
    ax_init_df.plot(x_init_df, y_init_df,lw=LINEWIDTH,c=COLS[n],marker="")
    ax_spec.plot(x_spec, y_spec,lw=LINEWIDTH,c=COLS[n],marker="",label=LABELS[n])

if(COMPARE):
    ax_spec.plot(pT_compare, spec_compare,lw=LINEWIDTH,label=r"ALICE$-$FluiduM",c="k",ls=LS_COMPARE,marker="")
    ax_spec.legend(fontsize=LEGENDSIZE,fancybox=True, framealpha=0.85,shadow=False)

ax_spec.set_yscale("log")

ax_init_f.set_ylabel(FIELD_YLABEL, fontsize=AXISLABELSIZE)
ax_init_df.set_ylabel(DFIELD_YLABEL, fontsize=AXISLABELSIZE)
ax_spec.set_ylabel(SPEC_YLABEL, fontsize=AXISLABELSIZE)

ax_init_df.set_xlabel(DFIELD_XLABEL, fontsize=AXISLABELSIZE)
ax_spec.set_xlabel(SPEC_XLABEL, fontsize=AXISLABELSIZE)

xticks = [0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2]
xticklabels = [r"$0$",r"$\pi/8$",r"$\pi/4$",r"$3\pi/8$",r"$\pi/2$"]
ax_init_f.set_xticks(xticks, xticklabels)
# ax_init_df.set_xticks(xticks, xticklabels)
ax_init_f.set_xticklabels(xticklabels,visible=False)

ax_spec.set_xlim(0,2)

yticks = ax_init_df.get_yticks()
yticklabels = ax_init_df.get_yticklabels()
yticklabels[-1] = ""
ax_init_df.set_yticks(yticks,yticklabels)

ax_init_f.tick_params(size=TICKSIZE,width=TICKWIDTH,labelsize=TICKLABELSIZE)
ax_init_df.tick_params(size=TICKSIZE,width=TICKWIDTH,labelsize=TICKLABELSIZE)
ax_spec.tick_params(size=TICKSIZE,width=TICKWIDTH,labelsize=TICKLABELSIZE)

ax_init_f.tick_params(bottom=True,left=True,top=True,right=True)
ax_init_df.tick_params(bottom=True,left=True,top=True,right=True)
ax_spec.tick_params(bottom=True,left=True,top=True,right=True)

ax_init_f.set_xlim(0,np.pi/2)
# ax_init_df.set_xlim(0,np.pi/2)

ax_init_f.grid(False,which="both")
ax_init_df.grid(False,which="both")
ax_spec.grid(False,which="both")

fig_init.tight_layout()
fig_spec.tight_layout()

if(SAVE):
    fig_init.savefig(TITLEINIT+".svg")
    fig_spec.savefig(TITLESPEC+".svg")
    print("SAVED!")
else:
    print("NOT SAVED!")

fig_init.show()
fig_spec.show()

# %%
#############################################################
######### PLOT COMP FIELDS WITH VARYING RELPHASE ############
#############################################################

# parentdir = "Data/spec_comp_constfield_m140_varphase_20250116_120129"
# parentdir = "Data/spec_comp_constfield_m140_varphase_20250117_181049"
parentdir = "Data/spec_comp_constfield_m140_varphase_20250117_181304"
# parentdir = "Data/spec_comp_constfield_m140_varphase_20250120_094422"
folders = sorted(glob.glob(parentdir+"/*/"))

TITLEINIT = parentdir.replace("Data/","Images/") + "_init"
TITLESPEC = parentdir.replace("Data/","Images/") + "_spec"

# SAVE = False
SAVE = True

FIGSIZE = (14,7)
LINEWIDTH = 1.5
TICKLABELSIZE = 20
AXISLABELSIZE = 25
TICKSIZE = 10
TICKWIDTH = 2

FIELD_XLABEL = r"$\alpha$"
DFIELD_XLABEL = r"$\alpha$"
SPEC_XLABEL = r"$p_T\ [GeV]$"

FIELDRE_YLABEL = r"$\Re\phi_{DCC}\ [GeV]$"
FIELDIM_YLABEL = r"$\Im\phi_{DCC}\ [GeV]$"
DFIELDRE_YLABEL = r"$\Re n^\mu\partial_\mu\phi_{DCC}\ [GeV]$"
DFIELDIM_YLABEL = r"$\Im n^\mu\partial_\mu\phi_{DCC}\ [GeV]$"
SPEC_YLABEL = r"$(2\pi p_T)^{-1}dN_{DCC}/(dp_Td\eta_p)\ [GeV^{-2}]$"
SPECANTI_YLABEL = r"$(2\pi p_T)^{-1}d\overline{N}_{DCC}/(dp_Td\eta_p)\ [GeV^{-2}]$"

CMAP = LinearSegmentedColormap.from_list("custom", [(0,0,1,0.6),(1,0,0,0.6)])
# CMAP = LinearSegmentedColormap.from_list("custom", [(0,0,1),(1,0,0)])
CMAP_LBWH = [0.025, 0.025, 0.05, 0.45]
CMAP_LABELSIZE = 20
CMAP_TICKSIZE = 15

phases = np.zeros(len(folders))
LC_ARRAY = phases
LC_LABEL = r"$\Delta\varphi$"

fig_init = plt.figure(figsize=FIGSIZE)
gs_init = gridspec.GridSpec(nrows=2, ncols=2, hspace=0,wspace=0)
ax_init_f_re, ax_init_df_re, ax_init_f_im, ax_init_df_im = fig_init.add_subplot(gs_init[0,0]), fig_init.add_subplot(gs_init[1,0]),fig_init.add_subplot(gs_init[0,1]), fig_init.add_subplot(gs_init[1,1])
ax_init_f_re.sharex(ax_init_df_re)
ax_init_f_im.sharex(ax_init_df_im)
ax_init_f_re.sharey(ax_init_f_im)
ax_init_df_re.sharey(ax_init_df_im)

fig_spec = plt.figure(figsize=FIGSIZE)
gs_spec = gridspec.GridSpec(ncols=2,nrows=1, hspace=0,wspace=0)
ax_spec, ax_spec_anti = fig_spec.add_subplot(gs_spec[0]),fig_spec.add_subplot(gs_spec[1])

lines_init_f_re= []
lines_init_df_re = []
lines_init_f_im= []
lines_init_df_im = []
lines_spec = []
lines_spec_anti = []

cols = [(t, 0, 1-t) for t in np.linspace(0,1,len(folders),endpoint=True)]
for (n,folder) in enumerate(folders):
    df_spec = pd.read_csv(folder+"spectr.txt",comment="#") 
    df_spec_anti = pd.read_csv(folder+"spectr_anti.txt",comment="#") 
    df_init_f = pd.read_csv(folder+"field0.txt",comment="#")
    df_init_df = pd.read_csv(folder+"field0_deriv.txt",comment="#") 

    x_init_f = df_init_f["alpha"].to_numpy()
    x_init_df = df_init_df["alpha"].to_numpy()
    x_spec = df_spec["pT"].to_numpy() 

    y_init_f_re = df_init_f["field0Re"].to_numpy()
    y_init_f_im = df_init_f["field0Im"].to_numpy()
    y_init_df_re = df_init_df["Dfield0Re"].to_numpy()
    y_init_df_im = df_init_df["Dfield0Im"].to_numpy()
    y_spec = df_spec["abs2Re"].to_numpy()
    y_spec_anti = df_spec_anti["abs2Re"].to_numpy()

    phases[n] = np.angle(y_init_df_re + 1j * y_init_df_im)[0]-np.angle(y_init_f_re + 1j * y_init_f_im)[0]

    col = cols[n]

    line_init_f_re = np.column_stack((x_init_f, y_init_f_re))
    line_init_f_im = np.column_stack((x_init_f, y_init_f_im))
    line_init_df_re = np.column_stack((x_init_df, y_init_df_re))
    line_init_df_im = np.column_stack((x_init_df, y_init_df_im))
    line_spec = np.column_stack((x_spec, y_spec))
    line_spec_anti = np.column_stack((x_spec, y_spec_anti))

    lines_spec.append(line_spec)
    lines_spec_anti.append(line_spec_anti)
    lines_init_f_re.append(line_init_f_re)
    lines_init_f_im.append(line_init_f_im)
    lines_init_df_re.append(line_init_df_re)
    lines_init_df_im.append(line_init_df_im)  
    
lc_init_f_re = LineCollection(lines_init_f_re,array=LC_ARRAY,cmap=CMAP,linewidths=LINEWIDTH)
lc_init_f_im = LineCollection(lines_init_f_im,array=LC_ARRAY,cmap=CMAP,linewidths=LINEWIDTH)
lc_init_df_re = LineCollection(lines_init_df_re,array=LC_ARRAY,cmap=CMAP,linewidths=LINEWIDTH)
lc_init_df_im = LineCollection(lines_init_df_im,array=LC_ARRAY,cmap=CMAP,linewidths=LINEWIDTH)
lc_spec = LineCollection(lines_spec,array=LC_ARRAY,cmap=CMAP,linewidths=LINEWIDTH)
lc_spec_anti = LineCollection(lines_spec_anti,array=LC_ARRAY,cmap=CMAP,linewidths=LINEWIDTH)

ax_init_f_re.add_collection(lc_init_f_re)
ax_init_f_im.add_collection(lc_init_f_im)
ax_init_df_re.add_collection(lc_init_df_re)
ax_init_df_im.add_collection(lc_init_df_im)
ax_spec.add_collection(lc_spec)
ax_spec_anti.add_collection(lc_spec_anti)

ax_init_f_re.autoscale_view()
ax_init_f_im.autoscale_view()
ax_init_df_re.autoscale_view()
ax_init_df_im.autoscale_view()
ax_spec.autoscale_view()
ax_spec_anti.autoscale_view()

cax = ax_spec_anti.inset_axes(CMAP_LBWH)
cbar = fig_spec.colorbar(lc_spec, cax=cax)
cbar.set_label(LC_LABEL, fontsize=CMAP_LABELSIZE)
cbar.ax.tick_params(labelsize=CMAP_TICKSIZE)
cbar.ax.set_yticks([0,np.pi/2,np.pi],[r"$0$",r"$\pi/2$",r"$\pi$"])

ax_spec.set_yscale("log")
ax_spec_anti.set_yscale("log")

ax_init_f_re.set_ylabel(FIELDRE_YLABEL, fontsize=AXISLABELSIZE)
ax_init_f_im.set_ylabel(FIELDIM_YLABEL, fontsize=AXISLABELSIZE)
ax_init_df_re.set_ylabel(DFIELDRE_YLABEL, fontsize=AXISLABELSIZE)
ax_init_df_im.set_ylabel(DFIELDIM_YLABEL, fontsize=AXISLABELSIZE)
ax_spec.set_ylabel(SPEC_YLABEL, fontsize=AXISLABELSIZE)
ax_spec_anti.set_ylabel(SPECANTI_YLABEL, fontsize=AXISLABELSIZE)

ax_init_f_im.yaxis.set_ticks_position("right")
ax_init_f_im.yaxis.set_label_position("right")
ax_init_df_im.yaxis.set_ticks_position("right")
ax_init_df_im.yaxis.set_label_position("right")

ax_spec_anti.yaxis.set_ticks_position("right")
ax_spec_anti.yaxis.set_label_position("right")

ax_init_df_re.set_xlabel(DFIELD_XLABEL, fontsize=AXISLABELSIZE)
ax_init_df_im.set_xlabel(DFIELD_XLABEL, fontsize=AXISLABELSIZE)
ax_spec.set_xlabel(SPEC_XLABEL, fontsize=AXISLABELSIZE)
ax_spec_anti.set_xlabel(SPEC_XLABEL, fontsize=AXISLABELSIZE)

xticks = [0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2]
xticklabels = [r"$0$",r"$\pi/8$",r"$\pi/4$",r"$3\pi/8$",r"$\pi/2$"]
ax_init_f_im.set_xticks(xticks, xticklabels,visible=False)
xticklabels[-1]=""
ax_init_f_re.set_xticks(xticks, xticklabels,visible=False)

yticks = ax_init_df_re.get_yticks()
yticklabels = ax_init_df_re.get_yticklabels()
yticklabels[-1] = ""
ax_init_df_re.set_yticks(yticks,yticklabels)

xticks = ax_spec.get_xticks()
xticklabels = ax_spec.get_xticklabels()
xticklabels[-1] = ""
ax_spec.set_xticks(xticks,xticklabels)

ax_init_f_re.tick_params(size=TICKSIZE,width=TICKWIDTH,labelsize=TICKLABELSIZE)
ax_init_f_im.tick_params(size=TICKSIZE,width=TICKWIDTH,labelsize=TICKLABELSIZE)
ax_init_df_re.tick_params(size=TICKSIZE,width=TICKWIDTH,labelsize=TICKLABELSIZE)
ax_init_df_im.tick_params(size=TICKSIZE,width=TICKWIDTH,labelsize=TICKLABELSIZE)
ax_spec.tick_params(size=TICKSIZE,width=TICKWIDTH,labelsize=TICKLABELSIZE)
ax_spec_anti.tick_params(size=TICKSIZE,width=TICKWIDTH,labelsize=TICKLABELSIZE)

ax_init_f_re.tick_params(bottom=False,left=True,top=True,right=False)
ax_init_f_im.tick_params(bottom=False,left=False,top=True,right=True)
ax_init_df_re.tick_params(bottom=True,left=True,top=True,right=False)
ax_init_df_im.tick_params(bottom=True,left=False,top=True,right=True)
ax_spec.tick_params(bottom=True,left=True,top=True,right=False)
ax_spec_anti.tick_params(bottom=True,left=False,top=True,right=True)

ax_init_f_re.set_xlim(0,np.pi/2)
ax_init_f_im.set_xlim(0,np.pi/2)
ax_spec.set_xlim(0,2)
ax_spec_anti.set_xlim(0,2)

ax_init_f_re.grid(False,which="both")
ax_init_f_im.grid(False,which="both")
ax_init_df_re.grid(False,which="both")
ax_init_df_im.grid(False,which="both")
ax_spec.grid(False,which="both")
ax_spec_anti.grid(False,which="both")

fig_init.tight_layout()
fig_spec.tight_layout()

if(SAVE):
    fig_init.savefig(TITLEINIT+".svg")
    fig_spec.savefig(TITLESPEC+".svg")
    print("SAVED!")
else:
    print("NOT SAVED!")

fig_init.show()
fig_spec.show()

# %%
#############################################################
######### PLOT PRIMARY VS DECAY SPECTRA ############
#############################################################

# parentdir = "Data/decay_real_constfield_m600_varPT"
parentdir = "Data/decay_real_constfield_varm"
folders = sorted(glob.glob(parentdir+"/*/"))

TITLE_PRIMESPEC = parentdir.replace("Data/","Images/")+"_primespec"
TITLE_DECAYSPEC = parentdir.replace("Data/","Images/")+"_decayspec"

# SAVE = False
SAVE = True

FIGSIZE = (7,7)
LINEWIDTH = 1.5
TICKLABELSIZE = 20
AXISLABELSIZE = 25 
TICKSIZE = 10
TICKWIDTH = 2

PRIMESPEC_XLABEL = r"$q_T\ [GeV]$"
DECAYSPEC_XLABEL = r"$p_T\ [GeV]$"

PRIMESPEC_YLABEL = r"$(2\pi q_T)^{-1}dN_a/(dq_Td\eta_q)\ [GeV^{-2}]$"
DECAYSPEC_YLABEL = r"$(2\pi p_T)^{-1}dN_b/(dp_Td\eta_p)\ [GeV^{-2}]$"

CMAP = LinearSegmentedColormap.from_list("custom", [(0,0,1,0.6),(1,0,0,0.6)])
# CMAP = LinearSegmentedColormap.from_list("custom", [(0,0,1),(1,0,0)])
CMAP_LBWH = [0.025, 0.025, 0.05, 0.45]
CMAP_LABELSIZE = 20
CMAP_TICKSIZE = 15

COMPARE = False
# COMPARE = True
LS_COMPARE = "-"

# LABELLC = False
LABELLC = True


LC_ARRAY = np.arange(len(folders))
masses = np.zeros(len(folders))
qTmaxs = np.zeros(len(folders))
if(LABELLC):
    LC_ARRAY = masses
    LC_LABEL = r"$m_a\ [GeV]$"
    # LC_ARRAY = qTmaxs
    # LC_LABEL = r"$q_{T,\text{max}}\ [GeV]$"

fig_primespec, ax_primespec = plt.subplots(figsize=FIGSIZE)
fig_decayspec, ax_decayspec = plt.subplots(figsize=FIGSIZE)

lines_primespec = []
lines_decayspec = []

cols = [(t, 0, 1-t) for t in np.linspace(0,1,len(folders),endpoint=True)]
for (n,folder) in enumerate(folders):
    with open(folder+"decayspec.txt") as file:
        mylines = file.readlines()
        masses[n] = float(mylines[3].replace("# ma:\t",""))
        qTmaxs[n] = float(mylines[2].replace("# qmax:\t",""))

    df_primespec = pd.read_csv(folder+"primespec.txt",comment="#") 
    df_decayspec = pd.read_csv(folder+"decayspec.txt",comment="#") 

    qT = df_primespec["q"].to_numpy()
    pT = df_decayspec["p"].to_numpy()

    primespec = df_primespec["primespecRe"].to_numpy()
    decayspec = df_decayspec["finalspecRe"].to_numpy()

    col = cols[n]

    line_primespec = np.column_stack((qT, primespec))
    line_decayspec = np.column_stack((pT, decayspec))

    lines_primespec.append(line_primespec)
    lines_decayspec.append(line_decayspec)
    
lc_primespec = LineCollection(lines_primespec,array=LC_ARRAY,cmap=CMAP,linewidths=LINEWIDTH)
lc_decayspec = LineCollection(lines_decayspec,array=LC_ARRAY,cmap=CMAP,linewidths=LINEWIDTH)

ax_primespec.add_collection(lc_primespec)
ax_decayspec.add_collection(lc_decayspec)

ax_primespec.autoscale_view()
ax_decayspec.autoscale_view()

if(COMPARE):
    ax_decayspec.plot(pT_compare, spec_compare,lw=2*LINEWIDTH,label=r"ALICE$-$FluiduM",c="k",ls=LS_COMPARE,marker="")

if(LABELLC):
    cax = ax_decayspec.inset_axes(CMAP_LBWH)
    cbar = fig_decayspec.colorbar(lc_decayspec, cax=cax)
    cbar.set_label(LC_LABEL, fontsize=CMAP_LABELSIZE)
    cbar.ax.tick_params(labelsize=CMAP_TICKSIZE)

ax_primespec.set_yscale("log")
ax_decayspec.set_yscale("log")

ax_primespec.set_xlabel(PRIMESPEC_XLABEL, fontsize=AXISLABELSIZE)
ax_decayspec.set_xlabel(DECAYSPEC_XLABEL, fontsize=AXISLABELSIZE)

ax_primespec.set_ylabel(PRIMESPEC_YLABEL, fontsize=AXISLABELSIZE)
ax_decayspec.set_ylabel(DECAYSPEC_YLABEL, fontsize=AXISLABELSIZE)

ax_primespec.tick_params(size=TICKSIZE,width=TICKWIDTH,labelsize=TICKLABELSIZE)
ax_decayspec.tick_params(size=TICKSIZE,width=TICKWIDTH,labelsize=TICKLABELSIZE)

ax_primespec.tick_params(bottom=True,left=True,top=True,right=True)
ax_decayspec.tick_params(bottom=True,left=True,top=True,right=True)

ax_primespec.grid(False,which="both")
ax_decayspec.grid(False,which="both")

fig_primespec.tight_layout()
fig_decayspec.tight_layout()

if(SAVE):
    fig_primespec.savefig(TITLE_PRIMESPEC+".svg")
    fig_decayspec.savefig(TITLE_DECAYSPEC+".svg")
    print("SAVED!")
else:
    print("NOT SAVED!")

fig_primespec.show()
fig_decayspec.show()

# %%
#############################################################
######### PLOT WEIGHTED SIGMA RESONANCE DECAY SPEC ############
#############################################################

parentdir = "Data/decay_real_constfield_varm"
folders = sorted(glob.glob(parentdir+"/*/"))

# SAVE = False
SAVE = True

TITLE_DECAYSPEC = parentdir.replace("Data/","Images/")+"_weighted_decayspec"
TITLE_FULLSPEC = parentdir.replace("Data/","Images/")+"_weighted_fullspec"
TITLE_WEIGHTS = parentdir.replace("Data/","Images/")+"_weighted_weights"
TITLE_TOTALSPEC = parentdir.replace("Data/","Images/")+"_weighted_totalspec"

FIGSIZE = (7,7)
LINEWIDTH = 1.5
TICKLABELSIZE = 20
AXISLABELSIZE = 25 
TICKSIZE = 10
TICKWIDTH = 2

DECAYSPEC_XLABEL = r"$p_T\ [GeV]$"
DECAYSPEC_YLABEL = r"$(2\pi p_T)^{-1}dN_{\pi^\pm}/(dp_Td\eta_p)\ [GeV^{-2}]$"

TOTALSPEC_XLABEL = r"$p_T\ [GeV]$"
TOTALSPEC_YLABEL = r"$(2\pi p_T)^{-1}dN/(dp_Td\eta_p)\ [GeV^{-2}]$"

FULLSPEC_YLABEL = r"$(2\pi p_T)^{-1}dN/(dp_Td\eta_p)\ [GeV^{-2}]$"
FULLSPEC_LEGEND = r"decay spectrum"+"\n"+r"from $\sigma\to\pi\pi$"

WEIGHTS_XLABEL = r"$k\ [GeV]$"
WEIGHTS_YLABEL = r"$S(k)\ [GeV^{-2}]$"
WEIGHTS_LEGEND = r"spectral density"+"\n"+r"of the $f_0(500)$"

CMAP = LinearSegmentedColormap.from_list("custom", [(0,0,1,0.6),(1,0,0,0.6)])
# CMAP = LinearSegmentedColormap.from_list("custom", [(0,0,1),(1,0,0)])
CMAP_LBWH = [0.025, 0.025, 0.05, 0.45]
CMAP_LABELSIZE = 20
CMAP_TICKSIZE = 15

# COMPARE = False
COMPARE = True
LS_COMPARE = "-"

LEGEND= r"$\sigma_{DCC}\to\pi\,\pi$"

# LABELLC = False
LABELLC = True

COL = (1,0,0)

###
def data_to_bins(datax, datay, bins_lower, bins_upper):
    bins_x = ( bins_upper + bins_lower ) / 2

    bins_all = [[] for i in range(len(bins_x))]

    for (x,y) in list(zip(datax, datay)):
        try:
            binidx = np.where((x >= bins_lower).astype(int) * (x <= bins_upper).astype(int))[0][0]
        except:
            continue
        bins_all[binidx].append(y)

    bins_y = []
    for bin_y in bins_all:
        bins_y.append(np.mean(bin_y))

    return np.array(bins_x), np.array(bins_y)
###

###
df_alice = pd.read_csv("./../../Mathematica/data/HEPData-ins1222333-v1-Table_1.csv",comment="#")
pTs_alice = df_alice["PT [GEV]"].to_numpy()[:41].astype(float)
spec_alice = df_alice["(1/Nev)*(1/(2*PI*PT))*D2(N)/DPT/DYRAP [GEV**-2]"].to_numpy()[:41].astype(float)
spec_alice_err = df_alice["sys +"].to_numpy()[:41].astype(float)+df_alice["stat +"].to_numpy()[:41].astype(float)
bins_lower = df_alice["PT [GEV] LOW"].to_numpy()[:41].astype(float)
bins_upper = df_alice["PT [GEV] HIGH"].to_numpy()[:41].astype(float)

df = pd.read_csv("./../../Mathematica/data/pionListBig.txt")
pTs_fluidum= np.array(df.keys()).astype(float)
spec_fluidum_pi0 = df.iloc[0].to_numpy()
spec_fluidum_piplus = df.iloc[1].to_numpy()

def interplog(datax, datay, newx):
    return np.exp(scipy.interpolate.interp1d(datax, np.log(datay))(newx))

def mylogpolyfit(datax, datay):
    popt = np.polyfit(datax,np.log(datay),10)
    return np.poly1d(popt)
###

LC_ARRAY = np.arange(len(folders))
masses = np.zeros(len(folders))
if(LABELLC):
    LC_ARRAY = masses
    LC_LABEL = r"$m_\sigma\ [GeV]$"

### EXTRACT MASSES FOR EACH SPECTRUM
for (n,folder) in enumerate(folders):
    with open(folder+"decayspec.txt") as file:
        mylines = file.readlines()
        masses[n] = float(mylines[3].replace("# ma:\t",""))

### COMPUTE SPECTRAL WEIGHTS OF SIGMA RESONANCE
mpi = 0.14

Mpole = 0.449
Gpole = 2* 0.275

msigma = np.sqrt(1/4 * (16 * mpi**2 + 
                        np.sqrt(16 * Gpole**2 * Mpole**2 + 
                                (-16 * mpi**2 - Gpole**2 + 4*Mpole**2)**2)))
Gam = np.sqrt(1/2 * (16 * mpi**2 + Gpole**2 - 
     4*  Mpole**2 + np.sqrt(16*Gpole**2 * Mpole**2 + (-16 * mpi**2 - Gpole**2 + 4*Mpole**2)**2)))

def Delta(s):
    return 1/(s-msigma**2+1j*Gam*np.sqrt(s-(2*mpi)**2))

def S(k):
    # return -1/np.pi*np.imag(Delta(k**2))*(k<=0.7)
    return -1/np.pi*np.imag(Delta(k**2))

weights = 2*masses*S(masses)*np.ptp(masses)/len(masses)
weights /= np.sum(weights)

fig_weights, ax_weights = plt.subplots(figsize=FIGSIZE)
ax_weights.plot([0.28,*masses],S(np.array([0.28,*masses])),c="b",lw=2*LINEWIDTH,marker="",label=WEIGHTS_LEGEND)
ax_weights.legend(fontsize=LEGENDSIZE,fancybox=True, framealpha=0.85,shadow=False)

ax_weights.set_xlim(0,np.max(masses))
ax_weights.set_xlabel(WEIGHTS_XLABEL, fontsize=AXISLABELSIZE)
ax_weights.set_ylabel(WEIGHTS_YLABEL, fontsize=AXISLABELSIZE)
ax_weights.tick_params(size=TICKSIZE,width=TICKWIDTH,labelsize=TICKLABELSIZE)
ax_weights.tick_params(bottom=True,left=True,top=True,right=True)
ax_weights.grid(False,which="both")
fig_weights.tight_layout()

### GO TO PLOTTING
fig_fullspec, ax_fullspec = plt.subplots(figsize=FIGSIZE)
fig_totalspec, ax_totalspec = plt.subplots(figsize=FIGSIZE)
fig_decayspec, ax_decayspec = plt.subplots(figsize=FIGSIZE)

lines_decayspec = []
xs_fullspec, ys_fullspec = np.zeros(shape=(2,1))

for (n,folder) in enumerate(folders):
    df_decayspec = pd.read_csv(folder+"decayspec.txt",comment="#")

    pT = df_decayspec["p"].to_numpy()
    decayspec = df_decayspec["finalspecRe"].to_numpy()

    line_decayspec = np.column_stack((pT, decayspec))
    lines_decayspec.append(line_decayspec)

    xs_fullspec = pT
    ys_fullspec = ys_fullspec + weights[n] * decayspec

lc_decayspec = LineCollection(lines_decayspec,array=LC_ARRAY,cmap=CMAP,linewidths=LINEWIDTH)
ax_decayspec.add_collection(lc_decayspec)
ax_decayspec.autoscale_view()

ax_fullspec.plot(xs_fullspec, ys_fullspec,lw=2*LINEWIDTH,c="b",marker="",label=FULLSPEC_LEGEND)
ax_fullspec.set_xlim(np.min(xs_fullspec),np.max(xs_fullspec))

ax_fullspec.set_yscale("log")
ax_decayspec.set_yscale("log")
ax_totalspec.set_yscale("log")

###
fluidum_loginterp = mylogpolyfit(pTs_fluidum, spec_fluidum_piplus)
y_fluidum_interp = np.exp(fluidum_loginterp(xs_fullspec))
ax_totalspec.fill_between(xs_fullspec,y_fluidum_interp,y_fluidum_interp+2*ys_fullspec,facecolor=(*COL,0.2),edgecolor=COL,lw=LINEWIDTH,label=LEGEND)
ax_totalspec.plot(xs_fullspec, y_fluidum_interp,marker="",label="Fluidum",lw=LINEWIDTH)
ax_totalspec.errorbar(pTs_alice, spec_alice, spec_alice_err,label="ALICE",c="b",fmt="o",markersize=MARKERSIZE/12,lw=LINEWIDTH)
ax_totalspec.set_xlim(0,0.7)
ax_totalspec.set_ylim(8e1,3e3)

ax_totalspec.legend(fontsize=LEGENDSIZE,fancybox=True, framealpha=0.85,shadow=False)

ax_totalspec.tick_params(size=0.7*TICKSIZE,width=0.7*TICKWIDTH,labelsize=0.7*TICKLABELSIZE,which="minor")
locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8))
ax_totalspec.yaxis.set_minor_locator(locmin)
ax_totalspec.yaxis.set_minor_formatter(matplotlib.ticker.LogFormatterSciNotation(base=10,labelOnlyBase=False,minor_thresholds=(5,2.5))) # all of the data spans 5 decades and we want to see all minor ticks if we zoom in on a region of 2.5 decades
###

if(COMPARE):
    ax_fullspec.plot(pT_compare, spec_compare,lw=2*LINEWIDTH,label=r"ALICE$-$FluiduM",c="k",ls=LS_COMPARE,marker="")
ax_fullspec.legend(fontsize=LEGENDSIZE,fancybox=True, framealpha=0.85,shadow=False)

if(LABELLC):
    cax = ax_decayspec.inset_axes(CMAP_LBWH)
    cbar = fig_decayspec.colorbar(lc_decayspec, cax=cax)
    cbar.set_label(LC_LABEL, fontsize=CMAP_LABELSIZE)
    cbar.ax.tick_params(labelsize=CMAP_TICKSIZE)


ax_fullspec.set_xlabel(DECAYSPEC_XLABEL, fontsize=AXISLABELSIZE)
ax_decayspec.set_xlabel(DECAYSPEC_XLABEL, fontsize=AXISLABELSIZE)
ax_totalspec.set_xlabel(TOTALSPEC_XLABEL, fontsize=AXISLABELSIZE)

ax_fullspec.set_ylabel(FULLSPEC_YLABEL, fontsize=AXISLABELSIZE)
ax_decayspec.set_ylabel(DECAYSPEC_YLABEL, fontsize=AXISLABELSIZE)
ax_totalspec.set_ylabel(TOTALSPEC_YLABEL, fontsize=AXISLABELSIZE)

ax_fullspec.tick_params(size=TICKSIZE,width=TICKWIDTH,labelsize=TICKLABELSIZE)
ax_decayspec.tick_params(size=TICKSIZE,width=TICKWIDTH,labelsize=TICKLABELSIZE)
ax_totalspec.tick_params(size=TICKSIZE,width=TICKWIDTH,labelsize=TICKLABELSIZE)

ax_fullspec.tick_params(bottom=True,left=True,top=True,right=True)
ax_decayspec.tick_params(bottom=True,left=True,top=True,right=True)
ax_totalspec.tick_params(bottom=True,left=True,top=True,right=True)

ax_fullspec.grid(False,which="both")
ax_decayspec.grid(False,which="both")
ax_totalspec.grid(False,which="both")

fig_fullspec.tight_layout()
fig_decayspec.tight_layout()
fig_totalspec.tight_layout()

if(SAVE):
    fig_weights.savefig(TITLE_WEIGHTS+".svg")
    fig_fullspec.savefig(TITLE_FULLSPEC+".svg")
    fig_decayspec.savefig(TITLE_DECAYSPEC+".svg")
    fig_totalspec.savefig(TITLE_TOTALSPEC+".svg")
    print("SAVED!")
else:
    print("NOT SAVED!")

fig_weights.show()
fig_fullspec.show()
fig_decayspec.show()
fig_totalspec.show()

#%%

#############################################################
######### PLOT SIMPLIFIED FREEZEOUT SURFACE ############
#############################################################

from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit
from scipy.integrate import quad
from scipy.optimize import minimize

TITLE = "Images/freezeout_simplified"

LABELEXACT = "exact freezeout geometry"
LABELSIMPLE = "simplified freezeout geometry"

SAVE = False
# SAVE = True

FIGSIZE = (7,7)
LINEWIDTH = 2.5
TICKLABELSIZE = 25
AXISLABELSIZE = 30 
TICKSIZE = 10
TICKWIDTH = 2

ARROW_WIDTH = 0.3
ARROW_LENGTH = 0.2
ARROW_OVERHANG = 0.1
rmax = 11
taumax = 15

df = pd.read_csv("./../../Mathematica/data/ExampleFreezeOutCorrected.csv")

alphas = df["alpha"].to_numpy()
taus = df["tau"].to_numpy()
rs = df["r"].to_numpy()

tau = CubicSpline(alphas, taus)
r = CubicSpline(alphas, rs)

tau0 = tau(np.pi/2)

def d(alpha, R, T):
    if(alpha <= np.arctan2(R,T)):
        return T/np.cos(alpha)
    if(np.arctan2(R,T) < alpha):
        return R/np.sin(alpha)
    
def Dd(alpha, R, T):
    if(alpha <= np.arctan2(R,T)):
        return np.tan(alpha) * T/np.cos(alpha)
    if(np.arctan2(R,T) < alpha):
        return -R/(np.sin(alpha)*np.tan(alpha))

def myfit_tau(alpha, R, T):
    return d(alpha, R, T) * np.cos(alpha)

def myfit_Dtau(alpha, R, T):
    # return Dd(alpha, R, T) * np.cos(alpha) - d(alpha, R, T) * np.sin(alpha)
    if(alpha <= np.arctan2(R,T)):
        return 0
    if(np.arctan2(R,T) < alpha):
        # return -R/(np.tan(alpha)**2) - R
        return -R/np.sin(alpha)**2

def myfit_r(alpha, R, T):
    x = d(alpha, R, T)
    return d(alpha, R, T) * np.sin(alpha)

def myfit_Dr(alpha, R, T):
    # return Dd(alpha, R, T) * np.sin(alpha) + d(alpha, R, T) * np.cos(alpha)
    if(alpha <= np.arctan2(R,T)):
        # return T * np.tan(alpha)**2 + T
        return T/np.cos(alpha)**2
    if(np.arctan2(R,T) < alpha):
        return 0

def local_cost(alpha, R, T):
    return (tau(alpha) - tau0 - myfit_tau(alpha, R, T))**2 + (r(alpha) - myfit_r(alpha, R, T))**2

def mycost(x):
    R, T = x
    return quad(local_cost, 0, np.pi/2,args=(R,T))[0]

# FIND OPTIMAL PARAMETERS FOR RECTANGULAR REGION
result = minimize(mycost, (9, 13))
Ropt, Topt = result.x

mytaus = np.array([myfit_tau(a,Ropt, Topt) + tau0 for a in alphas])
myrs = np.array([myfit_r(a,Ropt, Topt) for a in alphas])

myDtaus = np.array([myfit_Dtau(a,Ropt, Topt) for a in alphas])
myDrs = np.array([myfit_Dr(a,Ropt, Topt) for a in alphas])

numDtaus = (mytaus[1:] - mytaus[:-1])/(alphas[1]-alphas[0])
numDrs = (myrs[1:] - myrs[:-1])/(alphas[1]-alphas[0])

dataDtaus = df["Dtau"].to_numpy()
dataDrs = df["Dr"].to_numpy()

datanumDtaus = (taus[1:] - taus[:-1])/(alphas[1]-alphas[0])
datanumDrs = (rs[1:] - rs[:-1])/(alphas[1]-alphas[0])


# PLOT
# fig, ax = plt.subplots()
# ax.plot(alphas, taus,c="b",marker="",label=r"$\tau$")
# ax.plot(alphas, mytaus,c="b",ls="--",marker="")
# ax.plot(alphas, rs,c="r",marker="",label=r"$r$")
# ax.plot(alphas, myrs,c="r",ls="--",marker="")
# ax.legend()
# plt.show()

# fig, ax = plt.subplots()
# ax.plot(alphas[:-1], numDtaus,c="b",marker="",label=r"$D\tau$, finite diff")
# ax.plot(alphas, myDtaus,c="b",ls="--",marker="",label=r"$D\tau$, analytic")
# ax.plot(alphas[:-1], datanumDtaus,c="b",ls="--",marker="x",label=r"$D\tau$, data finite diff")
# ax.plot(alphas, dataDtaus,c="b",ls="-.",marker="",label=r"$D\tau$, data")

# ax.plot(alphas[:-1], numDrs,c="r",marker="",label=r"$Dr$, finite diff")
# ax.plot(alphas, myDrs,c="r",ls="--",marker="",label=r"$Dr$, analytic")
# ax.plot(alphas[:-1], datanumDrs,c="r",ls="--",marker="x",label=r"$Dr$, data finite diff")
# ax.plot(alphas, dataDrs,c="r",ls="-.",marker="",label=r"$Dr$, data")
# ax.legend()
# plt.show()

fig, ax = plt.subplots(figsize=FIGSIZE)
ax.plot(rs,taus,marker="",c="b")
# ax.plot(rs,taus,marker="",label=LABELEXACT,c="b")
# ax.plot(myrs, mytaus,marker="",label=LABELSIMPLE,c="r")
# ax.legend()
plt.show()

ax.hlines(tau0,0,rmax,color="k",ls="--",lw=LINEWIDTH)
ax.hlines(Topt,0,Ropt,color="b",lw=LINEWIDTH,ls="--")
ax.vlines(Ropt,tau0,Topt,color="b",lw=LINEWIDTH,ls="--")

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
ax.set_xticks([Ropt])
ax.set_yticks([tau0,Topt])
ax.set_xticklabels([r"$\mathcal{R}$"])
ax.set_yticklabels([r"$\tau_0$",r"$\mathcal{T}$"])
ax.xaxis.set_tick_params(width=2,length=7,direction='out')
ax.yaxis.set_tick_params(width=2,length=7,direction='out')
ax.tick_params(size=TICKSIZE,width=TICKWIDTH,labelsize=TICKLABELSIZE)
ax.xaxis.set_ticks_position("bottom")
ax.yaxis.set_ticks_position("left")

ax.set_xlabel(r"$r$",rotation=0,size=AXISLABELSIZE)
ax.set_ylabel(r"$\tau$",rotation=0,size=AXISLABELSIZE)

ax.xaxis.set_label_coords(0.95,0)
ax.yaxis.set_label_coords(-0.03,0.95)

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

if(SAVE):
    fig.savefig(TITLE+".svg")
    print("SAVED!")
else:
    print("NOT SAVED!")

fig.show()

#%%
#############################################################
#### COMPARE ONE SPECTRUM WITH ANALYTIC PARAMETRIZATION #####
#############################################################

import scipy.special as spc
from matplotlib.widgets import Button, Slider

SAVE = False
# SAVE = True

TITLESPEC = "Images/Analytic_Compare"

FIGSIZE = (7,7)
LINEWIDTH = 2
MARKERSIZE = 50
LEGENDSIZE = 20
TICKLABELSIZE = 20
AXISLABELSIZE = 30
TICKSIZE = 10
TICKWIDTH = 2

LABELEXACT = "numerically exact\nspectrum"
LABELANALYT = "approximated\n spectrum"

folder = "Data/spec_real_constfield_m140_20250116_115949/spec_20250121_103632_00/"

GeVtoIfm = 5.0677
fmtoIGeV = GeVtoIfm
tau0 = 0.4
Topt = 12.16
Ropt = 8.52

R0 = Ropt
T0 = Topt
Amax = 0.2
Bmax = 0.1
m0 = 0.14

# good fit values for Data/spec_real_constfield_m140_20250116_115949/spec_20250121_103632_00
# scale, A1, B1, A2, B2 = -1.74, 0.1854, -0.0224, 0.1875, 0.0453

def w(m,p):
    return np.sqrt(m**2 + p**2)

def JT(p,m,Tau,R, A,B):
    return 2 * (np.pi**2) * Tau * R * (fmtoIGeV)**2 / p * (
        B * (-spc.y0(Tau * w(m,p) * GeVtoIfm) + 1j * spc.j0(Tau * w(m,p) * GeVtoIfm)) + 
        A * (-spc.y1(Tau * w(m,p) * GeVtoIfm) + 1j * spc.j1(Tau * w(m,p) * GeVtoIfm)) * w(m,p)
    ) * spc.j1(R * p * GeVtoIfm)

def JR(p,m,Tau,R, A,B):
    return -2 * (np.pi**2) * R * (fmtoIGeV)**2 / w(m,p) * (
        B * spc.j0(R * p * GeVtoIfm) + 
        A * spc.j1(R * p * GeVtoIfm) * p
    ) * (
        Tau * (-spc.y1(Tau * w(m,p) * GeVtoIfm) + 1j * spc.j1(Tau * w(m,p) * GeVtoIfm)) -
        tau0 * (-spc.y1(tau0 * w(m,p) * GeVtoIfm) + 1j * spc.j1(tau0 * w(m,p) * GeVtoIfm))
    )

def specan(p,m, Tau, R, AT, BT, AR, BR):
    return 0.5 * 1/(2*np.pi)**3 * np.abs(JT(p,m, Tau, R, AT, BT) + JR(p,m, Tau, R, AR, BR))**2

fig_init = plt.figure(figsize=FIGSIZE)
gs = gridspec.GridSpec(nrows=2, ncols=1, hspace=0.0)
ax_init_f, ax_init_df = fig_init.add_subplot(gs[0]), fig_init.add_subplot(gs[1])
ax_init_f.sharex(ax_init_df)
fig_spec, ax_spec = plt.subplots(figsize=FIGSIZE)

df_spec = pd.read_csv(folder+"spectr.txt",comment="#")
df_init_f = pd.read_csv(folder+"field0.txt",comment="#")
df_init_df = pd.read_csv(folder+"field0_deriv.txt",comment="#") 

df_spec = pd.read_csv(folder+"spectr.txt",comment="#") 
df_init_f = pd.read_csv(folder+"field0.txt",comment="#")
df_init_df = pd.read_csv(folder+"field0_deriv.txt",comment="#") 

x_init_f = df_init_f["alpha"].to_numpy()
x_init_df = df_init_df["alpha"].to_numpy()
x_spec = df_spec["pT"].to_numpy() 

y_init_f = df_init_f["field0Re"].to_numpy()
y_init_df = df_init_df["Dfield0Re"].to_numpy()
y_spec = df_spec["abs2Re"].to_numpy()

ax_init_f.plot(x_init_f, y_init_f,lw=LINEWIDTH,marker="")
ax_init_df.plot(x_init_df, y_init_df,lw=LINEWIDTH,marker="")

ax_spec.plot(x_spec, y_spec, marker="",lw=LINEWIDTH,c="b",label=LABELEXACT)
line, = ax_spec.plot(x_spec, specan(x_spec,m0,Topt, Ropt, 0.08, 0.02, 0.08, 0.02), marker="",lw=LINEWIDTH,c="r",label=LABELANALYT)

ax_spec.set_yscale("log")
ax_spec.grid(False, which="both")
ax_spec.set_ylabel(r"$(2\pi p_T)^{-1}dN/(dp_Td\eta_p)\ [{GeV}^{-2}]$",fontsize=AXISLABELSIZE)
ax_spec.set_xlabel(r"$p_T\ [GeV]$",fontsize=AXISLABELSIZE)

ax_spec.legend(fontsize=LEGENDSIZE,fancybox=True, framealpha=0.85,shadow=False)

figsl = plt.figure()
ax1 = figsl.add_axes([0.1, 0.1, 0.03, 0.8])
sl1 = Slider(
    ax=ax1,
    valmin=0,
    valmax=20,
    label="Tau",
    valinit=T0,
    orientation="vertical"
)

ax2 = figsl.add_axes([0.2, 0.1, 0.03, 0.8])
sl2 = Slider(
    ax=ax2,
    valmin=0,
    valmax=10,
    label="R",
    valinit=R0,
    orientation="vertical"
)

ax3 = figsl.add_axes([0.3, 0.1, 0.03, 0.8])
sl3 = Slider(
    ax=ax3,
    valmin=-5,
    valmax=5,
    label="scale",
    valinit=0,
    orientation="vertical"
)

ax4 = figsl.add_axes([0.4, 0.1, 0.03, 0.8])
sl4 = Slider(
    ax=ax4,
    valmin=-Amax,
    valmax=Amax,
    label="A1",
    valinit=1,
    orientation="vertical"
)

ax5 = figsl.add_axes([0.5, 0.1, 0.03, 0.8])
sl5 = Slider(
    ax=ax5,
    valmin=-Bmax,
    valmax=Bmax,
    label="B1",
    valinit=0,
    orientation="vertical"
)

ax6 = figsl.add_axes([0.6, 0.1, 0.03, 0.8])
sl6 = Slider(
    ax=ax6,
    valmin=-Amax,
    valmax=Amax,
    label="A2",
    valinit=1,
    orientation="vertical"
)

ax7 = figsl.add_axes([0.7, 0.1, 0.03, 0.8])
sl7 = Slider(
    ax=ax7,
    valmin=-Bmax,
    valmax=Bmax,
    label="B2",
    valinit=0,
    orientation="vertical"
)

ax8 = figsl.add_axes([0.8, 0.1, 0.03, 0.8])
sl8 = Slider(
    ax=ax8,
    valmin=0,
    valmax=1,
    label="m",
    valinit=m0,
    orientation="vertical"
)

# The function to be called anytime a slider's value changes
def update(val):

    myspec = specan(x_spec,sl8.val, sl1.val, sl2.val,sl4.val, sl5.val, sl6.val, sl7.val)
    line.set_ydata(np.exp(sl3.val)*myspec)

    fig_spec.canvas.draw_idle()

sl1.on_changed(update)
sl2.on_changed(update)
sl3.on_changed(update)
sl4.on_changed(update)
sl5.on_changed(update)
sl6.on_changed(update)
sl7.on_changed(update)
sl8.on_changed(update)

fig_init.tight_layout()
fig_spec.tight_layout()

if(SAVE):
    fig_spec.savefig(TITLESPEC+".svg")
    print("SAVED!")
else:
    print("NOT SAVED!")

fig_spec.show()

# %%
#############################################################
###### PLOT DCC CONTRIBUTION TO THE TOTAL SPECTRUM #########
#############################################################

import matplotlib.ticker

# specfile = sorted(glob.glob("Data/spec_real_constfield_m140_20250116_115949/*/"))[0]
# specfile = sorted(glob.glob("Data/spec_real_consteps_m140_20250117_113521/*/"))[-5]
# specfile = sorted(glob.glob("Data/spec_real_taudep_m140_20250117_162422/*/"))[-5]
# specfile = sorted(glob.glob("Data/spec_real_constfield_m280_20250122_110401/*/"))[0]
# specfile = sorted(glob.glob("Data/spec_real_consteps_m280_20250122_110422/*/"))[-5]
# specfile = sorted(glob.glob("Data/spec_real_taudep_m280_20250122_110439/*/"))[-1]
# specfile = sorted(glob.glob("Data/spec_real_constfield_m420_20250122_111331/*/"))[0]
specfile = sorted(glob.glob("Data/spec_real_consteps_m420_20250122_111354/*/"))[-5]
# specfile = sorted(glob.glob("Data/spec_real_taudep_m420_20250122_111415/*/"))[-1]

SCALE = 2.5

SAVETITLE = "Images/"+specfile[:-1].replace("Data/","").replace("/spec","")+("_%.1f"%(SCALE)).replace(".","-")
print(SAVETITLE)

SAVE = False
# SAVE = True

FIGSIZE = (7,7)
LINEWIDTH = 2
TICKLABELSIZE = 20
AXISLABELSIZE = 25
TICKSIZE = 10
TICKWIDTH = 2

MARKERSIZE = 50
LEGENDSIZE = 20

COL = (1,0,0)

XLIM = (0,0.7)
YLIM = (8e1,4e4)

DFIELD_XLABEL = r"$\alpha$"
SPEC_XLABEL = r"$p_T\ [GeV]$"

FIELD_YLABEL = r"$\phi_{DCC}\ [GeV]$"
DFIELD_YLABEL = r"$n^\mu\partial_\mu\phi_{DCC}\ [GeV]$"
SPEC_YLABEL = r"$(2\pi p_T)^{-1}dN/(dp_Td\eta_p)\ [GeV^{-2}]$"

with open(specfile+"spectr.txt") as f:
    mass = float(f.readlines()[3].replace("# particle mass:\t","").replace("\n",""))

LEGEND = "DCC, "+r"$m={%.2f}\ GeV$"%mass
###
def data_to_bins(datax, datay, bins_lower, bins_upper):
    bins_x = ( bins_upper + bins_lower ) / 2

    bins_all = [[] for i in range(len(bins_x))]

    for (x,y) in list(zip(datax, datay)):
        try:
            binidx = np.where((x >= bins_lower).astype(int) * (x <= bins_upper).astype(int))[0][0]
        except:
            continue
        bins_all[binidx].append(y)

    bins_y = []
    for bin_y in bins_all:
        bins_y.append(np.mean(bin_y))

    return np.array(bins_x), np.array(bins_y)
###

df_alice = pd.read_csv("./../../Mathematica/data/HEPData-ins1222333-v1-Table_1.csv",comment="#")
pTs_alice = df_alice["PT [GEV]"].to_numpy()[:41].astype(float)
spec_alice = df_alice["(1/Nev)*(1/(2*PI*PT))*D2(N)/DPT/DYRAP [GEV**-2]"].to_numpy()[:41].astype(float)
spec_alice_err = df_alice["sys +"].to_numpy()[:41].astype(float)+df_alice["stat +"].to_numpy()[:41].astype(float)
bins_lower = df_alice["PT [GEV] LOW"].to_numpy()[:41].astype(float)
bins_upper = df_alice["PT [GEV] HIGH"].to_numpy()[:41].astype(float)

df = pd.read_csv("./../../Mathematica/data/pionListBig.txt")
pTs_fluidum= np.array(df.keys()).astype(float)
spec_fluidum_pi0 = df.iloc[0].to_numpy()
spec_fluidum_piplus = df.iloc[1].to_numpy()

def interplog(datax, datay, newx):
    return np.exp(scipy.interpolate.interp1d(datax, np.log(datay))(newx))

def mylogpolyfit(datax, datay):
    popt = np.polyfit(datax,np.log(datay),10)
    return np.poly1d(popt)

###
fig_init = plt.figure(figsize=FIGSIZE)
gs = gridspec.GridSpec(nrows=2, ncols=1, hspace=0.0)
ax_init_f, ax_init_df = fig_init.add_subplot(gs[0]), fig_init.add_subplot(gs[1])
ax_init_f.sharex(ax_init_df)
fig_spec, ax_spec = plt.subplots(figsize=FIGSIZE)

df_spec = pd.read_csv(specfile+"spectr.txt",comment="#") 
df_init_f = pd.read_csv(specfile+"field0.txt",comment="#")
df_init_df = pd.read_csv(specfile+"field0_deriv.txt",comment="#") 

x_init_f = df_init_f["alpha"].to_numpy()
x_init_df = df_init_df["alpha"].to_numpy()
x_spec = df_spec["pT"].to_numpy() 

y_init_f = df_init_f["field0Re"].to_numpy()
y_init_df = df_init_df["Dfield0Re"].to_numpy()
y_spec = df_spec["abs2Re"].to_numpy()

###
ax_init_f.plot(x_init_f, np.sqrt(SCALE) * y_init_f,lw=LINEWIDTH,marker="")
ax_init_df.plot(x_init_df, np.sqrt(SCALE) * y_init_df,lw=LINEWIDTH,marker="")

fluidum_loginterp = mylogpolyfit(pTs_fluidum, spec_fluidum_piplus)
y_fluidum_interp = np.exp(fluidum_loginterp(x_spec))
ax_spec.fill_between(x_spec,y_fluidum_interp,y_fluidum_interp+SCALE*y_spec,facecolor=(*COL,0.2),edgecolor=COL,lw=LINEWIDTH,label=LEGEND)
ax_spec.plot(x_spec, y_fluidum_interp,marker="",label=r"Fluid$\mathit{u}$m",lw=LINEWIDTH)

# ax_spec.plot(pTs_fluidum, spec_fluidum_piplus)
ax_spec.errorbar(pTs_alice, spec_alice, spec_alice_err,label="ALICE",c="b",fmt="o",markersize=MARKERSIZE/12,lw=LINEWIDTH)


ax_init_f.set_ylabel(FIELD_YLABEL, fontsize=AXISLABELSIZE)
ax_init_df.set_ylabel(DFIELD_YLABEL, fontsize=AXISLABELSIZE)
ax_spec.set_ylabel(SPEC_YLABEL, fontsize=AXISLABELSIZE)

ax_init_df.set_xlabel(DFIELD_XLABEL, fontsize=AXISLABELSIZE)
ax_spec.set_xlabel(SPEC_XLABEL, fontsize=AXISLABELSIZE)

xticks = [0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2]
xticklabels = [r"$0$",r"$\pi/8$",r"$\pi/4$",r"$3\pi/8$",r"$\pi/2$"]
ax_init_f.set_xticks(xticks, xticklabels)
# ax_init_df.set_xticks(xticks, xticklabels)
ax_init_f.set_xticklabels(xticklabels,visible=False)

ax_spec.set_xlim(*XLIM)
ax_spec.set_ylim(*YLIM)

ax_spec.set_yscale("log")
locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8))
ax_spec.yaxis.set_minor_locator(locmin)
ax_spec.yaxis.set_minor_formatter(matplotlib.ticker.LogFormatterSciNotation(base=10,labelOnlyBase=False,minor_thresholds=(5,2.5))) # all of the data spans 5 decades and we want to see all minor ticks if we zoom in on a region of 2.5 decades

yticks = ax_init_df.get_yticks()
yticklabels = ax_init_df.get_yticklabels()
yticklabels[-1] = ""
ax_init_df.set_yticks(yticks,yticklabels)

ax_init_f.tick_params(size=TICKSIZE,width=TICKWIDTH,labelsize=TICKLABELSIZE)
ax_init_df.tick_params(size=TICKSIZE,width=TICKWIDTH,labelsize=TICKLABELSIZE)
ax_spec.tick_params(size=TICKSIZE,width=TICKWIDTH,labelsize=TICKLABELSIZE)
ax_spec.tick_params(size=0.7*TICKSIZE,width=0.7*TICKWIDTH,labelsize=0.7*TICKLABELSIZE,which="minor")

ax_init_f.tick_params(bottom=True,left=True,top=True,right=True)
ax_init_df.tick_params(bottom=True,left=True,top=True,right=True)
ax_spec.tick_params(bottom=True,left=True,top=True,right=True)

ax_init_f.set_xlim(0,np.pi/2)
# ax_init_df.set_xlim(0,np.pi/2)

ax_init_f.grid(False,which="both")
ax_init_df.grid(False,which="both")
ax_spec.grid(False,which="both")

ax_spec.legend(fontsize=LEGENDSIZE,fancybox=True, framealpha=0.85,shadow=False)

fig_init.tight_layout()
fig_spec.tight_layout()

if(SAVE):
    fig_spec.savefig(SAVETITLE+".svg")
    print("SAVED!")
else:
    print("NOT SAVED!")

fig_init.show()
fig_spec.show()
# %%
