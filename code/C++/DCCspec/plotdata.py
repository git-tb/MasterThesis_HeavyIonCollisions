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
    print("SAVED!")
else:
    print("NOT SAVED!")

fig_spec.show()
fig_diff.show()

# %%
#############################################################
######### PLOT 1-PARAMETER FAMILY OF REAL FIELDS ############
#############################################################

# parentdir = "Data/decay_real_constfield_m600_varPT"
# parentdir = "Data/spec_real_consteps_m140_20250117_113521"
# parentdir = "Data/spec_real_consteps_varm_20250116_114850"
# parentdir = "Data/spec_real_taudep_m140_20250117_162422"
parentdir = "Data/spec_real_constfield_varm"
folders = sorted(glob.glob(parentdir+"/*/"))

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

# LABELBYMASS = False
LABELBYMASS = True

# COMPARE = False
COMPARE = True
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

SAVE = False
# SAVE = True

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
ax_init_f_re.set_xlim(0,np.pi/2)
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

SAVE = False
# SAVE = True

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

if(LABELBYMASS):
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

SAVE = False
# SAVE = True

FIGSIZE = (7,7)
LINEWIDTH = 1.5
TICKLABELSIZE = 20
AXISLABELSIZE = 25 
TICKSIZE = 10
TICKWIDTH = 2

DECAYSPEC_XLABEL = r"$p_T\ [GeV]$"
DECAYSPEC_YLABEL = r"$(2\pi p_T)^{-1}dN_{\pi^\pm}/(dp_Td\eta_p)\ [GeV^{-2}]$"

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

# LABELLC = False
LABELLC = True

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
    # return -1/np.pi*np.imag(Delta(k**2))*(k<=1)
    return -1/np.pi*np.imag(Delta(k**2))

weights = 2*masses*S(masses)*np.ptp(masses)/len(masses)
weights /= np.sum(weights)

fig_weights, ax_weights = plt.subplots(figsize=FIGSIZE)
ax_weights.plot([0.28,*masses],S(np.array([0.28,*masses])),c="b",lw=2*LINEWIDTH,marker="",label=WEIGHTS_LEGEND)
ax_weights.legend(fontsize=LEGENDSIZE,fancybox=True, framealpha=0.85,shadow=False)

ax_weights.set_xlabel(WEIGHTS_XLABEL, fontsize=AXISLABELSIZE)
ax_weights.set_ylabel(WEIGHTS_YLABEL, fontsize=AXISLABELSIZE)
ax_weights.tick_params(size=TICKSIZE,width=TICKWIDTH,labelsize=TICKLABELSIZE)
ax_weights.tick_params(bottom=True,left=True,top=True,right=True)
ax_weights.grid(False,which="both")
fig_weights.tight_layout()

### GO TO PLOTTING
fig_fullspec, ax_fullspec = plt.subplots(figsize=FIGSIZE)
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

if(COMPARE):
    ax_fullspec.plot(pT_compare, spec_compare,lw=2*LINEWIDTH,label=r"ALICE$-$FluiduM",c="k",ls=LS_COMPARE,marker="")
ax_fullspec.legend(fontsize=LEGENDSIZE,fancybox=True, framealpha=0.85,shadow=False)

if(LABELBYMASS):
    cax = ax_decayspec.inset_axes(CMAP_LBWH)
    cbar = fig_decayspec.colorbar(lc_decayspec, cax=cax)
    cbar.set_label(LC_LABEL, fontsize=CMAP_LABELSIZE)
    cbar.ax.tick_params(labelsize=CMAP_TICKSIZE)

ax_fullspec.set_yscale("log")
ax_decayspec.set_yscale("log")

ax_fullspec.set_xlabel(DECAYSPEC_XLABEL, fontsize=AXISLABELSIZE)
ax_decayspec.set_xlabel(DECAYSPEC_XLABEL, fontsize=AXISLABELSIZE)

ax_fullspec.set_ylabel(FULLSPEC_YLABEL, fontsize=AXISLABELSIZE)
ax_decayspec.set_ylabel(DECAYSPEC_YLABEL, fontsize=AXISLABELSIZE)

ax_fullspec.tick_params(size=TICKSIZE,width=TICKWIDTH,labelsize=TICKLABELSIZE)
ax_decayspec.tick_params(size=TICKSIZE,width=TICKWIDTH,labelsize=TICKLABELSIZE)

ax_fullspec.tick_params(bottom=True,left=True,top=True,right=True)
ax_decayspec.tick_params(bottom=True,left=True,top=True,right=True)

ax_fullspec.grid(False,which="both")
ax_decayspec.grid(False,which="both")

fig_fullspec.tight_layout()
fig_decayspec.tight_layout()

if(SAVE):
    print("SAVED!")
else:
    print("NOT SAVED!")

fig_weights.show()
fig_fullspec.show()
fig_decayspec.show()