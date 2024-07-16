import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import glob

### INTERPOLATED FUNCTIONS
# files_interp = glob.glob("data/*interp.txt")
# files_samp = glob.glob("data/*samp.txt")
# print(list(zip(files_interp,files_samp)))
# for (file_i, file_s) in list(zip(files_interp,files_samp)):
#     df_i = pd.read_csv(file_i,header=None,sep=";")
#     df_s = pd.read_csv(file_s,header=None,sep=";")

#     fig, ax = plt.subplots()
#     ax.scatter(*df_s.to_numpy().T,s=20,c="r",marker="+",label=file_s)
#     ax.plot(*df_i.to_numpy().T,label=file_i)
#     ax.legend()
# plt.show()

### INITIAL CONDITIONS & SPECTRUM
df_initial = pd.read_csv("data/field0.txt",header=None,sep=";")
df_initial_deriv = pd.read_csv("data/field0_deriv.txt",header=None,sep=";")
df_spec = pd.read_csv("data/spectr.txt",header=None,sep=";")
df_spec_anti = pd.read_csv("data/spectr_anti.txt",header=None,sep=";")

fig, ((ax0,ax1),(ax2,ax3)) = plt.subplots(ncols=2,nrows=2,figsize=(10,10))

ax0.plot(df_initial.to_numpy().T[0],df_initial.to_numpy().T[1],label="Real",lw=3)
ax0.plot(df_initial.to_numpy().T[0],df_initial.to_numpy().T[2],label="Imag")
ax1.plot(df_initial_deriv.to_numpy().T[0],df_initial_deriv.to_numpy().T[1],label="Real",lw=3)
ax1.plot(df_initial_deriv.to_numpy().T[0],df_initial_deriv.to_numpy().T[2],label="Imag")
ax2.plot(df_spec.to_numpy().T[0],df_spec.to_numpy().T[1],label="spectr")
ax3.plot(df_spec_anti.to_numpy().T[0],df_spec_anti.to_numpy().T[1],label="spectr anti")

ax0.legend()
ax0.grid()
ax1.legend()
ax1.grid()
ax2.set_yscale('log')
ax2.legend()
ax2.grid()
ax3.set_yscale('log')
ax3.legend()
ax3.grid()
plt.show()

# df = pd.read_csv("../../Mathematica/data//ExampleFreezeOut.csv")
# alphas = df["alpha"].to_numpy()
# taus = df["tau"].to_numpy()
# rs = df["r"].to_numpy()
# Dtaus = df["Dtau"].to_numpy()
# Drs = df["Dr"].to_numpy()
# urs = df["ur"].to_numpy()
# utaus = df["utau"].to_numpy()

# plt.plot(alphas,taus,label="tau")
# plt.legend()
# plt.show()
# plt.plot(alphas,rs,label="r")
# plt.legend()
# plt.show()
# plt.plot(alphas,Dtaus,label="Dtau")
# plt.legend()
# plt.show()
# plt.plot(alphas,Drs,label="Dr")
# plt.legend()
# plt.show()
# plt.plot(alphas,urs,label="ur")
# plt.legend()
# plt.show()
# plt.plot(alphas,utaus,label="utau")
# plt.legend()
# plt.show()