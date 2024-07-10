import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import glob

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

df_spec = pd.read_csv("data/spectr.txt",header=None,sep=";")

fig, ax = plt.subplots()
ax.set_yscale('log')
ax.scatter(*df_spec.to_numpy().T)
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