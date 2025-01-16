#%%
import subprocess
import pandas as pd
import numpy as np
from IPython import get_ipython
import glob
import os

#%%
###############################################################
###############################################################
# COMPUTE SPECTRA FROM LIST OF INITIAL DATA
###############################################################
###############################################################

pTmax = 2
NpT = 200

parentdir = "Data/init_real_constfield_m140_20250116_115949"
initfiles = sorted(glob.glob(parentdir+"/*.csv"))

# define masses by hand...
masses = 140 * np.ones(len(initfiles))

# ...or read from each file
for (i,file) in enumerate(initfiles):
    with open(file) as f:
        masses[i]=float(f.readlines()[2].replace("# mass:","").replace("(GeV)\n",""))

for (i,file) in enumerate(initfiles):
    with open(file) as f:
        result = subprocess.run(args=[
            "./bin/specV2",
            "--m=%f"%(masses[i]),
            "--pTmax=%f"%(pTmax),
            "--NpT=%d"%(NpT),
            "--epsabs=0",
            "--epsrel=1e-5",
            "--iter=10000",
            "--initpath=%s"%(file),
            "--parentdir=%s"%("Data")
    ])
    print(result)

newdir = parentdir.replace("init","spec")
lastspecs = glob.glob("Data/spec_????????_??????")
for spec in lastspecs:
    subprocess.run(args=["mv",spec,newdir])
# %%
