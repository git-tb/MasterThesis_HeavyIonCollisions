#%%
import subprocess
import pandas as pd
import numpy as np
from IPython import get_ipython
import glob
import os
import time

#%%
###############################################################
###############################################################
# COMPUTE SPECTRA FROM LIST OF INITIAL DATA
###############################################################
###############################################################

pTmax = 2
NpT = 200

# parentdir = "Data/init_real_constfield_m140_20250116_115949"
# parentdir = "Data/init_comp_constfield_m140_varphase_20250116_120129"
# parentdir = "Data/init_comp_constfield_m140_varphase_20250117_181049"
parentdir = "Data/init_comp_constfield_m140_varphase_20250117_181304"
# parentdir = "Data/init_real_consteps_m140_20250117_113521"
# parentdir = "Data/init_real_consteps_varm_20250116_114850"
# parentdir = "Data/init_real_taudep_m140_20250117_162422"
initfiles = sorted(glob.glob(parentdir+"/*.csv"))

# define masses by hand...
masses = 140 * np.ones(len(initfiles))

# ...or read from each file
for (i,file) in enumerate(initfiles):
    with open(file) as f:
        lines = f.readlines()
        j = [k for (k,x) in enumerate(["mass" in l for l in lines]) if x][0]
        masses[i]=float(lines[j].replace("# mass:","").replace("(GeV)\n",""))

for (i,file) in enumerate(initfiles):
    time.sleep(1)
    with open(file) as f:
        result = subprocess.run(args=[
            "./bin/specV2",
            "--m=%f"%(masses[i]),
            "--pTmax=%f"%(pTmax),
            "--NpT=%d"%(NpT),
            "--epsabs=0",
            "--epsrel=1e-5",
            "--iter=10000",
            "--parentdir=%s"%("Data"),
            "--initpath=%s"%(file)
    ])
    print(result)

newdir = parentdir.replace("init","spec")
subprocess.run(args=["mkdir",newdir])
lastspecs = glob.glob("Data/spec_????????_??????")
for spec in lastspecs:
    subprocess.run(args=["mv",spec,newdir])