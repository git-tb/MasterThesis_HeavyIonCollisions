#%%
import subprocess
import pandas as pd
import numpy as np
from IPython import get_ipython
import glob
import os
import time
import datetime

#%%
###############################################################
###############################################################
# COMPUTE SPECTRA FROM LIST OF INITIAL DATA
###############################################################
###############################################################

pTmax = 2
NpT = 200

# parentdir = "Data/init_real_constfield_m140_20250116_115949"
# parentdir = "Data/init_real_constfield_m280_20250122_110401"
# parentdir = "Data/init_real_constfield_m420_20250122_111331"
# parentdir = "Data/init_comp_constfield_m140_varphase_20250116_120129"
# parentdir = "Data/init_comp_constfield_m140_varphase_20250117_181049"
# parentdir = "Data/init_comp_constfield_m140_varphase_20250117_181304"
# parentdir = "Data/init_comp_constfield_m140_varphase_20250120_094422"
# parentdir = "Data/init_real_consteps_m140_20250117_113521"
# parentdir = "Data/init_real_consteps_m280_20250122_110422"
# parentdir = "Data/init_real_consteps_m420_20250122_111354"
# parentdir = "Data/init_real_consteps_varm_20250116_114850"
# parentdir = "Data/init_real_taudep_m140_20250117_162422"
# parentdir = "Data/init_real_taudep_m280_20250122_110439"
# parentdir = "Data/init_real_taudep_m420_20250122_111415"
initfiles = sorted(glob.glob(parentdir+"/*.csv"))

# define masses by hand...
masses = 0.28 * np.ones(len(initfiles))

# ...or read from each file
for (i,file) in enumerate(initfiles):
    with open(file) as f:
        lines = f.readlines()
        j = [k for (k,x) in enumerate(["mass" in l for l in lines]) if x][0]
        masses[i]=float(lines[j].replace("# mass:","").replace("(GeV)\n",""))

for (i,file) in enumerate(initfiles):
    # time.sleep(1)
    foldername = "spec_{:%Y%m%d_%H%M%S}_".format(datetime.datetime.now())+str(i).zfill(int(np.ceil(np.log10(len(initfiles)))))
    print("save to", foldername)
    result = subprocess.run(args=[
        "./bin/specV2",
        "--m=%f"%(masses[i]),
        "--pTmax=%f"%(pTmax),
        "--NpT=%d"%(NpT),
        "--epsabs=0",
        "--epsrel=1e-5",
        "--iter=10000",
        "--parentdir=%s"%("Data"),
        "--initpath=%s"%(file),
        "--foldername="+foldername
    ])
    print(result)

newdir = parentdir.replace("init","spec")
subprocess.run(args=["mkdir",newdir])
lastspecs = glob.glob("Data/spec_????????_??????_*")
for spec in lastspecs:
    subprocess.run(args=["mv",spec,newdir])

#%%
###############################################################
###############################################################
# COMPUTE SPECTRA FROM SAME INITIAL DATA OVER DIFFERENT PT RANGES
###############################################################
###############################################################

initfile = "Data/init_real_constfield_m600_20250121_105135/init00.csv"

pTmaxs = 1  * np.arange(1,11)
NpTs =   100 * np.arange(1,11)

mass = 0.6

for (i,(pTmax,NpT)) in enumerate(list(zip(pTmaxs,NpTs))):
    foldername = "spec_{:%Y%m%d_%H%M%S}_".format(datetime.datetime.now())+str(i).zfill(int(np.ceil(np.log10(len(NpTs)))))
    print("save to", foldername)
    result = subprocess.run(args=[
        "./bin/specV2",
        "--m=%f"%(mass),
        "--pTmax=%f"%(pTmax),
        "--NpT=%d"%(NpT),
        "--epsabs=0",
        "--epsrel=1e-5",
        "--iter=10000",
        "--parentdir=%s"%("Data"),
        "--initpath=%s"%(initfile),
        "--foldername="+foldername
    ])
    print(result)

newdir = "Data/spec_real_constfield_m%d_varPT"%(1000*mass)
subprocess.run(args=["mkdir",newdir])
lastspecs = glob.glob("Data/spec_????????_??????_*")
for spec in lastspecs:
    subprocess.run(args=["mv",spec,newdir])

#%%
###############################################################
###############################################################
# COMPUTE SPECTRA FROM SAME INITIAL DATA FOR DIFFERENT MASSES
###############################################################
###############################################################

initfile = "Data/init_real_constfield_m600_20250121_105135/init00.csv"

pTmax = 3
NpT = 300

masses = np.linspace(2*0.14,1.5,100)

for (i,mass) in enumerate(masses):
    foldername = "spec_{:%Y%m%d_%H%M%S}_".format(datetime.datetime.now())+str(i).zfill(int(np.ceil(np.log10(len(masses)))))
    print("save to", foldername)
    result = subprocess.run(args=[
        "./bin/specV2",
        "--m=%f"%(mass),
        "--pTmax=%f"%(pTmax),
        "--NpT=%d"%(NpT),
        "--epsabs=0",
        "--epsrel=1e-5",
        "--iter=10000",
        "--parentdir=%s"%("Data"),
        "--initpath=%s"%(initfile),
        "--foldername="+foldername
    ])
    print(result)

newdir = "Data/spec_real_constfield_varm"
subprocess.run(args=["mkdir",newdir])
lastspecs = glob.glob("Data/spec_????????_??????_*")
for spec in lastspecs:
    subprocess.run(args=["mv",spec,newdir])

#%%
###############################################################
###############################################################
# COMPUTE DECAY SPECTRA FROM DIFFERENT PT RANGES OF SAME SPECTRUM
###############################################################
###############################################################

parentdir = "Data/spec_real_constfield_m600_varPT"
initfiles = sorted(glob.glob(parentdir+"/**/spectr.txt"))

ma, mb, mc = 0.6, 0.14, 0.14
pTmax = 2
NpT = 200
B = 1
Q = 1

for (i,file) in enumerate(initfiles):
    foldername = "decay_{:%Y%m%d_%H%M%S}_".format(datetime.datetime.now())+str(i).zfill(int(np.ceil(np.log10(len(initfiles)))))
    print("save to", foldername)
    result = subprocess.run(args=[
        "./bin/decayV2",
        "--ma=%f"%(ma),
        "--mb=%f"%(mb),
        "--mc=%f"%(mc),
        "--pTmax=%f"%(pTmax),
        "--NpT=%d"%(NpT),
        "--epsabs=0",
        "--epsrel=1e-5",
        "--iter=100000",
        "--primespecpath=%s"%(file),
        "--parentdir=%s"%("Data"),
        "--foldername="+foldername,
        "--B=%f"%(B),
        "--Q=%f"%(Q)
    ])
    print(result)

newdir = "Data/decay_real_constfield_m%d_varPT"%(1000*mass)
subprocess.run(args=["mkdir",newdir])
lastspecs = glob.glob("Data/decay_????????_??????_*")
for spec in lastspecs:
    subprocess.run(args=["mv",spec,newdir])

#%%
###############################################################
###############################################################
# COMPUTE DECAY SPECTRA FROM SPECTRA OF DIFFERENT MASSES
###############################################################
###############################################################

parentdir = "Data/spec_real_constfield_varm"
initfiles = sorted(glob.glob(parentdir+"/**/spectr.txt"))

ma, mb, mc = 0.6, 0.14, 0.14
pTmax = 2
NpT = 200
B = 1
Q = 1

for (i,file) in enumerate(initfiles):
    with open(file) as f:
        ma = float(f.readlines()[3].replace("# particle mass:\t",""))
    foldername = "decay_{:%Y%m%d_%H%M%S}_".format(datetime.datetime.now())+str(i).zfill(int(np.ceil(np.log10(len(initfiles)))))
    print("save to", foldername)
    print(ma)
    result = subprocess.run(args=[
        "./bin/decayV2",
        "--ma=%f"%(ma),
        "--mb=%f"%(mb),
        "--mc=%f"%(mc),
        "--pTmax=%f"%(pTmax),
        "--NpT=%d"%(NpT),
        "--epsabs=0",
        "--epsrel=1e-5",
        "--iter=100000",
        "--primespecpath=%s"%(file),
        "--parentdir=%s"%("Data"),
        "--foldername="+foldername,
        "--B=%f"%(B),
        "--Q=%f"%(Q)
    ])
    print(result)

newdir = "Data/decay_real_constfield_varm"
subprocess.run(args=["mkdir",newdir])
lastspecs = glob.glob("Data/decay_????????_??????_*")
for spec in lastspecs:
    subprocess.run(args=["mv",spec,newdir])