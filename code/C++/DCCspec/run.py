import subprocess

#%%

pTmax = 1
NpT = 100
initpath = "data/init_20240725_151035/initialfields_pi0.csv"

subprocess.run(args=[
    "./bin/spec",
    "--pTmax=%".format(pTmax),
    "--NpT=%".format(NpT),
    "--initpath="+initpath,
])