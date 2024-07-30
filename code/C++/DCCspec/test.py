import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

dfp = pd.read_csv("data/decayspec_20240729_224916/primespec_interp.txt",comment="#")
plt.plot(dfp["q"].to_numpy(), dfp["primespecRe"].to_numpy())
plt.gca().set_yscale("log")
plt.show()

df = pd.read_csv("data/decayspec_20240729_224916/decayspec.txt",comment="#")
spec = df["finalspecRe"].to_numpy()
ps = df["p"].to_numpy()

plt.plot(ps,spec)
plt.gca().set_yscale("log")
plt.show()