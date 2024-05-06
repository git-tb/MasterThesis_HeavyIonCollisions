from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

plt.style.use("mplstyles/myclassic_white.mplstyle")

df = pd.read_csv("../Mathematica/data/ExampleFreezeOut.csv",names=["v","???","alpha","tau","r","tauDeriv","rDeriv"])
print(df)

# ============================================
# ===============PARAMETERS===================
# ============================================
msigma = 1              # sigma mass
epsilon = 0.160054      # energy density at freeze out
fpi = 1                 # pion decay constant = sigma vev
theta0 = 0              # integration constant in complex phase


# ============================================
alpha = df["alpha"].to_numpy()
v = df["v"].to_numpy()
tau = df["tau"].to_numpy()
r = df["r"].to_numpy()
tauDeriv = df["tauDeriv"].to_numpy()
rDeriv = df["rDeriv"].to_numpy()

dalpha = alpha[1:] - alpha[:-1]
utau = np.cosh(np.arctanh(v))
ur = np.sinh(np.arctanh(v))

Npoints = alpha.shape[0]

# ============================================
chi2 = (-msigma**2 + np.sqrt(12*epsilon/fpi + 2*msigma**2))/6
theta = theta0 + np.sqrt(chi2)*np.cumsum(dalpha*(tauDeriv*utau - rDeriv*ur)[:-1])
rho = np.sqrt(fpi * (2*chi2 + msigma**2)/(msigma**2) )

# ============================================
# =================PLOT=======================
# ============================================
plt.scatter(np.arange(Npoints-1),np.cos(theta))
print("rho={rho:.5f}".format(rho=rho))
plt.show()
