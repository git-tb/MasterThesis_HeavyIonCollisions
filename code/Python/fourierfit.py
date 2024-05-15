import numpy as np
from matplotlib import pyplot as plt

"""
    Consider a (1+1)-dimensional problem with coordinate
"""

Npoints = 100
alphas = np.linspace(0,np.pi/2,Npoints)
taus = np.sin(alphas)
rs = np.cos(alphas)

pmax = 100
Ncoeffs = Npoints
ps = np.linspace(0,pmax,Ncoeffs)
coeffsRE, coeffsIM = np.ones((2,Ncoeffs))

def phi(x,coeffsRE, coeffsIM):
    return np.exp(1j*1)


plt.plot(np.real(phi()))
