import numpy as np
from matplotlib import pyplot as plt
import scipy.optimize
import scipy.special

#########################################################################################################
#########################################################################################################
# region PROBLEM PARAMETERS
m=1
def omega(p):
    return np.sqrt(p**2+m**2)

Npoints = 100
xs = np.linspace(0,10,Npoints)

pmax = 50
Ncoeffs = 4*Npoints
ps = np.linspace(0,pmax,Ncoeffs)
dp = ps[1]-ps[0]
coeffs = 0.2*(np.random.random(2*Ncoeffs)*2-1)
coeffs = 0*coeffs
coeffs[50] = 1
# endregion
#########################################################################################################
#########################################################################################################

#########################################################################################################
#########################################################################################################
# region OPTIMIZATION TARGET
phitarget = np.exp(1j*xs)*np.exp(-(xs-4)**2)

def phi(t,x,ps,coeffs):
    tw = t*omega(ps)
    px = x*ps
    return dp*np.sum(
        (coeffs[::2] + 1j*coeffs[1::2])*np.exp(1j*(tw-px))
    )

def phixs(t,xs,ps,coeffs):
    tw = t*omega(ps)
    px = np.outer(xs,ps)
    expmat = np.exp(1j*(tw-px))
    coeffsmat = coeffs[::2] + 1j*coeffs[1::2]
    return dp*np.matmul(expmat,coeffsmat)

def targetfun(coeffs):
    t = 0
    return(np.sum(np.abs(
        phitarget-phixs(0,xs,ps,coeffs))
        ))
# endregion
#########################################################################################################
#########################################################################################################

def plot(coeffs):
    # field = np.array([phi(0,x,ps,coeffs) for x in xs])
    field = phixs(0,xs,ps,coeffs)
    plt.plot(xs, np.real(field),c="b",ls="-")
    plt.plot(xs, np.imag(field),c="b",ls="-.")
    plt.plot(xs, np.real(phitarget),c="r",ls="-")
    plt.plot(xs, np.imag(phitarget),c="r",ls="-.")
    plt.show()

#########################################################################################################
#########################################################################################################
# region RANDOM OPTIMIZER
def randomoptimizer(coeffs,maxiter,recordx=None,recordf=None):
    for i in range(maxiter):
        ridx = int(np.random.random()*Ncoeffs)
        ramp = (2*np.random.random() - 1)

        oldtarget = targetfun(coeffs)
        coeffs[ridx] += ramp
        newtarget = targetfun(coeffs)

        if newtarget > oldtarget:
            coeffs[ridx] -= ramp
            newtarget = oldtarget
        else:
            pass
        print(i+1,"/",maxiter)

        if(recordx):
            recordx.append(coeffs)
        if not(recordf == None):
            recordf.append(newtarget)

    return coeffs
# endregion
#########################################################################################################
#########################################################################################################

# recordf = []
# coeffs_opt = randomoptimizer(coeffs,1000,recordf=recordf)
# plt.plot(recordf)
# plt.show()

optimizerecordx, optimizerecordf = [], []
# def callbackf(intermediate_result: scipy.optimize.OptimizeResult):
#     print(type(intermediate_result))
#     optimizerecordf.append(intermediate_result.fun)
def callbackf(coeffs):
    optimizerecordf.append(targetfun(coeffs))
# result = scipy.optimize.minimize(targetfun, coeffs, method="SLSQP", options={"maxiter":10,"disp":True},callback=callbackf)
result = scipy.optimize.minimize(targetfun, coeffs, method="SLSQP", options={"maxiter":10,"disp":True},callback=callbackf)
plt.plot(optimizerecordf)
plt.show()
plot(result.x)
plt.plot(result.x)
plt.show()