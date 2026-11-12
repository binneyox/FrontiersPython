import numpy as np
from  scipy.integrate import solve_ivp as solve
from scipy.special import erf
import matplotlib.pyplot as plt

rho0=1
PI=3.14159
beta=1

def rhoI(phi) :
    return rho0*np.exp(-beta*phi)

def derivs(r,w) :
    dMdr=4*PI*r*r*rhoI(w[1])
    dPhidr=w[0]/(r*r)
    return (dMdr,dPhidr)

rstart=1e-3
rmax=100
rs=(rstart,rmax)
Mstart=4./3.*PI*rstart*rstart*rstart*rho0
Phistart=.5*Mstart/rstart
w=(Mstart,Phistart)
acc=(1e-4,1e-4)
soln=solve(derivs,rs,w,method="RK45",rtol=1e-4,atol=acc,dense_output=True,t_eval=None)
n=soln.t.size-1
rs=[]
rhos=[]
x=rstart
for x in np.arange(10*rstart,rmax,.1) :
     rs.append(np.log10(x))
     rhos.append(np.log10(rhoI(soln.sol(x)[1])))
plt.plot(rs,rhos)
plt.ylim(-5,.1)
plt.show()
 