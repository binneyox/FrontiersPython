import numpy as np
from  scipy.integrate import solve_ivp as solve
from scipy.special import erf
import matplotlib.pyplot as plt

rtpi=np.sqrt(3.14159)
rho0=1

def rhoW(phi) :
    if phi<=0 : return 0
    xm=np.sqrt(phi)
    return rhoK(phi)-2/15.*phi*phi*xm/rho0

def rhoK(phi) :
    if phi<=0 : return 0
    xm=np.sqrt(phi)
    return (np.exp(phi)*erf(xm)-2*xm/rtpi*(1+2*phi/3))/rho0
#    return (np.exp(phi)*.5*(-xm*np.exp(-xm*xm) + erf(xm)*.5*rtpi)-phi*xm/3.)/rho0

def derivs(r,w) :
    dMdr=r*r*rhoK(w[1])
    dPhidr=-9*w[0]/(r*r)
    return (dMdr,dPhidr)

def phizero(r,w) :
    return w[1]

for phi0 in np.arange(3,9.1,3) :
    print(phi0)
    rho0=1
    rho0=rhoW(phi0)
    rmin=1e-2
    rmax=200
    w=[rmin*rmin*rmin/3,phi0]
    rs=(rmin,rmax)
    acc=(1e-4,1e-4)
    soln=solve(derivs,rs,w,method="RK45",events=phizero,rtol=1e-4,atol=acc,dense_output=True,t_eval=None)
    n=soln.t.size-1
    rs=[]
    rhos=[]
    x=rmin
    for x in np.arange(rmin,soln.t_events[0],.1) :
        rs.append(np.log10(x))
        rhos.append(np.log10(rhoW(soln.sol(x)[1])))
    plt.plot(rs,rhos)
    plt.ylim(-9,.1)
    print(soln.t_events[0])

plt.xlabel("log r")
plt.ylabel("log(rho/rho(0))")
plt.show()
