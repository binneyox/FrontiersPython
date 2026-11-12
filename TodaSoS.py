import numpy as np
from  scipy.integrate import solve_ivp as solve
import matplotlib.pyplot as plt
twort=2*np.sqrt(3)

def E2(x,v) :
    return .5*(v[0]*v[0]+v[1]*v[1]) + TodaPhi(x)

def E1(w) :
    x=(w[0],w[1])
    v=(w[2],w[3])
    return E2(x,v)

def TodaPhi(x) :
    return (np.exp(2*x[1])*np.cosh(twort*x[0]) + .5 *np.exp(-4*x[1]))/12.-0.125

def TodaForce(x) :
    return (-np.exp(2*x[1])*np.sinh(twort*x[0])*twort/12., -(np.exp(2*x[1])*np.cosh(twort*x[0]) - np.exp(-4*x[1]))/6.)
    
def Todaderivs(t,w) :
    x=(w[0],w[1])
    f=TodaForce(x)
    return (w[2],w[3],f[0],f[1])

def xcross(t,y) :
    return y[0]

def add_orbit(E,x) :
    v0=np.sqrt(2*(E-TodaPhi(x)))
    theta=-.3
    v=(np.cos(theta)*v0, np.sin(theta)*v0)
    w=(x[0],x[1],v[0],v[1])
    print (E1(w))
    xcross.direction=1
    ts=(0,1000)
    acc=(1e-10,1e-10,1e-10,1e-10)
    soln=solve(Todaderivs,ts,w,method="RK45",events=xcross,rtol=1e-5,atol=acc,dense_output=True,t_eval=None)
    n=soln.t.size-1
    w=(soln.y[0,n],soln.y[1,n],soln.y[2,n],soln.y[3,n])
    print (E1(w))
    print (soln.t.size,soln.y.shape)
    ys=[]
    pys=[]
    nCross=len(soln.t_events[0])
    print(nCross," crosses")
    for i in range (0,nCross-1,1) :
            ys.append(soln.sol(soln.t_events[0][i])[1])
            pys.append(soln.sol(soln.t_events[0][i])[3])
    plt.scatter(ys,pys,s=1)

def main() :
    E=.09
    x=(0,.2)
    plt.xlim(-.5,.5)
    plt.ylim(-0.5,.5)
    plt.xlabel('y')
    plt.ylabel('py')
    norb=10
    for i in range (0, norb, 1) :
        x=(0,-.3+.5*i/(norb-1))
        add_orbit(E,x)
    plt.show()

main()
