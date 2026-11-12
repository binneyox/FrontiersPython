import numpy as np
from  scipy.integrate import solve_ivp as solve
import matplotlib.pyplot as plt

q=.9
Rc2=.1
Lz2=30

def E2(x,v) :
    return .5*(v[0]*v[0]+v[1]*v[1])+logPhi(x)

def E1(w) :
    x=(w[0],w[1])
    v=(w[2],w[3])
    return E2(x,v)

def logPhi(x) :
    R2=x[0]*x[0]+x[1]*x[1]/q+Rc2
    return np.log(R2)+.5*Lz2/(x[0]*x[0])

def HHPhi(x) :
    return .5*(x[0]*x[0] + y[0]*y[0] + 2*x[0]*x[0]*y[1] - 2*y[0]*y[0]*y[0]/3)

def Force(x) :
    R2=x[0]*x[0]+x[1]*x[1]/q+Rc2
    return (-2*x[0]/R2+Lz2/(x[0]*x[0]*x[0]),-2*x[1]/(q*R2))

def derivs(t,w) :
    x=(w[0],w[1])
    v=(w[2],w[3])
    f=Force(x)
    return (w[2],w[3],f[0],f[1])

def zcross(t,y) :
    return y[1]

def main() :
#    zcross.direction=1
    x=(5,.3)
    f=Force((x[0],0))
    print("Lz2 ",x[0]*x[0]*(-f[0]*x[0]))
    v=(.2,.1)
    w=(x[0],x[1],v[0],v[1])
    print (E1(w))
    ts=(0,2000)
    acc=(1e-5,1e-5,1e-5,1e-5)
    soln=solve(derivs,ts,w,method="RK45",events=zcross,rtol=1e-5,atol=acc,dense_output=True,t_eval=None)
    n=soln.t.size-1
    w=(soln.y[0,n],soln.y[1,n],soln.y[2,n],soln.y[3,n])
    print (E1(w))
    print (soln.t.size,soln.y.shape)
    Rs=[]
    pRs=[]
    nCross=len(soln.t_events[0])
    print(nCross," crosses")
    for i in range (0,nCross-1,1) :
        if( soln.sol(soln.t_events[0][i])[3]>0 ) :
            Rs.append(soln.sol(soln.t_events[0][i])[0])
            pRs.append(soln.sol(soln.t_events[0][i])[2])
    plt.xlim(2.5,6)
    plt.ylim(-.7,.7)
    plt.xlabel('R')
    plt.ylabel('z')
    plt.plot(soln.y[0],soln.y[1])
    plt.show()

    plt.scatter(Rs,pRs,s=10)
    plt.xlabel('R')
    plt.ylabel('pR')
    plt.show()

main()

