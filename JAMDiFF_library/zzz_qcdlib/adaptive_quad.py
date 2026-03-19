import numpy as np
from numba import jit

xg4,wg4=np.polynomial.legendre.leggauss(4)
xg8,wg8=np.polynomial.legendre.leggauss(8)

@jit(nopython=True, cache=False)
def get_gauss_knots(a,b,xg):
    x0=0.5*(b-a)*xg+0.5*(a+b)
    jac=0.5*(b-a)
    return x0,jac

@jit(nopython=True, cache=False)
def quad(integrand,args,a,b,ngrid=3,eps_rel=1e-3,eps_abs=1e-3):
    x0,jac=get_gauss_knots(a,b,xg4)
    integrand_vals=np.array([integrand(_,*args) for _ in x0])
    int1=np.sum(wg4*jac*integrand_vals)
    
    grid=np.linspace(a,b,ngrid)
    int2=0
    for i in range(ngrid-1):
        x0,jac=get_gauss_knots(grid[i],grid[i+1],xg4)
        integrand_vals=np.array([integrand(_,*args) for _ in x0])
        int2+=np.sum(wg4*jac*integrand_vals)

    rel_err=np.abs(int2-int1)/int2
    abs_err=np.abs(int2-int1)
    
    if rel_err<eps_rel or abs_err<eps_abs:
        #print(a,b,int2)
        return int2
    else:
        int3=0
        for i in range(ngrid-1):
            int3+=quad(integrand,args,grid[i],grid[i+1],ngrid,eps_rel,eps_abs)
        return int3