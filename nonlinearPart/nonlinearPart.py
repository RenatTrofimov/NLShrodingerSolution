from scipy import integrate 
from numpy import sin, cos, exp, pi, sqrt
from scipy.special import gamma, factorial
from numba import njit, prange
import numpy as np

def get(U):
    return GG(U)

def GG(U):
    mm = np.float64(7)
    kk = np.float64(9)
    W0 = np.float64(4.3e-12)
    j = np.float64(1/10.6e-15)
    _FF = np.zeros(int(mm*kk))
    _A1 = np.zeros(int(mm*kk))
    _Gp = np.zeros(int(mm*kk))
    def gammaPart(q, U):
        r = np.arange(5)
        up = ((-1)**r)*(q**(2*r+1))
        down = 2**(2*r+1)*gamma(r+1)*gamma(r+2)
        return np.sum(up*np.absolute(U)**(2*r)/down)

    def b(i):
        return cos(pi*i/mm)
    def ee(x, i): 
        return W0 * sqrt(1.0 + 4.0*cos(x)*b(i)+4.0*b(i)*b(i))
    def A1(i,k):
        if _A1[int(i*kk + k)] == 0.0:
            _A1[int(i*kk + k)] = integrate.quad(lambda x: ee(x, i)*cos(k*x), -pi, pi)[0]/pi
        return _A1[int(i*kk + k)] 
    def tetta(x,i,k):    
        return exp(-(j*A1(i,0)/2+np.sum([j*A1(i,_k)*cos(_k*x) for _k in np.arange(1, kk)])))/(1+ exp(-(j*A1(i,0)/2+np.sum([j*A1(i,_k)*cos(_k*x) for _k in np.arange(1, kk)]))))

    def FF(i,k):
        def func(i,k):
            return integrate.quad(lambda x :tetta(x,i,k)*cos(k*x),-pi,pi)[0]
        def func2(i,k):    
            return integrate.quad(lambda x :tetta(x,i,k),-pi,pi)[0]
        if _FF[int(i*kk + k)] == 0.0:
            _FF[int(i*kk + k)] = -k*A1(i,k)*func(i,k)/np.sum([func2(_i,k) for _i in np.arange(1, mm)])/W0

        return _FF[int(i*kk + k)]
    _sum = 0
    for _k in np.arange(kk):
        _sum += gammaPart(_k, U)*np.sum([FF(_i,_k) for _i in np.arange(1, mm)])
        #print(_k, "-", _sum)
    return _sum

