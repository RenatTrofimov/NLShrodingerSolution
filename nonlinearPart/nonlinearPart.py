from scipy import integrate 
from numpy import sin, cos, exp, pi, sqrt
from scipy.special import gamma, factorial
from numba import njit, prange
import numpy as np
def getNonLinearPart(U):
    return effectiveEq(U, 10, 10)

def b_sq(s,q, max_q, max_s):
    tempeture = 1.0
    gamma_0 = 2.7*1.6e-9
    q_s = np.arange(max_q)
    def alpha_sq(s,q):
        dx = 3*1.42/2
        plank = 1.6e-25
        m = max_s
        a = dx/plank
        def eps(p_x):
            return gamma_0*sqrt(1+4*cos(a*p_x)*cos(pi*s/m)+4*cos(pi*s/m)**2)*cos(p_x*q*dx/plank)
        
        return dx/(pi*plank) * integrate.quad(eps, -pi*plank/dx, pi*plank/dx)[0]
    def up(r): 
        return cos(q*r)*down(r)
    def down(r):
        return exp(-np.sum(cos(r*q_s))/tempeture)
    
    return -q * (alpha_sq(s,q)/gamma_0) * (integrate.quad(up, -pi, pi)[0]) / (integrate.quad(down, -pi, pi))[0]

def sumPart(q, U):
    r = np.arange(10)
    up = ((-1)**r)*(q**(2*r+1))*np.absolute(U)**(2*r)
    down = 2**(2*r+1)*gamma(r+1)*gamma(r+2)
    return np.sum(up/down)

def effectiveEq(U, q_max, s_max):
    q = np.arange(q_max)
    s = np.arange(s_max)
    answ = 0
    for _q in q:
        for _s in s:
            answ+=sumPart(_q,U)*b_sq(_s, _q, q_max, s_max)
    return answ

mm = 7
kk = 10
W0 = 4.3e-12
j = 1/10.6e-15
def b(i):
    return cos(pi*i/mm)
def ee(x, i): 
    return W0 * (1 + 4*cos(x)*b(i)+4*b(i)*b(i))**0.5
def A1(i,k):
    return integrate.quad(lambda x: ee(x, i)*cos(k*x), -pi, pi)[0]/pi

def tetta(x,i,k):    
    def up():
        return exp(-(j*A1(i,0)/2+np.sum([j*A1(i,_k)*cos(_k*x) for _k in np.arange(1, kk)])))
    return up()/(1+up())
def FF(i,k):
    def func(i,k):
        return integrate.quad(lambda x :tetta(x,i,k)*cos(k*x),-pi,pi)[0]
    return -(k*A1(i,k)/W0)*func(i,k)/np.sum([func(_i,k) for _i in np.arange(1, mm)])

def GG():
    Sum = 0
    for _i in prange(mm):
        for _k in prange(kk):
            Sum += FF(_i,_k)
    print(Sum)