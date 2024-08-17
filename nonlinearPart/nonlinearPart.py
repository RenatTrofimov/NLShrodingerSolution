from scipy import integrate 
from numpy import sin, cos, exp, pi
from scipy.special import gamma, factorial
import numpy as np
def get(U):
    effectiveEq(U, 10, 10)

def alpha_sq(s,q):
    return .01

def b_sq(s,q, max_q):
    tempeture = 1.0
    gamma_0 = 2.7
    q_s = np.arange(max_q, dtype = np.float)
    
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
            answ+=sumPart(_q,U)*b_sq(_s, _q, q_max)
    return answ
    
