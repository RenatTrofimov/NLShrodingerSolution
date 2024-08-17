from scipy import integrate 
from numpy import sin, cos, exp, pi
import numpy as np
def get():
    1.

def alpha_sq(s,q):
    return 1

def b_sq(s,q, max_q):
    tempeture = 1.0
    gamma_0 = 2.7
    q_s = np.arange(max_q, dtype = np.float)
    
    def up(r): 
        return cos(q*r)*down(r)
    def down(r):
        return exp(-np.sum(cos(r*q_s))/tempeture)
    
    return -q * (alpha_sq(s,q)/gamma_0) * (integrate.quad(up, -pi, pi)[0]) / (integrate.quad(down, -pi, pi))[0]
