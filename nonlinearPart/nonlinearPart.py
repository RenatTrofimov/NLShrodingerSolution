from scipy import integrate 
from numpy import sin, cos, exp, pi
import numpy as np
def get():
    1.

def alpha_sq(s,q):
    return s*q

def b_sq(s,q):
    tempeture = 1.0
    gamma_0 = 1.0
    q_s = np.arange(1,10, dtype=float64)
    return -q * (alpha_sq(s,q)/gamma_0) * (integrate.quad(lambda r: cos(q*r)*exp(-np.sum(cos(q_s*r)*alpha_sq(s,q_s))/tempeture), -pi, pi)) / (integrate.quad(lambda r: exp(-np.sum(cos(q_s*r)*alpha_sq(s,q_s))/tempeture), -pi, pi))
    

result = integrate.quad(lambda x: sin(x), 0, 4.5)
print (result)