from TM.thomasMethod import thomas
from nonlinearPart.nonlinearPart import get
from matplotlib import pyplot as plt
from scipy.linalg import solve_banded
import numpy as np
# Main Program
#step 1: Create Coefficient and LHS matrices store (a,b,c,d) seperately
hx = 0.1
ht = 0.01
x = np.arange(0,5,hx)
Xc = np.exp(-(x-2)**2)*np.exp(1j*(x-2))
Xn = np.insert(Xc[:-1], 0, 0)
Xp = np.append(Xc[1:], 0)

a = np.zeros(x.shape[0])
c = np.zeros(x.shape[0])
a[1:x.shape[0]-1] = a[1:x.shape[0]-1] + 1/(2*hx**2)
c[:x.shape[0]-2] = c[:x.shape[0]-2] + 1/(2*hx**2)
nnl = np.array([get(i) for i in Xc])
b = Xc*(2j/ht + hx**(-2))-(Xn+Xp)/(2*hx**2) + Xc*nnl

ab = np.array([a, np.ones(x.shape[0])*(2j/hx - hx**(-2)), c])

Xc = solve_banded((1, 1), ab, b)

plt.plot(x, np.abs(Xc)**2)
plt.show()