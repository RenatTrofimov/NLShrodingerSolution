from TM.thomasMethod import thomas
from nonlinearPart.nonlinearPart import get
import numpy as np
# Main Program
#step 1: Create Coefficient and LHS matrices store (a,b,c,d) seperately
len = 5

hx = 0.1
ht = 0.1
x = np.arange(0,6,hx)
Xc = np.exp(-(x-3)**2)
Xn = np.insert(Xc[:-1], 0, 0)
Xp = np.append(Xc[1:], 0)


z = np.complex128(0.1)
y = np.complex128(0.1)
C = Xc*0.1
c0 = np.complex128(0.1)
w0 = np.complex128(0.1) 
A = np.complex128(1.0)
k = np.complex128(1.0)
B = 2.*1j*k*c0/w0
D = Xc*0.1

a = (0.5*A/hx**2)*np.ones(Xc.shape[0], dtype = np.complex128)
b = (B/ht - A/hx**2)*np.ones(Xc.shape[0], dtype = np.complex128)
c = (0.5*A/hx**2)*np.ones(Xc.shape[0], dtype = np.complex128)
d = Xc*(B/ht + A/hx**2) + C - 0.5*A*(Xn - Xp)

# Call thomas function
X = thomas(a,b,c,d, np.dtype(np.complex128))
# Print x
print('The values of x are: ',x)
from matplotlib import pyplot as plt
plt.plot(np.abs(X))
plt.show()