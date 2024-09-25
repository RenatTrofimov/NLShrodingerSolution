from TM.thomasMethod import thomas
from nonlinearPart.nonlinearPart import get
from matplotlib import pyplot as plt
import numpy as np
# Main Program
#step 1: Create Coefficient and LHS matrices store (a,b,c,d) seperately
len = 5

hx = 0.06
ht = 0.01
x = np.arange(0,6,hx)
Xc = np.exp(-(x-3)**2)
Xn = np.insert(Xc[:-1], 0, 0)
Xp = np.append(Xc[1:], 0)
plt.plot(x, np.abs(Xc))
plt.show()
C = np.fromiter((get(_x) for _x in Xc), Xc.dtype, count=Xc.shape[0])
w0 = np.complex64(1e14) 
w = np.complex64(5e14) 
eps = np.complex64(4.0) 
A = np.complex64(1.0)

B = 2.*1j*w*np.sqrt(eps)/w0


a0 = np.ones(Xc.shape[0], dtype = np.complex64)
b0 = np.ones(Xc.shape[0], dtype = np.complex64)
c0 = np.ones(Xc.shape[0], dtype = np.complex64)

a = (0.5*A/hx**2)*a0
b = (B/ht - A/hx**2)*a0
c = (0.5*A/hx**2)*c0
d = Xc*(B/ht + A/hx**2) - C*Xc - 0.5*A*(Xn + Xp)

X = thomas(a,b,c,d, np.dtype(np.complex64))
print("run")
for i in range(10):
    
    
    Xn = np.insert(X[:-1], 0, 0)
    Xp = np.append(X[1:], 0)
    
    d = X*(B/ht + A/hx**2) - np.fromiter((get(_x) for _x in X), X.dtype, count=X.shape[0])*Xc - 0.5*A*(Xn + Xp)
    
    X = thomas(a,b,c,d, np.dtype(np.complex64))

plt.plot(x, np.abs(X))
plt.plot(x, np.abs(Xc))
plt.show()