from TM.thomasMethod import thomas
#from nonlinearPart.nonlinearPart import get
import numpy as np
from matplotlib import pyplot as plt

# Main Program
#step 1: Create Coefficient and LHS matrices store (a,b,c,d) seperately
len = 100

hx = 6/len
ht = 0.1
lamda = 1.0

x = np.arange(0,6,hx)
Xc = np.exp(-(x-3)**2)
Xn = np.insert(Xc[:-1], 0, 0)
Xp = np.append(Xc[1:], 0)

Vp = np.zeros(len)
Vn = np.absolute(Xc)

a = 1/hx**2*np.ones(len)
b = (1j/ht - 1/hx**2 - lamda*Vn/2)
c = 1/hx**2*np.ones(len)
d = -(1j/ht - 1/hx**2 - lamda*Vn/2)*Xc - Xn/hx**2 - Xp/hx**2



for i in range(500):
    Xc = thomas(a,b,c,d, np.dtype(np.complex128))
    
    Vp = np.copy(Vn)
    
    Xn = np.insert(Xc[:-1], 0, 0)
    Xp = np.append(Xc[1:], 0)
    
    Vn = Vp - 2*np.absolute(Xc)
    
    a = 1/hx**2*np.ones(len)
    b = (1j/ht - 1/hx**2 - lamda*Vn/2)
    c = 1/hx**2*np.ones(len)
    d = -(1j/ht - 1/hx**2 - lamda*Vn/2)*Xc - Xn/hx**2 - Xp/hx**2
    





plt.plot(Xc)
plt.show()

""""
a = D*np.ones(len-1, dtype = np.complex128)
b = (1. - 2.*D)*np.ones(len, dtype = np.complex128)
c = D*np.ones(len-1, dtype = np.complex128)
d = C*z/B

# Call thomas function
x = thomas(a,b,c,d, np.dtype(np.complex128))
# Print x
print('The values of x are: ',x)

plt.plot(x, Xp)
plt.show()
"""