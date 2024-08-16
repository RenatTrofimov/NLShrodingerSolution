from TM.thomasMethod import thomas
import numpy as np
# Main Program
#step 1: Create Coefficient and LHS matrices store (a,b,c,d) seperately
len = 5
z = np.complex128(0.1)
y = np.complex128(0.1)
C = np.ones(len, dtype = np.complex128)
c0 = np.complex128(0.1)
w0 = np.complex128(0.1) 
A = np.complex128(1.0)
k = np.complex128(1.0)
B = 2.*1j*k*c0/w0
D = A*z/B*y**2

a = D*np.ones(len-1, dtype = np.complex128)
b = (1. - 2.*D)*np.ones(len, dtype = np.complex128)
c = D*np.ones(len-1, dtype = np.complex128)
d = C*z/B

# Call thomas function
x = thomas(a,b,c,d, np.dtype(np.complex128))
# Print x
print('The values of x are: ',x)