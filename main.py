from TM.thomasMethod import thomas
import numpy as np
# Main Program
#step 1: Create Coefficient and LHS matrices store (a,b,c,d) seperately
z = 0.1
y = 0.1
C = 1.0
c0 = 0.1
w0 = 0.1 
A = 1.0
k = 1.0
B = 2.*1j*k*c0/w0
D = A*z/B*y**2

a = D
b = (1. - D)
c = D
d = C*z/B

# Call thomas function
x = thomas(a,b,c,d)
# Print x
print('The values of x are: ',x)