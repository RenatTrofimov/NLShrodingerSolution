from TM.thomasMethod import thomas
#from nonlinearPart.nonlinearPart import get
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.linalg import solve_banded
# Main Program
#step 1: Create Coefficient and LHS matrices store (a,b,c,d) seperately
len = 250

hx = 10/len
ht = 0.0001
lamda = 1.0

x = np.arange(0,10,hx)
k=1*np.pi
Xc = np.exp(-(x-3)**2)*np.exp(1j*k*(x-3))
Xn = np.insert(Xc[:-1], 0, 0)
Xp = np.append(Xc[1:], 0)

#plt.plot(Xc)
#plt.plot(Xn)
#plt.plot(Xp)
#plt.show()
Vp = np.zeros(len)
Vn = np.absolute(Xc)*np.absolute(Xc)
k = 0.5/hx**2
a = k*np.ones(len)
b = (1j/ht - 1/hx**2 - lamda*Vn/2)
c = k*np.ones(len)
d = (1j/ht + 1/hx**2 + lamda*Vn/2)*Xc - Xn*k - Xp*k

A = np.zeros(len)
C = np.zeros(len)
A[1:] = a[1:]
C[:-1] = a[:-1]
ab = np.array([A,b,C])

fig, ax = plt.subplots()

line, = ax.plot(x, Xc)
def update_cos(frame, line, x):
    global a,b,c,d, Xc, Xn, Xp,Vn,Vp, hx, ht
    global A,ab,C
    Xc = thomas(a,b,c,d, np.dtype(np.complex128))
    Xc1 = solve_banded((1, 1), ab, d)
    Vp = np.copy(Vn)
    
    Xn = np.insert(Xc[:-1], 0, 0)
    Xp = np.append(Xc[1:], 0)
    
    Vn = -Vp + 2*np.absolute(Xc)*np.absolute(Xc)
    k = 0.5/hx**2
    a = k*np.ones(len)
    b = (1j/ht - 1/hx**2 - lamda*Vn/2)
    c = k*np.ones(len)
    d = (1j/ht + 1/hx**2 + lamda*Vn/2)*Xc - Xn*k - Xp*k


    A = np.zeros(len)
    C = np.zeros(len)
    A[1:] = a[1:]
    C[:-1] = a[:-1]
    ab = np.array([A,b,C])
    print(sum(np.absolute(Xc)*np.absolute(Xc)) - sum(np.absolute(Xc1)*np.absolute(Xc1)))
    line.set_ydata( np.absolute(Xc1)**2 )
    return [line]
phasa = np.arange(0, 100, 0.01)
animation = FuncAnimation(
    fig,                # фигура, где отображается анимация
    func=update_cos,    # функция обновления текущего кадра
    frames=phasa,       # параметр, меняющийся от кадра к кадру
    fargs=(line, x),    # дополнительные параметры для функции update_cos
    interval=1,       # задержка между кадрами в мс
    blit=True,          # использовать ли двойную буферизацию
    repeat=False)   

plt.show()

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