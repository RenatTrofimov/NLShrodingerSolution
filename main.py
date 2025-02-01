from TM.thomasMethod import thomas
from nonlinearPart.nonlinearPart import get
from matplotlib import pyplot as plt
from scipy.linalg import solve_banded
from joblib import Parallel, delayed
import numpy as np
import multiprocessing

def execute(values):
	global nnl, Xc
	for i in range(values[0],values[1]):
		nnl[i] = get(Xc[i])
		print(i)


def process(i):
    return get(Xc[i])

# Main Program
#step 1: Create Coefficient and LHS matrices store (a,b,c,d) seperately
hx = 0.001
ht = 0.001
x = np.arange(0,5,hx, dtype=np.complex128)
Xc = np.exp(-1*(x-2)**2)*np.exp(1j*(x-2))
X0 = np.exp(-1*(x-2)**2)*np.exp(1j*(x-2))
Xn = np.insert(Xc[:-1], 0, 0)
Xp = np.append(Xc[1:], 0)
a = np.zeros(x.shape[0], dtype=np.complex128)
c = np.zeros(x.shape[0], dtype=np.complex128)
a[2:x.shape[0]-1] = a[2:x.shape[0]-1] + 1/(2*hx**2)
c[:x.shape[0]-2] = c[:x.shape[0]-2] + 1/(2*hx**2)
ab = np.ones(x.shape[0], dtype=np.complex128)
ab[1:-1] = ab[1:-1]*(10*2j/hx - hx**(-2))
ab = np.array([a, ab, c], dtype=np.complex128)
nnl = np.zeros(x.shape[0])
pool_obj = multiprocessing.Pool()
pool_indexes = np.zeros((16,2), dtype=np.int32)
for index in range(16):
	pool_indexes[index][0] = int(x.shape[0]/16)*index
	pool_indexes[index][1] = int(x.shape[0]/16)*(index + 1) - 1
pool_indexes[0][0] = 0
pool_indexes[15][1] = x.shape[0] - 1
print(pool_indexes)
for j in range(10):
	for i in range(100):
	
		#nnl = np.array([get(i) for i in Xc])
		pool_obj.map(execute,pool_indexes)
		#nnl = Parallel(n_jobs=24)(delayed(process)(i) for i in range(pool_indexes))
		b = Xc*(10*2j/ht + 1/(hx**2))-(Xn+Xp)/(2*hx**2) #+ Xc*nnl
		b[0] = b[-1] = 0
	
		Xc = solve_banded((1, 1), ab, b)
		Xn = np.insert(Xc[:-1], 0, 0)
		Xp = np.append(Xc[1:], 0)

	
	plt.plot(x, np.abs(Xc**2))
	print(j)
pool_obj.close()
plt.plot(x, np.abs(X0**2))

print(np.abs(np.sum(Xc))**2 - np.abs(np.sum(X0))**2)
plt.show()