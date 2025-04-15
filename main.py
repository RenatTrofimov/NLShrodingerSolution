import numpy as np
from numpy import sin, cos, pi, sqrt, exp
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.special import gamma
from joblib import Parallel, delayed
import multiprocessing
import time

class NonlinearPart():
	def __init__(self):
		self.times = 0
		self.mm = np.float64(7)
		self.kk = np.float64(9)
		self.W0 = np.float64(4.3e-12)
		self.j = np.float64(1/10.6e-15)
		self.r = 5
		self._r = np.arange(self.r)
		self._sumFF = np.empty(int(self.mm*self.kk), dtype=np.complex128)
		self._A1 = np.empty(int(self.mm*self.kk))
		self._cachedGammaUp = np.empty((int(self.kk), self.r), dtype=np.complex128)
		self._cachedGammaDown = np.empty((self.r), dtype=np.complex128)
		self._cachedGammaUp[:,:] = np.nan
		self._cachedGammaDown[:] = np.nan
		self._sumFF[:] = np.nan
		self._A1[:] = np.nan
		self.GG(1)

	def getNNP(self, inArray, outArray):
		#outArray = np.array([self.GG(element) for element in inArray])
		outArray = np.array(Parallel(n_jobs=8)(delayed(self.chuck)(element) for element in np.array_split(inArray, 8))).ravel()
		#outArray = np.array([self.chuck(element) for element in np.array_split(inArray, 16)]).ravel()

	def chuck(self, inArray):
		start = time.time()
		a = np.array([self.GG(element) for element in inArray])
		end = time.time()
		#print(end-start)
		return a

		#return np.array([self.GG(element) for element in inArray])
	def get(self, U):
		start = time.time()
		a = self.GG(U)
		end = time.time()
		print(end-start)
		return a
		#return self.GG(U)

	def gammaPart(self, q, U):
		if np.isnan(np.sum(self._cachedGammaUp[int(q)])):
			self._cachedGammaUp[int(q)] = ((-1)**self._r) * (q**(2*self._r + 1))
		if np.isnan(np.sum(self._cachedGammaDown)):
			self._cachedGammaDown = 2**(2*self._r + 1) * gamma(self._r + 1) * gamma(self._r + 2)
		return np.sum(self._cachedGammaUp[int(q)] * np.absolute(U)**(2*self._r) / self._cachedGammaDown)
	def b(self, i):
			return cos(pi*i/self.mm)
	def ee(self,  x, i):
			return self.W0 * sqrt(1.0 + 4.0*cos(x)*self.b(i)+4.0*self.b(i)*self.b(i))
	def A1(self,  i,k):
			if np.isnan(self._A1[int(i*self.kk + k)]):
				self._A1[int(i*self.kk + k)] = integrate.quad(lambda x: self.ee(x, i)*cos(k*x), -pi, pi)[0]/pi
			return self._A1[int(i*self.kk + k)]
	def tetta(self, x,i,k):
			return exp(-(self.j*self.A1(i,0)/2+np.sum([self.j*self.A1(i,_k)*cos(_k*x) for _k in np.arange(1, int(self.kk))])))/(1+ exp(-(self.j*self.A1(i,0)/2+np.sum([self.j*self.A1(i,_k)*cos(_k*x) for _k in np.arange(1, int(self.kk))]))))
	def FF(self, i,k):
			def func(i,k):
					return integrate.quad(lambda x :self.tetta(x,i,k)*cos(k*x),-pi,pi)[0]
			def func2(i,k):
					return integrate.quad(lambda x :self.tetta(x,i,k),-pi,pi)[0]
			return -k*self.A1(i,k)*func(i,k)/np.sum([func2(_i,k) for _i in np.arange(1, self.mm)])/self.W0
	def GG(self, U):
		_sum = 0
		if not self.times:
			for _k in np.arange(self.kk):
				if np.isnan(self._sumFF[int(_k)]):
					self._sumFF[int(_k)] = np.sum([self.FF(_i,_k) for _i in np.arange(1, self.mm)])
				_sum += self.gammaPart(_k, U)*self._sumFF[int(_k)]
				
		#print(_sum)
		else:
			for _k in np.arange(self.kk):
				_sum += (np.sum(self._cachedGammaUp[int(_k)] * np.absolute(U)**(2*self._r) / self._cachedGammaDown))*self._sumFF[int(_k)]
				#print(_k, "-", self._sumFF[int(_k)])
		self.times+=1
		return _sum


def Besse():
	M = 1000
	xEnd = 10
	tEnd = 5           
	dt = 0.001         
	dx = 2*xEnd / M
	jj = 1j           
	r = dt / (2 * dx * dx)
	lamda = 1
	nn = 1
	Nt = tEnd / dt

	x = np.zeros(M)
	Resh = np.zeros(M, dtype=np.complex128)
	U0 = np.zeros(M, dtype=np.complex128)
	U1 = np.zeros(M, dtype=np.complex128)
	V0 = np.zeros(M)
	V1 = np.zeros(M)
	alfa_plus = np.zeros(M, dtype=np.complex128)
	alfa_minus = np.zeros(M, dtype=np.complex128)

	# Initialize A_plus and A_minus as sparse matrices (more efficient)
	A_plus = np.zeros((M, M), dtype=np.complex128)
	A_minus = np.zeros((M, M), dtype=np.complex128)

	for i in range(M):
		x[i] = -xEnd + dx * (i + 1)
		U0[i] = np.exp(-(x[i]**2)) 

	# Set up the A_plus and A_minus matrices
	for i in range(1, M):
		A_plus[i, i-1] = 1
		A_minus[i, i-1] = -1

	for i in range(M-1):
		A_plus[i, i+1] = 1
		A_minus[i, i+1] = -1

	A_plus[0, M-1] = 1
	A_minus[0, M-1] = -1
	A_plus[M-1, 0] = 1
	A_minus[M-1, 0] = -1
	plt.plot(x, np.real(U0))
	nlp = NonlinearPart()
	while nn < Nt:
		for i in range(M):
			V0[i] = np.abs(U0[i])**2
			V1[i] = -V0[i] + 2 * np.abs(U0[i])**2
		nlp.getNNP(np.copy(V1), V1)
		for i in range(M):
			alfa_plus[i] = jj/r - 2 - dx**2 * V1[i]
			alfa_minus[i] = -jj/r - 2 - dx**2 * V1[i]
			A_plus[i, i] = alfa_plus[i]
			A_minus[i, i] = -alfa_minus[i]

		B = np.dot(A_minus, U0)
		U1 = np.linalg.solve(A_plus, B)  # More efficient than explicit inverse

		U0 = U1.copy()
		nn += 1
		print(nn)

	for i in range(M):
		Resh[i] = np.exp(jj * (x[i] - (lamda + 1) * dt * Nt))

	plt.plot(x, np.real(U1))
	plt.show()

	# Call the function
Besse()