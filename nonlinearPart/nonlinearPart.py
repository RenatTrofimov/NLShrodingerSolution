from scipy import integrate 
from numpy import sin, cos, exp, pi, sqrt
from scipy.special import gamma, factorial
from numba import njit, prange
from joblib import Parallel, delayed
import numpy as np
import time
import multiprocessing


class NonlinearPart():
	mm = np.float64(7)
	kk = np.float64(9)
	W0 = np.float64(4.3e-12)
	j = np.float64(1/10.6e-15)
	_FF = np.zeros(int(mm*kk))
	_A1 = np.zeros(int(mm*kk))
	_Gp = np.zeros(int(mm*kk))
	_sumFF = np.zeros(int(mm*kk))
	_cachedGammaUp = np.empty((int(kk), 5))
	_cachedGammaDown = np.empty((5))
	r = np.arange(5)
	def getNNP(self, inArray, outArray):
		
		outArray = Parallel(n_jobs=512)(delayed(self.GG)(element) for element in inArray)
		
	def get(self, U):
		start = time.time()
		a = self.GG(U)
		end = time.time()
		print(end-start)
		return a
		#return self.GG(U)

	def gammaPart(self, q, U):
        	if np.sum(self._cachedGammaUp[int(q)]):
        		self._cachedGammaUp[int(q)] = ((-1)**self.r)*(q**(2*self.r+1))
	        if np.sum(self._cachedGammaDown):
        		self._cachedGammaDown = 2**(2*self.r+1)*gamma(self.r+1)*gamma(self.r+2)
        	return np.sum(self._cachedGammaUp[int(q)]*np.absolute(U)**(2*self.r)/self._cachedGammaDown)

	def b(self, i):
        	return cos(pi*i/self.mm)
	def ee(self,  x, i): 
        	return self.W0 * sqrt(1.0 + 4.0*cos(x)*self.b(i)+4.0*self.b(i)*self.b(i))
	def A1(self,  i,k):
        	if self._A1[int(i*self.kk + k)] == 0.0:
        	    self._A1[int(i*self.kk + k)] = integrate.quad(lambda x: self.ee(x, i)*cos(k*x), -pi, pi)[0]/pi
        	return self._A1[int(i*self.kk + k)] 
	def tetta(self, x,i,k):    
        	return exp(-(self.j*self.A1(i,0)/2+np.sum([self.j*self.A1(i,_k)*cos(_k*x) for _k in np.arange(1, self.kk)])))/(1+ exp(-(self.j*self.A1(i,0)/2+np.sum([self.j*self.A1(i,_k)*cos(_k*x) for _k in np.arange(1, self.kk)]))))

	def FF(self, i,k):
        	def func(i,k):
            		return integrate.quad(lambda x :self.tetta(x,i,k)*cos(k*x),-pi,pi)[0]
        	def func2(i,k):    
            		return integrate.quad(lambda x :self.tetta(x,i,k),-pi,pi)[0]
        	if self._FF[int(i*self.kk + k)] == 0.0:
           		self._FF[int(i*self.kk + k)] = -k*self.A1(i,k)*func(i,k)/np.sum([func2(_i,k) for _i in np.arange(1, self.mm)])/self.W0

	        return self._FF[int(i*self.kk + k)]

	def GG(self, U):
		_sum = 0
		for _k in np.arange(self.kk):
			if self._sumFF[int(_k)] == 0.0:
				self._sumFF[int(_k)] = np.sum([self.FF(_i,_k) for _i in np.arange(1, self.mm)])
			_sum += self.gammaPart(_k, U)*self._sumFF[int(_k)]
			#print(_k, "-", self._sumFF[int(_k)])
		#print("+")
		return _sum

#a = NonlinearPart()
#for i in range(1000):
	#print(a.get(0.001*(i+i*1j)))
