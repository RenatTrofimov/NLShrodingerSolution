import numpy as np
from numpy import sin, cos, pi, sqrt, exp
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.special import gamma
from joblib import Parallel, delayed
import multiprocessing
import time

class NNs_NonlinearPart():
	def energy_spectrum(self, ksi, s):
			return self.gamma0*sqrt(1 + 4*cos(ksi)*cos(pi*s/self.m) +4*cos(pi*s/self.m)*cos(pi*s/self.m) )
	def __init__(self):
		#CNT Constants and parameters:
		self.el = -4.8e-10
		self.h1 = 1.055e-27
		self.k_B = 1.38e-16
		self.gamma0 = (2.7)*1.6e-12
		self.b_CNT = 0.142e-7
		self.a_CNT = (3/2)*self.b_CNT
		self.m = 7
		self.T = 77
		self.rel_perm=4
		self.el_concentration=1e18
		self.omega0=2*abs(self.el)*self.a_CNT*sqrt(pi*self.el_concentration*self.gamma0)/self.h1


		#beam wave-vector: -------------------------------------------------------
		self.omega=4e14 #%1e14
		self.kappa=2*sqrt(self.rel_perm)*self.omega/self.omega0 # wave vector
		self.l_max=5
		self.N_garm = 9
		self.g=0.25
		self.E0 = 0.10e7 #% V/cm
		self.E0 = 1*(self.E0)/300 #% SGS(E)
		self.A0 = self.E0*abs(self.el)*self.a_CNT/(self.h1*self.omega) 
		self.steps = 10000
		#Definitions of the arrays: ----------------------------------------------
		self.Integral_1     = np.zeros((self.N_garm, self.m))
		self.delta          = np.zeros((self.N_garm, self.m))
		self.delta_relative = np.zeros((self.N_garm, self.m))
		self.Integral_1_0   = np.zeros(self.m)
		self.delta_0	  	= np.zeros(self.m)
		self.delta_0_relative  = np.zeros(self.m)
		self.Integral_2     = np.zeros((self.N_garm, self.m))
		self.Summa_v_integralah = np.zeros(self.m)
		self.Integral_3   = np.zeros(self.m)
		self.F = np.zeros((self.N_garm, self.m))
		self.G = np.zeros((self.N_garm, self.m))

		self.Coeff = np.zeros(self.l_max)
    
		for l in np.arange (1 , self.l_max):
			self.Coeff[l] = ((-1)**l)*l/( gamma(l)*gamma(l+1)*(2**(2*l))) 
		self.int1()
	
	def int1(self):
		#% Integral_1[r,s]: --------------------------------------------------------
		for r in np.arange( 1, self.N_garm):
			for s in np.arange( 1, self.m):
				self.Integral_1[r, s] = 0
				for k in np.arange( 1 , 2*self.steps):
					p = -pi + (pi/self.steps)*k
					self.Integral_1[r, s] = self.Integral_1[r, s] + 1*(pi/self.steps)*self.energy_spectrum(p, s)*cos(r*p)
				#end
				#% Coefficients in the Fourier expansion: delta[r,s]: --------------
				self.delta[r, s] = (1/pi)* self.Integral_1[r, s]
				self.delta_relative[r,s] = self.delta[r, s]/self.gamma0
				#% -----------------------------------------------------------------
			#end
		#end
		#% -------------------------------------------------------------------------
		#% Integral_1[0,s], when r = 0: --------------------------------------------
		for s in np.arange( 1, self.m):
			self.Integral_1_0[s] = 0
			for k in np.arange( 1 , 2*self.steps):
				p = -pi + (pi/self.steps)*k
				self.Integral_1_0[s] = self.Integral_1_0[s] + 1*(pi/self.steps)*self.energy_spectrum(p, s)*cos(0*p)
			#end
		#end
		#% -------------------------------------------------------------------------
		#% delta[0,s], when r = 0: -------------------------------------------------
		for s in np.arange( 1, self.m):
			self.delta_0[s] = (1/pi)* self.Integral_1_0[s]
			self.delta_0_relative[s] = self.delta_0[s]/self.gamma0
		#end
		#% -------------------------------------------------------------------------
		#% Integral_2[r,s]: --------------------------------------------------------
		for r in np.arange( 1, self.N_garm):
			for s in np.arange( 1, self.m):
				self.Integral_2[r, s] = 0
				for k in np.arange( 1 , 2*self.steps):
					p = -pi + (pi/self.steps)*k
					self.Summa_v_integralah[s] = 0
					for r2 in np.arange( 1, self.N_garm):
						self.Summa_v_integralah[s] = self.Summa_v_integralah[s] + (self.delta[r2, s]/(self.k_B*self.T))*cos(r2*p)
					#end
					self.Integral_2[r, s] = self.Integral_2[r, s] + 1.0*(pi/self.steps)*cos(r*p)/(1.0 + exp( (self.delta_0[s]/(2*self.k_B*self.T)) + self.Summa_v_integralah[s]))
					
				#end
			#end
		#end
		#% -------------------------------------------------------------------------
		#% Integral_3[r,s]: --------------------------------------------------------
		for s in np.arange( 1, self.m):
			self.Integral_3[s] = 0
			for k in np.arange( 1 , 2*self.steps):
				p = -pi + (pi/self.steps)*k
				self.Summa_v_integralah[s] = 0
				for r2 in np.arange( 1, self.N_garm):
					self.Summa_v_integralah[s] =  self.Summa_v_integralah[s] + (self.delta[r2, s]/(self.k_B*self.T))*cos(r2*p)
				#end
				self.Integral_3[s] = self.Integral_3[s] + 1*(pi/self.steps)/(1 + exp( (self.delta_0[s]/(2*self.k_B*self.T)) + self.Summa_v_integralah[s] ))
				
			#end
		#end
		#% -------------------------------------------------------------------------
		#% self.Summa_v_znamenatele: ----------------------------------------------------
		self.Summa_v_znamenatele = 0
		for s in np.arange( 1, self.m):
			self.Summa_v_znamenatele = self.Summa_v_znamenatele + self.Integral_3[s]
		#end
		#% -------------------------------------------------------------------------

		#% F[r,s]: -----------------------------------------------------------------
		for r in np.arange( 1, self.N_garm):
			for s in np.arange( 1, self.m):
				self.F[r, s] = -r*(self.delta[r, s]/self.gamma0)*(self.Integral_2[r, s]/self.Summa_v_znamenatele)
			#end
		#end
		#% -------------------------------------------------------------------------

		#% G[r]: -------------------------------------------------------------------
		for r in np.arange( 1, self.N_garm):
			self.G[r] = 0
			for s in np.arange( 1, self.m):
				self.G[r] = self.G[r] + self.F[r, s]
			#end
		#end
	
	def nnl(self, inV, outV):
		for i in np.arange(self.M):
			for r in np.arange(1, self.N_garm):
				for l in np.arange(1, self.N_garm):
					outV[i] = outV[i] + self.G[r] *  r**(2*l) * inV[i]**l * self.Coeff[l]
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
	def chuck(self, inArray):
		return np.array([self.GG(element) for element in inArray])
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
	nlp = NNs_NonlinearPart()
	while nn < Nt:
		for i in range(M):
			V0[i] = np.abs(U0[i])**2
			V1[i] = -V0[i] + 2 * np.abs(U0[i])**2
		nlp.nnl(np.copy(V1), V1)
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