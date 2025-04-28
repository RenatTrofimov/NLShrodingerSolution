import numpy as np
from numpy import sin, cos, pi, sqrt, exp
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.special import gamma
from joblib import Parallel, delayed
import multiprocessing
import time

import numpy as np
import matplotlib.pyplot as plt
from math import factorial, sqrt, pi, cos, exp

def Besse_CNT():
    # Simulation parameters
    M = 1000
    tEnd = 2
    xEnd = 10
    dt = 0.001
    dx = 2 * xEnd / M
    r = dt / (2 * dx * dx)
    lamda = 1
    nn = 1
    Nt = int(tEnd / dt)
    jj = 1j  # imaginary unit in Python

    # Physical constants
    el = -4.8e-10
    h1 = 1.055e-27
    k_B = 1.38e-16  # Boltzmann constant (CGS)
    
    # CNT Constants and parameters
    gamma0 = 2.7 * 1.6e-12
    b_CNT = 0.142e-7
    a_CNT = (3/2) * b_CNT
    m = 7
    T = 77  # Temperature
    rel_perm = 4
    el_concentration = 1e18
    omega0 = 2 * abs(el) * a_CNT * sqrt(pi * el_concentration * gamma0) / h1
    
    omega = 4e14
    kappa = 2 * sqrt(rel_perm) * omega / omega0  # wave vector
    l_max = 5
    
    def energy_spectrum(ksi, s):
        return gamma0 * sqrt(1 + 4 * cos(ksi) * cos(pi * s / m) + 4 * cos(pi * s / m)**2)
    
    N_garm = 9
    g = 0.25
    E0 = 0.10e7  # V/cm
    E0 = 1 * (E0) / 300  # SGS(E)
    A0 = E0 * abs(el) * a_CNT / (h1 * omega)
    
    Integral_1 = np.zeros((N_garm, m))
    delta = np.zeros((N_garm, m))
    delta_relative = np.zeros((N_garm, m))
    Integral_1_0 = np.zeros(m)
    delta_0 = np.zeros(m)
    delta_0_relative = np.zeros(m)
    Integral_2 = np.zeros((N_garm, m))
    Summa_v_integralah = np.zeros(m)
    Integral_3 = np.zeros(m)
    F = np.zeros((N_garm, m))
    G = np.zeros(N_garm)
    
    steps = 10000
    
    # Integral_1[r,s]
    for r in range(N_garm):
        for s in range(m):
            Integral_1[r, s] = 0
            for k in range(1, 2*steps + 1):
                p = -pi + (pi/steps) * k
                Integral_1[r, s] += 1 * (pi/steps) * energy_spectrum(p, s) * cos((r+1) * p)
            
            # Coefficients in the Fourier expansion: delta[r,s]
            delta[r, s] = (1/pi) * Integral_1[r, s]
            delta_relative[r, s] = delta[r, s] / gamma0
    
    # Integral_1[0,s], when r = 0
    for s in range(m):
        Integral_1_0[s] = 0
        for k in range(1, 2*steps + 1):
            p = -pi + (pi/steps) * k
            Integral_1_0[s] += 1 * (pi/steps) * energy_spectrum(p, s) * cos(0 * p)
        
        delta_0[s] = (1/pi) * Integral_1_0[s]
        delta_0_relative[s] = delta_0[s] / gamma0
    
    # Integral_2[r,s]
    for r in range(N_garm):
        for s in range(m):
            Integral_2[r, s] = 0
            for k in range(1, 2*steps + 1):
                p = -pi + (pi/steps) * k
                Summa_v_integralah[s] = 0
                for r2 in range(N_garm):
                    Summa_v_integralah[s] += (delta[r2, s]/(k_B*T)) * cos((r2+1)*p)
                
                Integral_2[r, s] += 1 * (pi/steps) * cos((r+1)*p) / (1 + exp((delta_0[s]/(2*k_B*T)) + Summa_v_integralah[s]))
    
    # Integral_3[r,s]
    for s in range(m):
        Integral_3[s] = 0
        for k in range(1, 2*steps + 1):
            p = -pi + (pi/steps) * k
            Summa_v_integralah[s] = 0
            for r2 in range(N_garm):
                Summa_v_integralah[s] += (delta[r2, s]/(k_B*T)) * cos((r2+1)*p)
            
            Integral_3[s] += 1 * (pi/steps) / (1 + exp((delta_0[s]/(2*k_B*T)) + Summa_v_integralah[s]))
    
    # Summa_v_znamenatele
    Summa_v_znamenatele = np.sum(Integral_3)
    
    # F[r,s]
    for r in range(N_garm):
        for s in range(m):
            F[r, s] = - (r+1) * (delta[r, s]/gamma0) * (Integral_2[r, s]/Summa_v_znamenatele)
    
    # G[r]
    for r in range(N_garm):
        G[r] = np.sum(F[r, :])
    
    x = np.zeros(M)
    Resh = np.zeros(M)
    U0 = np.zeros(M, dtype=complex)
    U1 = np.zeros(M, dtype=complex)
    V0 = np.zeros(M)
    V1 = np.zeros(M)
    nonlin = np.zeros(M)
    alfa_plus = np.zeros(M, dtype=complex)
    alfa_minus = np.zeros(M, dtype=complex)
    A_plus = np.zeros((M, M), dtype=complex)
    A_minus = np.zeros((M, M), dtype=complex)
    
    for i in range(M):
        x[i] = 0 + dx * (i+1)
        U0[i] = A0 * exp(-(x[i] - xEnd)**2 / g)
    
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
    
    Coeff = np.zeros(l_max)
    for l in range(l_max):
        Coeff[l] = ((-1)**(l+1)) * (l+1) / (factorial(l+1) * factorial(l+2) * (2**(2*(l+1))))
    
    while nn < Nt:
        for i in range(M):
            V0[i] = abs(U0[i])**2
            V1[i] = -V0[i] + 2 * abs(U0[i])**2
            
            nonlin[i] = 0
            for r in range(N_garm):
                for l in range(l_max):
                    nonlin[i] += G[r] * (r+1)**(2*(l+1)) * V1[i]**(l+1) * Coeff[l]
            
            alfa_plus[i] = 2 * kappa * jj / r - 2 - dx**2 * nonlin[i]
            alfa_minus[i] = -2 * kappa * jj / r - 2 - dx**2 * nonlin[i]
            A_plus[i, i] = alfa_plus[i]
            A_minus[i, i] = -alfa_minus[i]
        
        B = np.dot(A_minus, U0)
        U1 = np.linalg.solve(A_plus, B)
        
        U0 = U1.copy()
        nn += 1
    
    plt.plot(x, abs(U1)/A0)
    plt.show()

Besse_CNT()

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