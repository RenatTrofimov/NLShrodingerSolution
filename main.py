import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags, eye
from scipy.sparse.linalg import inv

<<<<<<< HEAD
import pdb

import mpmath
import numpy as np
import matplotlib.pyplot as plt
from math import factorial, pi, cos

def Besse_CNT():
	# Инициализация mpmath с высокой точностью
	mpmath.mp.dps = 10  # Количество значащих цифр
	
	# Конвертация констант в mpmath формате
	def mpf(x):
		return mpmath.mpf(str(x))
	
	# Simulation parameters
	M = 1000
	tEnd = mpf(20)
	xEnd = mpf(20)
	dt = mpf('0.001')
	dx = mpf(2)*xEnd/mpf(M)
	r = dt/(mpf(2)*dx**2)
	lamda = mpf(1)
	nn = 1
	Nt = tEnd/dt
	jj = mpmath.mpc(1j)  # мнимая единица

	# Physical constants
	el = mpf('-4.8e-10')
	h1 = mpf('1.055e-27')
	k_B = mpf('1.38e-16')  # Boltzmann constant (CGS)
	
	# CNT Constants and parameters
	gamma0 = mpf('2.7') * mpf('1.6e-12')
	b_CNT = mpf('0.142e-7')
	a_CNT = (mpf(3)/mpf(2)) * b_CNT
	m = 7
	T = mpf(77)  # Temperature
	rel_perm = mpf(4)
	el_concentration = mpf('1e18')
	omega0 = mpf(2)*mpmath.fabs(el)*a_CNT*mpmath.sqrt(pi*el_concentration*gamma0)/h1
	
	# Beam wave-vector
	omega = mpf('4e14')
	kappa = mpf(2)*mpmath.sqrt(rel_perm)*omega/omega0
	l_max = 5
	
	# Energy spectrum function
	def energy_spectrum(ksi, s):
		return gamma0 * mpmath.sqrt(1 + 4*mpmath.cos(ksi)*mpmath.cos(pi*s/mpf(m)) + 
								4*mpmath.cos(pi*s/mpf(m))**2)
	
	N_garm = 9
	g = mpf('0.25')
	E0 = mpf('0.10e7')  # V/cm
	E0 = mpf(1)*(E0)/mpf(300)
	A0 = E0*mpmath.fabs(el)*a_CNT/(h1*omega)
	
	# Initialize arrays
	Integral_1 = mpmath.matrix(N_garm, m)
	delta = mpmath.matrix(N_garm, m)
	delta_relative = mpmath.matrix(N_garm, m)
	Integral_1_0 = mpmath.matrix(1, m)
	delta_0 = mpmath.matrix(1, m)
	delta_0_relative = mpmath.matrix(1, m)
	Integral_2 = mpmath.matrix(N_garm, m)
	Summa_v_integralah = mpmath.matrix(1, m)
	Integral_3 = mpmath.matrix(m, 1)
	F = mpmath.matrix(N_garm, m)
	G = mpmath.matrix(N_garm, 1)
	
	steps = 1000
	
	# Integral calculations with mpmath
	for r in range(N_garm):
		for s in range(m):
			Integral_1[r,s] = mpf(0)
			for k in range(1, 2*steps + 1):
				p = -mpmath.pi + (mpmath.pi/steps)*k
				Integral_1[r,s] += (mpmath.pi/steps)*energy_spectrum(p,s)*mpmath.cos((r+1)*p)
			
			delta[r,s] = (1/mpmath.pi)*Integral_1[r,s]
			delta_relative[r,s] = delta[r,s]/gamma0
	
	for s in range(m):
		Integral_1_0[0,s] = mpf(0)
		for k in range(1, 2*steps + 1):
			p = -mpmath.pi + (mpmath.pi/steps)*k
			Integral_1_0[0,s] += (mpmath.pi/steps)*energy_spectrum(p,s)
		
		delta_0[0,s] = (1/mpmath.pi)*Integral_1_0[0,s]
		delta_0_relative[0,s] = delta_0[0,s]/gamma0
	
	for r in range(N_garm):
		for s in range(m):
			Integral_2[r,s] = mpf(0)
			for k in range(1, 2*steps + 1):
				p = -mpmath.pi + (mpmath.pi/steps)*k
				Summa_v_integralah[0,s] = mpf(0)
				for r2 in range(N_garm):
					Summa_v_integralah[0,s] += (delta[r2,s]/(k_B*T))*mpmath.cos((r2+1)*p)
				
				exponent = (delta_0[0,s]/(2*k_B*T)) + Summa_v_integralah[0,s]
				term = 1 + mpmath.exp(exponent)
				Integral_2[r,s] += (mpmath.pi/steps)*mpmath.cos((r+1)*p)/term
	print( "i2" )
	for s in range(m):
		Integral_3[s,0] = mpf(0)
		for k in range(1, 2*steps + 1):
			p = -mpmath.pi + (mpmath.pi/steps)*k
			Summa_v_integralah[0,s] = mpf(0)
			for r2 in range(N_garm):
				Summa_v_integralah[0,s] += (delta[r2,s]/(k_B*T))*mpmath.cos((r2+1)*p)
			
			exponent = (delta_0[0,s]/(2*k_B*T)) + Summa_v_integralah[0,s]
			term = 1 + mpmath.exp(exponent)
			Integral_3[s,0] += (mpmath.pi/steps)/term
	
	Summa_v_znamenatele = mpmath.fsum(Integral_3[:,0])
	print( "i3" )
	for r in range(N_garm):
		for s in range(m):
			F[r,s] = -(r+1)*(delta[r,s]/gamma0)*(Integral_2[r,s]/Summa_v_znamenatele)
	print( "F[r,0]" )
	for r in range(N_garm):
		G[r,0] = mpmath.fsum(F[r,:])
		
	print( "G[r,0]" )
	# Initialize arrays for simulation
	x = np.zeros(M)
	Resh = np.zeros(M)
	U0 = np.zeros(M, dtype=complex)
	U1 = np.zeros(M, dtype=complex)
	V0 = np.zeros(M)
	V1 = np.zeros(M)
	nonlin = np.zeros(M)
	nonlin2 = np.zeros(M)
	alfa_plus = np.zeros(M, dtype=complex)
	alfa_minus = np.zeros(M, dtype=complex)
	A_plus = np.zeros((M, M), dtype=complex)
	A_minus = np.zeros((M, M), dtype=complex)
	
	# Initial conditions
	for i in range(M):
		x[i] = 0 + dx * (i+1)
		U0[i] = A0 * exp(-(x[i] - xEnd)**2 / g)
	
	# Setup A_plus and A_minus matrices
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
	
	# Coefficients for nonlinear term
	Coeff = [mpf(0) for _ in range(l_max)]
	for l in range(l_max):
		Coeff[l] = ((-1)**(l+1))*(l+1)/(mpmath.factorial(l+1)*mpmath.factorial(l+2)*(2**(2*(l+1))))
	Coeff = np.array([float(str(i)) for i in Coeff ])
	
	G = np.array([float(str(i)) for i in G ])
	A_shtr = 0
	G1 = np.array([G[r] * np.cos((r + 1) * A_shtr) for r in np.arange(G.shape[0])])
	G2 = np.array([G[r] * np.sin((r + 1) * A_shtr) for r in np.arange(G.shape[0])])
	
	II = np.eye(M, dtype=complex)
	# Main simulation loop
	while nn < Nt:
		
		for i in range(M):
			V0[i] = abs(U0[i])**2
			V1[i] = -V0[i] + 2 * abs(U0[i])**2
			
			nonlin[i] = 0
			nonlin2[i] = 0
			for r in range(N_garm):
				for l in range(l_max):
					nonlin[i] += G1[r] * (r)**(2*l+1) * V1[i]**(l) * Coeff[l]
					nonlin2[i] +=  G2[r] * r**(2*l) * V1[i]**(l) * Coeff[l]*((l+1)-r**2*V1[i]/2/(l+2))
			r = dt/(mpf(2)*dx**2)
			alfa_plus[i] = 2 * kappa * jj / r - 2 - dx**2 * nonlin[i]
			alfa_minus[i] = -2 * kappa * jj / r - 2 - dx**2 * nonlin[i]
			A_plus[i, i] = alfa_plus[i]
			A_minus[i, i] = -alfa_minus[i]
			II[i,i] = 2*dx*dx*nonlin2[i]
		print( nn )
		B = np.dot(A_minus, U0) + float(str(2*dx*dx))*nonlin2
		U1 = scipy.sparse.linalg.spsolve(scipy.sparse.csr_matrix(A_plus), B)
		U0 = U1.copy()
		nn += 1
		np.save(f"D:/A = 0/{dt*nn}.npy", U0 )
		if not nn%250:
			plt.plot(x, abs(U1)/A0)
	
	# Convert results to float for plotting
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
		pdb.set_trace()
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
=======
def Besse_CNT_Acoustic():
	# Simulation parameters
	M = 4000
	tEnd = 25
	xEnd = 15
	dt = 0.008
	dx = 2 * xEnd / M
	rr = dt / (2 * dx**2)
>>>>>>> 9d154b12b692bde4edc37a7ea00f6131fcfd59bd
	lamda = 1
	nn = 1
	Nt = int(tEnd / dt)
	j = 1j  # imaginary unit

	el = -4.8e-10  # electron charge (CGS)
	h1 = 1.055e-27  # reduced Planck constant (CGS)
	k_B = 1.38e-16  # Boltzmann constant (CGS)
	
	# CNT Constants and parameters
	gamma0 = 2.7 * 1.6e-12
	b_CNT = 0.142e-7
	a_CNT = (3/2) * b_CNT
	m = 7
	T = 77  # Temperature (K)
	rel_perm = 4
	el_concentration = 1e18
	omega0 = 2 * abs(el) * a_CNT * np.sqrt(np.pi * el_concentration * gamma0) / h1

	# Beam wave-vector
	omega = 1e14
	kappa = 2 * np.sqrt(rel_perm) * omega / omega0  # wave vector
	l_max = 5
	energy_spectrum = lambda ksi, s: gamma0 * np.sqrt(1 + 4 * np.cos(ksi) * np.cos(np.pi * s / m) + 
									4 * np.cos(np.pi * s / m)**2)
	N_garm = 9
	g = 0.25
	E0 = 1e4  # V/cm
	E0 = E0 / 300  # Convert to SGS units
	A0 = E0 * abs(el) * a_CNT / (h1 * omega)
	A_shtr = 0.1  # A0 * 0.05

	# Initialize arrays
	Integral_1 = np.zeros((N_garm, m), dtype=np.longdouble)  # или np.float64
	delta = np.zeros((N_garm, m), dtype=np.longdouble)
	delta_relative = np.zeros((N_garm, m), dtype=np.longdouble)
	Integral_1_0 = np.zeros(m, dtype=np.longdouble)
	delta_0 = np.zeros(m, dtype=np.longdouble)
	delta_0_relative = np.zeros(m, dtype=np.longdouble)
	Integral_2 = np.zeros((N_garm, m), dtype=np.longdouble)
	Summa_v_integralah = np.zeros(m, dtype=np.longdouble)
	Integral_3 = np.zeros(m, dtype=np.longdouble)
	F = np.zeros((N_garm, m), dtype=np.longdouble)
	G = np.zeros(N_garm, dtype=np.longdouble)
	G1 = np.zeros(N_garm, dtype=np.longdouble)
	G2 = np.zeros(N_garm, dtype=np.longdouble)

	# Number of steps for integral calculation
	steps = 10000
	p_values = np.linspace(-np.pi, np.pi, 2 * steps, endpoint=False)

	# Compute Integral_1[r, s] and delta[r, s]
	for r in range(N_garm):
		for s in range(m):
			energy = energy_spectrum(p_values, s)
			Integral_1[r, s] = np.sum(energy * np.cos((r + 1) * p_values)) * (np.pi / steps)
			delta[r, s] = (1 / np.pi) * Integral_1[r, s]
			delta_relative[r, s] = delta[r, s] / gamma0

	# Compute Integral_1_0[s] and delta_0[s]
	for s in range(m):
		energy = energy_spectrum(p_values, s)
		Integral_1_0[s] = np.sum(energy) * (np.pi / steps)
		delta_0[s] = (1 / np.pi) * Integral_1_0[s]
		delta_0_relative[s] = delta_0[s] / gamma0
	
	def safe_exp_division(x):
		# Сначала обрабатываем крайние случаи
		result = np.zeros_like(x)
		mask_pos = x > 100
		mask_neg = x < -100
		mask_mid = ~(mask_pos | mask_neg)
	
		result[mask_pos] = 0.0
		result[mask_neg] = 1.0
		# Только для средних значений вычисляем exp
		x_mid = x[mask_mid]
		result[mask_mid] = 1 / (1 + np.exp(x_mid))
	
		return result

	# Compute Integral_2[r, s]
	for r in range(N_garm):
		for s in range(m):
			sum_in_exp = np.zeros(2 * steps)
			for r2 in range(N_garm):
				sum_in_exp += (delta[r2, s] / (k_B * T)) * np.cos((r2 + 1) * p_values)
			exponent = (delta_0[s] / (2 * k_B * T)) + sum_in_exp
			integrand = np.cos((r + 1) * p_values) * safe_exp_division(exponent)
			Integral_2[r, s] = np.sum(integrand) * (np.pi / steps)

	# Compute Integral_3[s]
	for s in range(m):
		sum_in_exp = np.zeros(2 * steps)
		for r2 in range(N_garm):
			sum_in_exp += (delta[r2, s] / (k_B * T)) * np.cos((r2 + 1) * p_values)
		exponent = (delta_0[s] / (2 * k_B * T)) + sum_in_exp
		integrand = safe_exp_division(exponent)
		Integral_3[s] = np.sum(integrand) * (np.pi / steps)

	Summa_v_znamenatele = np.sum(Integral_3)
	# 1. Добавим проверку знаменателя с небольшой константой EPS
	EPS = 1e-15  # Малое число для защиты от деления на 0
	# Compute F[r, s]
	for r in range(N_garm):
		for s in range(m):
			# Проверяем знаменатель
			safe_denominator = Summa_v_znamenatele if abs(Summa_v_znamenatele) > EPS else EPS
			
			# Проверяем числитель
			numerator = Integral_2[r, s]
			if not np.isfinite(numerator):
				numerator = 0.0  # Заменяем NaN/inf на 0
			
			# Безопасное вычисление
			term1 = delta[r, s] / gamma0
			term2 = numerator / safe_denominator
			F[r, s] = - (r + 1) * term1 * term2
		
			# Дополнительная проверка результата
			if not np.isfinite(F[r, s]):
				F[r, s] = 0.0

	# Compute G[r], G1[r], G2[r]
	for r in range(N_garm):
		G[r] = np.sum(F[r, :])
		G1[r] = G[r] * np.cos((r + 1) * A_shtr)
		G2[r] = G[r] * np.sin((r + 1) * A_shtr)

	# Coefficients for nonlinear terms
	Coeff = np.zeros(l_max)
	from scipy.special import factorial
	for l in range(l_max):
		Coeff[l] = ((-1)**l) * (l + 1) / (factorial(l + 1) * factorial(l + 2) * (2**(2 * l)))

	# Initialize spatial grid and solution arrays
	x = np.linspace(0, 2 * xEnd, M, endpoint=False)
	U0 = A0 * np.exp(-(x - xEnd)**2 / g)
	V0 = np.abs(U0)**2
	U1 = np.zeros(M, dtype=complex)
	V1 = np.zeros(M)

	# Construct sparse matrices for periodic boundary conditions
	diagonals = [[1] * (M - 1), [-1] * (M - 1)]
	offsets = [1, -1]
	A_plus = diags(diagonals, offsets, shape=(M, M), dtype=complex)
	A_minus = diags(diagonals, offsets, shape=(M, M), dtype=complex)
	A_plus += diags([1, 1], [-(M - 1), M - 1], shape=(M, M))  # Periodic BC
	A_minus += diags([-1, -1], [-(M - 1), M - 1], shape=(M, M))  # Periodic BC

	II = eye(M, dtype=complex)
	from scipy.sparse.linalg import spsolve
	# Time-stepping loop
	# Предварительные вычисления (вынесены за цикл)
	r_values = np.arange(1, N_garm+1)  # [1, 2, ..., N_garm]
	l_values = np.arange(l_max)
	r_powers = r_values[:, np.newaxis] ** (2*l_values + 1)  # (N_garm × l_max)
	r_powers_2 = r_values[:, np.newaxis] ** (2*l_values)    # (N_garm × l_max)

	for nn in range(1, Nt + 1):
		V1 = -V0 + 2 * np.abs(U0)**2
		
		# Векторизованное вычисление nonlin1 и nonlin2
		V1_expanded = V1[:, np.newaxis, np.newaxis]  # (M × 1 × 1)
		V1_pows = V1_expanded ** l_values             # (M × 1 × l_max)
		
		# Вычисление nonlin1 (векторизованная версия)
		term1 = G1[:, np.newaxis] * r_powers * Coeff  # (N_garm × l_max)
		nonlin1 = np.sum(term1 * V1_pows, axis=(1,2))  # Сумма по r и l
		
		# Вычисление nonlin2 (векторизованная версия)
		term2_part = (l_values + 2) - (r_values[:, np.newaxis]**2 * V1_expanded) / (2 * (l_values + 3))
		term2 = G2[:, np.newaxis] * r_powers_2 * Coeff * term2_part
		nonlin2 = np.sum(term2 * V1_pows, axis=(1,2))
		
		# Обновление матриц
		alfa_plus = 2 * kappa * j / rr - 2 - dx**2 * nonlin1
		alfa_minus = -2 * kappa * j / rr - 2 - dx**2 * nonlin1
		
		A_plus.setdiag(alfa_plus)
		A_minus.setdiag(-alfa_minus)
		II.setdiag(2 * dx**2 * nonlin2)
		
		# Решение системы
		B = A_minus.dot(U0) + II.dot(U0)
		U1 = spsolve(A_plus.tocsc(), B)
		
		# Вывод прогресса
		if nn % 125 == 0:
			print(f"Iteration {nn}, U1[1000] = {U1[1000]}")
		
		# Обновление переменных
		U0, V0 = U1.copy(), V1.copy()
		
		# Визуализация
		if (nn) % 625 == 0:
			plt.figure(figsize=(10, 6))
			plt.plot(x, np.abs(U1) / A0)
			plt.title(f"Wave amplitude at iteration {nn}")
			plt.xlabel("Position")
			plt.ylabel("Normalized amplitude")
			plt.grid(True)
			plt.show()

print("Simulation completed!")

if __name__ == "__main__":
	Besse_CNT_Acoustic()