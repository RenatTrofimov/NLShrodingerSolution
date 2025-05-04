import numpy as np
from numpy import sin, cos, pi, sqrt, exp
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.special import gamma
from joblib import Parallel, delayed
import multiprocessing
import time
import scipy
import numpy as np
import matplotlib.pyplot as plt
from math import factorial, sqrt, pi, cos, exp
import cupy as cp
import pdb
import cupy as cp
import mpmath
import numpy as np
import matplotlib.pyplot as plt
from math import factorial, pi, cos
import cupyx
def custom_factorial(n):
    return cp.exp(cupyx.scipy.special.gammaln(n + 1))

def alfaP(x, kappa, r, dx):
	return 2 * kappa * 1j / r - 2 - dx**2 * x

def alfaM(x, kappa, r, dx):
	return -2 * kappa * 1j / r - 2 - dx**2 * x

def log_gamma_approx(n):
	# Аппроксимация логарифма гамма-функции (формула Стирлинга)
	return (n - 0.5) * cp.log(n) - n + 0.5 * cp.log(2 * cp.pi)

def nlp(x, G, r, l):
	# Вычисление аппроксимации gammaln
	gammaln_l2 = log_gamma_approx(l + 2)
	gammaln_l3 = log_gamma_approx(l + 3)
	
	return cp.sum(
		G * (r + 1) ** (2 * (l + 1)) *
		x ** (l + 1) * ((-1) ** (l + 1)) * (l + 1) /
		(cp.exp(gammaln_l2) * cp.exp(gammaln_l3) * (2 ** (2 * (l + 1))))
	)

def Besse_CNT():
	# Инициализация mpmath с высокой точностью
	mpmath.mp.dps = 10  # Количество значащих цифр
	
	# Конвертация констант в mpmath формате
	def mpf(x):
		return mpmath.mpf(str(x))
	
	# Simulation parameters
	M = 1000
	tEnd = mpf(2)
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
	
	steps = 1
	
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
	# Main simulation loop
	gV1 = cp.array(V1).astype(cp.complex64)
	gU0 = cp.array(U0).astype(cp.complex64)
	gU1 = cp.array(U0).astype(cp.complex64)
	gCoeff = cp.array(Coeff)
	gG = cp.array(G)
	gA_plus = cp.array(A_plus).astype(cp.complex64)
	gA_minus = cp.array(A_minus).astype(cp.complex64)
	kernel_code = r'''
	#include <cupy/complex.cuh>

	extern "C" __global__
	void complex_multiply(
		const complex<float>* V1,
		const float* Coeff,
		const float* G,
		complex<float>* A_plus,
		complex<float>* A_minus,
		const float dx,
		const float kappa,
		const float _r,
		const int N_garm,
		const int l_max,
		const int n
	) {
		int idx = min(blockIdx.x * blockDim.x + threadIdx.x, n - 1);
		
		complex<float> a_plus_val(0.0f, 0.0f);
		
		for(int r = 0; r < N_garm; r++) {
			for(int l = 0; l < l_max; l++) {
				float temp1 = powf((float)(r + 1), (float)(2 * (l + 1)));
				
				complex<float> temp2(1.0f, 0.0f);
				for(int p = 0; p < (2 * (l + 1)); p++) {
					temp2 *= V1[idx];
				}
				
				a_plus_val += G[r] * temp1 * temp2 * Coeff[l];
			}
		}
		
		float dx_sq = dx * dx;
		complex<float> kappa_term(0.0f, 2.0f * kappa / _r);
		
		A_plus[idx * n + idx] = kappa_term - complex<float>(2.0f, 0.0f) - complex<float>(dx_sq, 0.0f) * a_plus_val;
		A_minus[idx * n + idx] = -1.0f * (-kappa_term - complex<float>(2.0f, 0.0f) - complex<float>(dx_sq, 0.0f) * a_plus_val);
	}
	'''
	complex_mult_kernel = cp.RawKernel(kernel_code, 'complex_multiply')

	threads_per_block = 256
	blocks_per_grid = (M + threads_per_block - 1) // threads_per_block
	from cupyx.scipy.sparse.linalg import spsolve
	import cupyx.scipy.sparse as cusp
	while nn < Nt:
		gV1 = cp.absolute(gU0)**2
		complex_mult_kernel(
			(blocks_per_grid,), 
			(threads_per_block,), 
			(gV1, gCoeff, gG, gA_plus, gA_minus, float(str(dx)), float(str(kappa)), float(str(r)), N_garm, l_max, M)
		)
		print( nn )
		gU0 = cp.dot(gA_minus, gU0)
		#gU1 = cupyx.scipy.sparse(gA_plus, gU0)
		gU1 = spsolve(cusp.csr_matrix(gA_plus), gU0) 
		gV1 = cp.absolute(gU1)**2
		complex_mult_kernel(
			(blocks_per_grid,), 
			(threads_per_block,), 
			(gV1, gCoeff, gG, gA_plus, gA_minus, float(str(dx)), float(str(kappa)), float(str(r)), N_garm, l_max, M)
		)
		gU1 = cp.dot(gA_minus, gU1)
		#gU0 = cp.linalg.solve(gA_plus, gU1)
		gU0 = spsolve(cusp.csr_matrix(gA_plus), gU1)
		nn += 1
	
	# Convert results to float for plotting
	plt.plot(x, abs(cp.asarray(gU1))/A0)
	plt.show()

Besse_CNT()
