import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags, eye
from scipy.sparse.linalg import inv

import pdb

import mpmath
import numpy as np
import matplotlib.pyplot as plt
from math import factorial, pi, cos, exp
import scipy
def Besse_CNT():
	# Инициализация mpmath с высокой точностью
	mpmath.mp.dps = 10  # Количество значащих цифр
	
	# Конвертация констант в mpmath формате
	def mpf(x):
		return mpmath.mpf(str(x))
	
	# Simulation parameters
	M = 2000
	tEnd = mpf(20)
	xEnd = mpf(10)
	dt = mpf('0.008')
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
	E0 = mpf('10e6')  # V/cm
	E0 = mpf(1)*(E0)/mpf(300)
	A0 = E0*mpmath.fabs(el)*a_CNT/(h1*omega)
	print(A0)
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
	A_shtr = 0.1
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
		if not nn%625:
			plt.plot(x, abs(U1)/A0)
	
	# Convert results to float for plotting
	plt.plot(x, abs(U1)/A0)
	#np.save(f"D:/A = 0/{dt*nn}.npy", U0 )
	plt.show()

Besse_CNT()