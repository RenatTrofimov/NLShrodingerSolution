import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags, eye
from scipy.sparse.linalg import inv

def Besse_CNT_Acoustic():
	# Simulation parameters
	M = 4000
	tEnd = 25
	xEnd = 15
	dt = 0.008
	dx = 2 * xEnd / M
	rr = dt / (2 * dx**2)
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