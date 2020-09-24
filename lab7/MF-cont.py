"""
Esse codigo simula o comportamento do motor em MF
"""
# Importando bibliotecas
from scipy.integrate import odeint
from numpy import zeros, linspace, pi, ones
from random import uniform
import matplotlib.pyplot as plt


# Definindo dinamica da planta
def din_motor(y, t, u):
    omega = y
    # Dinamica do motor
    omega_p = -20*3.1415/pi*2*(0.5 + (4 - 4) )*omega + 104.7198*3.1415/pi*(1 + (2 - 2) )*u

    return omega_p


# Parametros de simulacao
Ta = 2e-3  # intervalo de amostragem
Tsim = 2
kend = int(Tsim/Ta)

# scopes
u = zeros(kend)  # % de duty cycle
e = zeros(kend)  # erro de rastreamento
omega = zeros(kend)  # velocidade real do motor
omega_med = zeros(kend)  # velocidade medida do motor
# Referencia
omega_ref =  200.0*ones(kend)  # rad. Note que theta_ref eh um vetor (ou lista)

a1_barra =  1.009
a2_barra = -0.9881
b_barra = -0.9972


# loops
for k in range(kend-1):
    tempo_atual = k*Ta
    if tempo_atual >= 1:
        # Calculo do erro de rastreamento
        e[k] = omega_ref[k] - omega_med[k]   # Utilize a velocidade medida omega_med no calculo doe rro

    # Calculando a entrada
    u[k] =   -b_barra*u[k-1] + a1_barra*e[k] + a2_barra*e[k-1]

    # Limitando acao de controle
    u[k] = min(u[k], 100)  # sup
    u[k] = max(u[k], -100)  # inf

    # Evoluindo a din. da planta
    sol = odeint(din_motor, omega[k], [0.0, Ta], args=(u[k],))
    omega[k+1] = sol[:, 0][-1]
    omega_med[k + 1] = omega[k+1] + uniform(-5, 5)  # o termo da direita eh o ruido


# Plotando resultados
fig1 = plt.figure()
plt.plot(linspace(0, kend, kend)*Ta, omega_ref, lw=2, label='\omega_ref (rad/s)')
plt.plot(linspace(0, kend, kend)*Ta, omega_med, lw=2, label='\omega (rad/s)')
plt.plot(linspace(0, kend, kend-1)*Ta, u[0:-1], lw=2, label='u (%)')
plt.xlabel('Tempo (s)')
plt.legend()
plt.grid(True)
plt.show()