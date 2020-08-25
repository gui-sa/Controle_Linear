"""
Esse codigo simula o comportamento do motor em MA
"""
# Importando bibliotecas
from scipy.integrate import odeint
from numpy import zeros, linspace
from random import uniform
import matplotlib.pyplot as plt


# Definindo dinamica da planta
def din_motor(y, t, u, am):
    Km = 1000  # ganho do motor
    omega = y
    # Dinamica do motor
    omega_p = -am*omega + Km*u

    return omega_p


# Parametros de simulacao
am = uniform(0.8, 1.2)*20  # polo do motor
Ta = 10e-3  # intervalo de amostragem
Tsim = 2
kend = int(Tsim/Ta)

# scopes
u = zeros(kend)  # % de duty cycle
omega = zeros(kend)  # velocidade do motor

# loops
for k in range(kend-1):
    # Amplitude de entrada da planta apos 1s
    if k*Ta >= 1:
        u[k] = 25

    # Limitando acao de controle
    u[k] = min(u[k], 100)  # sup
    u[k] = max(u[k], 0)  # inf

    # Evoluindo a din. da planta
    sol = odeint(din_motor, omega[k], [0.0, Ta], args=(u[k], am))
    omega[k+1] = sol[:, 0][-1] + uniform(-15, 15)  # o termo da direita eh o ruido


# Plotando resultados
fig1 = plt.figure()
plt.plot(linspace(0, kend, kend)*Ta, omega, lw=2, label='\omega (rpm)')
plt.plot(linspace(0, kend, kend-1)*Ta, u[0:-1], lw=2, label='u (%)')
plt.xlabel('Tempo (s)')
plt.legend()
plt.grid(True)
plt.show()




