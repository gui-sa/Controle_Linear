# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 17:19:29 2020

@author: Luiz Renato
"""

"""
Esse codigo simula o comportamento do motor em MA
"""
# Importando bibliotecas
from scipy.integrate import odeint
from numpy import zeros, linspace, pi, mean
from random import uniform
import matplotlib.pyplot as plt


# Definindo dinamica da planta
def din_motor(y, t, u):
    omega = y#velocidade angular
    # Dinamica do motor
    omega_p = -20*3.1415/pi*2*(0.5 + (4 - 4) )*omega + 104.7198*3.1415/pi*(1 + (2 - 2) )*u

    return omega_p#aceleracao


# Parametros de simulacao
Ta = 1e-4  # intervalo de amostragem
Tsim = 2#tempo de simulacao
kend = int(Tsim/Ta)#quantidade de passos

# scopes
u = zeros(kend)  # % de duty cycle
omega = zeros(kend)  # velocidade real do motor
omega_med = zeros(kend)  # velocidade medida do motor

# loops
for k in range(kend-1):
    # Amplitude de entrada da planta apos 1s
    if k*Ta >= 1:
        u[k] = 80#80% de duty cycle

    # Limitando acao de controle->limite a tensao do motor
    u[k] = min(u[k], 100)  # sup
    u[k] = max(u[k], -100)  # inf

    # Evoluindo a din. da planta
    sol = odeint(din_motor, omega[k], [0.0, Ta], args=(u[k],))
    omega[k+1] = sol[:, 0][-1]
    omega_med[k + 1] = omega[k+1] + uniform(-5, 5)  # o termo da direita eh o ruido


# Plotando resultados
fig1 = plt.figure()
plt.plot(linspace(0, kend, kend)*Ta, omega_med, lw=2, label='\omega (rpm)')
plt.plot(linspace(0, kend, kend-1)*Ta, u[0:-1], lw=2, label='u (%)')
plt.xlabel('Tempo (s)')
plt.legend()
plt.grid(True)
plt.show()