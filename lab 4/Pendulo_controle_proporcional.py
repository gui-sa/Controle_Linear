# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 18:04:07 2020

@author: Luiz Renato
"""

"""
Esse codigo simula o comportamento pendulo invertido
"""
# Importando bibliotecas
from scipy.integrate import odeint
from numpy import zeros, ones, arange, pi, sin
from random import uniform
import matplotlib.pyplot as plt


# Definindo dinamica da planta
def din_pend(y, t, tau):
    # Parametros da planta
    m = 0.007  # kg
    L = 0.0475  # m
    g = 9.81  # m/s^2

    x1, x2 = y
    # Dinamica do pendulo
    x2p = 1/(m*L**2)*tau - g/L*sin(x1)
    return [x2, x2p]


# Parametros de simulacao
Ta = 1e-3  # intervalo de amostragem
Tsim = 1.5  # tempo de simulacao
kend = int(Tsim/Ta)

# scopes
tau = zeros(kend)  # % de duty cycle
thetap = zeros(kend)  # velocidade angular do motor
theta = zeros(kend)  # posicao angular real do motor
e = zeros(kend)  # erro de rastreamento

# Referencia
theta_ref = 0*pi/180*ones(kend)  # rad. Em relacao a theta = pi

# Condicao inicial
theta[0] = 30*pi/180  # rad
kp=0.9
# Loop
for k in range(kend-1):
    # Calculo do erro de rastreamento
    e[k] = theta_ref[k]-theta[k]

    # Calculo da acao de controle
    tau[k] = e[k]*kp

    # Limitando acao de controle
    tau[k] = min(tau[k], 100)  # sup
    tau[k] = max(tau[k], -100)  # inf

    # Evoluindo a din. da planta
    x0 = [theta[k] + pi, thetap[k]]   # condicao inicial
    sol = odeint(din_pend, x0, [0.0, Ta], args=(tau[k],))
    theta[k + 1] = sol[:, 0][-1] - pi
    thetap[k+1] = sol[:, 1][-1]


# Plotando resultados
fig1 = plt.figure()
plt.plot(arange(0, Tsim, Ta), theta*180/pi, lw=2, label=r'\theta_{med} (deg)')
plt.plot(arange(0, Tsim, Ta), theta_ref*180/pi*ones(kend), lw=2, label=r'\theta_{ref} (deg)')
plt.xlabel('Tempo (s)')
plt.legend()
plt.grid(True)

fig2 = plt.figure()
plt.plot(arange(0, Tsim-Ta, Ta), tau[0:-1], lw=2)
plt.xlabel('Tempo (s)')
plt.ylabel(r'\tau (N)')
plt.grid(True)

plt.show()