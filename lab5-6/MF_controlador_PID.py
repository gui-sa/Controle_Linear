# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 16:52:37 2020

@author: Luiz Renato
"""

"""
Esse codigo simula o comportamento do motor em MA
"""
# Importando bibliotecas
from scipy.integrate import odeint
from numpy import zeros, ones, arange, pi
from random import uniform
import matplotlib.pyplot as plt


# Definindo dinamica da planta
def din_motor(y, t, u):
    x1, x2, x3 = y
    # Dinamica do motor
    x3p = -180*x3 - 20*160*x2 + 104.7198*160*u - 5000*x1
    return [x2, x3, x3p]


# Parametros de simulacao
Ta = 1e-4  # intervalo de amostragem
Tsim = 5  # tempo de simulacao
kend = int(Tsim/Ta)

# scopes
u = zeros(kend)  # controle total
u_p = zeros(kend)  # acao proporcional
u_i = zeros(kend)  # acao derivativa
u_d = zeros(kend)  # acao integral
omega_p = zeros(kend)  # aceleracao angular do motor
omega = zeros(kend)  # velocidade angular do motor
theta = zeros(kend)  # posicao angular real do motor
theta_med = zeros(kend)  # posicao angular medida do motor
e = zeros(kend)  # erro de rastreamento
e_i = zeros(kend)  # integral do erro
e_d = zeros(kend)  # derivada do erro

#constantes do PID
td=0.0139
ti=0.0555
kp=20.268
ki=kp/ti
kd=td*kp

# Referencia
theta_ref = 120*pi/180*ones(kend)  # rad. Note que theta_ref eh um vetor (ou lista)

# loop
for k in range(kend-1):
    # Aplicando entrada apos 1s (k*Ta = tempo atual)
    if k * Ta >= 1:
        e[k] = theta_ref[k]-theta_med[k]  # utilize theta_med (medido) no calculo do erro

    # Calculo da integral e da derivada do erro
    e_d[k] = (e[k] - e[k - 1]) / Ta  # derivada do erro
    e_i[k] = e_i[k - 1] + e[k] * Ta  # integral do erro
    # Calculando as acoes de controle e o controle total
    u_p[k] = e[k]*kp
    u_d[k] = e_d[k]*kd
    #if (u[k-1]==-100 or u[k-1]==100):
     #   u_i[k] =0
    #else:
    u_i[k] = e_i[k]*ki

    u[k] = u_p[k]+u_d[k]+u_i[k]

    # Limitando acao de controle
    u[k] = min(u[k], 100)  # sup
    u[k] = max(u[k], -100)  # inf

    # Evoluindo a din. da planta
    sol = odeint(din_motor, [theta[k], omega[k], omega_p[k]], [0.0, Ta], args=(u[k],))
    theta[k + 1] = sol[:, 0][-1]
    theta_med[k + 1] = theta[k + 1] + uniform(-0.05, 0.05)  # o termo da direita eh o ruido
    omega[k+1] = sol[:, 1][-1]
    omega_p[k + 1] = sol[:, 2][-1]


# Plotando resultados
fig1 = plt.figure()
plt.plot(arange(0, Tsim, Ta), theta_med*180/pi, lw=2, label=r'\theta_{med} (deg)')
plt.plot(arange(0, Tsim, Ta), theta_ref*180/pi*ones(kend), lw=2, label=r'\theta_{ref} (deg)')
plt.xlabel('Tempo (s)')
plt.legend()
plt.grid(True)

fig2 = plt.figure()
plt.plot(arange(0, Tsim-Ta, Ta), u[0:-1], lw=2, label=r'u')
plt.xlabel('Tempo (s)')
plt.legend()
plt.grid(True)

plt.show()