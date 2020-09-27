# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 19:10:45 2020

@author: Luiz Renato
"""

"""
Esse codigo simula o comportamento do motor em MF
"""
# Importando bibliotecas
from scipy.integrate import odeint
from numpy import zeros, ones, linspace, pi, mean
from random import uniform
import matplotlib.pyplot as plt


# Definindo dinamica da planta
def din_motor(y, t, u):
    x1, x2 = y#x1->theta x->velocidade 
    # Dinamica do motor
    x2p = -20*3.1415/pi*2*(0.5 + (4 - 4) )*x2 + 104.7198*3.1415/pi*(1 + (2 - 2) )*u#aceleracao
    return [x2, x2p]


# Parametros de simulacao
Ta = 1e-4  # intervalo de amostragem
Tsim = 6 # tempo de simulacao
kend = int(Tsim/Ta)#numero de passos

# scopes
u = zeros(kend)  # % de duty cycle
omega = zeros(kend)  # velocidade angular do motor
theta = zeros(kend)  # posicao angular real do motor
theta_med = zeros(kend)  # posicao angular medida do motor
e = zeros(kend)  # erro de rastreamento
#kp=12.26
kp=100
alfa=0.0

# Referencia
theta_ref = zeros(kend)
# loop
for k in range(kend-1):
    # Calculo do erro de rastreamento
    e[k] = theta_ref[k] -theta[k]  # utilize theta_med

    # Aplicando entrada apos 1s (k*Ta = tempo atual)
    if k*Ta >= 1:
        #theta_ref[k]=1.4*k*Ta
        theta_ref[k]=120*pi/180
        e[k] = theta_ref[k] -theta[k]
        u[k] =  e[k]*kp-alfa*(theta[k]-theta[k-1])/Ta

    # Limitando acao de controle
    u[k] = min(u[k], 100)  # sup
    u[k] = max(u[k], -100)  # inf

    # Evoluindo a din. da planta
    sol = odeint(din_motor, [theta[k], omega[k]], [0.0, Ta], args=(u[k],))
    theta[k + 1] = sol[:, 0][-1]
    theta_med[k + 1] = theta[k + 1] + uniform(-0.05, 0.05)  # o termo da direita eh o ruido
    omega[k+1] = sol[:, 1][-1]



# Plotando resultados
fig1 = plt.figure()
plt.plot(linspace(0, kend, kend)*Ta, theta_med*180/pi, lw=2, label=r'\theta_{med} (deg)')
plt.plot(linspace(0, kend, kend)*Ta, theta_ref*180/pi*ones(kend), lw=2, label=r'\theta_{ref} (deg)')
plt.xlabel('Tempo (s)')
plt.legend()
plt.grid(True)

plt.show()