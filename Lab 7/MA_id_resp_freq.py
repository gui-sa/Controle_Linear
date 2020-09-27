"""
Esse codigo simula o comportamento do motor em MA
"""
# Importando bibliotecas
from scipy.integrate import odeint
from numpy import zeros, linspace, pi, sin
from random import uniform
import matplotlib.pyplot as plt
import tqdm

# Definindo dinamica da planta
def din_motor(y, t, u):
    omega = y
    # Dinamica do motor
    omega_p = -20*3.1415/pi*2*(0.5 + (4 - 4) )*omega + 104.7198*3.1415/pi*(1 + (2 - 2) )*u

    return omega_p


# Parametros de simulacao
Ta = 2e-3  # intervalo de amostragem
Tsim = 1
kend = int(Tsim/Ta)

# scopes
u = zeros(kend)  # % de duty cycle
omega = zeros(kend)  # velocidade real do motor
omega_med = zeros(kend)  # velocidade medida do motor

# loops
for k in tqdm.tqdm(range(kend-1)):
    tempo_atual = k*Ta

    # Calculando a entrada
    freq = 1000#freqencua
    u[k] = 80*sin(freq*k*Ta) #Amplitude 80 e freq de 0.1 rad/s

    # Limitando acao de controle
    u[k] = min(u[k], 100)  # sup
    u[k] = max(u[k], -100)  # inf

    # Evoluindo a din. da planta
    sol = odeint(din_motor, omega[k], [0.0, Ta], args=(u[k],))
    omega[k+1] = sol[:, 0][-1]
    omega_med[k + 1] = omega[k+1] + uniform(-5, 5)  # o termo da direita eh o ruido


#%% Plotando resultados
fig1 = plt.figure()
plt.plot(linspace(0, kend, kend)*Ta, omega_med, lw=2, label='\omega (rad/s)')
plt.plot(linspace(0, kend, kend-1)*Ta, u[0:-1], lw=2, label='u (%)')
plt.xlabel('Tempo (s)')
plt.legend()
plt.grid(True)
plt.show()