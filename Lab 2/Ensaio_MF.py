"""
Esse codigo simula o comportamento do motor em MA
"""
# Importando bibliotecas
from scipy.integrate import odeint
from numpy import zeros, linspace
from random import uniform
import matplotlib.pyplot as plt



# Definindo dinamica da planta
def din_motor(y, t, u, am):#dinamica do motor que depende de velocidade, tempo, duty cycle, polos e ganho
    Km = 1000  # ganho do motor, sobre a tensao
    omega = y #velocidade do motor
    # Dinamica do motor
    omega_p = -am*omega + Km*u #aceleracao = -polos*velocidade + ganho*duty cycle
    
    return omega_p #aceleracao motor


# Parametros de simulacao
am = uniform(0.8, 1.2)*20  # polo do motor, quantidade de ima na carcaça, essa quantidade seria variada
Ta = 10e-3  # intervalo de amostragem
Tsim = 10 #tempo final
kend = int(Tsim/Ta) # quantidade de iteracoes

# scopes-> cria vetores nulos
u = zeros(kend)  # % de duty cycle
omega = zeros(kend)  # velocidade do motor
e = zeros(kend)  # erro de rastreamento
int_e = zeros(kend)  # integral do erro

# Referencia
omega_ref = 3000  # velocidade de referencia (rpm)

#Constantes
Kp=0.15 #constante proporcional
Ki=0.1 #constante de integracao
# loop percorrendo o tempo do experimento
for k in range(kend-1):
    # Calculo do erro de rastreamento
    e[k] = omega_ref-omega[k]  #erro de rastramento nesse instante (rpm)
    int_e[k]=int_e[k-1]+Ta*e[k] #calculo da integral para o instante k
    
    # Amplitude de entrada da planta apos 1s
    if k*Ta >= 1:#se o tempo for maior que 1 segundo (resposta ao degraum duty cycle em 1s)
        u[k] = e[k]*Kp + int_e[k]*Ki #atribua o valor do duty cycle

    # Limitando acao de controle->limita o duty cycle entre 0 e 100
    u[k] = min(u[k], 100)  # maximo 100
    u[k] = max(u[k], 0)  # minimo 0


    # Evoluindo a dinamica da planta
    sol = odeint(din_motor, omega[k], [0.0, Ta], args=(u[k], am))#resolve a EDO da dinamica do motor, encontrando a velocidade do motor (a cada for atualiza todos os parametros, menos o tempo, utiliza incremento de tempos, usa um PVI do instante anterior)
    omega[k+1] = sol[:, 0][-1] + uniform(-15, 15)  #velocidade = velocidade encontrada na EDO + ruido


# Plotando resultados
fig1 = plt.figure() #cria ou ativa uma figura
plt.plot(linspace(0, kend, kend)*Ta, omega, lw=2, label='\omega (rpm)')# plota a velocidade x tempo
plt.plot(linspace(0, kend, kend-1)*Ta, u[0:-1], lw=2, label='u (%)')#plota o duty cycle x tempo
plt.xlabel('Tempo (s)') #legenda do x
plt.legend()
plt.grid(True)
plt.show()