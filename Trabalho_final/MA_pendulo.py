#Este codigo vem para simular o comportamento de um aeropendulo
#Materia: controle linear

# Importando bibliotecas
from scipy.integrate import odeint
import numpy as np
from random import uniform
import matplotlib.pyplot as plt
import tqdm

#%% Parametros do Sistema
Lh = 0.32 #metros
Kh = 2.12829*(10 ** -5)#N/(rad/s)²
I = 0.0264 #kgm²
m = 0.3182 # kg
g = 9.81 #m/s²
b = 0.006856 #(rad/s)^(-1)

#%% Dinamica do aeropendulo

# Definindo dinamica da planta
#Esta funcao segue a documentacao da funcao odeint : https://www.youtube.com/watch?v=PfXJWa4TrXY (17-09-2020)
def din_aeropendulo(y, t, v):#dinamica do aeropendulo y = theta_pp /  t= tempo discreto /v = rotacao do motor
    thetaf,theta_pf = y #y é uma lista com [theta,theta_p]
    theta_pf_e_ppf = [ theta_pf , ((Lh*Kh/I)*(v ** 2) - (Lh*m*g/I)*np.sin(thetaf) - (b/I)*theta_pf) ]#aceleracao angular
    return theta_pf_e_ppf #velocidade angular e aceleracao angular


#%%======= Parametros de simulacao
# Parametros de simulacao
Ta = 1e-3  # intervalo de amostragem
Tsim = 100 #tempo final
kend = int(Tsim/Ta) # quantidade de iteracoes

# scopes
#Equilibrio
theta_eq=np.zeros(kend)+2*np.pi/9 #40 graus [rad] theta_eq utilizado para converter theta em phi
omega_eq=np.zeros(kend)+np.sqrt((m*g*np.sin(theta_eq))/Kh) #[rad/s] omega_eq utilizado para converter omega para v
#velocidades
omega =np.zeros(kend)  # velocidade de rotacao do motor real
v =omega-omega_eq  # velocidade de rotacao do motor linearizada
v_p=np.zeros(kend)#acao de controle proporcinal
v_i=np.zeros(kend)#acao de controle integral
v_d=np.zeros(kend)#acao de controle derivativo
#angulos
theta = np.zeros(kend)# Angulo do pendulo com a estrutura real [rad]
phi = theta-theta_eq# Angulo do pendulo com o ponto de equilibrio [rad]
theta_med = np.zeros(kend) # Angulo do pendulo com a estrutura
phi_med=theta_med-theta_eq# Angulo do pendulo linearizado

theta_p = np.zeros(kend) # velocidade de aeropendulo real
theta_pp = np.zeros(kend) # aceleracao do aeropendulo real
e = np.zeros(kend)  #Erro do sistema
e_d = np.zeros(kend)  #Erro do sistema
e_i = np.zeros(kend)  #Erro do sistema
phi_ref = np.zeros(kend) # Angulacao requerida ou referenciada linearizada
theta_ref=np.zeros(kend)#Angulo de referencia
#Referencias
theta_ref[:]=0 #0 [rad]
phi_ref[:]=theta_ref[:]-theta_eq # [rad]


#%%=======================Loop percorrendo o tempo do experimento
for k in tqdm.tqdm((range(kend-1))):

    if(Ta*k >= 1):#REFERENCIA
        omega[k]=omega_eq[k]
        
    
    #ATUADOR
    omega[k] = min(omega[k], 375)  # maximo 375 rad/s
    omega[k] = max(omega[k], 0)  # minimo 0
    
    #SISTEMA NAO LINEAR
    sol = odeint(din_aeropendulo, [theta[k], theta_p[k]], [ Ta*k, Ta*(k+1) ], args= (omega[k-150],))#args recebe a velocidade real do motor (ou seja 150 ms atrasado)
    theta[k+1]=sol[1,0]
    theta_p[k+1] = sol[1,1]
    theta_med[k+1] = theta[k+1] + uniform(-0.02, 0.02)
#%%======================Plotando resultado:

plt.figure()
plt.plot(np.linspace(0, Tsim,kend) , theta*180.0/np.pi, lw =2.0 , color="k", label = "Posicao aeropendulo [º]")
omega[-1]=omega[-2]
plt.plot(np.linspace(0, Tsim,kend), omega , lw = 2.0, color = "b", label = "Velocidade rotacao das helices [rad/s]")
plt.xlabel("Tempo [s]")
plt.ylabel("Y")
plt.title("Simulacao do aeropendulo")
plt.legend()
plt.grid()
plt.show()
