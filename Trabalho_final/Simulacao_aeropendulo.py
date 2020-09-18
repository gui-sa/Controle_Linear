#Este codigo vem para simular o comportamento de um aeropendulo
#Materia: controle linear

# Importando bibliotecas
from scipy.integrate import odeint
import numpy as np
from random import uniform
import matplotlib.pyplot as plt
import tqdm

#%% Dinamica do aeropendulo

# Definindo dinamica da planta
#Esta funcao segue a documentacao da funcao odeint : https://www.youtube.com/watch?v=PfXJWa4TrXY (17-09-2020)
def din_aeropendulo(y, t, v):#dinamica do aeropendulo y = theta_pp /  t= tempo discreto /v = rotacao do motor
    theta,theta_p = y #y é uma lista com [theta,theta_p]
    Lh = 0.32 #metros
    Kh = 2.12829*(10 ** -5)#N/(rad/s)²
    I = 0.0264 #kgm²
    m = 0.3182 # kg
    g = 9.81 #m/s²
    b = 0.006856 #(rad/s)^(-1)
    omega = theta_p #velocidade do aeropendulo
    omega_p = [ omega , ((Lh*Kh/I)*(v ** 2) - (Lh*m*g/I)*np.sin(theta) - (b/I)*omega) ]
    return omega_p #omega_p = theta_pp



#%%======= Parametros de simulacao
# Parametros de simulacao
Ta = 1e-3  # intervalo de amostragem
Tsim = 50 #tempo final
kend = int(Tsim/Ta) # quantidade de iteracoes

# scopes
v =np.zeros(kend)  # velocidade de rotacao do motor
theta = np.zeros(kend) # Angulo do pendulo com a estrutura
theta_med = np.zeros(kend) # Angulo do pendulo com a estrutura
theta_p = np.zeros(kend) # valocidade de aeropendulo
theta_pp = np.zeros(kend) # aceleracao do aeropendulo
erro = np.zeros(kend)  #Erro do sistema
theta_ref = np.zeros(kend) # Angulacao requerida ou referenciada


#Parametros do controlador
Kp = 20000 #Só pra testar mesmo

#%%=======================Loop percorrendo o tempo do experimento
for k in tqdm.tqdm((range(kend-1))):

    if(Ta*k >= 10):
        theta_ref[k] = 40 #graus
    
    theta_ref[k] = theta_ref[k]*np.pi/180  #converte para radianos
    
    erro[k] = theta_ref[k] - theta_med[k]
    v[k]= Kp*erro[k]
        
    v[k] = min(v[k], 375)  # maximo 375 rad/s
    v[k] = max(v[k], 0)  # minimo 0
    
    sol = odeint(din_aeropendulo, [theta[k] , theta_p[k]], [ 0.0, Ta ], args= (v[k],))
    theta[k+1],theta_p[k+1] = sol[1,:] 
    theta_med[k+1] = theta[k+1] + uniform(-0.01, 0.01)
#%%======================Plotando resultado:

plt.figure()
plt.plot(np.linspace(0, Tsim,kend) , theta*180.0/np.pi, lw =2.0 , color="k", label = "Posicao aeropendulo")
plt.plot(np.linspace(0, Tsim,kend) , theta_ref*180.0/np.pi, lw =2.0 , color="r", label = "Referencia-angulo")
plt.xlabel("Tempo [s]")
plt.ylabel("Angulo [graus]")
plt.title("Simulacao do aeropendulo")
plt.legend()
plt.show()

plt.Figure()
plt.plot(np.linspace(0, Tsim,kend), v , lw = 2.0, color = "b", label = "Velocidade rotacao das helices")
plt.xlabel("Tempo [s]")
plt.ylabel("Velociade de rotação [rad/s]")
plt.title("Simulacao do aeropendulo - velocidade de rotação")
plt.legend()
plt.show()
    
    
    
    