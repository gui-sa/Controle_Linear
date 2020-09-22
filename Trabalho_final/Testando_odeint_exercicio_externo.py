#Este codigo vem para simular o comportamento de um aeropendulo
#Materia: controle linear

# Importando bibliotecas
from scipy.integrate import odeint
import numpy as np
from random import uniform
import matplotlib.pyplot as plt
import tqdm

#%% EDO

def EDO_2ordem(y, t):#dinamica do aeropendulo y = theta_pp /  t= tempo discreto /v = rotacao do motor
    return  [ y[1] , -2*y[1] - 2*y[0] + 2*np.cos(2*t)]
#y′′+2y′+2y=cos(2x),y(0)=0,y′(0)=0
#z≡y′⇒z′+2z+2y=cos(2x),z(0)=y(0)=0.


#%%======= Parametros de simulacao
# Parametros de simulacao
Ta = 1e-3  # intervalo de amostragem = 1 ms
Tsim = 10 #tempo final
kend = int(Tsim/Ta) # quantidade de iteracoes


# scopes
y = np.zeros(kend) # Y
y_p = np.zeros(kend) # derivada de y no tempo
y_pp = np.zeros(kend) # segunda derivada de y no tempo


#%%=======================Loop percorrendo o tempo do experimento
#for k in range(kend-1):
for k in tqdm.tqdm((range(kend-1))):
    
    
    sol = odeint(EDO_2ordem, [y[k] , y_p[k]], [ k*Ta , (k + 1)*Ta ])
    y[k+1],y_p[k+1] = sol[1,:] 
#%%======================Plotando resultado:

plt.figure()
plt.plot(np.linspace(0, Tsim,kend) , y, lw =2.0 , color="k", label = "EDO")
plt.xlabel("X")
plt.ylabel("Y")
plt.title("Resolucao da EDO")
plt.legend()
plt.show()

    
    
    
    