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
#angulos
theta = np.zeros(kend)# Angulo do pendulo com a estrutura real [rad]
phi = theta-theta_eq# Angulo do pendulo com o ponto de equilibrio [rad]
theta_med = np.zeros(kend) # Angulo do pendulo com a estrutura
phi_med=theta_med-theta_eq# Angulo do pendulo linearizado
phi_estimado = np.zeros(kend) # Angulo do pendulo com a estrutura

theta_p = np.zeros(kend) # velocidade de aeropendulo real
theta_pp = np.zeros(kend) # aceleracao do aeropendulo real
erro = np.zeros(kend)  #Erro do sistema
phi_ref = np.zeros(kend) # Angulacao requerida ou referenciada linearizada
theta_ref=np.zeros(kend)#Angulo de referencia
#Referencias
theta_ref[:]=0 #0 [rad]
phi_ref[:]=theta_ref[:]-theta_eq # [rad]

#Parametros do controlador 
#K = np.array([-21.52366751,   2.96947505]) #obtidos por ackerman e tr 1 segundo
K = np.array([[0.86434184, 6.19977139]])

#Ke = np.array([[202.19694102],
#                [4474.20068109]])#obtidos por ackerman e tr 1 segundo

Ke =0# np.array([[ 20.23242662, 141.64814536]])
Ke = np.transpose(Ke)

A = np.array([[0,1],
              [-28.9847261,-0.25969]])

B = np.array([[0],
              [1]])

C = np.array([0.15841992,0])

x_til = np.zeros((kend,2,1))
y_til = np.zeros(kend)
x_til_ponto = np.zeros((kend,2,1))

rss=np.zeros(kend)
xss=np.zeros((kend,2,1))
uss=np.zeros(kend)
erro=np.zeros(kend)

Nx=np.array([[6.31233749],
                 [0]])
Nu=182.96137317

#%%=======================Loop percorrendo o tempo do experimento
for k in tqdm.tqdm((range(kend-1))):
    
    if(Ta*k >= 0):#REFERENCIA
        theta_ref[k] = 40*np.pi/180#angulo de referencial real [rad] 
        phi_ref[k] = theta_ref[k]-theta_eq[k] #angulo para o controlador[rad]
        rss[k]=phi_ref[k]
    if(Ta*k >= 20):
        theta_ref[k] = 30*np.pi/180#angulo de referencial real [rad] 
        phi_ref[k] = theta_ref[k]-theta_eq[k] #angulo para o controlador[rad]
        rss[k]=phi_ref[k]
    if(Ta*k >= 40):
        theta_ref[k] = 30*np.pi/180#angulo de referencial real [rad] 
        phi_ref[k] = theta_ref[k]-theta_eq[k] #angulo para o controlador[rad]
        rss[k]=phi_ref[k]
    if(Ta*k >= 44):
        theta_ref[k] = 34*np.pi/180#angulo de referencial real [rad] 
        phi_ref[k] = theta_ref[k]-theta_eq[k] #angulo para o controlador[rad]
        rss[k]=phi_ref[k]
    if(Ta*k >= 48):
        theta_ref[k] = 38*np.pi/180#angulo de referencial real [rad] 
        phi_ref[k] = theta_ref[k]-theta_eq[k] #angulo para o controlador[rad]
        rss[k]=phi_ref[k]
    if(Ta*k >= 52):
        theta_ref[k] = 42*np.pi/180#angulo de referencial real [rad] 
        phi_ref[k] = theta_ref[k]-theta_eq[k] #angulo para o controlador[rad]
        rss[k]=phi_ref[k]
    if(Ta*k >= 56):
        theta_ref[k] = 46*np.pi/180#angulo de referencial real [rad] 
        phi_ref[k] = theta_ref[k]-theta_eq[k] #angulo para o controlador[rad]
        rss[k]=phi_ref[k]              
    if(Ta*k >= 60):
        theta_ref[k] = 50*np.pi/180#angulo de referencial real [rad] 
        phi_ref[k] = theta_ref[k]-theta_eq[k] #angulo para o controlador[rad]
        rss[k]=phi_ref[k]
    #CONTROLADOR
    
    uss[k]=Nu*rss[k]
    xss[k]=np.dot(Nx,rss[k])
    
    
    v[k] = uss[k] - np.dot(K,(x_til[k]-xss[k])) #Realimentação de estados
    
    #v[k]=np.dot(-K,x_til[k])
    omega[k]=v[k]+omega_eq[k]#rotacao controlador->linear
        
    
    #ATUADOR
    omega[k] = min(omega[k], 375)  # maximo 375 rad/s
    omega[k] = max(omega[k], 0)  # minimo 0
    
    #SISTEMA NAO LINEAR
    sol = odeint(din_aeropendulo, [theta[k], theta_p[k]], [ Ta*k, Ta*(k+1) ], args= (omega[k],))#args recebe a velocidade real do motor (ou seja 150 ms atrasado)
    theta[k+1]=sol[1,0]
    phi[k+1]=theta[k+1]-theta_eq[k+1]
    theta_p[k+1] = sol[1,1]
    theta_med[k+1] = theta[k+1]# + uniform(-0.02, 0.02)
    phi_med[k+1]=theta_med[k+1]-theta_eq[k+1]


    #Observador:
    y_til[k]  = np.dot(C,x_til[k]) #Saida estimada
    erro[k] = phi_med[k] - y_til[k] 
    x_til_ponto[k] = np.dot(A, x_til[k]) + np.dot(B,v[k]) + np.dot(Ke, erro[k]) #x_til_ponto estimado
    
    x_til[k+1] = x_til[k] + x_til_ponto[k]*Ta #integracao numerica por retangulo de Euler


#%%======================Plotando resultado:

plt.figure()
plt.plot(np.linspace(0, Tsim,kend) , theta*180.0/np.pi, lw =2.0 , color="k", label = "Posicao aeropendulo")
theta_ref[-1]=theta_ref[-2]
plt.plot(np.linspace(0, Tsim,kend) , theta_ref*180.0/np.pi, lw =2.0 , color="r", label = "Referencia-angulo")
plt.xlabel("Tempo [s]")
plt.ylabel("Angulo [graus]")
plt.title("Simulacao do aeropendulo")
plt.grid()
plt.legend()
plt.show()

plt.figure()
omega[-1]=omega[-2]
plt.plot(np.linspace(0, Tsim,kend), omega , lw = 2.0, color = "b", label = "Velocidade rotacao das helices")
plt.xlabel("Tempo [s]")
plt.ylabel("Velociade de rotação [rad/s]")
plt.title("Simulacao do aeropendulo - velocidade de rotação - ação de controle")
plt.grid()
plt.legend()
plt.show()