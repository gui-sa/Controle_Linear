#O objetivo deste calculo é obter os valores do compensador baseado no requisito de projeto:
 
import numpy as np
import matplotlib.pyplot as plt
from control import tf, feedback, bode_plot, pzmap, pade, margin, rlocus

    
#%% Planta
s = tf("s")

Lh = 0.32 #metros
Kh = 2.12829*(10 ** -5)#N/(rad/s)²
I = 0.0264 #kgm²
m = 0.3182 # kg
g = 9.81 #m/s²
b = 0.006856 #(rad/s)^(-1)
Td = 0.15 #segundos
G = ((2*Lh*Kh/I)*np.sqrt((m*g*np.sin(2*np.pi/9))/Kh))/(s*s + (b/I)*s + (Lh*m*g*np.cos(2*np.pi/9.0)/I))#Funcao transferencia da nossa planta em MA
ordem_pade = 3
TfTD = tf(pade(Td,ordem_pade)[0],pade(Td,ordem_pade)[1])
G1 = G*TfTD
#%% Requisitos de projeto

MS = 10 #%

# gm,pm,wf,wc = margin(G)
# print("\nO sistema G(s) sem atraso possui GM =" + str(gm) +", PM = " + str(gm) + ", frequencia de corte (Wc) = " + str(wc) + ", e frequencia de fase(Wf) = " + str(wf))


#%%  funcoes para obtencao dos parametros

def csi_sis(MS):
    return np.sqrt((np.log(MS/100) ** 2)/(np.pi ** 2 + np.log(MS/100) ** 2))# coeficiente de amortecimento


def tr_sis(csi, wn):
    return (np.pi - np.arccos(csi))/(wn * np.sqrt(1 - csi ** 2 )) #Tempo de subida 

def tss_sis(csi, wn):
    return 4.0/(csi*wn)#tempo de estabelecimento

def wn_by_csi_tss(csi, tss):
    return 4.0/(csi*tss)#


#%% Comecando o calculo do compensador
