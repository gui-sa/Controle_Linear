#O objetivo deste calculo é obter os valores do compensador baseado no requisito de projeto:
 
import numpy as np
import matplotlib.pyplot as plt
from control import tf,rlocus , acker, obsv, ctrb, pole, pade, feedback, sisotool, step_response , bode_plot, gangof4, margin, pzmap, nyquist_plot, ss2tf, tf2ss, canonical_form

    
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


teste = tf2ss(G)  #Trabalhando com G, e deixando os polos dominantes em relacao aos polos do PADE
A = np.flip(teste.A)
B = np.flip(teste.B)
C = np.flip(teste.C)
D = np.flip(teste.D)
polos = np.linalg.eig(A)[0]

control_matrix = ctrb(A, B)

control_bool = np.linalg.det(control_matrix)
if (control_bool != 0):
    control_bool = True
    
observability_matrix = obsv(A, C)

observability_bool = np.linalg.det(control_matrix)
if (observability_bool != 0):
    observability_bool = True

# K = acker(A, B, [])
# Ke = acker(A, C, [])


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

def wn_by_csi_tr(csi,tr):
    return (np.pi - np.arccos(csi))/(tr * np.sqrt(1 - csi ** 2 )) #Tempo de subida 

def polos_by_csi_wn(csi, wn):
    return [-csi*wn +  wn*np.sqrt(csi**2 -1 +0j) , -csi*wn -wn*np.sqrt(csi**2 - 1 + 0j) ]
#%% Comecando o calculo do compensador
csi = 0.5911 #Trabalhando com MS = 10%
polos_desejados = polos_by_csi_wn(csi, wn_by_csi_tr(csi, 1))


K = acker(A, B, polos_desejados)#Obtendo a matriz linha do regulador K 
Ke = acker(np.transpose(A), np.transpose(C), np.array(polos_desejados)*10)#Obtendo matriz coluna Ke do estimador. Os polos do estimador devem ser 10x maiores (mais rapidos)






