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
#ordem_pade = 3
#TfTD = tf(pade(Td,ordem_pade)[0],pade(Td,ordem_pade)[1])
#G1 = G*TfTD
#%% Requisitos de projeto


#rlocus(G1)
#sisotool(G1)

teste = tf2ss(G)  #Trabalhando com G, e deixando os polos dominantes em relacao aos polos do PADE
A = np.flip(teste.A)
B = np.flip(teste.B)
C = np.flip(teste.C)
D = np.flip(teste.D)


ABCD=np.array([[0, 1, 0],
            [-28.9847261,  -0.25969697, 1],
            [0.15841992, 0, 0]])


polos = np.linalg.eig(A)[0]

control_matrix = ctrb(A, B)

control_bool = np.linalg.det(control_matrix)
if (control_bool != 0):
    control_bool = True
    
observability_matrix = obsv(A, C)

observability_bool = np.linalg.det(control_matrix)
if (observability_bool != 0):
    observability_bool = True


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
csi = csi_sis(10) #Trabalhando com MS = 10%
polos_desejados = polos_by_csi_wn(csi, wn_by_csi_tr(csi, 0.5))

#polos_desejados = [(-0.1298+5.3822j),(-0.1298-5.3822j)]

K = acker(A, B, polos_desejados)#Obtendo a matriz linha do regulador K 
Ke = acker(np.transpose(A), np.transpose(C),np.real(np.array(polos_desejados))*5)#Obtendo matriz coluna Ke do estimador. Os polos do estimador devem ser 10x maiores (mais rapidos)

Inv=np.linalg.inv(ABCD)

N=np.dot(Inv,np.array([[0],
                       [0],
                       [1]]))



polos_controladora_checkup=np.linalg.eig(A-np.dot(B,K))[0]
polos_observer_checkup = np.linalg.eig(A-np.dot(np.transpose(Ke),C))[0]

