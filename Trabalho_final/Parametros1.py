#O objetivo deste calculo é obter os valores do compensador baseado no requisito de projeto:
 
import numpy as np
import matplotlib.pyplot as plt
from control import tf, feedback, bode_plot, pzmap, pade

    
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
tss = 20 #segundos

bode_plot([G1,G],  dB =True ,  omega_num=30000)

polos = pzmap(G , Plot = False)

print ("\n polos de G = " + str(polos[0]) + "\n\n")
#Estes sao os valores dos polos



#%%  funcoes para obtencao dos parametros

def csi_sis(MS):
    return np.sqrt((np.log(MS/100) ** 2)/(np.pi ** 2 + np.log(MS/100) ** 2))# coeficiente de amortecimento


def tr_sis(csi, wn):
    return (np.pi - np.arccos(csi))/(wn * np.sqrt(1 - csi ** 2 )) #Tempo de subida 

def tss_sis(csi, wn):
    return 4.0/(csi*wn)#tempo de estabelecimento

def wn_by_csi_tss(csi, tss):
    return 4.0/(csi*tss)#

#%% 

csi = csi_sis(10) #este csi é o minimo para alcançar o maximo de 10% de maximo sobressinal
csi = np.linspace(csi,0.99,150)#Criei um vetor de csis maiores que o de cima 

wn = wn_by_csi_tss(csi, tss)
tr = tr_sis(csi, wn)
menor_tr = np.min(tr)
print("\nMenor tempo de subida= " + str(menor_tr) + " para csi = " + str(csi_sis(10))  + " e wn = " + str(wn_by_csi_tss(csi_sis(10), tss)) + "\n\n")

plt.figure()
plt.plot(wn, tr, lw= 2.0, color = "k",label = "Relacao wn com tempo de subida")
plt.plot(csi, tr, lw= 2.0, color = "b", label = "Relacao csi com tempo de subida")
plt.legend()
plt.grid()
plt.ylabel("Tempo de subida [s]")
plt.title("Obtendo o minimo tempo de subida possivel")
plt.show()

validacao_1 = csi*wn  #Para bater com tss, csi*wn é um constante e vale 4.0/tss ou seja 0.2
#Concluindo: temos que usar o csi no limiar do MS  [csi = csi_sis(10)] ou seja 0,5911550337988974


#%%



