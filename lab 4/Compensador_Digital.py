# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 21:40:09 2020

@author: Luiz Renato
"""

"""
Codigo para obter o equivalente digital para C(s)
"""
import control.matlab as control

# Parametros calculados do compensador
K = 0.062
a = 19.54
b = 92.14

# Definindo FT do compensador
C = K*control.tf([1, a], [1, b])

# Intervalo de amostragem
Ta = 1e-3  # s.

# Obtendo o equivalente de discreto pelo metodo de Tustin
Cd = control.c2d(C, Ta, 'tustin')
Cd_zpk = control.tf2zpk(Cd.num[0][0], Cd.den[0][0])
# Extraindo parametros
abarra = -Cd_zpk[0]
bbarra = -Cd_zpk[1]
Kbarra = Cd_zpk[2]

# Imprimindo valores
print(f'abarra = {abarra}')
print(f'bbarra = {bbarra}')
print(f'Kbarra = {Kbarra}')