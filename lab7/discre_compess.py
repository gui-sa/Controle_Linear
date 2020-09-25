"""
Codigo para obter o equivalente digital para C(s)
"""
import control.matlab as control

# Parametros calculados do compensador
a = 10.5
b = 1.42

# Definindo FT do compensador
C = control.tf([1, a], [1, b])

# Intervalo de amostragem
Ta = 2e-3  # s.

# Obtendo o equivalente de discreto pelo metodo de Tustin
Cd = control.c2d(C, Ta, 'tustin')

# Imprimindo valores
print(Cd)

