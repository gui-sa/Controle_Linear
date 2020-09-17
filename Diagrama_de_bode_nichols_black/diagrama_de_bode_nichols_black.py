import numpy as np
import control.matlab as ml
import matplotlib.pyplot as plt
from scipy import signal

#num = np.array([0, 100])
#den1 = np.polymul(np.array([1, 0]), np.array([1, 2]))
#den = np.polymul(np.array([1, 1]), np.array([1, 10]))
#G = ml.tf([1], [1, 1])
num = [100]
den = [1, 11, 10]

sys = signal.TransferFunction(num, den)
G = ml.tf(num, den)

print(sys)
#print(ml.pole(sys))
#print(ml.zero(sys))

fig = plt.figure(0)
plt.xlabel('Frequência [rad/s]')
plt.ylabel('Magnitude [dB]')
#magnitude, fase, w = ml.bode(G, dB = True)
w, mag, phase = signal.bode(sys)
plt.semilogx(w, mag)
plt.grid(True, which = 'both')
#plt.plot()

plt.figure(1)
plt.xlabel('Frequência [rad/s]')
plt.ylabel('Fase [°]')
plt.semilogx(w, phase)
plt.grid(True, which = 'both')


plt.figure(2)
ml.nichols(G)
plt.plot()

plt.show()