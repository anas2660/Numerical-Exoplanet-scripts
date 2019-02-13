import matplotlib.pyplot as plt
import numpy as np


eulerdata = np.genfromtxt('eulerint.csv',
                     delimiter=',')
rungekuttadata = np.genfromtxt('rungekuttaint.csv',
                     delimiter=',')

plt.figure(1)
plt.plot(eulerdata[:, 0], eulerdata[:, 1], 'r')
plt.plot(eulerdata[:, 0], eulerdata[:, 2], 'g')
plt.grid()
plt.ylim([-0.5, 1.0])
plt.xlabel('$\\xi$')
plt.ylabel('$\\theta$')
plt.title('Euler integration')

plt.figure(2)
plt.plot(rungekuttadata[:, 0], rungekuttadata[:, 1], 'r')
plt.plot(rungekuttadata[:, 0], rungekuttadata[:, 2], 'g')
plt.grid()
plt.ylim([-0.5, 1.0])
plt.xlabel('$\\xi$')
plt.ylabel('$\\theta$')
plt.title('Runge-kutta integration')

plt.show()


