import matplotlib.pyplot as plt
from numpy import genfromtxt
eulerdata = genfromtxt('eulerint.csv', delimiter=',')
rungekuttadata = genfromtxt('rungekuttaint.csv', delimiter=',')
def plotsetup(n):
    plt.figure(n)
    plt.grid()
    plt.ylim([-0.5, 1.0])
    plt.xlabel('$\\xi$')
    plt.ylabel('$\\theta$')
plotsetup(1)
plt.plot(eulerdata[:, 0], eulerdata[:, 1], 'r')
plt.plot(eulerdata[:, 0], eulerdata[:, 2], 'g')
plt.title('Euler integration')
plotsetup(2)
plt.plot(rungekuttadata[:, 0], rungekuttadata[:, 1], 'r')
plt.plot(rungekuttadata[:, 0], rungekuttadata[:, 2], 'g')
plt.title('Runge-kutta integration')
plt.show()
