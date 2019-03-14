import matplotlib.pyplot as plt
from numpy import log10, genfromtxt
ERR = genfromtxt('error.csv', delimiter=',')
plt.plot(log10(ERR[5:, 0]), log10(ERR[5:, 1]), 'r', label="Euler")
plt.plot(log10(ERR[5:, 0]), log10(ERR[5:, 2]), 'g', label="RK2")
plt.plot(log10(ERR[5:, 0]), log10(ERR[5:, 3]), 'g', label="RK4")
plt.grid()
plt.xlabel('$log_{10}(Steps)$')
plt.ylabel('$log_{10}(Error)$')
plt.title('Error convergence')
plt.legend()
plt.savefig("out/error.pdf", bbox_inches='tight', pad_inches=0)
plt.show()
