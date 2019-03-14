import matplotlib.pyplot as plt
from numpy import log, genfromtxt
ERR = genfromtxt('error.csv', delimiter=',')
plt.plot(log(ERR[5:, 0]), log(ERR[5:, 1]), 'r', label="Euler")
plt.plot(log(ERR[5:, 0]), log(ERR[5:, 2]), 'b', label="RK2")
plt.plot(log(ERR[5:, 0]), log(ERR[5:, 3]), 'g', label="RK4")
plt.grid()
plt.xlabel('$ln(Steps)$')
plt.ylabel('$ln(Error)$')
plt.title('Error convergence')
plt.legend()
plt.savefig("out/error.pdf", bbox_inches='tight', pad_inches=0)
plt.show()
