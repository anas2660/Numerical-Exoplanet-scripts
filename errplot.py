import matplotlib.pyplot as plt
import numpy as np

errdata = np.genfromtxt('error.csv',
                        delimiter=',')
cutoff = 10
plt.plot(np.log(errdata[cutoff:, 0]),
         np.log(errdata[cutoff:, 1]), 'r', label="Euler")
plt.plot(np.log(errdata[cutoff:, 0]),
         np.log(errdata[cutoff:, 2]), 'g', label="RK4")
plt.grid()
plt.xlabel('$ln(Steps)$')
plt.ylabel('$ln(Error)$')
plt.title('Error convergence')
plt.legend()
plt.savefig("out/error.pdf", bbox_inches='tight', pad_inches=0)
plt.show()
