import matplotlib.pyplot as plt
import numpy as np


errdata = np.genfromtxt('error.csv',
                        delimiter=',')
cutoff = 10
plt.plot(np.log(errdata[cutoff:, 0]), np.log(errdata[cutoff:, 1]), 'r', label="Euler")
plt.plot(np.log(errdata[cutoff:, 0]), np.log(errdata[cutoff:, 2]), 'g', label="RK4")
#print(errdata[0, 2])
plt.grid()
plt.xlabel('$Steps$')
plt.ylabel('$Error$')
plt.title('Error convergence')
plt.legend()
plt.show()

