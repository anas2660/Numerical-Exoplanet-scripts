import matplotlib.pyplot as plt
import numpy as np


data = np.genfromtxt('eulerint.csv',
                     delimiter=',')

plt.plot(data[:,0], data[:,1], 'r')
plt.plot(data[:,0], data[:,2], 'g')
plt.grid()
plt.ylim([-0.5,1.0])
plt.show()


