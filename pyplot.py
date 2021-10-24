import numpy as np
import matplotlib.pyplot as plt

file = 'state_ana.txt'

field = np.loadtxt(file)

plt.pcolor(field)
plt.show()
