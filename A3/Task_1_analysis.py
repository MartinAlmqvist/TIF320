import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from ase.io.trajectory import Trajectory
from ase.io import read, write

data = np.genfromtxt('logfile.txt', dtype=float)


# Plotting
x = data[:, 0]
temp = data[:, 4]
E_tot = data[:, 1]
E_pot =data[:, 2]
E_kin = data[:, 3]

print(len(x))

plt.plot(x, temp)
plt.show()
