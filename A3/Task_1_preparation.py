import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from ase.io.trajectory import Trajectory
from ase.io import read, write

def set_plot_style():
    sns.set_context("paper", font_scale=2)
    sns.set_style("darkgrid")
    sns.set_palette("deep")
    sns.set(font='sans-serif')

    fig_width = 8
    fig_height = 6

    plt.rcParams['figure.figsize'] = (fig_width, fig_height)
    plt.rcParams['figure.dpi'] = 100
    plt.rcParams['lines.linewidth'] = 1
    plt.rcParams['lines.markersize'] = 4
    plt.rcParams['font.size'] = 2.5 * fig_width
    plt.rcParams['xtick.labelsize'] = 1.8 * fig_height
    plt.rcParams['ytick.labelsize'] = 1.8 * fig_height
    plt.grid(True)

set_plot_style()

data = np.genfromtxt('cluster24.txt', dtype=float)


# Plotting
x = data[:, 0]
temp = data[:, 4]
E_tot = data[:, 1]
E_pot =data[:, 2]
E_kin = data[:, 3]

print(len(temp))
length = 21
bath_size = int(len(temp)/length)
print(bath_size)


averages = []

""" for i in range(len(temp)):
    average = np.mean(temp[i * bath_size: i * bath_size + bath_size])
    #average = np.mean(temp[:i])
    averages.append(average) """

for i in range(length):
    average = np.mean(temp[i * bath_size: i * bath_size + bath_size])
    #average = np.mean(temp[:i])
    averages.append(average)    

# Plot the averages
plt.plot(range(length), averages, marker='o', linestyle='-', color='b')
plt.xlabel('Batch number')
plt.ylabel('Batch average Temperature')
plt.title('Average Temperature over Batch')
plt.savefig('Figures/Task_1_Batch.png')
plt.show()

""" traj = Trajectory('cluster24.traj')
index = int(len(traj)*(4/6))
atoms = traj[index]

write("Task1.xyz",atoms) """


plt.plot(x, temp)
plt.xlabel('Time [ps]')
plt.ylabel('Temperature [K]')
plt.title('Temperature Convergence over Time')
plt.savefig('Figures/Task_1_temp.png')
plt.show()

plt.plot(x, E_tot, label = 'Total Energy')
plt.plot(x, E_pot, label = 'Potential Energy')
#plt.plot(x, E_kin, label = 'Kinetic Energy')
plt.xlabel('Time [ps]')
plt.ylabel('Energy [eV]')
plt.title('Energy Convergence over Time')
plt.legend()
plt.savefig('Figures/Task_1_energy_zoom.png')
plt.show()  