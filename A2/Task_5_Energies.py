import matplotlib.pyplot as plt
from ase.db import connect
from ase.db import connect
from ase.io import read, write
import numpy as np
import seaborn as sns
import vesta

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
# Set the plot style
set_plot_style()
fig_width = 7

# Load the trajectory file
trajectory_file = 'A2_7/all_candidates.traj'
atoms = read(trajectory_file, index=':')

energies = []
for i in range(len(atoms)):
    energy = atoms[i].get_total_energy()
    energies.append(energy)
    if i==0:
        print(energy)
    if i==60:
        print(energy)    
    

plt.plot(range(len(atoms)), energies, linewidth=2.0) 
plt.tight_layout() 
plt.title('Change in Energy during Trajectory - 7 atoms', fontsize=2 * fig_width)
plt.xlabel('Configuration', fontsize=2 * fig_width)
plt.ylabel('Energy [eV]', fontsize=2 * fig_width)
plt.tick_params(axis='both', which='both', labelsize=1.8 * fig_width)
plt.tight_layout()
plt.savefig('Figures/Task_5_Energy_7.png')
plt.show()   