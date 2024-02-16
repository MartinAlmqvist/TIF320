import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.signal import savgol_filter
import seaborn as sns
from ase.io import read

size_factor = 3

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
    plt.rcParams['font.size'] = size_factor * fig_width
    plt.rcParams['xtick.labelsize'] = size_factor * fig_height  
    plt.rcParams['ytick.labelsize'] = size_factor * fig_height  
    plt.grid(True)

# Set the plot style
set_plot_style()
fig_width = 6

# Load trajectory file (.traj)
# Update the trajectory file path for your water simulation
trajectory_file = 'Data/cluster24.traj'
atoms = read(trajectory_file)
oxygen_atoms = [atom for atom in atoms if atom.symbol == 'O']

# Compute all distances between oxygen atoms using the Minimum Image Convention (mic)
oxygen_distances = atoms.get_all_distances(mic=True)
oxygen_distances = oxygen_distances.flatten()
print(oxygen_distances)

bins = 100
counts_good, bin_edges = np.histogram(oxygen_distances, bins, range=(0, 6))

# Length of one bin
dr_good = bin_edges[2] - bin_edges[1]

# Number of atoms and volume
N = len(oxygen_atoms)
V = atoms.get_volume()
print(N)
print(V)
total_time = 6.0000  # ps
time_step = 0.0005  # ps
timesteps = int(total_time / time_step)

# Calculate the bins
gr_good = np.zeros((bins))
r_array_good = np.zeros((bins))
for i in range(1, bins):
    r_array_good[i] = bin_edges[i]
    rho_local_good = counts_good[i] * 2
    # Use the average density for liquid water
    average_density = N / V  # Average density
    
    gr_good[i] = (rho_local_good / average_density)/timesteps
    

plt.plot(r_array_good, gr_good, label=f'raw $g(r)$', linewidth=2.0)
plt.title('Radial distribution function for liquid water', fontsize=size_factor * fig_width)
plt.ylabel('$g(r)$', fontsize=size_factor * fig_width)
plt.xlabel('r [Å]', fontsize=size_factor * fig_width)
plt.legend(fontsize=size_factor * fig_width, loc='upper right')
plt.show()
