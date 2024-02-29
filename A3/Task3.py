from ase import Atoms
from ase.io import read
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.signal import savgol_filter
import seaborn as sns
from ase.io.trajectory import Trajectory
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
    plt.rcParams['font.size'] = 2.5 * fig_width
    plt.rcParams['xtick.labelsize'] = 1.8 * fig_height  
    plt.rcParams['ytick.labelsize'] = 1.8 * fig_height  
    plt.grid(True)
# Set the plot style
    
set_plot_style()
fig_width = 6


distances = np.loadtxt('distances_task3_with_duplicates.txt', delimiter='\t')
traj = Trajectory("Data/cluster24.traj")

box_lengths = traj[-1].get_cell_lengths_and_angles()[:3]  # Assuming the box doesn't change over the trajectory
num_bins = 400
rdf, bins = np.histogram(distances, bins=num_bins, density=False)

bin_centers = 0.5 * (bins[1:] + bins[:-1])
dr = bins[1] - bins[0]
box_volume = np.prod(box_lengths)
num_particles = 24
density = num_particles/(box_volume)
num_frames = 4000

rdf_normalized = rdf / (4 * np.pi * bin_centers**2 * (dr) *density*num_frames*num_particles)


range_indices = np.where((bin_centers >= 2.9) & (bin_centers <= 4.4))[0]
min_index = range_indices[np.argmin(rdf_normalized[range_indices])]
rdf_integral_long = np.trapz(rdf_normalized[:min_index] * bin_centers[:min_index]**2, bin_centers[:min_index])
print(min_index)

# Calculate coordination number using the correct formula (integral of RDF up to first minimum)
coordination_number_short = 4 * np.pi * rdf_integral_long*density
print(coordination_number_short)

# Plotting
plt.plot(bin_centers, rdf_normalized, linewidth  =2)
plt.axvline(x=bin_centers[min_index], linestyle = 'dashed', color = 'black', label =f'First minimum = {bin_centers[min_index]:.2f} Å')
plt.xlabel('r [Å]', fontsize=2 * fig_width)
plt.ylabel('Amplitude', fontsize=2 * fig_width)
plt.title('RDF Between Oxygen Atoms in Water', fontsize=2 * fig_width)
plt.tight_layout()
plt.legend(fontsize = 1.8*fig_width)
plt.savefig('Figures/Task3_RDF.png')
plt.show()    


