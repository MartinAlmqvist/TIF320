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

def Get_RDF(distances, traj):
    
    box_lengths = traj[-1].get_cell_lengths_and_angles()[:3]  # Assuming the box doesn't change over the trajectory
    box_volume = np.prod(box_lengths)

    num_bins = 500
    rdf, bins = np.histogram(distances, bins=num_bins, density=False)
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    dr = bins[1] - bins[0]
    
    num_particles = 25
    density = num_particles/box_volume
    num_frames = len(distances)/num_particles

    rdf_normalized = rdf / (4 * np.pi * bin_centers**2 * dr * density* num_frames)

    return rdf_normalized, bin_centers, density

distances_short = np.loadtxt('distances_task2_short.txt', delimiter='\t')
traj_short = Trajectory("Data/A3_Task1_Mars.traj")

distances_long = np.loadtxt('distances_task2_long.txt', delimiter='\t')
traj_long = Trajectory("Data/NaCluster24.traj")

rdf_short, bins_short, density_short = Get_RDF(distances_short, traj_short)
rdf_long, bins_long, density_long = Get_RDF(distances_long, traj_long)


range_indices = np.where((bins_short >= 2.7) & (bins_short <= 3.5))[0]
min_index_short = range_indices[np.argmin(rdf_short[range_indices])]
rdf_integral_short = np.trapz(rdf_short[:min_index_short],  bins_short[:min_index_short])


range_indices = np.where((bins_long >= 2.7) & (bins_long <= 3.5))[0]
min_index_long = range_indices[np.argmin(rdf_long[range_indices])]
rdf_integral_long = np.trapz(rdf_long[:min_index_long],  bins_long[:min_index_long])

# Print the integral value
print(f"Integral of RDF short: {rdf_integral_short}")
print(f"Integral of RDF long: {rdf_integral_long}")

min_index_short = min_index_long
rdf_integral_short = np.trapz(rdf_short[:min_index_short] * bins_short[:min_index_short]**2, bins_short[:min_index_short])
rdf_integral_long = np.trapz(rdf_long[:min_index_long] * bins_long[:min_index_long]**2, bins_long[:min_index_long])

# Calculate coordination number using the correct formula (integral of RDF up to first minimum)
coordination_number_short = 4 * np.pi * rdf_integral_short*density_short
coordination_number_long = 4 * np.pi * rdf_integral_long*density_long

print(coordination_number_short, coordination_number_long)


# Plotting
plt.plot(bins_short, rdf_short, label ='2 ps run')
plt.plot(bins_long, rdf_long, label = '7 ps run')
plt.axvline(x=bins_long[min_index_long], linestyle = 'dashed', color = 'black', label =f'First minimum = {bins_long[min_index_long]:.2f}')
plt.xlabel('r [Å]', fontsize=2 * fig_width)
plt.ylabel('Amplitude', fontsize=2 * fig_width)
plt.title('RDF Between Na ion and Oxygen Atoms', fontsize=2 * fig_width)
plt.legend(fontsize=1.8 * fig_width)
plt.tight_layout()
plt.savefig('Figures/Task2_RDF.png')
plt.show()    
