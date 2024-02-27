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
    plt.rcParams['font.size'] = size_factor * fig_width
    plt.rcParams['xtick.labelsize'] = size_factor * fig_height
    plt.rcParams['ytick.labelsize'] = size_factor * fig_height
    plt.grid(True)

# Set the plot style
set_plot_style()
fig_width = 6

distances = np.loadtxt('distances_task2.txt', delimiter='\t')
traj = Trajectory("Data/A3_Task1_Mars.traj")

box_lengths = traj[-1].get_cell_lengths_and_angles()[:3]  # Assuming the box doesn't change over the trajectory
num_bins = 100
rdf, bins = np.histogram(distances, bins=num_bins, density=True, range=(2,6))
bin_centers = 0.5 * (bins[1:] + bins[:-1])
dr = bins[1] - bins[0]
box_volume = np.prod(box_lengths)
density = 25/box_volume

rdf_normalized = rdf / (4 * np.pi * bin_centers**2 * (dr) *density )


bin_centers = 0.5 * (bins[1:] + bin[int(len(bins)*(6-2)*3.1)])

rdf_integral = np.trapz(rdf_normalized, bin_centers)

# Print the integral value
print(f"Integral of RDF: {rdf_integral}")


# Plotting
plt.plot(bin_centers, rdf_normalized)
plt.xlabel('Distance')
plt.ylabel('Amplitude')
plt.title('RDF Between Oxygen Atoms in Water')
plt.show()    
