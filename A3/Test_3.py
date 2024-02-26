import numpy as np
import matplotlib.pyplot as plt
from ase.io.trajectory import Trajectory
from itertools import combinations

# Load the trajectory
traj = Trajectory("A3/Data/cluster24.traj")
unique_elements = set()

for frame in traj:
    # Get the chemical symbols for atoms in the frame
    symbols = frame.get_chemical_symbols()
    
    # Add unique symbols to the set
    unique_elements.update(symbols)

unique_elements = sorted(unique_elements)
print("Unique elements in the trajectory:", unique_elements)

# Specify the two elements you want to compare
element_1 = unique_elements[0]
element_2 = unique_elements[1]

# Extract indices of atoms of the two specified elements
indices_H = [i for i, atom in enumerate(traj[0]) if atom.symbol == element_1]
indices_O = [i for i, atom in enumerate(traj[0]) if atom.symbol == element_2]

# Initialize list to store distances
distances = []

# Iterate through frames in the trajectory
for frame in traj:
    # Get positions of atoms in the frame
    positions = frame.get_positions()

    # Extract positions of atoms belonging to the specified elements
    positions_H = positions[indices_H]
    positions_O = positions[indices_O]

    # Calculate distances between unique pairs of atoms
    for pos1, pos2 in combinations(positions_O, 2):
        distance_frame = np.linalg.norm(pos1 - pos2, axis=-1, keepdims=True)
        distances.append(distance_frame)

# Flatten distances list
distances = np.array(distances)

# Compute the volume of the system
# Assuming the system is cubic, you can calculate the volume from the box lengths
box_lengths = traj[-1].get_cell_lengths_and_angles()[:3]  # Assuming the box doesn't change over the trajectory
volume = np.prod(box_lengths)

density = 24/volume 



# Calculate the radial distribution function
num_bins = 100
rdf, bins = np.histogram(distances, bins=num_bins, density=True, range=(2,6))
bin_centers = 0.5 * (bins[1:] + bins[:-1])

dr = bins[1] - bins[0]
r_max = 6
r_vec = np.linspace(0.001,r_max,num_bins)
dominator = 4/3*np.pi * r_vec**3* density

# Normalize RDF by volume
rdf /= dominator

# Plotting
plt.plot(bin_centers, rdf)
plt.xlabel('Distance')
plt.ylabel('Radial Distribution Function')
plt.title(f'Radial Distribution Function between {element_2} atoms')
plt.savefig("loooooooooooool.png")


# # Convert distances to numpy array
# distances = np.array(distances)

# # Define bins for RDF calculation
# r_min = 0.0
# r_max = np.max(distances)
# number_of_bins = 50
# bin_edges = np.linspace(r_min, r_max, number_of_bins + 1)

# # Initialize RDF array
# rdf = np.zeros(number_of_bins)

# # Iterate through bins
# for i in range(number_of_bins):
#     bin_start = bin_edges[i]
#     bin_end = bin_edges[i + 1]

#     # Count distances falling within the current bin
#     count = np.sum((distances >= bin_start) & (distances < bin_end))
#     rdf[i] = count

# # Normalize RDF
# bin_volume = (4 / 3) * np.pi * (bin_edges[1:]**3 - bin_edges[:-1]**3)
# rdf /= bin_volume

# # Plot the RDF
# plt.plot(bin_edges[:-1], rdf, label='RDF')
# plt.xlabel('Distance')
# plt.ylabel('RDF')
# plt.title('Radial Distribution Function (RDF) between {} and {}'.format(element_1, element_2))
# plt.legend()
# plt.savefig("Task2_RDF.png")
# plt.show()
