from ase import Atoms
from ase.io import read
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.signal import savgol_filter
import seaborn as sns
from ase.io.trajectory import Trajectory

trajectory_file = 'Data/cluster24.traj'
print(len(trajectory_file))
all_distances = [] 
for i in range(8000, 12000):
    atoms = read(trajectory_file, index = i)  # look at fram i of the trajectory

    oxygen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'O']

    # Extract positions of only oxygen atoms
    oxygen_positions = atoms.positions[oxygen_indices]

    # set up a lattice with the oxygens in order to be able to call on the get_all_distances
    oxygen_symbols = ['O'] * len(oxygen_positions)
    oxygen_atoms = Atoms(symbols=oxygen_symbols, positions=oxygen_positions, cell=atoms.cell, pbc=atoms.pbc)

    # Calculate distances between oxygen atoms
    distances_between_oxygen = oxygen_atoms.get_all_distances(mic=True)
    lower_triangle = np.tril(distances_between_oxygen, k=-1)

    # Remove all duplicates
    all_distances.append(lower_triangle.flatten())
    

all_distances = np.array(all_distances)
print(len(all_distances))

nonzero_indices = np.nonzero(all_distances)

# Use the indices to get the non-zero elements
distances = all_distances[nonzero_indices]

print(len(distances))

np.savetxt('distances_task3.txt', distances, delimiter='\t')
