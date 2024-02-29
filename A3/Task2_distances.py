from ase import Atoms
from ase.io import read
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.signal import savgol_filter
import seaborn as sns
from ase.io.trajectory import Trajectory

trajectory_file = 'Data/A3_Task1_Mars.traj'
#trajectory_file = 'Data/NaCluster24.traj'
traj = Trajectory(trajectory_file)

unique_elements = set()

""" for frame in traj:
   # Get the chemical symbols for atoms in the frame
   symbols = frame.get_chemical_symbols()
   
   # Add unique symbols to the set
   unique_elements.update(symbols)

unique_elements = sorted(unique_elements)
print("Unique elements in the trajectory:", unique_elements)
element_1 = unique_elements[0]
element_2 = unique_elements[1]
print(element_1, element_2)
 """

all_distances = []
total_frames = len(traj)
start_index = int(0 * total_frames)
print(total_frames)
all_distances = []

for i in range(start_index, total_frames):
    frame = traj[i]
    atomic_numbers = frame.get_atomic_numbers()
    oxygen_indices = np.where(atomic_numbers == 8)[0]
    Na_index = np.where(atomic_numbers == 11)[0]

    # Extract positions of only oxygen atoms
    distances = frame.get_distances(Na_index, oxygen_indices, mic = True)
    all_distances.append(distances)

print(np.shape(all_distances)) 

all_distances = np.array(all_distances).flatten()
np.savetxt('distances_task2_short.txt', all_distances, delimiter='\t')