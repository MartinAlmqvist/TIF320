from ase.io import read
from ase import Atoms
from ase.io.trajectory import Trajectory

# Load the trajectory
traj = Trajectory("Data/cluster24.traj")

# Read the first frame from the trajectory
atoms = traj[0]

# Calculate distances between atoms
atom_distance = atoms.get_all_distances()

print(atom_distance)
