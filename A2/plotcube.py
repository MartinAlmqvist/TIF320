
from ase import Atoms
from ase.optimize import GPMin
from ase.io import write
from ase.db import connect
import numpy as np
from ase.units import Bohr


import vesta

# def visualize_with_vesta(cube_files):
#     # Create a VESTA session
#     session = vesta.Session()
    
#     # Open each cube file in VESTA
#     for cube_file in cube_files:
#         session.load_file(cube_file)

# List of cube files to visualize
cube_files = ['A2/A2_6/Cubefiles/na_atoms_wavefunction_0.cube', 'A2/A2_6/Cubefiles/na_atoms_wavefunction_1.cube']  # Add more cube files as needed
local_client = vesta.LocalClient()
# vesta.show(cube_files[0])

# vesta.

# Visualize cube files in VESTA
#visualize_with_vesta(cube_files)
