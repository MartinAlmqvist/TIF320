import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms
from gpaw import GPAW , PW
from ase.io import read, write
from ase.optimize import GPMin
from ase.build import bulk

#Part one
x = np.linspace(3.5, 4.5, 10)  # Generate 10 points from 3.5 to 4.5

calc = GPAW(xc='PBE',
            mode=PW(450),
            kpts=(12, 12, 12),
            txt='calculation.txt')

# Open file in append mode ("a") to accumulate results
import numpy as np
from gpaw import GPAW, PW, FermiDirac
from ase.build import bulk

x = np.linspace(3.5, 4.5, 10)  # Generate 10 points from 3.5 to 4.5

calc = GPAW(xc='PBE',
            mode=PW(450),
            kpts=(12, 12, 12),
            txt='calculation.txt')

# Open files in append mode ("a") to accumulate results
with open("total_energy_gold.txt", "a") as au,\
     open("total_energy_pt.txt", "a") as pt,\
     open("total_energy_rh.txt", "a") as rh:
    for indx in x:
        Au_bulk = bulk('Au', 'fcc', indx)
        Au_bulk.set_calculator(calc)
        Pt_bulk = bulk('Pt', 'fcc', indx)
        Pt_bulk.set_calculator(calc)
        Rh_bulk = bulk('Rh', 'fcc', indx)
        Rh_bulk.set_calculator(calc)
        total_energy_au = Au_bulk.get_total_energy()
        total_energy_pt = Pt_bulk.get_total_energy()
        total_energy_rh = Rh_bulk.get_total_energy()
        au.write(f'Total energy: {total_energy_au} eV, A = {indx}\n')
        pt.write(f'Total energy: {total_energy_pt} eV, A = {indx}\n')
        rh.write(f'Total energy: {total_energy_rh} eV, A = {indx}\n')

















# #set upp bulk gold, platinum, and rhodium and compute the lattice parameters


# # Define lattice constant
# a = 4.0

# # # Define lattice vectors for a simple cubic lattice
# # cell = [[a, 0, 0], [0, a, 0], [0, 0, a]]

# # Create a simple cubic lattice with 6 Na atoms
# #positions = [(0, 0, 0), (a, 0, 0), (0, a, 0), (a/2, a/2, a/2), (a/2, a/2, -a/2), (a, a, 0)]
# positions_gold = [(0, 0, 0),(0, 0, a),(0, a, 0),  (0, a, a),(a/2, a/2, a/2)]

# symbols = 'Au5'
# Au_atoms = Atoms(symbols='Au5', positions=positions_gold)
# Au_atoms.center(a)
# write("Gold_Before.xyz",Au_atoms)
# # Set up the structure in GPAW
# Au_atoms.set_calculator(calc)

# # Save the electron density to a cube file
# # calc.write('electron_density.cube', data=calc.get_all_electron_density())

# # Relax
# dyn = GPMin(Au_atoms, trajectory='relax_ref.traj', logfile='relax_ref.log')
# dyn.run(fmax=0.02, steps=100)
# write("Gold_After.xyz",Au_atoms)





#Part two
#n the second part, you will use the energies you have calculated and model the kinetics of the
#reaction over the three different catalysts