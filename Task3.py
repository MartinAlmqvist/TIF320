from ase import Atoms
from gpaw import GPAW, FermiDirac
from ase.optimize import GPMin
from ase.io import read, write
import numpy as np

# Define lattice constant
a = 4.5

# Define lattice vectors for a simple cubic lattice
# cell = [[a, 0, 0], [0, a, 0], [0, 0, a]]

# Create a simple cubic lattice with 6 Na atoms
positions = [(0, 0, 0), (a, 0, 0), (0, a, 0), (0, 0, a), (a, a, 0), (a, 0, a)]


symbols = 'Na6'
na_atoms = Atoms(symbols='Na6', positions=positions)
na_atoms.center(5)

write("Before.xyz",na_atoms)
# Code from ga.py
calc = GPAW(nbands=10, #Number of electronic bands
            h=0.25, #Grid spacing [Å]
            txt='out.txt',
            occupations=FermiDirac(0.05),
            setups={'Na': '1'},
            mode='lcao',
            basis='dzp')

# Set up the structure in GPAW
na_atoms.set_calculator(calc)

# Relax
dyn = GPMin(na_atoms, trajectory='relax_ref.traj', logfile='relax_ref.log')
dyn.run(fmax=0.02, steps=100)

# Get the total energy of the relaxed structure
total_energy = na_atoms.get_potential_energy()
write("After.xyz",na_atoms)

# Save the wavefunction in a .gpw file
calc.write('na_atoms_wavefunction.gpw')

# Save the total energy in a .txt file
#with open('total_energy.txt', 'w') as f:
#    f.write(f'Total energy: {total_energy} eV') 

f = open("total_energy.txt", "w")
f.write(f'Total energy: {total_energy} eV')
f.close()
