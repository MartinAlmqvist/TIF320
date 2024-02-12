from ase.db import connect
from ase.io import read, write
from ase import Atoms
from gpaw import GPAW, FermiDirac
from ase.optimize import GPMin
from ase.io import read, write
import numpy as np

try:
    db6 = connect("A2_6/gadb.db")
    
    # Attempt to retrieve atoms with ID 74
    entry = db6.get(id=193)
    entry2 = db6.get(id=24)
    
    if entry is not None:
        atoms6_lowest = entry.toatoms() #id 193
        atom6_second_lowest  = entry2.toatoms() # id 24
    else:
        print("No matching entry found in the database for ID 193 for A6.")
    
except Exception as e:
    print("An error occurred:", e)

# save the lattice to an xyz file
#write("Task8_6atoms_lowest_Before.xyz",atoms6_lowest)
#write("Task8_6atoms_second_lowest_Before.xyz",atoms6_lowest)

#calculate initial energy before relaxation
initial_energy_lowest = atoms6_lowest.get_potential_energy()
initial_energy_second_lowest = atom6_second_lowest.get_potential_energy()

# Code from ga.py
calc = GPAW(nbands=10, #Number of electronic bands
            h=0.25, #Grid spacing [Å]
            txt='out.txt',
            occupations=FermiDirac(0.05),
            setups={'Na': '1'},
            mode='lcao',
            basis='dzp')

# Set up the structure in GPAW
atoms6_lowest.set_calculator(calc)
atom6_second_lowest.set_calculator(calc)

# Relax
dyn_lowest = GPMin(atoms6_lowest, trajectory='relax_ref_lowest.traj', logfile='relax_ref_lowest.log')
dyn_lowest.run(fmax=0.02, steps=100)

dyn_second_lowest = GPMin(atom6_second_lowest, trajectory='relax_ref_second_lowest.traj', logfile='relax_ref_second_lowest.log')
dyn_second_lowest.run(fmax=0.02, steps=100)

# Get the total energy of the relaxed structure
final_energy_lowest = atoms6_lowest.get_potential_energy()
final_energy_second_lowest = atom6_second_lowest.get_potential_energy()

# save the relaxed lattice to an xyz file
#write("Task8_6atoms_lowest_After.xyz",atoms6_lowest)
#write("Task8_6atoms_second_lowest_After.xyz",atom6_second_lowest)

# Save the wavefunction in a .gpw file
#calc.write('na_atoms_wavefunction.gpw')

with open("energy_before_after.txt", "w") as f:
    f.write(f'Initial energy lowest: {initial_energy_lowest} eV\n')
    f.write(f'Final energy lowest: {final_energy_lowest} eV\n')
    f.write(f'Initial energy second lowest: {final_energy_lowest} eV\n')
    f.write(f'Final energy second lowest: {final_energy_second_lowest} eV\n')
