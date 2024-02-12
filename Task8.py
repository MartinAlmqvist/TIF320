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
        write("atoms6_minimum.xyz", atoms6_lowest)
        write("atom6_24.xyz", atom6_second_lowest)
        print("Atoms written to atoms6_minimum.xyz")
    else:
        print("No matching entry found in the database for ID 193 for A6.")
    
except Exception as e:
    print("An error occurred:", e)


write("Task8_6atoms_Before.xyz",atoms6_lowest)
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

# Save the electron density to a cube file
#calc.write('electron_density.cube', data=calc.get_all_electron_density())

# Relax
dyn = GPMin(atoms6_lowest, trajectory='relax_ref.traj', logfile='relax_ref.log')
dyn.run(fmax=0.02, steps=100)

# Get the total energy of the relaxed structure
total_energy = atoms6_lowest.get_potential_energy()
write("Task8_6atoms_After.xyz",atoms6_lowest)

# Save the wavefunction in a .gpw file
#calc.write('na_atoms_wavefunction.gpw')

f = open("total_energy.txt", "w")
f.write(f'Total energy: {total_energy} eV')
f.close()    
