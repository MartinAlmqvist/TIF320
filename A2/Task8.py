# scp Task8.py lovakes@vera1.c3se.chalmers.se:/cephyr/users/lovakes/Vera/tif320-computational-materials-and-molecular-physics/Na-clusters-GA-search/TIF320/Task8
from ase.db import connect
from ase.io import write
from ase import Atoms
from gpaw import GPAW, FermiDirac, PW
from ase.optimize import GPMin
import numpy as np

directory = "Results/5_2/"

try:
    db6 = connect("A2_6/gadb.db")
    
    entry = db6.get(id=193)
    entry2 = db6.get(id=50)
    
    if entry and entry2 is not None:
        atoms6_lowest = entry.toatoms()
        atom6_second_lowest  = entry2.toatoms() 
    else:
        print("No matching entry found in the database for ID 193 and 50 for A6.")
    
except Exception as e:
    print("An error occurred:", e)

# calculate initial energy before relaxation
initial_energy_lowest = atoms6_lowest.get_total_energy()
initial_energy_second_lowest = atom6_second_lowest.get_total_energy()

# Code from ga.py
calc = GPAW(nbands=10,  # Number of electronic bands
            h=0.25,  # Grid spacing [Å]
            txt='out.txt',
            occupations=FermiDirac(0.05),
            setups={'Na': '1'},
            mode='pw',
            #basis='tzp',
            xc='PBE',
            )

# Set up the structure in GPAW
atoms6_lowest.set_calculator(calc)
atom6_second_lowest.set_calculator(calc)

# Relax
dyn_lowest = GPMin(atoms6_lowest, trajectory=f'{directory}relax_ref_lowest.traj', logfile=f'{directory}relax_ref_lowest.log')
dyn_lowest.run(fmax=0.01, steps=100)

dyn_second_lowest = GPMin(atom6_second_lowest, trajectory=f'{directory}relax_ref_second_lowest.traj', logfile=f'{directory}relax_ref_second_lowest.log')
dyn_second_lowest.run(fmax=0.01, steps=100)

# Get the total energy of the relaxed structure
final_energy_lowest = atoms6_lowest.get_total_energy()
final_energy_second_lowest = atom6_second_lowest.get_total_energy()

# save the energy difference
diff_lowest = initial_energy_lowest - final_energy_lowest
diff_second_lowest = initial_energy_second_lowest - final_energy_second_lowest

# save the relaxed lattice to an xyz file
write(f'{directory}Task8_6atoms_lowest_After.xyz', atoms6_lowest)
write(f'{directory}Task8_6atoms_second_lowest_After.xyz', atom6_second_lowest)

with open(f'{directory}Task_8_5_2.txt', "w") as f:
    f.write(f'Initial energy lowest: {initial_energy_lowest} eV\n')
    f.write(f'Final energy lowest: {final_energy_lowest} eV \n')
    f.write(f'Energy diff: {diff_lowest} eV\n \n')
    f.write(f'Initial energy second lowest: {initial_energy_second_lowest} eV\n')
    f.write(f'Final energy second lowest: {final_energy_second_lowest} eV\n')
    f.write(f'Energy diff: {diff_second_lowest} eV\n \n')
