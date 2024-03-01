import numpy as np
from ase.build import molecule
from gpaw import GPAW, PW
from ase import Atoms
from ase.optimize import GPMin

# Define function to get GPAW calculator
def get_calc(element):
    calc = GPAW(xc='PBE',  
                mode=PW(450), 
                kpts=(1, 1, 1),  # gamma-point spinpol=False,  # non-spin-polarized calculation
                spinpol=True,
                txt=f'{element}_Task3.txt')
    return calc

# Define the size of the box
box_size = 12.0  # Å

# Create O2 and CO molecules
O2 = molecule("O2",cell=(box_size, box_size, box_size))
CO = molecule("CO",cell=(box_size, box_size, box_size))
O2.center()
CO.center()

# Get GPAW calculators
calc_O2 = get_calc("O2")
calc_CO = get_calc("CO")

# Set calculators for the molecules
O2.set_calculator(calc_O2)
CO.set_calculator(calc_CO)

# Get energies
E_O2_pot = O2.get_potential_energy()
E_CO_pot = CO.get_potential_energy()
E_O2 = O2.get_total_energy()
E_CO = CO.get_total_energy()

# Write energies to files
with open("Energy_O2.txt", "a") as O2_file,\
     open("Energy_CO.txt", "a") as CO_file:
    O2_file.write(f'Total energy of O2: {E_O2} eV, potential energy: {E_O2_pot} eV\n')
    CO_file.write(f'Total energy of CO: {E_CO} eV, potential energy: {E_CO_pot} eV\n')
