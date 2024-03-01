# import numpy as np
# from ase.build import molecule
# from ase.vibrations import Vibrations
# from gpaw import GPAW, PW
# from ase import Atoms
# from ase.optimize import BFGS
# from ase.io import write


# # Define function to get GPAW calculator
# def get_calc(element):
#     calc = GPAW(xc='PBE',  
#                 mode=PW(450), 
#                 kpts=(1, 1, 1),  # gamma-point
#                 spinpol=True,  # non-spin-polarized calculation
#                 txt=f'{element}_Task5.txt')
#     return calc

# # Define the size of the box
# box_size = 12.0  # Å

# # Create O2 and CO molecules
# O2 = molecule("O2", cell=(box_size, box_size, box_size))
# CO = molecule("CO", cell=(box_size, box_size, box_size))
# O2.center()
# CO.center()

# # Get GPAW calculators
# calc_O2 = get_calc("O2")
# calc_CO = get_calc("CO")

# # Set calculators for the molecules
# O2.set_calculator(calc_O2)
# CO.set_calculator(calc_CO)

# # Optimize the structures
# BFGS(O2).run(fmax=0.01)
# BFGS(CO).run(fmax=0.01)

# # Perform vibrational analysis
# vib_O2 = Vibrations(O2)
# vib_CO = Vibrations(CO)
# vib_O2.run()
# vib_CO.run()

# # Write summary files
# vib_O2.summary(log='O2_vibrations_summary.txt')
# vib_CO.summary(log='CO_vibrations_summary.txt')

# # Write vibrational modes to files
# vib_O2.write_mode(-1)
# vib_CO.write_mode(-1)

# # Save the vibrational mode files
# write("O2_modes.xyz", vib_O2)
# write("CO_modes.xyz", vib_CO)
import numpy as np
from ase.build import molecule
from gpaw import GPAW, PW
from ase.vibrations import Vibrations

# Define function to get GPAW calculator
def get_calc(element):
    calc = GPAW(xc='PBE',  
                mode=PW(450), 
                kpts=(1, 1, 1),  # gamma-point
                spinpol=False,  # non-spin-polarized calculation
                txt=f'{element}_Task5.txt')
    return calc

# Create gas-phase molecule
molecule_name = "O2"  # Change to the molecule you want to analyze
gas_molecule = molecule(molecule_name)

# Get GPAW calculator
calc_gas = get_calc(molecule_name)

# Set calculator for the gas-phase molecule
gas_molecule.set_calculator(calc_gas)

# Optimize the gas-phase molecule
gas_molecule.calc.set(fmax=0.01)
dyn = gas_molecule.calc.get_dynamics(fmax=0.01)
dyn.run(0.1)  # Run the dynamics for a short time to relax the structure

# Perform vibrational analysis
vib_gas = Vibrations(gas_molecule)
vib_gas.run()

# Access additional functionalities
# Fold frequencies and intensities within a given range
folded_freq, folded_intensity = vib_gas.fold(start=800.0, end=4000.0, npts=None, width=4.0, type='Gaussian', normalize=False)

# Get vibration energies in eV
energies = vib_gas.get_energies()

# Get vibration frequencies in cm^-1
frequencies = vib_gas.get_frequencies()


# Get specific vibration mode
mode_number = 0  # Change to the mode number you want to retrieve
mode = vib_gas.get_mode(mode_number)


