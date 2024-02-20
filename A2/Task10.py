from ase import Atoms
from gpaw import GPAW, FermiDirac, restart
from ase.optimize import GPMin
from ase.io import write
from ase.db import connect
import numpy as np
from ase.units import Bohr

try:
    db = connect("A2_6/gadb.db")
    entry = db.get(id=193)
    
    # Check if the entry exists
    if entry is not None:
        atoms = entry.toatoms()
        
        # Set up GPAW calculator
        calc = GPAW(nbands=10,  # Number of electronic bands
                    h=0.25,  # Grid spacing [Å]
                    txt='out.txt',
                    occupations=FermiDirac(0.05),
                    setups={'Na': '1'},
                    mode='lcao',
                    basis='dzp')
        atoms.set_calculator(calc)
        
        # Relax the structure
        dyn = GPMin(atoms, trajectory='relaxed_structure.traj', logfile='relax.log')
        dyn.run(fmax=0.02, steps=100)
        write("relaxed_structure.xyz", atoms)
        
        # Save the wavefunction in a .gpw file and then load that state
        calc.write('na_atoms_wavefunction.gpw')
        # atoms, calc = restart('na_atoms_wavefunction.gpw')
        
        # Loop over all wavefunctions and write their cube files
        nbands = calc.get_number_of_bands()
        for band in range(nbands):
            wf = calc.get_pseudo_wave_function(band=band)
            fname = f'na_atoms_wavefunction_{band}.cube'
            print('writing wf', band, 'to file', fname)
            write(fname, atoms, data=wf * Bohr**1.5)
    else:
        print("No matching entry found in the database for ID 193.")
except Exception as e:
    print("An error occurred:", e)
