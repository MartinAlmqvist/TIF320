from ase import Atoms
from gpaw import GPAW, FermiDirac
from ase.optimize import GPMin

# Create a simple cubic lattice with 6 Na atoms
a = 3.5  # Lattice constant
positions = [(0, 0, 0), (a, 0, 0), (0, a, 0), (0, 0, a), (a, a, 0), (a, 0, a)]
symbols = ['Na'] * 6
na_atoms = Atoms(symbols=symbols, positions=positions)


# Code from ga.py
calc = GPAW(nbands=10,
            h=0.25,
            txt='out.txt',
            occupations=FermiDirac(0.05),
            setups={'Na': '1'},
            mode='lcao',
            basis='dzp')


# set up the structure in GPAW
na_atoms.set_calculator(calc)

# Relax
dyn = GPMin(na_atoms, trajectory='relax_ref.traj', logfile='relax_ref.log')
dyn.run(fmax=0.02, steps=100)

# Get the total energy of the relaxed structure
total_energy = na_atoms.get_potential_energy()

# Save the wavefunction in a .gpw file
calc.write('na_atoms_wavefunction.gpw')

with open('total_energy.txt', 'w') as f:
    f.write(f'Total energy: {total_energy} eV')
