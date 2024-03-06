from ase.build import molecule
from ase.calculators.emt import EMT
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo
from ase.io import read

atoms = read('E_O2_pot_relaxed_structure.xyz')
atoms.calc = EMT()
BFGS(atoms).run(fmax=0.01)

vib = Vibrations(atoms)
vib.run()
vib.summary(log='E_O2_EMT_summary.txt')
vib.write_mode(-1)
vib_energies = vib.get_energies()
potentialenergy = atoms.get_potential_energy()

thermo = IdealGasThermo(vib_energies=vib_energies,
                        potentialenergy=potentialenergy,
                        atoms=atoms,
                        geometry='linear',
                        symmetrynumber=2, spin=0)

Entropy = thermo.get_entropy(temperature=300, pressure=101325.)
with open(f"Energy_O2_task5.txt", "a") as elementfile:
    elementfile.write(str(Entropy) + '\n')

atoms = read('E_CO_pot_relaxed_structure.xyz')
atoms.calc = EMT()
BFGS(atoms).run(fmax=0.01)

vib = Vibrations(atoms)
vib.run()
vib.summary(log='E_CO_EMT_summary.txt')
vib.write_mode(-1)
vib_energies = vib.get_energies()
potentialenergy = atoms.get_potential_energy()

thermo = IdealGasThermo(vib_energies=vib_energies,
                        potentialenergy=potentialenergy,
                        atoms=atoms,
                        geometry='linear',
                        symmetrynumber=2, spin=0)

Entropy = thermo.get_entropy(temperature=300, pressure=101325.)
with open(f"Energy_CO_task5.txt", "a") as elementfile:
    elementfile.write(str(Entropy) + '\n')

