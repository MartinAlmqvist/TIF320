import numpy as np
from ase.build import fcc111
from ase.optimize import GPMin
from ase.io import write
from gpaw import GPAW, PW

# https://wiki.fysik.dtu.dk/ase/ase/build/surface.html
lattice_constants = {'Au': 4.167, 'Pt': 4.056, 'Rh': 3.833}

def get_calc(element):
    calc = GPAW(xc='PBE',
                mode=PW(450),
                kpts=(4, 4, 1),
                txt=f'{element}_Task3.txt')
    return calc

elements = ["Au", "Pt", "Rh"]
for element in elements:
    with open(f"total_energy_{element}.txt", "a") as elementfile:
        # Create slab
        surface = fcc111(element, (3, 3, 3), a=lattice_constants[element], vacuum=6.0)
        surface.set_calculator(get_calc(element))
        write(f"{element}_before_structure.xyz", surface)
        
        # Relax slab
        dyn = GPMin(surface, trajectory=f'{element}_slab_relax_ref.traj', logfile=f'{element}_slab_relax_ref.log')
        dyn.run(fmax=0.01, steps=100)
        E_surface = surface.get_potential_energy()
        write(f"{element}_relaxed_structure.xyz", surface)
        
        # Calculate surface energy
        A = surface.get_volume() / (2 * surface.cell[2, 2])  # Surface area of the slab
        surface_energy = E_surface / A  # Surface energy per unit area
        
        # Write surface energy to file
        elementfile.write(f'Surface energy of {element}: {surface_energy:.4f} eV/Angstrom^2\n')
