# https://wiki.fysik.dtu.dk/ase/ase/build/surface.html
#https://wiki.fysik.dtu.dk/ase/tutorials/db/db.html
from gpaw import GPAW, PW
from ase.build import fcc111, add_adsorbate, molecule
from ase.optimize import GPMin
from ase.io import write

lattice_constants = {'Au': 4.08, 'Pt': 3.969, 'Rh': 3.840}

def get_calc(element):
    calc = GPAW(xc='PBE',
                mode=PW(450),
                kpts=(4, 4, 1),
                txt=f'{element}/Task6.txt')
    return calc

surface_elements = ["Au", "Pt", "Rh"]
adsorbate_elements = ["CO","O"]
positions = ["fcc", 'ontop', 'bridge','hcp']

for element in surface_elements:
    with open(f"{element}/total_energy_{element}.txt", "a") as elementfile:
        for adsorbate in adsorbate_elements:
            for position in positions:
                surface = fcc111(element, (3, 3, 3), a=lattice_constants[element], vacuum=6.0)
                calc = get_calc(element)
                surface.set_calculator(calc)
                # Add the adsorbate
                if adsorbate == adsorbate_elements[0]:  # Check if adsorbate is 'CO'
                    adsorbate = molecule(adsorbate)
                add_adsorbate(surface, adsorbate, height=2.0, position=position)
                # Relax Structure
                dyn = GPMin(surface, trajectory=f'{element}/slab_relax_{adsorbate}_{position}.traj',
                            logfile=f'{element}/slab_relax_{adsorbate}_{position}.log')
                dyn.run(fmax=0.1)  # Relax until forces are below 0.1 eV/Å
                E_surface_adsorbate = surface.get_potential_energy()
                write(f"{element}/{adsorbate}_{position}_relaxed_structure.xyz", surface)
                elementfile.write(f'{element} and {position}, {E_surface_adsorbate} \n')
                calc.write(f'{element}_{position}_surface_relaxed.gpaw')
