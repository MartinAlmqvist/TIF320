from ase.db import connect
from ase.io import read, write
from ase import Atoms
from gpaw import GPAW, FermiDirac
from ase.optimize import GPMin
import numpy as np
from ase.units import fs, kB
from ase.io.trajectory import Trajectory
from ase.md.npt import NPT


atoms = read("Task_1_our.xyz")
calc = GPAW(
    mode = 'lcao',
    xc = 'PBE',
    basis = 'dzp',
    symmetry = {'point_group': False }, # Turn off point - group symmetry
    charge = 1, # Charged system
    txt = 'output.gpaw-out' # Redirects calculator output to this file !
    )

atoms.set_calculator(calc)

dyn = NPT( # Some MD method
    atoms,
    temperature_K = 350 ,
    timestep = 0.5* fs , # 200 fs equals 2ps of runtime is 10 steps of the MD simulation is made
    ttime = 20*fs , # Don ’t forget the fs !
    externalstress = 0 , # We don ’t use the barostat , but this needs to be set anyway !
    logfile = 'mdOutput.log',  # Outputs temperature ( and more ) to file at each timestep
    pfactor = None 
    )

trajectory = Trajectory('A3_Task1.traj', 'w' , atoms)
dyn.attach(trajectory.write , interval =1) # Write the current positions etc . to file each timestep

dyn.run(4000) # Run 10 steps of MD simulation

write('A3_Task1_After.xyz', atoms)
calc.write('atoms_wavefunction.gpw')
