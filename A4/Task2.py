# DFT energies of a single atom for each material
gold_atom_energy = -0.158 
platinum_atom_energy = -0.704  
rhodium_atom_energy = -1.218  

#from Task 1
N = 1
gold_tot_energy = -3.1322 /N 
platinum_tot_energy = -6.434/N  
rhodium_tot_energy = -7.307  /N

# Calculate cohesive energies for each material
gold_cohesive_energy = -gold_tot_energy + gold_atom_energy
platinum_cohesive_energy = -platinum_tot_energy + platinum_atom_energy
rhodium_cohesive_energy = -rhodium_tot_energy + rhodium_atom_energy

# Output cohesive energies
print("Cohesive energy of Gold:", gold_cohesive_energy, "eV")
print("Cohesive energy of Platinum:", platinum_cohesive_energy, "eV")
print("Cohesive energy of Rhodium:", rhodium_cohesive_energy, "eV")


