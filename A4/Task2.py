# DFT energies of a single atom for each material
gold_atom_energy = -0.158 
platinum_atom_energy = -0.704  
rhodium_atom_energy = -1.218  

#from Task 1
N = 4
gold_min_energy = -3.146 /N 
platinum_min_energy = -5.109/N  
rhodium_min_energy = -6.152  /N

# Calculate cohesive energies for each material
gold_cohesive_energy = gold_min_energy - gold_atom_energy
platinum_cohesive_energy = platinum_min_energy - platinum_atom_energy
rhodium_cohesive_energy = rhodium_min_energy - rhodium_atom_energy

# Output cohesive energies
print("Cohesive energy of Gold:", gold_cohesive_energy, "eV")
print("Cohesive energy of Platinum:", platinum_cohesive_energy, "eV")
print("Cohesive energy of Rhodium:", rhodium_cohesive_energy, "eV")


