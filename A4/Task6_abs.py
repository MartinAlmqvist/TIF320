
#from Task 1
gold_tot_energy = -3.1322
platinum_tot_energy = -6.434
rhodium_tot_energy = -7.307

#Task4
gold_atom_energy = -14.1942 
platinum_atom_energy = -8.7277  
rhodium_atom_energy = -1.218  

#Task6
Au_O_fcc= -83.27843459687264 
Au_O_ontop= -81.9046172313769 
Au_O_bridge= -83.27555078273981 
Au_O_hcp= -81.9046172313769 
positions = ["fcc", 'ontop', 'bridge','hcp']

def Eads(Au_fcc):
    return Au_fcc - gold_atom_energy- gold_tot_energy


def get_lowest(Eads_list):
    value = min(Eads_list)
    name_ind = Eads_list.index(value)
    name = positions[name_ind]
    return value, name


Eads_au_O_fcc = Eads(Au_O_fcc)
Eads_au_O_ontop = Eads(Au_O_ontop)
Eads_au_O_bridge = Eads(Au_O_bridge)
Eads_au_O_hcp = Eads(Au_O_hcp)

Eads_au_O_list = [Eads_au_O_fcc,Eads_au_O_ontop,Eads_au_O_bridge,Eads_au_O_hcp]
E_ads_au_O_opt, E_ads_O_struct = get_lowest(Eads_au_O_list)
E_ads_au_CO_opt, E_ads_CO_struct = get_lowest(Eads_au_O_list)

print(E_ads_au_O_opt,E_ads_O_struct)

E_ads_CO_opt = 0

# E_ads_O_opt = Eads_au_fcc
# E_ads_CO_opt = E_ads_CO_fcc
E_a = -0.3*(E_ads_au_O_opt + E_ads_CO_opt) + 0.22
print(f'E_a for atom O with structure {E_ads_O_struct} and CO {E_ads_CO_struct} CO : {E_a}')

# print(f'Ea = {E_a}')