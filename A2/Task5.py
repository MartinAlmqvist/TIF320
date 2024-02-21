from ase.db import connect
from ase.io import read, write

try:
    db6 = connect("A2_6/gadb.db")
    
    # Attempt to retrieve atoms with ID 74
    entry = db6.get(id=193)
    entry2 = db6.get(id=50)
    
    if entry and entry2 is not None:
        atoms_lowest = entry.toatoms() #for 7 atoms
        atom6_seclowest  = entry2.toatoms()
        write("atoms6_minimum.xyz", atoms_lowest)
        write("atom6_second_min.xyz", atom6_seclowest)
        print("Atoms written to atoms6_minimum.xyz")
    else:
        print("No matching entry found in the database for ID 193 for A6.")
    
except Exception as e:
    print("An error occurred:", e)

try:
    db7 = connect("A2_7/gadb.db")
    
    # Attempt to retrieve atoms with ID 74
    entry = db7.get(id=74)
    entry2 = db7.get(id=98)
    if entry is not None:
        atoms7_lowest = entry.toatoms() #for 7 atoms
        atoms7_seclowest = entry2.toatoms()
        write("atoms7_minimum.xyz", atoms7_lowest)
        write("atoms7_second_min.xyz", atoms7_seclowest)
        print("Atoms written to atoms7_minimum.xyz")
    else:
        print("No matching entry found in the database for ID 74 for A7.")
    
except Exception as e:
    print("An error occurred:", e)

try:
    db8 = connect("A2_8/gadb.db")
    
    # Attempt to retrieve atoms with ID 74
    entry_lowest = db8.get(id=119)
    entry_second_lowest = db8.get(id=34)
    
    if entry_lowest and entry_second_lowest is not None:
        atoms8 = entry_lowest.toatoms() 
        atoms8_2 = entry_second_lowest.toatoms() 
        write("atoms8_minimum.xyz", atoms8)
        write("atoms8_second_minimum.xyz", atoms8_2)
        print("Atoms written to atoms8_minimum.xyz")
    else:
        print("No matching entry found in the database for ID 119 for A8.")
    
except Exception as e:
    print("An error occurred:", e)

