from ase.db import connect
from ase.io import read, write

try:
    db6 = connect("A2_6/gadb.db")
    
    # Attempt to retrieve atoms with ID 74
    entry = db6.get(id=193)
    entry2 = db6.get(id=24)
    
    if entry is not None:
        atoms6 = entry.toatoms() #for 7 atoms
        atom6_24  = entry2.toatoms()
        write("atoms6_minimum.xyz", atoms6)
        write("atom6_24.xyz", atom6_24)
        print("Atoms written to atoms6_minimum.xyz")
    else:
        print("No matching entry found in the database for ID 193 for A6.")
    
except Exception as e:
    print("An error occurred:", e)



try:
    db7 = connect("A2_7/gadb.db")
    
    # Attempt to retrieve atoms with ID 74
    entry = db7.get(id=74)
    entry2 = db7.get(id=50)
    if entry is not None:
        atoms7 = entry.toatoms() #for 7 atoms
        atoms7_50 = entry.toatoms()
        write("atoms7_minimum.xyz", atoms7)
        write("atoms7_50.xyz", atoms7_50)
        print("Atoms written to atoms7_minimum.xyz")
    else:
        print("No matching entry found in the database for ID 74 for A7.")
    
except Exception as e:
    print("An error occurred:", e)


try:
    db8 = connect("A2_8/gadb.db")
    
    # Attempt to retrieve atoms with ID 74
    entry = db8.get(id=119)
    
    if entry is not None:
        atoms8 = entry.toatoms() #for 7 atoms
        write("atoms8_minimum.xyz", atoms8)
        print("Atoms written to atoms8_minimum.xyz")
    else:
        print("No matching entry found in the database for ID 119 for A8.")
    
except Exception as e:
    print("An error occurred:", e)

