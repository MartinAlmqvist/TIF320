from ase.db import connect

db = connect("gadb.db")
atoms = db.get("id =1").toatoms()