import numpy as np
import matplotlib.pyplot as plt
from gpaw import GPAW , PW

#Part one
calc = GPAW ( xc = ' PBE ',
    mode = PW (450),
    kpts =(12 , 12 , 12), #Suitable cutoff energy and number of k-points for the bulks are 450 eV and (12, 12, 12), respectively. 
    txt = ' calculation.txt ')









#Part two
#n the second part, you will use the energies you have calculated and model the kinetics of the
#reaction over the three different catalysts