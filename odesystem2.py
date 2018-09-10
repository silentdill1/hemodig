from enzymes2 import enzymes
from enzymes2 import plas1
from enzymes2 import fal2
from enzymes2 import hdp
from foodvacuole import get_volume
from foodvacuole import get_hb_abundance_change
from foodvacuole import get_volume_change
from datainterpretation import names

import numpy as np

def get_void_array():
    voidarray = []
    for i in range(0, len(names)):
        voidarray.append(0)
    return voidarray

UCF_MM_TO_M = 0.001  # unit conversion from mM to M
timeGrid = np.linspace(2, 47, 101)
initialAbundances = get_void_array()
# [Hb, Fpp, Hz, 20, 22, ... 146, 1, 2 ... 20, 22 ... 70]
# number = peptide length in AAs, second part for fpp free peptides
for indexValuePair in names:
    if indexValuePair[1] == '10wFpp':
        iPwF = indexValuePair[0] # Startindex for peptides with fpp
    elif indexValuePair[1] == '1':
        iPwoF = indexValuePair[0] # Startindex for peptides without fpp





def derivative(abundances, t):
    concentrations = []
    foodVacuoleVolume = get_volume(t)
    for abundance in abundances:
        concentrations.append(abundance/foodVacuoleVolume)
    eci = det_relative_enzyme_concentration_index(t)
    for enzyme in enzymes:



def det_relative_enzyme_concentration_index(t):
    if t % 2 <= 1:
        return int(t/2)-1  # index for protein abundance starts at 0 for 2hpi
    else:
        return int(t/2)


def det_absolute_concentration_of_possible_substrates(enzyme, concentrations):
    acos = 0  # absolute concentration of substrates in M
    if enzyme.indices[0] != 'None':  # operates on species wFpp
        minIndex = enzyme.indices[0][0]
        maxIndex = enzyme.indices[0][1]
        for i in range(minIndex, maxIndex+1):
            acos += concentrations[i]
    if enzyme.indices[1] != 'None':  # operates on species woFpp
        minIndex = enzyme.indices[1][0]
        maxIndex = enzyme.indices[1][1]
        for i in range(minIndex, maxIndex+1):
            acos += concentrations[i]
    return acos


