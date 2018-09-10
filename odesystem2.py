from enzymes2 import enzymes
from enzymes2 import plas1
from enzymes2 import fal2
from enzymes2 import hdp
from foodvacuole import get_volume
from foodvacuole import get_hb_abundance_change
from datainterpretation import names

import numpy as np

def get_void_array():
    void_array = []
    for i in range(0, len(names)):
        void_array.append(0)
    return void_array

timeGrid = np.linspace(2, 47, 101)
initialAbundances = get_void_array()
# [Hb, Fpp, Hz, 20, 22, ... 146, 1, 2 ... 20, 22 ... 70]
# number = peptide length in AAs, second part for peptides without Ferriprotoporphyrin IX


def derivative(abundances, t):
    concentrations = []
    abundance_changes = get_void_array()
    food_vacuole_volume = get_volume(t)
    for abundance in abundances:
        concentrations.append(abundance/food_vacuole_volume)
    eci = det_relative_enzyme_concentration_index(t)

# hb uptake
    abundance_changes[0] = get_hb_abundance_change(t)

# plasmepsin I equations
    abundance_changes[0] = -plas1.kCatKm * plas1.abundance[eci] * concentrations[0] * food_vacuole_volume # hb
    abundance_changes[names['146wFpp']] = 2 * plas1.kCatKm * plas1.abundance[eci] * concentrations[0] * food_vacuole_volume
    # abundance change of peptide chain with length of 146 AS containing fpp
    abundance_changes[names['108wFpp']] = 2 * plas1.kCatKm * plas1.abundance[eci] * concentrations[0] * food_vacuole_volume
    abundance_changes[names['33']] = 2 * plas1.kCatKm * plas1.abundance[eci] * concentrations[0] * food_vacuole_volume

# falcipain II equations

# heme detoxification protein equations
    abundance_changes[1] = -hdp.kCatKm * hdp.abundance[eci] * concentrations[1] * food_vacuole_volume # fpp degradation
    abundance_changes[2] = hdp.kCatKm * hdp.abundance[eci] * concentrations[1] * food_vacuole_volume # hz production

# other enzymes
    for enzyme in enzymes:



def det_relative_enzyme_concentration_index(t):
    if t % 2 <= 1:
        return int(t/2)-1  # index for protein abundance starts at 0 for 2hpi
    else:
        return int(t/2)


def det_absolute_concentration_of_possible_substrates(enzyme, concentrations):
    acos = 0  # absolute concentration of substrates in M
    if enzyme.indices[0] != 'None':  # operates on species wFpp
        min_index = enzyme.indices[0][0]
        max_index = enzyme.indices[0][1]
        for i in range(min_index, max_index+1):
            acos += concentrations[i]
    if enzyme.indices[1] != 'None':  # operates on species woFpp
        min_index = enzyme.indices[1][0]
        max_index = enzyme.indices[1][1]
        for i in range(min_index, max_index+1):
            acos += concentrations[i]
    return acos


