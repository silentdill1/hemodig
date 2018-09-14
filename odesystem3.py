import numpy as np
from datainterpretation2 import names, lengths
from foodvacuole import get_volume, get_hb_abundance_change
from peptides import Peptides, Peptide
from enzymes3 import plas1, fal2, hdp, enzymes


def get_relative_enzyme_concentration_index(t):
    if t % 2 <= 1:
        return int(t/2)-1  # index for protein abundance starts at 0 for 2hpi
    else:
        return int(t/2)


def get_absolute_concentration_of_possible_substrates(enzyme, concentrations):
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


def get_available_enzyme_concentration(concentration_of_substrate, absolute_concentration_of_substrates):
    # returns factor from 1 to 0 determining the ratio of enzyme concentration
    # most likely to be involved in reaction with given substrate
    return concentration_of_substrate / absolute_concentration_of_substrates


def get_void_array():
    void_array = []
    for i in range(0, len(names)):
        void_array.append(0)
    return void_array


def get_fragments(sequence, cut_index):
    list_segment1 = []
    list_segment2 = []
    i = 0
    while i <= cut_index:
        list_segment1.append(sequence[i])
        i += 1
    while i != len(sequence):
        list_segment2.append(sequence[i])
        i += 1
    return [tuple(list_segment1), tuple(list_segment2)]


timeGrid = np.linspace(2, 47, 101)
initialAbundances = get_void_array()
currentPeptideFragments = Peptides(108)
# [Hb, Fpp, Hz, 20, 22, ... 146, 1, 2 ... 20, 22 ... 70]
# number = peptide length in AAs


def derivative(abundances, t):
    concentrations = []
    abundance_changes = get_void_array()
    food_vacuole_volume = get_volume(t)
    for abundance in abundances:
        concentrations.append(abundance/food_vacuole_volume)
    eci = get_relative_enzyme_concentration_index(t)

# hb uptake
    abundance_changes[names['Hb']] += get_hb_abundance_change(t)

# plasmepsin equations
    abundance_changes[names['Hb']] += -plas1.kCatKm * plas1.abundance[eci] * concentrations[0] * food_vacuole_volume
    # TODO: integrate hb beta chain
    abundance_changes[names['108']] += 2 * plas1.kCatKm * plas1.abundance[eci] * concentrations[0] * food_vacuole_volume
    abundance_changes[names['33']] += 2 * plas1.kCatKm * plas1.abundance[eci] * concentrations[0] * food_vacuole_volume
    segments = get_fragments( , 33)
    currentPeptideFragments.peptides.append(Peptide())
    # heme detoxification protein equations
    abundance_changes[names['Fpp']] += -hdp.kCatKm * concentrations[1] * food_vacuole_volume
    # fpp degradation
    abundance_changes[names['Hz']] += hdp.kCatKm * concentrations[1] * food_vacuole_volume
    # hz production

# other enzymes