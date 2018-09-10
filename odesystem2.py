from enzymes2 import enzymes
from enzymes2 import plas1
from enzymes2 import fal2
from enzymes2 import hdp
from foodvacuole import get_volume
from foodvacuole import get_hb_abundance_change
from datainterpretation import names
from datainterpretation import lengths

import numpy as np


def get_void_array():
    void_array = []
    for i in range(0, len(names)):
        void_array.append(0)
    return void_array


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


def get_number_of_cuts(enzyme, peptide_length):
    number_of_cuts = 0
    current_length = peptide_length
    while (current_length - enzyme.stepSize) >= 0:
        current_length -= enzyme.stepSize
        number_of_cuts += 1
    return number_of_cuts


def get_available_enzyme_concentration(concentration_of_substrate, absolute_concentration_of_substrates):
    # returns factor from 1 to 0 determining the ratio of enzyme concentration
    # most likely to be involved in reaction with given substrate
    return concentration_of_substrate / absolute_concentration_of_substrates


timeGrid = np.linspace(2, 47, 101)
initialAbundances = get_void_array()
# [Hb, Fpp, Hz, 20, 22, ... 146, 1, 2 ... 20, 22 ... 70]
# number = peptide length in AAs, second part for peptides without Ferriprotoporphyrin IX (fpp)


def derivative(abundances, t):
    concentrations = []
    abundance_changes = get_void_array()
    food_vacuole_volume = get_volume(t)
    for abundance in abundances:
        concentrations.append(abundance/food_vacuole_volume)
    eci = get_relative_enzyme_concentration_index(t)

# hb uptake
    abundance_changes[names['Hb']] += get_hb_abundance_change(t)

# plasmepsin I equations
    abundance_changes[names['Hb']] += -plas1.kCatKm * plas1.abundance[eci] * concentrations[0] * food_vacuole_volume
    abundance_changes[names['146wFpp']] += 2 * plas1.kCatKm * plas1.abundance[eci] * concentrations[0] * food_vacuole_volume
    # abundance change of peptide chain with length of 146 AS containing fpp
    abundance_changes[names['108wFpp']] += 2 * plas1.kCatKm * plas1.abundance[eci] * concentrations[0] * food_vacuole_volume
    # TODO: eigentlich 33 AS
    abundance_changes[names['32']] += 2 * plas1.kCatKm * plas1.abundance[eci] * concentrations[0] * food_vacuole_volume

# heme detoxification protein equations
    abundance_changes[names['Fpp']] += -hdp.kCatKm * concentrations[1] * food_vacuole_volume
    # fpp degradation
    abundance_changes[names['Hz']] += hdp.kCatKm * concentrations[1] * food_vacuole_volume
    # hz production

# other enzymes
    for enzyme in enzymes:
        absolute_concentration_of_possible_substrates = get_absolute_concentration_of_possible_substrates(enzyme, concentrations)
        if enzyme.indices[0] != 'None':  # operates on sequences with fpp
            start_index = enzyme.indices[0][0]
            end_index = enzyme.indices[0][1]
            max_enzyme_abundance = enzyme.MAX_ENZYME_ABUNDANCE
            k_cat_k_m = enzyme.kCatKm
            for i in range(start_index, end_index+1):
                if concentrations[i] != 0:  # check if substrate has concentration
                    p = get_available_enzyme_concentration(concentrations[i],
                                                           absolute_concentration_of_possible_substrates)
                    # TODO: get abundance data according to eci if possible
                    abundance_changes[i] += -p * max_enzyme_abundance * k_cat_k_m * concentrations[i]
                    # decay of substrate
                    if enzyme is fal2:
                        abundance_changes[1] += p * max_enzyme_abundance * k_cat_k_m * concentrations[i]
                        # falcipain II releases fpp from substrate
                    length_of_substrate = lengths[i]
                    n = get_number_of_cuts(enzyme, length_of_substrate)
                    length_of_fragment1 = length_of_substrate
                    # length of one of the products of current cut
                    if enzyme is fal2:  # the fpp is no longer part of any peptide, both are added to non fpp array
                        while (length_of_fragment1 - enzyme.stepSize) > 0:
                            length_of_fragment1 -= enzyme.stepSize
                            length_of_fragment2 = length_of_substrate - length_of_fragment1
                            index_of_fragment1 = names[str(length_of_fragment1)]
                            index_of_fragment2 = names[str(length_of_fragment2)]
                            abundance_change_of_fragments = p / n * max_enzyme_abundance * k_cat_k_m * concentrations[i]
                            abundance_changes[index_of_fragment1] += abundance_change_of_fragments
                            abundance_changes[index_of_fragment2] += abundance_change_of_fragments
                    else:
                        while (length_of_fragment1 - enzyme.stepSize) > 0:  # longer fragment maintains fpp (if-part)
                            length_of_fragment1 -= enzyme.stepSize
                            length_of_fragment2 = length_of_substrate - length_of_fragment1
                            if length_of_fragment1 >= length_of_fragment2:
                                index_of_fragment1 = names[str(length_of_fragment1)+'wFpp']
                                index_of_fragment2 = names[str(length_of_fragment2)]
                            else:
                                index_of_fragment1 = names[str(length_of_fragment1)]
                                index_of_fragment2 = names[str(length_of_fragment2)+'wFpp']
                            abundance_change_of_fragments = p / n * max_enzyme_abundance * k_cat_k_m * concentrations[i]
                            abundance_changes[index_of_fragment1] += abundance_change_of_fragments
                            abundance_changes[index_of_fragment2] += abundance_change_of_fragments
        if enzyme.indices[1] != 'None':  # operates on sequences without fpp
            start_index = enzyme.indices[1][0]
            end_index = enzyme.indices[1][1]
            for i in range(start_index, end_index+1):
                if concentrations[i] != 0:  # check if substrate has concentration
                    max_enzyme_abundance = enzyme.MAX_ENZYME_ABUNDANCE
                    k_cat_k_m = enzyme.kCatKm
                    p = get_available_enzyme_concentration(concentrations[i],
                                                           absolute_concentration_of_possible_substrates)
                    # TODO: get abundance data according to eci if possible
                    abundance_changes[i] += -p * max_enzyme_abundance * k_cat_k_m * concentrations[i]
                    # decay of substrate
                    length_of_substrate = lengths[i]
                    if enzyme.endBased:
                        n = 1  # end based enzymes cut only once at N terminus of protein
                    else:
                        n = get_number_of_cuts(enzyme, length_of_substrate)
                    length_of_fragment1 = length_of_substrate
                    # length of one of the products of current cut
                    while (length_of_fragment1 - enzyme.stepSize) > 0:
                        length_of_fragment1 -= enzyme.stepSize
                        length_of_fragment2 = length_of_substrate - length_of_fragment1
                        index_of_fragment1 = names[str(length_of_fragment1)]
                        index_of_fragment2 = names[str(length_of_fragment2)]
                        abundance_change_of_fragments = p / n * max_enzyme_abundance * k_cat_k_m * concentrations[i]
                        abundance_changes[index_of_fragment1] += abundance_change_of_fragments
                        abundance_changes[index_of_fragment2] += abundance_change_of_fragments
                        if enzyme.endBased:
                            break  # simulates only one cut for end based enzymes

    return abundance_changes










