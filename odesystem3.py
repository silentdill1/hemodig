import numpy as np
from datainterpretation2 import names, lengths
from foodvacuole import get_volume, get_hb_abundance_change
from peptides import Peptides, Peptide
from enzymes3 import plas1, hdp, fal2, enzymes
from dataimport import alphaHbChain, betaHbChain
from copy import deepcopy


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


def get_void_array():
    void_array = []
    for i in range(0, len(names)):
        void_array.append(0)
    return void_array


timeGrid = np.linspace(2, 47, 101)
initialAbundances = get_void_array()
currentPeptideFragments = Peptides(108)
# [Hb, Fpp, Hz, 1, ... 108]
# number = peptide length in AAs


def get_fragments(sequence, cut_index):
    """ returns list of tuples of AAs in order fragment before cut index, fragment after cut index """
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


def check_peptide_existence(sequence, length, contains_fpp, new_peptide_fragments):
    """ returns peptide index in currentPeptideFragments.peptides[lengthIndex] (lengthIndex = length-1)
        or -1 if not found"""
    index = -1
    for peptideIndex in range(0, len(new_peptide_fragments.peptides[length-1])):
        peptide = new_peptide_fragments.peptides[length-1][peptideIndex]
        if peptide.sequence == sequence and peptide.containsFpp == contains_fpp:
            index = peptideIndex
    return index


def longer_fragment_tuple(length_of_fragment1, length_of_fragment2):
    if length_of_fragment1 >= length_of_fragment2:
        answers = (True, False)
        return answers
    else:
        answers = (False, True)
        return answers


def add_peptide_fragments(peptide, number_of_peptides_for_length, total_abundance_change, new_peptide_fragments):
    """

    :param peptide: peptide, peptide to be cleaved (cleavage sites have been set previously)
    :param number_of_peptides_for_length: postive integer, number of different peptides for length
    :param total_abundance_change: postive integer, change in abundance for all peptides of current length
    :return: void, modifies currentPeptideFragments
    """
    peptide.abundance -= total_abundance_change / number_of_peptides_for_length
    for cleavageSite in peptide.cleavageSites:
        LiN_for_current_cleavage_site = cleavageSite[1]
        segments = get_fragments(peptide.sequence, cleavageSite[0])
        length_of_fragment1 = cleavageSite[0]+1
        length_of_fragment2 = peptide.length - length_of_fragment1
        abundance_change_of_fragments = LiN_for_current_cleavage_site / peptide.sumOfLiNsForAllCleavageSites * total_abundance_change

        if peptide.containsFpp:
            fpp_distribution = longer_fragment_tuple(length_of_fragment1, length_of_fragment2)
        else:
            fpp_distribution = (False, False)

        fragment1_index_in_peptides_of_length = check_peptide_existence(segments[0], length_of_fragment1, new_peptide_fragments)
        if fragment1_index_in_peptides_of_length == -1:  # create new peptide
            fragment1 = Peptide(segments[0], length_of_fragment1, abundance_change_of_fragments, fpp_distribution[0])
            new_peptide_fragments.peptides[length_of_fragment1].append(fragment1)
        else:  # add abundance change to old peptide
            new_peptide_fragments.peptides[length_of_fragment1][fragment1_index_in_peptides_of_length].abundance += abundance_change_of_fragments

        fragment2_index_in_peptides_of_length = check_peptide_existence(segments[1], length_of_fragment2, new_peptide_fragments)
        if fragment2_index_in_peptides_of_length == -1:
            fragment2 = Peptide(segments[1], length_of_fragment2, abundance_change_of_fragments, fpp_distribution[1])
            new_peptide_fragments.peptides[length_of_fragment1].append(fragment2)
        else:
            new_peptide_fragments.peptides[length_of_fragment2][fragment2_index_in_peptides_of_length].abundance += abundance_change_of_fragments


def derivative(abundances, t):
    concentrations = []
    abundance_changes = get_void_array()
    new_peptide_fragments = deepcopy(currentPeptideFragments)
    food_vacuole_volume = get_volume(t)
    for abundance in abundances:
        concentrations.append(abundance/food_vacuole_volume)
    eci = get_relative_enzyme_concentration_index(t)

# hb uptake
    abundance_changes[names['Hb']] = get_hb_abundance_change(t)

# plasmepsin equations
    abundance_changes[names['Hb']] += -plas1.kCatKm * plas1.abundance[eci] * concentrations[0] * food_vacuole_volume
    # TODO: integrate hb beta chain
    abundance_changes[names['108']] = 2 * plas1.kCatKm * plas1.abundance[eci] * concentrations[0] * food_vacuole_volume
    abundance_changes[names['33']] = 2 * plas1.kCatKm * plas1.abundance[eci] * concentrations[0] * food_vacuole_volume
    segments = get_fragments(alphaHbChain, 32)

    index1 = check_peptide_existence(segments[0], 33, False, new_peptide_fragments)
    if index1 == -1:
        peptide33 = Peptide(segments[0], 33, abundance_changes[names['33']], False)
        new_peptide_fragments.peptides[33].append(peptide33)
    else:
        new_peptide_fragments.peptides[33][index1].abundance += abundance_changes[names['33']]

    index2 = check_peptide_existence(segments[1], 108, True, new_peptide_fragments)
    if index2 == -1:
        peptide108 = Peptide(segments[0], 108, abundance_changes[names['108']], True)
        new_peptide_fragments.peptides[108].append(peptide108)
    else:
        new_peptide_fragments.peptides[108][index2].abundance += abundance_changes[names['108']]

# heme detoxification protein equations
    abundance_changes[names['Fpp']] += -hdp.kCatKm * concentrations[1] * food_vacuole_volume
    # fpp degradation
    abundance_changes[names['Hz']] += hdp.kCatKm * concentrations[1] * food_vacuole_volume
    # hz production

# other enzymes
    for enzyme in enzymes:
        absolute_concentration_of_possible_substrates = get_absolute_concentration_of_possible_substrates(enzyme, concentrations)
        start_index = names[str(enzyme.minLength)]
        end_index = names[str(enzyme.maxLength)]
        max_enzyme_abundance = enzyme.maxAbundance
        k_cat_k_m = enzyme.kCatKm
        if enzyme.abundance:
            relative_enzyme_abundance = enzyme.abundance[eci]
        else:
            relative_enzyme_abundance = 1
        for i in range(start_index, end_index+1):  # i = index of substrate in ode solver array
            if concentrations[i] != 0:  # check if substrate has concentration
                p = concentrations[i] / absolute_concentration_of_possible_substrates  # likelihood of substrate being used
                conversion = p * max_enzyme_abundance * relative_enzyme_abundance * k_cat_k_m * concentrations[i]
                abundance_changes[i] += -conversion
                # decay of substrate
                increase_of_products = conversion
                peptides_for_length = currentPeptideFragments.peptides[lengths[str(i)]]  # takes peptides from CURRENT state
                number_of_peptides_for_length = len(peptides_for_length)
                for peptide in peptides_for_length:
                    if enzyme is fal2:  # fal2 removes fpp
                        peptide.containsFpp = False
                        abundance_changes[names['Fpp']] += increase_of_products / number_of_peptides_for_length
                    if enzyme.isExopeptidase and enzyme.aminoPeptidaseIndex != -1:
                        enzyme.first_aa_cleavage(peptide)
                        add_peptide_fragments(peptide, number_of_peptides_for_length, increase_of_products, new_peptide_fragments)
                        #  changes peptides in NEW state
                    else:
                        enzyme.configure_cleavage_sites(peptide)
        currentPeptideFragments = new_peptide_fragments  # actualization
        return abundance_changes
