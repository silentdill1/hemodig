import numpy as np
from datainterpretation2 import names, lengths
from foodvacuole import get_volume, get_hb_abundance_change
from peptides import Peptides, Peptide
from enzymes3 import plas1, hdp, fal2, enzymes
from dataimport import alphaHbChain, betaHbChain


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


def check_peptide_existence(sequence, length):
    """ returns peptide index in currentPeptideFragments.peptides[lengthIndex] (lengthIndex = length-1)
        or -1 if not found"""
    index = -1
    for peptideIndex in range(0, len(currentPeptideFragments.peptides[length-1])):
        peptide = currentPeptideFragments.peptides[length-1][peptideIndex]
        if peptide.sequence == sequence:
            index = peptideIndex
    return index

def add_peptide_fragments(peptide, number_of_peptides_for_length, total_abundance_change):
    """

    :param peptide: peptide, peptide to be cleaved (cleavage sites have been set previously)
    :param number_of_peptides_for_length: postive integer, number of different peptides for length
    :param total_abundance_change: postive integer, change in abundance for all peptides of current length
    :return: void, modifies currentPeptideFragments
    """
    peptide.abundance -= total_abundance_change / number_of_peptides_for_length
    for cleavageSite in peptide.cleavageSites:
        segments = get_fragments(peptide.sequence, cleavageSite[0])
        length_of_fragment1 = cleavageSite[0]+1
        length_of_fragment2 = peptide.lengthIndex
        fragment1_index_in_peptides_of_length = check_peptide_existence(segments[0], length_index_of_fragment1)
        if fragment1_index_in_peptides_of_length == -1:
            currentPeptideFragments.peptides[length_index_of_fragment1]


def derivative(abundances, t):
    concentrations = []
    abundance_changes = get_void_array()
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
    index = check_peptide_existence(segments[0], 32)
    if index == -1
        currentPeptideFragments.peptides[32].append(Peptide(segments[0], 33, ))

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
                peptides_for_length = currentPeptideFragments.peptides[lengths[str(i)]]
                number_of_peptides_for_length = len(peptides_for_length)
                for peptide in peptides_for_length:
                    if enzyme.isExopeptidase and enzyme.aminoPeptidaseIndex != -1:
                        enzyme.first_aa_cleavage(peptide)
                        add_peptide_fragments(peptide, number_of_peptides_for_length, increase_of_products)
                    else:
                        enzyme.configure_cleavage_sites(peptide)
