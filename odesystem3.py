import numpy as np
from datainterpretation2 import names, lengths
from foodvacuole import get_volume, get_hb_abundance_change
from peptides import Peptide, PeptideChange
from enzymes3 import plas1, hdp, fal2, enzymes
from dataimport import alphaHbChain, betaHbChain


def get_relative_enzyme_concentration_index(t):
    """

    :param t: timepoint
    :return: nearest index in abundance list for timepoint
    """
    if t % 2 <= 1:
        return int(t/2)-1  # index for protein abundance starts at 0 for 2hpi
    else:
        return int(t/2)


def get_absolute_concentration_of_possible_substrates(enzyme, concentrations):
    """

    :param enzyme: enzyme, used for to determine possible peptide lengths
    :param concentrations: list of current concentrations (calculated based on abundances and food vacuole volume for timepoint)
    :return: float, sum of concentrations of possible substrates in fmol
    """
    acos = 0  # absolute concentration of substrates in M
    min_index = names[str(enzyme.minLength)]
    max_index = names[str(enzyme.maxLength)]
    for i in range(min_index, max_index+1):
        acos += concentrations[i]
    return acos


def get_void_array():
    """

    :return: void list for ode solver (used for abundances and abundance changes
    """
    void_array = []
    for i in range(0, len(names)):
        void_array.append(0)
    return void_array


timeGrid = np.linspace(2, 47, 101)
initialAbundances = get_void_array()
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


def get_peptide_index(peptide_change, length, peptides_list):
    """ returns peptide index in current_peptide_fragments.peptides[length] (length of peptide in AAs)
        or -1 if not found"""
    index = -1
    for peptideIndex in range(0, len(peptides_list[length])):
        peptide = peptides_list[length][peptideIndex]
        if peptide.sequence == peptide_change.sequence and peptide.containsFpp == peptide_change.containsFpp:
            index = peptideIndex
    return index


def update_peptides_list(peptides_list, peptide_changes_list, time_step):
    for length in range(1, len(peptides_list)):
        for peptide_change in peptide_changes_list[length]:
            peptide_index = get_peptide_index(peptide_change, length, peptides_list)
            if peptide_index == -1:  # add new Peptide
                new_peptide_abundance = peptide_change.abundanceChange * time_step
                new_peptide = Peptide(peptide_change.sequence, peptide_change.length, new_peptide_abundance, peptide_change.containsFpp)
                peptides_list[length].append(new_peptide)
            else:
                peptides_list[length][peptide_index].abundance += peptide_change.abundanceChange * time_step


def longer_fragment_tuple(length_of_fragment1, length_of_fragment2):
    """
    (True, False) if length_of_fragment1 > length_of_fragment2
    etc.
    """
    if length_of_fragment1 >= length_of_fragment2:
        answers = (True, False)
        return answers
    else:
        answers = (False, True)
        return answers


def add_peptide_fragment_changes(peptide, number_of_peptides_for_length, total_abundance_change, peptide_abundance_changes_list):
    """
    adds resulting peptide fragment changes (dn/dt) of each cleavage site for given peptide to peptide_abundance_changes_list
    :param peptide: peptide, peptide to be cleaved (cleavage sites have been set previously)
    :param number_of_peptides_for_length: positive integer, number of different peptides for length
    :param total_abundance_change: positive integer, change in abundance for all peptides of current length
    :return: void, modifies currentPeptideFragments
    """
    substrate_peptide_abundance_change = PeptideChange(peptide.sequence, peptide.length, total_abundance_change / number_of_peptides_for_length, False)
    peptide_abundance_changes_list[substrate_peptide_abundance_change.length].append(substrate_peptide_abundance_change)
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

        fragment1_change = PeptideChange(segments[0], length_of_fragment1, abundance_change_of_fragments, fpp_distribution[0])
        peptide_abundance_changes_list[length_of_fragment1].append(fragment1_change)

        fragment2_change = PeptideChange(segments[1], length_of_fragment2, abundance_change_of_fragments, fpp_distribution[1])
        peptide_abundance_changes_list[length_of_fragment1].append(fragment2_change)


def derivative(abundances, t, current_peptide_fragments):
    """
    derivative function for ode solver
    """
# update of current_peptide_fragments.peptides_list
    if t != 0:
        time_step = t - current_peptide_fragments.lastTimePoint
        update_peptides_list(current_peptide_fragments.peptidesList, current_peptide_fragments.peptideChangesList, time_step)
    
    concentrations = []
    abundance_changes = get_void_array()
    current_peptide_fragments.reset_peptide_changes_list()
    peptide_abundance_changes_list = current_peptide_fragments.peptideChangesList  # holds abundance changes for given timestep
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

    peptide33_change = PeptideChange(segments[0], 33, abundance_changes[names['33']], False)
    peptide_abundance_changes_list[33].append(peptide33_change)

    peptide108 = PeptideChange(segments[1], 108, abundance_changes[names['108']], True)
    peptide_abundance_changes_list[108].append(peptide108)

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
                abundance_changes[i] -= conversion
                # decay of substrate
                increase_of_products = conversion
                peptides_for_length = current_peptide_fragments.peptidesList[lengths[str(i)]]  # takes peptides from CURRENT state
                number_of_peptides_for_length = len(peptides_for_length)
                for peptide in peptides_for_length:
                    if enzyme is fal2:  # fal2 removes fpp
                        if peptide.containsFpp:
                            peptide.containsFpp = False
                            abundance_changes[names['Fpp']] += increase_of_products / number_of_peptides_for_length
                    if enzyme.isExopeptidase and enzyme.aminoPeptidaseIndex != -1:
                        enzyme.first_aa_cleavage(peptide)
                        add_peptide_fragment_changes(peptide, number_of_peptides_for_length, increase_of_products, peptide_abundance_changes_list)
                    else:
                        enzyme.configure_cleavage_sites(peptide)
                        add_peptide_fragment_changes(peptide, number_of_peptides_for_length, increase_of_products, peptide_abundance_changes_list)
        current_peptide_fragments.lastTimePoint = t
        return abundance_changes
