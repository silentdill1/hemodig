import peptides
import numpy as np
import sys


class Peptidase:
    def __init__(self, name, pf_id_old, is_exopeptidase, cleavage_motifs, k_cat_k_m, max_enzyme_abundance,
                 min_length, max_length, start_aa_index=-1, end_aa_index=-1, amino_peptidase_index=-1):
        self.name = name
        self.pfidOld = pf_id_old
        self.isExopeptidase = is_exopeptidase
        self.cleavageMotifs = cleavage_motifs
        # tuple of tuples containing amino acid in one letter code and corresponding change
        # in likelihood number (LiN), e.g. ('R', 3) means if 'R' is found in the vicinity of
        # current cleavage site cleavage likelihood increases to from o/n to (o+3)/(n+3)
        # o - sum of LiNs of previously found amino acids for current cleavage site
        # n - sum of LiNs of previously found amino acids for all cleavage sites
        self.kCatKm = k_cat_k_m  # in 1/(s*M)
        self.maxEnzymeAbundance = max_enzyme_abundance  # in fmol
        self.minLength = min_length  # minimum length of peptide that could be cleaved
        self.maxLength = max_length
        self.startAAIndex = start_aa_index  # index of first amino acid that could be cleavage center for exopeptidase
        self.endAAIndex = end_aa_index  # index of last amino acid that could be cleavage center for exopeptidase
        self.aminoPeptidaseIndex = amino_peptidase_index

    def configure_cleavage_sites(self, peptide):  # for endopeptidases, falcilysin and dpap
        peptide.cleavageSites = []  # resetting attributes
        peptide.sumOfLiNsForAllCleavageSites = 0

        available_amino_acid_indices = []
        # indices of amino acids in peptide who are not yet part of a cleavage site
        available_motif_indices = []
        # indices of amino acids in peptide who are not yet part of a cleavage site

        if self.isExopeptidase:  # falcilysin operates only on N terminus of protein
            if self.startAAIndex == -1:
                sys.exit('No startAAIndex for exopeptidase')
            for i in range(self.startAAIndex, self.endAAIndex+1):
                available_amino_acid_indices.append(i)
        else:
            for i in range(0, peptide.length):
                available_amino_acid_indices.append(i)

        found_cleavage_site = False
        # is used only for falcilysin and dpap exonucleases,
        # because their cleavage likelihood depends on amino acids in vicinity,
        # changed to true after first cleavage site is found and terminates search

        for i in range(0, len(self.cleavageMotifs)):
                available_motif_indices.append(i)

        while len(available_motif_indices) != 0 and not found_cleavage_site:
            index_of_motif_index = np.random.randint(0, len(available_motif_indices))  # pick first AS for cleavage site randomly
            # TODO: randomized changes between current motif (currently looks for 'A' throughout the whole peptide chain,
            # TODO: until it takes new motif for first AA of cleavage site)
            motif_index = available_motif_indices[index_of_motif_index]
            motif = self.cleavageMotifs[motif_index]  # motif = tuple of AA Code and LiN, e.g. ('X', 4)
            index_of_amino_acid_index = 0

            while index_of_amino_acid_index < len(available_amino_acid_indices) and not found_cleavage_site:
                amino_acid_index = available_amino_acid_indices[index_of_amino_acid_index]
                amino_acid = peptide.sequence[amino_acid_index]
                if amino_acid_index != (peptide.length-1) and amino_acid == motif[0]:  # no cleavage after last AA
                    cleavage_site_index = len(peptide.cleavageSites)
                    # save index in cleavage site array for future addition of LiNs for AAs in vicinity
                    peptide.cleavageSites.append([amino_acid_index, motif[1]])  # create new cleavage site
                    # motif[1] is the LiN for current cleavage motif
                    LiN = motif[1]
                    available_amino_acid_indices.remove(amino_acid_index)
                    # current amino acid is no longer available as cleavage option,
                    # since a cleavage site there is already taken into account
                    for i in range(-3, 4):  # check nearby amino acids for changes in cleavage likelihood
                        if i == 0:
                            continue  # no need to check set for removed amino acid
                        else:
                            index_of_nearby_amino_acid = amino_acid_index + i

                        if index_of_nearby_amino_acid in available_amino_acid_indices:
                            for motif2 in self.cleavageMotifs:  # check all motifs again, therefore different variable name
                                if peptide.sequence[index_of_nearby_amino_acid] == motif2[0]:
                                    peptide.cleavageSites[cleavage_site_index][1] += motif2[1]  # add LiN of motif
                                    # cleavageSites = list of lists of position after which cut happens and LiN, e.g. [112, 5]
                                    LiN += motif2[1]
                                    available_amino_acid_indices.remove(index_of_nearby_amino_acid)  # amino acid no longer available for cleavage site
                                    if i < 0:  # removing AAs index before current index -> moves and would skip one index
                                        index_of_amino_acid_index -= 1
                    peptide.sumOfLiNsForAllCleavageSites += LiN
                    if self.isExopeptidase:
                        found_cleavage_site = True
                else:
                    index_of_amino_acid_index += 1

            available_motif_indices.remove(motif_index)

    def first_aa_cleavage(self, peptide):  # for single aa cleavage
        peptide.cleavageSites = []  # resetting attributes
        peptide.sumOfLiNsForAllCleavageSites = 0

        if self.aminoPeptidaseIndex == -1:
            sys.exit('No AminopeptidaseIndex given')
        else:
            api = self.aminoPeptidaseIndex  # index of amino acid to be checked for motif compatibility

        for motif in self.cleavageMotifs:
            if peptide.sequence[api] == motif[0]:
                peptide.cleavageSites.append([0, motif[1]])
                peptide.sumOfLiNsForAllCleavageSites = 1


def add_motif_tuples(motif_tuple1, motif_tuple2):
    list_sum_of_tuples = []
    doubles = []
    for motif_index1 in range(0, len(motif_tuple1)):
        for motif_index2 in range(0, len(motif_tuple2)):
            motif1 = motif_tuple1[motif_index1]
            motif2 = motif_tuple2[motif_index2]
            if motif1[0] == motif2[0]:
                doubles.append([motif_index1, motif_index2])
    for i in range(0, len(motif_tuple1)):
        no_double = True
        for pair in doubles:
            if i == pair[0]:
                no_double = False
        if no_double:
            list_sum_of_tuples.append(motif_tuple1[i])

    for i in range(0, len(motif_tuple2)):
            no_double = True
            for pair in doubles:
                if i == pair[1]:
                    no_double = False
            if no_double:
                    list_sum_of_tuples.append(motif_tuple2[i])
    for pair in doubles:
        list_sum_of_tuples.append((motif_tuple1[pair[0]][0], max(motif_tuple1[pair[0]][1], motif_tuple2[pair[1]][1])))

    return tuple(list_sum_of_tuples)


hydrophobicAminoAcidMotifs = (('A', 1), ('G', 1), ('H', 1), ('I', 1), ('L', 1), ('M', 1), ('F', 1),
                              ('P', 1), ('V', 1))
chargedAminoAcidMotifs = (('D', -10), ('E', -10), ('K', -10))
iPepMotifs = (('F', 2), ('I', 2))
lPepMotifs = (('R', 2), ('L', 2), ('F', 2))
flnMotifs = (('E', 2), ('M', 2), ('H', 2), ('S', 3), ('F', 2))
proApMotif = [('P', 1)]


plas2 = Peptidase('Plasmepsin II', 'PF14_0077', False, hydrophobicAminoAcidMotifs, 500 * 10**3, 30 * 10**(-6), 146, 80)
fln = Peptidase('Falcilysin', '', True, hydrophobicAminoAcidMotifs, 148 * 10**3, 80 * 10**(-6), 25, 8, 0, 8)
proAp = Peptidase('Prolin Aminopeptidase', '', True, tuple(proApMotif), 148 * 10**3, 80 * 10**(-6), 4, 2, amino_peptidase_index=0)
testSequence = ('A', 'X', 'X', 'A', 'L', 'X', 'X', 'X', 'X', 'A', 'T', 'L', 'F', 'L', 'L', 'X', 'X', 'A', 'T', 'L', 'F')
pep = peptides.Peptide(testSequence, len(testSequence))
print()





