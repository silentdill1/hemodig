import peptides
import numpy as np


class Enzyme:
    def __init__(self, name, pf_id_old, cleavage_motifs, k_cat_k_m, max_enzyme_abundance):
        self.name = name
        self.pfidOld = pf_id_old
        self.cleavageMotifs = cleavage_motifs
        # tuple of tuples containing amino acid in one letter code and corresponding change
        # in likelihood number (LiN), e.g. ('R', 3) means if 'R' is found in the vicinity of
        # current cleavage site cleavage likelihood increases to from o/n to (o+3)/(n+3)
        # o - sum of LiNs of previously found amino acids for current cleavage site
        # n - sum of LiNs of previously found amino acids for all cleavage sites
        self.kCatKm = k_cat_k_m  # in 1/(s*M)
        self.maxEnzymeAbundance = max_enzyme_abundance  # in fmol

    def configure_cleavage_sites(self, peptide):
        available_amino_acid_indices = []
        # indices of amino acids in peptide who are not yet part of a cleavage site
        available_motif_indices = []
        # indices of amino acids in peptide who are not yet part of a cleavage site
        for i in range(0, peptide.length):
            available_amino_acid_indices.append(i)
        for i in range(0, len(self.cleavageMotifs)):
                available_motif_indices.append(i)
        while len(available_motif_indices) != 0:
            index_of_motif_index = np.random.randint(0, len(available_motif_indices))  # pick first AS for cleavage site randomly
            # TODO: randomized changes between current motif (currently looks for 'A' throughout the whole peptide chain,
            # TODO: until it takes new motif for first AA of cleavage site)
            motif_index = available_motif_indices[index_of_motif_index]
            motif = self.cleavageMotifs[motif_index]  # motif = tuple of AA Code and LiN, e.g. ('X', 4)
            index_of_amino_acid_index = 0

            while index_of_amino_acid_index < len(available_amino_acid_indices):
                amino_acid_index = available_amino_acid_indices[index_of_amino_acid_index]
                amino_acid = peptide.sequence[amino_acid_index]
                if amino_acid_index != (peptide.length-1) and amino_acid == motif[0]:  # no cleavage after last AA
                    cleavage_site_index = len(peptide.cleavageSites)
                    # save index in cleavage site array for future addition of LiNs for AAs in vicinity
                    peptide.cleavageSites.append([amino_acid_index, motif[1]])  # create new cleavage site
                    # motif[1] is the LiN for current cleavage motif
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
                                    available_amino_acid_indices.remove(index_of_nearby_amino_acid)  # amino acid no longer available for cleavage site
                                    if i < 0:  # removing AAs index before current index -> moves and would skip one index
                                        index_of_amino_acid_index -= 1
                else:
                    index_of_amino_acid_index += 1

            available_motif_indices.remove(motif_index)


hydrophobicAminoAcidMotifs = (('A', 1), ('G', 1), ('H', 1), ('I', 1), ('L', 1), ('M', 1), ('F', 1),
                              ('P', 1), ('V', 1))
plas2 = Enzyme('Plasmepsin II', 'PF14_0077', hydrophobicAminoAcidMotifs, 500 * 10**3, 30 * 10**(-6))
testSequence = ('A', 'X', 'X', 'A', 'L', 'X', 'X', 'X', 'X', 'A', 'T', 'L', 'F', 'L', 'L', 'X', 'X', 'A', 'T', 'L', 'F')
pep = peptides.Peptide(testSequence, len(testSequence))
plas2.configure_cleavage_sites(pep)
print(pep.cleavageSites)






