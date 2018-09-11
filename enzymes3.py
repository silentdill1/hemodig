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
        available_amino_acids = set()
        # indices of amino acids in peptide who are not yet part of a cleavage site
        for i in range(0, peptide.length):
            available_amino_acids.add(i)
        for motif in self.cleavageMotifs:  # motif = tuple of AA Code and LiN, e.g. ('X', 4)
            for index_of_amino_acid in range(0, len(peptide.sequence)):
                amino_acid = peptide.sequence[index_of_amino_acid]
                if amino_acid == motif[0]:
                    cleavage_site_index = len(peptide.cleavageSites)
                    # save index in cleavage site array for future addition of LiNs for AAs in vicinity
                    peptide.cleavageSites.append([index_of_amino_acid, motif[1]])  # create new cleavage site
                    # motif[1] is the LiN for current cleavage motif
                    available_amino_acids.remove(index_of_amino_acid)
                    # current amino acid is no longer available as cleavage option,
                    # since a cleavage site there is already taken into account
                    for i in range(-3, 4):  # check nearby amino acids for changes in cleavage likelihood
                        if i == 0:
                            continue  # no need to check set for removed amino acid
                        else:
                            index_of_nearby_amino_acid = index_of_amino_acid + i

                        if index_of_nearby_amino_acid in available_amino_acids:
                            for motif2 in self.cleavageMotifs:  # check all motifs again, therefore different variable name
                                if peptide.sequence[index_of_nearby_amino_acid] == motif2[0]:
                                    peptide.cleavageSites[cleavage_site_index][1] += motif2[1]  # add LiN of motif
                                    # cleavageSites = list of lists of position after which cut happens and LiN, e.g. [112, 5]
                                    available_amino_acids.remove(index_of_nearby_amino_acid)  # amino acid no longer available for cleavage site


plas2 = Enzyme('Plasmepsin II', )





