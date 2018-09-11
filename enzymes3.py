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
        self.sumOfLiNsForAllCleavageSites = 0  # see n above, has to be set to 0 upon enzyme change
        self.kCatKm = k_cat_k_m
        self.maxEnzymeAbundance = max_enzyme_abundance
