class Peptide:
    def __init__(self, sequence, length, abundance, contains_fpp):
        self.sequence = sequence  # tuple of amino acids one letter code
        self.length = length  # length of peptide = index in peptides list
        self.abundance = abundance  # in fmol
        self.containsFpp = contains_fpp
        self.cleavageSites = []
        # list of lists containing indices for possible cleavage and
        # corresponding likelihood number (LiN) (see enzymes3 for documentation),
        # has to be updated according to current enzymes specificity
        self.sumOfLiNsForAllCleavageSites = 0
        # sum of LiNs for all cleavage sites for current enzyme,
        # has to be reset to 0 upon enzyme change


class Peptides:
    def __init__(self, number_of_possible_chain_lengths):
        self.peptides = []
        # list of lists containing currently existing peptides for chain length (index), 0 is not used
        # ascending
        for i in range(0, number_of_possible_chain_lengths+1):
            self.peptides.append([])
