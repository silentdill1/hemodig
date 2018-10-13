

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


class PeptideChange:
    """ represents abundance change (dn/dt) for last time point of ode solver for given sequence to be saved in peptideChangesList
        of currentPeptideFragments"""
    def __init__(self, sequence, length, abundance_change, contains_fpp):
        self.sequence = sequence
        self.length = length
        self.abundanceChange = abundance_change
        self.containsFpp = contains_fpp


class Peptides:
    def __init__(self, number_of_possible_chain_lengths):
        self.peptidesList = []
        self.numberOfPossibleChainLengths = number_of_possible_chain_lengths
        # list of lists containing currently existing peptides for chain length (index), 0 is not used
        # ascending
        for i in range(0, self.numberOfPossibleChainLengths+1):
            self.peptidesList.append([])
        self.numberOfTimePoints = 0
        self.peptideChangesList = []  # holds dn/dt for last timepoint
        self.reset_peptide_changes_list()

    def reset_peptide_changes_list(self):
        self.peptideChangesList = []
        for i in range(0, self.numberOfPossibleChainLengths+1):
            self.peptideChangesList.append([])

