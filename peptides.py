class Peptide:
    def __init__(self, sequence, length):
        self.sequence = sequence  # tuple of amino acids one letter code
        self.length = length
        self.cleavageSites = []
        # list of lists containing indices for possible cleavage and
        # corresponding likelihood number (see enzymes3 for documentation),
        # has to be updated according to current enzymes specificity


class Peptides:
    def __init__(self, number_of_possible_chain_lengths):
        self.peptides = []
        # list of lists containing currently existing peptides for chain length index+1
        # ascending
        for i in range(0, number_of_possible_chain_lengths):
            self.peptides.append([])
