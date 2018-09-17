# uses pandas data frame object
import pandas as pd
import sys


def det_data_frame_line(data_frame, enzyme):
    line = 0
    for i in range(data_frame.shape[0]):
        if data_frame.at[i, 'PlasmoDB ID \nor NCBI \nAccession \nNumber'] == enzyme.pfidOld:
            line = i
    return line


def get_abundance(data_frame, enzyme):

    values = []
    values_non_log = []         # reversed log2 transformed values
    values_converted = []       # scaled in respect to maximum

    for i in range(2, 50, 2):
        values.append(data_frame.at[enzyme.line, str(i)+'hpi'])

    for value in values:
        values_non_log.append(2**value)

    max_value = 0
    for value in values_non_log:
        if max_value <= value:
            max_value = value

    for value in values_non_log:
        values_converted.append(value/max_value)

    return values_converted


dfInit = pd.read_excel('enzyme initialization data.xlsx', header=0)  # data frame for initialization values
dfProt = pd.read_excel('proteome_plas.xlsx', sep='\t', header=1)  # data frame for protein abundance


def initialize_parameters(enzyme):
    initialization_parameters = []
    row = -1
    for row_index in range(0, dfInit.shape[0]):
        if enzyme.name == dfInit.at[row_index, 'Name']:
            row = row_index
    if row is -1:
        sys.exit('Enzyme not found')
    else:
        for i in range(1, dfInit.shape[1]):
            initialization_parameters.append(dfInit.at[row, dfInit.keys()[i]])
        enzyme.pfidOld = initialization_parameters[0]
        enzyme.kCatKm = initialization_parameters[1]
        enzyme.maxAbundance = initialization_parameters[2]
        if enzyme.pfidOld:
            enzyme.line = det_data_frame_line(dfProt, enzyme)
            enzyme.abundance = get_abundance(dfProt, enzyme)


AAThreeLetterCode = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H',
                     'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
                     'TYR': 'Y', 'VAL': 'V'}


def convert_sequence(aa_list):
    sequence_list = []
    for aa in aa_list:
        sequence_list.append(AAThreeLetterCode[aa])
    return tuple(sequence_list)


def import_aa_list(data_frame, name, length):
    aa_list = []
    for i in range(0, length):
        aa_list.append(data_frame.at[i, name])
    return aa_list


dfHbChains = pd.read_excel('hb_chains_sequence.xlsx', header=0)
alphaHbChain = convert_sequence(import_aa_list(dfHbChains, 'alpha Hb chain', 141))
betaHbChain = convert_sequence(import_aa_list(dfHbChains, 'beta Hb chain', 146))





