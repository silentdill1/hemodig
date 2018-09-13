# uses pandas dataframe object
import pandas as pd
import sys


def det_dataframe_line(dataframe, enzyme):
    line = 0
    for i in range(dataframe.shape[0]):
        if dataframe.at[i, 'PlasmoDB ID \nor NCBI \nAccession \nNumber'] == enzyme.pfidOld:
            line = i
    return line


def get_abundance(dataframe, enzyme):

    values = []
    values_non_log = []         # reversed log2 transformed values
    values_converted = []       # scaled in respect to maximum

    for i in range(2, 50, 2):
        values.append(dataframe.at[enzyme.line, str(i)+'hpi'])

    for value in values:
        values_non_log.append(2**value)

    maxValue = 0
    for value in values_non_log:
        if maxValue <= value:
            maxValue = value

    for value in values_non_log:
        values_converted.append(value/maxValue)

    return values_converted


df = pd.read_excel('enzyme initialization data.xlsx', header=0)


def initialize_parameters(enzyme):
    initialization_parameters = []
    row = -1
    for row_index in range(0, df.shape[0]):
        if enzyme.name == df.at[row_index, 'Name']:
            row = row_index
    if row is -1:
        sys.exit('Enzyme not found')
    else:
        for i in range(1, df.shape[1]):
            initialization_parameters.append(df.at[row, df.keys()[i]])
        enzyme.pfidOld = initialization_parameters[0]
        enzyme.kCatKm = initialization_parameters[1]
        enzyme.maxAbundance = initialization_parameters[2]






