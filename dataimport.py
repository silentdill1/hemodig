# uses pandas dataframe object


def plas_det_dataframe_line(dataframe, enzyme):
    line = 0
    for i in range(dataframe.shape[0]):
        if dataframe.at[i, 'PlasmoDB ID \nor NCBI \nAccession \nNumber'] == enzyme.pfidOld:
            line = i
    return line


def plas_get_abundance(dataframe, enzyme):

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





