import sys

names = {'Hb': 0, 'Fpp': 1, 'Hz': 2}
currentIndex = 3
for i in range(20, 148, 2):
    names[str(i)+'wFpp'] = currentIndex
    currentIndex += 1
for i in range(1, 21):
    names[str(i)] = currentIndex
    currentIndex += 1
for i in range(22, 72, 2):
    names[str(i)] = str(i)+'wFpp'
    currentIndex += 1
print(names)

def get_index_list(enzyme):
    if type is 0:  # initial peptide peptidases
        wFppIndices = [names[str(enzyme.minSize)+'wFpp'], names[str(enzyme.maxSize)+'wFpp']]
        woFppIndices = 'None'

    elif type is 1:  # long peptide peptidases
        wFppIndices = [names[str(enzyme.minSize)+'wFpp'], names[str(enzyme.maxSize)+'wFpp']]
        woFppIndices = [names[str(enzyme.minSize)], names['70']]

    elif (type is 2) or (type is 3):  # dipeptidyl aminopeptidase
        wFppIndices = 'None'
        woFppIndices = [names[str(enzyme.minSize)], names[str(enzyme.maxSize)]]

    else:
        sys.exit('No valid enzyme type')

    return [wFppIndices, woFppIndices]
