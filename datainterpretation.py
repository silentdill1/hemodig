import sys

names = {'Hb': 0, 'Fpp': 1, 'Hz': 2}  # lists protein names and corresponding indices of array for ode solver
lengths = {}  # lists indices of array for ode solver and corresponding peptide length
currentIndex = 3
for i in range(20, 148, 2):
    names[str(i)+'wFpp'] = currentIndex
    lengths[currentIndex] = i
    currentIndex += 1
for i in range(1, 21):
    names[str(i)] = currentIndex
    lengths[currentIndex] = i
    currentIndex += 1
for i in range(22, 72, 2):
    names[str(i)] = currentIndex
    lengths[currentIndex] = i
    currentIndex += 1


def get_index_list(enzyme):
    # returns list containing lists with start and stop indices for peptides with and without fpp (type)
    # or String 'None' if not operating on corresponding type of species
    category = enzyme.category
    print(category)
    if category is 0:  # initial peptide peptidases
        w_fpp_indices = [names[str(enzyme.minSize)+'wFpp'], names[str(enzyme.maxSize)+'wFpp']]
        wo_fpp_indices = 'None'

    elif category is 1:  # long peptide peptidases
        w_fpp_indices = [names[str(enzyme.minSize)+'wFpp'], names[str(enzyme.maxSize)+'wFpp']]
        wo_fpp_indices = [names[str(enzyme.minSize)], names['70']]

    elif (category is 2) or (category is 3):  # dipeptidyl aminopeptidase
        w_fpp_indices = 'None'
        wo_fpp_indices = [names[str(enzyme.minSize)], names[str(enzyme.maxSize)]]

    else:
        sys.exit('No valid enzyme category')

    return [w_fpp_indices, wo_fpp_indices]
