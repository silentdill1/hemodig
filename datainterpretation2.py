names = {'Hb': 0, 'Fpp': 1, 'Hz': 2}  # lists protein names and corresponding indices of array for ode solver
lengths = {}  # lists indices of array for ode solver and corresponding peptide length
currentIndex = 3
for i in range(1, 109):
    names[str(i)] = currentIndex  # keys are always strings for distinguishing purposes
    lengths[str(currentIndex)] = i
    currentIndex += 1
