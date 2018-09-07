names = [[0, 'Hb'], [1, 'Fpp'], [2, 'Hz']]
currentIndex = 3
for i in range(10, 148, 2):
    names.append([currentIndex, str(i)+'wFpp'])
    currentIndex += 1
for i in range(1, 21):
    names.append([currentIndex, str(i)])
    currentIndex += 1
for i in range(22, 72, 2):
    names.append([currentIndex, str(i)])
    currentIndex += 1
print(names)

