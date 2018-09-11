import matplotlib.pyplot as plt
from scipy.integrate import odeint
import odesystem2
from expdataimport import expValues
from datainterpretation import names

n = odeint(odesystem2.derivative, odesystem2.initialAbundances, odesystem2.timeGrid)

hbAbundance = []
hzAbundance = []
ipAbundance = []  # wFpp
pwFppAbundance = []  # peptides with fpp
lpAbundance = []
spAbundance = []
dpAbundance = []
asAbundance = []

currentIndex = 0
for values in n:
    hbAbundance.append(values[names['Hb']])
    hzAbundance.append(values[names['Hz']])

    ipAbundance.append(0)
    for i in range(names['80wFpp'], names['146wFpp']+1):
        ipAbundance[currentIndex] += values[i]

    pwFppAbundance.append(ipAbundance[currentIndex])
    for i in range(names['10wFpp'], names['80wFpp']+1):
        pwFppAbundance[currentIndex] += values[i]

    lpAbundance.append(0)
    for i in range(names['20wFpp'], names['80wFpp']+1):
        lpAbundance[currentIndex] += values[i]
    for i in range(names['20'], names['70']):
        lpAbundance[currentIndex] += values[i]

    spAbundance.append(0)
    for i in range(names['4'], names['20']+1):
        lpAbundance[currentIndex] += values[i]

    dpAbundance.append(values[names['2']])
    asAbundance.append(values[names['1']])
    currentIndex += 1

Abundances = [hbAbundance, hzAbundance, ipAbundance, pwFppAbundance, lpAbundance,
              spAbundance, dpAbundance]
AbundanceNames = ['hb', 'hz', 'ip', 'wfpp', 'lp', 'sp', 'dp']
fig = plt.figure()
plot = fig.add_subplot(111)
'''
for i in range(0, len(Abundances)):
    plot.plot(odesystem2.timeGrid, Abundances[i], label=AbundanceNames[i])
plot.plot(expValues[0], expValues[1], label="exp hz values")
'''
plot.plot(odesystem2.timeGrid, asAbundance, label='amino acids')
plot.set_ylabel("n [fmol]")
plot.set_xlabel("t [h]")
plot.legend(fontsize='small')
plt.show()
fig.savefig('plot_amino_acids')
