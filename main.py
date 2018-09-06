import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from foodvacuole import get_volume
import odesystem
from expdataimport import expValues

n = odeint(odesystem.derivative, odesystem.initialAbundances, odesystem.timeGrid)

hbAbundance = []
fppAbundance = []
hzAbundance = []
s5Abundance = []
s6Abundance = []

for values in n:
    UCF_TO_FMOL = 0.001
    hbAbundance.append(values[0])
    fppAbundance.append(values[5])
    hzAbundance.append(values[6])
    s5Abundance.append(values[8])
    s6Abundance.append(values[9])
Abundances = [hbAbundance, fppAbundance, hzAbundance, s5Abundance, s6Abundance]
AbundanceNames = ['hb', 'fpp', 'hz', 's5', 's6']
fig = plt.figure()
plot = fig.add_subplot(111)
for i in range(0, len(Abundances)):
    plot.plot(odesystem.timeGrid, Abundances[i], label=AbundanceNames[i])
plot.plot(expValues[0], expValues[1], label="exp hz values")
plot.set_ylabel("n [fmol]")
plot.set_xlabel("t [h]")
plot.legend(fontsize='small')
'''
plot2 = fig.add_subplot(222)
timegrid = np.linspace(0, 40, 100)
vollist = []

for time in timegrid:
    vollist.append(get_volume(time))

plot2.plot(timegrid, vollist)
'''
plt.show()
