import matplotlib.pyplot as plt
from scipy.integrate import odeint

import odesystem
from expdataimport import expValues

n = odeint(odesystem.derivative, odesystem.initialAbundances, odesystem.timeGrid)

hbAbundance = []
s0Abundance = []
s1Abundance = []
s2Abundance = []
s3Abundance = []
fppAbundance = []
hzAbundance = []
s4Abundance = []

for values in n:
    UCF_TO_FMOL = 0.001
    hbAbundance.append(values[0])
    s0Abundance.append(values[1])
    s1Abundance.append(values[2])
    s2Abundance.append(values[3])
    s3Abundance.append(values[4])
    fppAbundance.append(values[5])
    hzAbundance.append(values[6])
    s4Abundance.append(values[7])

Abundances = [hbAbundance, s0Abundance, s1Abundance, s2Abundance, s3Abundance,
              fppAbundance, hzAbundance, s4Abundance]
AbundanceNames = ['hb', 's0', 's1', 's2', 's3', 'fpp', 'hz', 's4']
fig = plt.figure()
plot = fig.add_subplot(111)
'''
for i in range(0, len(Abundances)-1):
    plot.plot(odesystem.timeGrid, Abundances[i], label=AbundanceNames[i])
plot.plot(expValues[0], expValues[1], label="exp hz values")
'''
plot.plot(odesystem.timeGrid, Abundances[5], label=AbundanceNames[5])
plot.set_ylabel("n [fmol]")
plot.set_xlabel("t [h]")
plot.legend()
plt.show()
