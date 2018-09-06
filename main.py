import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from foodvacuole import get_volume
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
for i in range(0, len(Abundances)-1):
    plot.plot(odesystem.timeGrid, Abundances[i], label=AbundanceNames[i])
plot.plot(expValues[0], expValues[1], label="exp hz values")
plot.set_ylabel("n [fmol]")
plot.set_xlabel("t [h]")
plot.legend()

plot2 = fig.add_subplot(222)
timegrid = np.linspace(0, 40, 100)
vollist = []

for time in timegrid:
    vollist.append(get_volume(time))

plot2.plot(timegrid, vollist)
plt.show()
