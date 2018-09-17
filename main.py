import matplotlib.pyplot as plt
from scipy.integrate import odeint
import odesystem3
# import odesystem2
from expdataimport import expValues
from datainterpretation2 import names
from peptides import Peptides
import pandas as pd

currentPeptideFragments = Peptides(108)
n = odeint(odesystem3.derivative, odesystem3.initialAbundances, odesystem3.timeGrid, args=(currentPeptideFragments,))

hbAbundance = []
hzAbundance = []
ipAbundance = []
lpAbundance = []
spAbundance = []
dpAbundance = []
asAbundance = []

for values in n:
    asAbundance.append(values[names['1']])

Abundances = [hbAbundance, hzAbundance, asAbundance]
AbundanceNames = ['hb', 'hz', 'AAs']
fig = plt.figure()
plot = fig.add_subplot(111)
plot.plot(odesystem3.timeGrid, asAbundance)
'''
for i in range(0, len(Abundances)):
    plot.plot(odesystem3.timeGrid, Abundances[i], label=AbundanceNames[i])
plot.plot(expValues[0], expValues[1], label="exp hz values")
'''
plot.set_ylabel("n [fmol]")
plot.set_xlabel("t [h]")
plot.legend(fontsize='small')
plt.show()
fig.savefig('plot_model3')
