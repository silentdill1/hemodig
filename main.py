import matplotlib.pyplot as plt
from scipy.integrate import odeint
import odesystem2
from expdataimport import expValues
from datainterpretation import names

n = odeint(odesystem2.derivative, odesystem2.initialAbundances, odesystem2.timeGrid)

hbAbundance = []
hzAbundance = []
asAbundance = []

for values in n:
    hbAbundance.append(values[names['Hb']])
    hzAbundance.append(values[names['Hz']])
    asAbundance.append(values[names['1']]/100)

Abundances = [hbAbundance, hzAbundance, asAbundance]
AbundanceNames = ['hb', 'hz', 'as']
fig = plt.figure()
plot = fig.add_subplot(111)
for i in range(0, len(Abundances)):
    plot.plot(odesystem2.timeGrid, Abundances[i], label=AbundanceNames[i])
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
