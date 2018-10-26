import matplotlib.pyplot as plt
from changed_ivp import solve_ivp  # operates on changed version that passes ode solver timesteps to function
import odesystem3
import odesystem2
from scipy.integrate import odeint
from expdataimport import expValues
from datainterpretation2 import names
# from datainterpretation import names
import pickle
'''
y = odeint(odesystem2.derivative, odesystem2.initialAbundances, odesystem2.timeGrid)
with open('data_odesys2.pickle', 'wb') as f:
    pickle.dump(y, f)


y = pickle.load(open('data_odesys2.pickle', 'rb'))
# odesys2:
hbAbundance = []
hzAbundance = []
ipAbundance = []
lpAbundance = []
spAbundance = []
dpAbundance = []
asAbundance = []
print(names)
for timepoint in y:
    hzAbundance.append(timepoint[names['Hz']])
    asAbundance.append(timepoint[names['1']]/7)
    ipAbundance.append(0)
    for i in range(108, 148, 2):
        ipAbundance[len(ipAbundance)-1]+=timepoint[names[str(i)+'wFpp']]
    lpAbundance.append(0)
    for i in range(40, 108, 2):
        lpAbundance[len(lpAbundance) - 1] += timepoint[names[str(i) + 'wFpp']]
    for i in range(40, 72, 2):
        lpAbundance[len(lpAbundance) - 1] += timepoint[names[str(i)]]
    spAbundance.append(0)
    for i in range(10, 40, 2):
        spAbundance[len(spAbundance) - 1] += timepoint[names[str(i) + 'wFpp']]
    for i in range(20, 40, 2):
        spAbundance[len(spAbundance) - 1] += timepoint[names[str(i)]]
    for i in range(3, 20, 2):
        spAbundance[len(spAbundance) - 1] += timepoint[names[str(i)]]
    dpAbundance.append(timepoint[names['2']])

Abundances = [ipAbundance, lpAbundance, spAbundance, dpAbundance]
AbundanceNames = ['Lange Peptide', 'Peptide mittlerer Länge', 'Kurze Peptide', 'Dipeptide']

fig = plt.figure()
plot = fig.add_subplot(111)
for i in range(0, len(Abundances)):
    plot.plot(odesystem2.timeGrid, Abundances[i], label=AbundanceNames[i])
# plot.plot(odesystem2.timeGrid, hzAbundance, label='Theoretical Hz values')
# plot.plot(expValues[0], expValues[1], '.', label="Experimental Hz values")
plot.plot(odesystem2.timeGrid, asAbundance, label='Aminosäuren (skaliert)')
plot.plot(odesystem2.timeGrid, hzAbundance, label='Hemozoin')

plot.set_ylabel("n [fmol]")
plot.set_xlabel("t [h]")
plot.legend(fontsize='small')
fig.savefig('mama.png')

'''
# odesys3:
solution = solve_ivp(odesystem3.derivative_for_ode_solver, [12, 30], odesystem3.initialAbundances, odesystem3.update_function_for_ode_solver)
with open('data.pickle', 'wb') as f:
    pickle.dump(solution, f)
    
solution = pickle.load(open('data.pickle', 'rb'))

AbundanceNames = ['hb', 'hz', 'AAs']
fig = plt.figure()
plot = fig.add_subplot(111)
# plot.plot(solution.t, solution.y[names['1']], label='AAs')
plot.plot(solution.t, solution.y[names['Hb']], label='Hb')


plot.plot(expValues[0], expValues[1], label="exp hz values")
plot.set_ylabel("n [fmol]")
plot.set_xlabel("t [h]")
plot.legend(fontsize='small')
plt.show()
fig.savefig('plot_model3')
