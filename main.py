import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp  # operates on changed version that passes ode solver timesteps to function
import odesystem3
# import odesystem2
from expdataimport import expValues
from datainterpretation2 import names
import pickle


solution = solve_ivp(odesystem3.wrapper, [2, 50], odesystem3.initialAbundances, t_eval=odesystem3.timeGrid)
with open('data.pickle', 'wb') as f:
    pickle.dump(solution, f)

hbAbundance = []
hzAbundance = []
ipAbundance = []
lpAbundance = []
spAbundance = []
dpAbundance = []
asAbundance = []


Abundances = [hbAbundance, hzAbundance, asAbundance]
AbundanceNames = ['hb', 'hz', 'AAs']
fig = plt.figure()
plot = fig.add_subplot(111)
plot.plot(odesystem3.timeGrid, solution.y[names['1']], label='AAs')
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
