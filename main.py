import matplotlib.pyplot as plt
from changed_ivp import solve_ivp  # operates on changed version that passes ode solver timesteps to function
import odesystem3
from expdataimport import expValues
from datainterpretation2 import names
import pickle


solution = solve_ivp(odesystem3.derivative_for_ode_solver, [12, 30], odesystem3.initialAbundances, odesystem3.update_function_for_ode_solver)
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
plot.plot(solution.t, solution.y[names['1']], label='AAs')
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
