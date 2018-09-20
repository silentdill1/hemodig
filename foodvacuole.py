from math import exp
from hostcell import get_hb_concentration
import matplotlib.pyplot as plt
import numpy as np

INITIAL_VOLUME = 0.5  # in fL
# TODO: concentrating factor bestimmen (Verhältnis Volumenzunahme Parasit / Vacuole)
CONCENTRATING_FACTOR = 10  # vesicle volume lost through osmosis
UCF_TO_FMOL = 10**(-3)


def get_volume(t):
    volume = 2.6+(0.16-2.6)/(1+exp((t-20.8)/1.5))
    return volume  # in fl


def get_hb_abundance_change(t):
    hb_abundance_change = get_volume_change(t) * get_hb_concentration(t) * CONCENTRATING_FACTOR * UCF_TO_FMOL
    return hb_abundance_change  # in fmol


def get_volume_change(t):
    # volumechange = -(0.16-2.6)/((1+exp((t-20.8)/1.5)**2)*1/1.5*exp((t-20.8)/1.5)) vertippt?
    volume_change = 1.62667 * exp(0.66667 * (-20.8 + t)) / (1 + exp(0.6667 * (-20.8 + t))) ** 2  # mathematica
    # dV/dt(t)
    return volume_change

'''
hbChange = []
timeGrid = np.linspace(2, 50, 101)
for time in timeGrid:
    hbChange.append(get_hb_abundance_change(time))
fig = plt.figure()
plot = fig.add_subplot(111)
plot.plot(timeGrid, hbChange)
fig.savefig('bla.png')
plt.interactive(False)
plt.show()
'''