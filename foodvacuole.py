from math import exp
from hostcell import get_hb_concentration

INITIAL_VOLUME = 0.5  # in fL
# TODO: concentrating factor bestimmen (Verhältnis Volumenzunahme Parasit / Vacuole)
CONCENTRATING_FACTOR = 70  # vesicle volume lost through osmosis

def get_volume(t):
    volume = 2.6+(0.16-2.6)/(1+exp((t-20.8)/1.5))
    return volume


def get_hb_abundance_change(t):
    hbabundancechange = get_volume_change(t) * get_hb_concentration(t) * CONCENTRATING_FACTOR  # in 10^-3 fmol
    return hbabundancechange

def get_volume_change(t):
    # volumechange = -(0.16-2.6)/((1+exp((t-20.8)/1.5)**2)*1/1.5*exp((t-20.8)/1.5)) vertippt?
    volumechange = 1.62667 * exp(0.66667 * (-20.8 + t)) / (1 + exp(0.6667 * (-20.8 + t))) ** 2  # mathematica
    # dV/dt(t)
    return volumechange