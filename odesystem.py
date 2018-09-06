from foodvacuole import get_volume
from foodvacuole import get_hb_abundance_change
from foodvacuole import get_volume_change
from enzymes import plasmepsine
from enzymes import falcipain2
from enzymes import hdp

import numpy as np

UCF_PER_S_TO_PER_H = 3600  # unit conversion from 1/s to 1/h
UCF_MM_TO_M = 0.001  # unit conversion from M to mM
initialAbundances = [0, 0, 0, 0, 0, 0, 0, 0]
timeGrid = np.linspace(2, 47, 101)


# TODO: modellieren der möglichen Mehrfachspaltung durch die Peptidasen -> Substratkonkurrenz
# n vector of abundances in fmol in order hb, s0, s1, s2, s3, fpp, hz, s4
def derivative(n, t):
    y = []   # y vector of concentrations in M
    foodVacuoleVolume = get_volume(t)    # volume of foodvacuole
    for abundance in n:
        y.append(abundance/foodVacuoleVolume)
    eci = det_enzyme_concentration_index(t)
    
    plas1Kcat = plasmepsine[0].kCat * UCF_PER_S_TO_PER_H  # in 1/h
    plas1Km = plasmepsine[0].km * UCF_MM_TO_M  # in M
    plas1MEC = plasmepsine[2].MAX_ENZYME_ABUNDANCE/foodVacuoleVolume
    # maximum enzyme concentration in M
    plas1abund = plasmepsine[0].abundance[eci]  # relative abundance scaled to maximum
    
    plas2Kcat = plasmepsine[1].kCat * UCF_PER_S_TO_PER_H
    plas2Km = plasmepsine[1].km * UCF_MM_TO_M
    plas2MEC = plasmepsine[2].MAX_ENZYME_ABUNDANCE/foodVacuoleVolume
    plas2abund = plasmepsine[1].abundance[eci]
    
    plas4Kcat = plasmepsine[2].kCat * UCF_PER_S_TO_PER_H
    plas4Km = plasmepsine[2].km * UCF_MM_TO_M
    plas4MEC = plasmepsine[2].MAX_ENZYME_ABUNDANCE/foodVacuoleVolume
    plas4abund = plasmepsine[2].abundance[eci]

    falcipain2EC = falcipain2.MAX_ENZYME_ABUNDANCE/foodVacuoleVolume

    hdpEC = hdp.MAX_ENZYME_ABUNDANCE/foodVacuoleVolume

    # abundance changes in fmol/h
    dhbdt = (get_hb_abundance_change(t)/foodVacuoleVolume * UCF_MM_TO_M
             - plas1Kcat * plas1abund * plas1MEC * y[0]/(plas1Km+y[0])) * foodVacuoleVolume
    ds0dt = (plas1Kcat * plas1MEC * plas1abund * y[0]/(plas1Km+y[0])
             - plas2Kcat * plas2MEC * plas2abund * y[1]/(plas2Km+y[1])) * foodVacuoleVolume
    ds1dt = (plas2Kcat * plas2MEC * plas2abund * y[1]/(plas2Km+y[1])
             - plas4Kcat * get_occupancy_ratio(t) * plas4MEC * plas4abund * y[2]/(plas4Km+y[2])) * foodVacuoleVolume
    ds2dt = (plas4Kcat * get_occupancy_ratio(t) * plas4MEC * plas4abund * y[2]/(plas4Km+y[2])
             - plas4Kcat * get_occupancy_ratio(t) * plas4MEC * plas4abund * y[3]/(plas4Km+y[3])) * foodVacuoleVolume
    ds3dt = (plas4Kcat * get_occupancy_ratio(t) * plas4MEC * plas4abund * y[3]/(plas4Km+y[3])
             - falcipain2.kCat_Km * falcipain2EC * y[4]) * foodVacuoleVolume
    ds4dt = falcipain2.kCat_Km * falcipain2EC * y[4] * foodVacuoleVolume
    dfppdt = 4 * (falcipain2.kCat_Km * falcipain2EC * y[4] - hdp.kCat_Km * hdpEC * y[5]) * foodVacuoleVolume
    dhzdt = 0.5 * (hdp.kCat_Km * hdpEC * y[5] * foodVacuoleVolume)
    return [dhbdt, ds0dt, ds1dt, ds2dt, ds3dt, dfppdt, dhzdt, ds4dt]


def det_enzyme_concentration_index(t):
    if t % 2 <= 1:
        return int(t/2)-1  # index for protein abundance starts at 0 for 2hpi
    else:
        return int(t/2)


def det_enzyme_abundance_indices(t):
    return [int(t % 2), int((t % 2) + 1)]


def det_enzyme_abundance(timepoints, enzyme, t):  # linearization of intervals between datapoints
    e1 = enzyme.abundance[timepoints[0]]
    e2 = enzyme.abundance[timepoints[1]]
    return e1 + (e2-e1) * t

def get_maximum_volume_change():
    maxvolchange = 0
    for t in timeGrid:
        volchange = get_volume_change(t)
        if volchange > maxvolchange:
            maxvolchange = volchange
    return maxvolchange


MAXIMUM_VOLUME_CHANGE = get_maximum_volume_change()  # maximum volume change of the food vacuole


# TODO: sinnvolle Implementierung des occupancy factors
# models the occupancy of plasmepsin iv with the lysis of membrane proteins of HSVs
def get_occupancy_ratio(t):
    return 1 # -get_volume_change(t)/MAXIMUM_VOLUME_CHANGE



