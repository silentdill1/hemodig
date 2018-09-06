import pandas as pd
from foodvacuole import get_volume
import matplotlib.pyplot as plt

df = pd.read_excel('hemozoin exp.xls', sep='\t', header=0)
unsortedPoints = []
for i in range(0, df.shape[0]):
    t = df.at[i, 't (h)']
    unsortedPoints.append([t, df.at[i, 'n (fmol)']])

def get_timepoint(point):
    return point[0]


unsortedPoints.sort(key=get_timepoint)
timepoints = []
hzAbundances = []
for point in unsortedPoints:
    timepoints.append(point[0])
    hzAbundances.append(point[1])
expValues = [timepoints, hzAbundances]

