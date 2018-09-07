import pandas as pd
import dataimport


class Enzyme(object):
    def __init__(self, name, pfidold, kcatkm, maxenzymeabundance, maxsize, minsize, stepsize, endbased):

            self.name = name
            self.pfidOld = pfidold  # pfid of reference enzyme used for abundance
            self.kCatKm = kcatkm  # in 1/(M*s)
            self.line = 0  # line in proteome table
            self.abundance = []  # relative protein concentration in 2hpi resolution
            self.MAX_ENZYME_ABUNDANCE = maxenzymeabundance  # in fmol
            self.maxSize = maxsize  # maximum protein length for cleavage
            self.minSize = minsize  # minimum protein length for cleavage
            self.stepSize = stepsize  # interval for cleavage sites
            self.endBased = endbased  # works only on the n terminus of the protein


dfProt = pd.read_excel('proteome_plas.xlsx', sep='\t', header=1)  # dataframe for protein abundance
plas1 = Enzyme('Plasmepsin I', 'PF14_0076', 600 * 10**3, 20 * 10**(-6), 0, 0, 0, False)
plas1.line = dataimport.det_dataframe_line(dfProt, plas1)     # 2 lines in dataframe only first is used
plas1.abundance = dataimport.get_abundance(dfProt, plas1)
ipep = Enzyme('Initial Peptidases', 'PF14_0077', 500 * 10**3, 30 * 10**(-6), 146, 80, 20)
# plas1 used for abundance values (theory = initial concentration spikes of plas4 are for membrane degradation)
ipep.line = dataimport.det_dataframe_line(dfProt, ipep)
ipep.abundance = dataimport.get_abundance(dfProt, ipep)
fal2 = Enzyme('Falcipain II', 'PF11_0165',  1511 * 10**3, 10 * 10**(-6), 80, 20, 10, False)
fal2.line = dataimport.det_dataframe_line(dfProt, fal2)
fal2.abundance = dataimport.get_abundance(dfProt, fal2)
hdp = Enzyme('Heme Detoxification Protein', 'PF14_0446', 4179 * 10**3, 10 * 10**(-6), 0, 0, 0, False)
lppep = Enzyme('Long Peptide Peptidases', 'PF11_0161', 0.4 * 10**3, )
lppep.line = dataimport.det_dataframe_line(dfProt, lppep)
lppep.abundance = dataimport.get_abundance(dfProt, lppep)
# hap is used for abundance values (no data on fal3)
fln = Enzyme('Falcilysin', 'PF13_0322', 0.000317, 0.0469, 10 * 10**(-6))
fln.line = dataimport.det_dataframe_line(dfProt, falcilysin)
fln.abundance = dataimport.get_abundance(dfProt, falcilysin)

