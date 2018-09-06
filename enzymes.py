import pandas as pd
import dataimport


class Enzyme(object):
    def __init__(self, name, pfidold, km, kcat, maxenzymeabundance):

            self.name = name
            self.pfidOld = pfidold
            self.km = km  # in mM
            self.kCat = kcat  # in 1/s
            self.line = 0  # line in proteome table
            self.abundance = []  # relative protein concentration in 2hpi resolution
            self.MAX_ENZYME_ABUNDANCE = maxenzymeabundance  # in fmol


class MassActionEnzyme(object):
    def __init__(self, name, kcatkm, maxenzymeabundance):
        self.name = name
        self.kCat_Km = kcatkm
        self.MAX_ENZYME_ABUNDANCE = maxenzymeabundance  # in fmol


dfPlas = pd.read_excel('proteome_plas.xlsx', sep='\t', header=1)  # dataframe for plasmepsin protein abundance
plas1 = Enzyme('Plasmepsin I', 'PF14_0076', 0.01, 6, 20 * 10**(-6))
# name, PFID (old), Km in mM, kCat in 1/s, maximum enzyme abundance in fmol
plas1.line = dataimport.plas_det_dataframe_line(dfPlas, plas1)     # 2 lines in dataframe only first is used
plas1.abundance = dataimport.plas_get_abundance(dfPlas, plas1)
plas2 = Enzyme('Plasmepsin II', 'PF14_0077', 0.04, 5, 15 * 10**(-6))
plas2.line = dataimport.plas_det_dataframe_line(dfPlas, plas2)
plas2.abundance = dataimport.plas_get_abundance(dfPlas, plas2)
plas4 = Enzyme('Plasmepsin IV', 'PF14_0075', 0.02, 15, 20 * 10**(-6))
plas4.line = dataimport.plas_det_dataframe_line(dfPlas, plas4)
plas4.abundance = dataimport.plas_get_abundance(dfPlas, plas4)
plasmepsine = [plas1, plas2, plas4]
falcipain2 = MassActionEnzyme('Falcipain II', 1511 * 10**3, 10 * 10**(-6))
# name, kcat/km in 1/(h*M), maximum enzyme abundance in fmol
hdp = MassActionEnzyme("Heme Detoxification Protein", 4179 * 10**3, 10 * 10**(-6))


