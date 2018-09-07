import pandas as pd
import dataimport
import matplotlib.pyplot as plt

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
    def __init__(self, name, pfidold, kcatkm, maxenzymeabundance):
        self.name = name
        self.pfidOld = pfidold
        self.kCat_Km = kcatkm
        self.MAX_ENZYME_ABUNDANCE = maxenzymeabundance  # in fmol
        self.line = 0  # line in proteome table
        self.abundance = []  # relative protein concentration in 2hpi resolution



dfProt = pd.read_excel('proteome_plas.xlsx', sep='\t', header=1)  # dataframe for protein abundance
plas1 = Enzyme('Plasmepsin I', 'PF14_0076', 0.01, 6, 20 * 10**(-6))
# name, PFID (old), Km in mM, kCat in 1/s, maximum enzyme abundance in fmol
plas1.line = dataimport.det_dataframe_line(dfProt, plas1)     # 2 lines in dataframe only first is used
plas1.abundance = dataimport.get_abundance(dfProt, plas1)
plas2 = Enzyme('Plasmepsin II', 'PF14_0077', 0.04, 5, 15 * 10**(-6))
plas2.line = dataimport.det_dataframe_line(dfProt, plas2)
plas2.abundance = dataimport.get_abundance(dfProt, plas2)
plas4 = Enzyme('Plasmepsin IV', 'PF14_0075', 0.02, 15, 20 * 10**(-6))
plas4.line = dataimport.det_dataframe_line(dfProt, plas4)
plas4.abundance = dataimport.get_abundance(dfProt, plas4)
plasmepsine = [plas1, plas2, plas4]
falcipain2 = MassActionEnzyme('Falcipain II', 'PF11_0165',  1511 * 10**3, 10 * 10**(-6))  # 45000 aus anderem datensatz
# name, pfidold, kcat/km in 1/(h*M), maximum enzyme abundance in fmol
falcipain2.line = dataimport.det_dataframe_line(dfProt, falcipain2)
falcipain2.abundance = dataimport.get_abundance(dfProt, falcipain2)
hdp = MassActionEnzyme('Heme Detoxification Protein', 'PF14_0446', 4179 * 10**3, 10 * 10**(-6))
falcipain3 = Enzyme('Falcipain III', 'PF11_0162', 0.067, 0.022, 10 ** 10**(-6))
falcipaine = [falcipain2, falcipain3]
hap = Enzyme('Histo Aspartic Protease', 'PF11_0161', 0.0034, 0.0016, 10 * 10**(-6))
hap.line = dataimport.det_dataframe_line(dfProt, hap)
hap.abundance = dataimport.get_abundance(dfProt, hap)

falcilysin = Enzyme('Falcilysin', 'PF13_0322', 0.000317, 0.0469, 10 * 10**(-6))
falcilysin.line = dataimport.det_dataframe_line(dfProt, falcilysin)
falcilysin.abundance = dataimport.get_abundance(dfProt, falcilysin)

dpap = Enzyme('Dipeptidyl aminopeptidase', 'PF11_0174', 1.4, 98, 10 * 10**(-6))
aaap = Enzyme('Aminoacylproline aminopeptidase', 'PF14_0517', 0.97, 16, 10 * 10**(-6))
alap = Enzyme('Alanyl aminopeptidase', 'MAL13P1.56', 0.322, 8.815, 10 * 10**(-6))
alap.line = dataimport.det_dataframe_line(dfProt, alap)
alap.abundance = dataimport.get_abundance(dfProt, alap)

fig = plt.figure()
plot = fig.add_subplot(111)
plot.plot(range(2, 50, 2), hap.abundance, label='hap')
plot.plot(range(2, 50, 2), falcilysin.abundance, label='fln')
plot.legend()
plt.show()


