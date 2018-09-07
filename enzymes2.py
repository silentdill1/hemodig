import pandas as pd
import dataimport
import

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
iPep = Enzyme('Initial Peptidases', 'PF14_0077', 500 * 10**3, 30 * 10**(-6), 146, 80, 20)
# plas1 used for abundance values (theory = initial concentration spikes of plas4 are for membrane degradation)
iPep.line = dataimport.det_dataframe_line(dfProt, iPep)
iPep.abundance = dataimport.get_abundance(dfProt, iPep)
fal2 = Enzyme('Falcipain II', 'PF11_0165',  1511 * 10**3, 10 * 10**(-6), 80, 20, 10, False)
fal2.line = dataimport.det_dataframe_line(dfProt, fal2)
fal2.abundance = dataimport.get_abundance(dfProt, fal2)
hdp = Enzyme('Heme Detoxification Protein', 'PF14_0446', 4179 * 10**3, 10 * 10**(-6), 0, 0, 0, False)
lpPep = Enzyme('Long Peptide Peptidases', 'PF11_0161', 0.4 * 10**3, 20 * 10**(-6), 80, 20, 10, False)
lpPep.line = dataimport.det_dataframe_line(dfProt, lpPep)
lpPep.abundance = dataimport.get_abundance(dfProt, lpPep)
# hap is used for abundance values (no data on fal3)
fln = Enzyme('Falcilysin', 'PF13_0322', 0.000317, 0.0469, 10 * 10**(-6), 80, 20, 10, False)
fln.line = dataimport.det_dataframe_line(dfProt, fln)
fln.abundance = dataimport.get_abundance(dfProt, fln)
dpAp = Enzyme('DiPeptidyl aminopeptidase', 'PF11_0174', 70 * 10**3, 50 * 10**(-6), 20, 4, 2, True)
LAMAp = Enzyme('Leucyl, Aspartyl, Methionyl Aminopeptidasen', '', 1.5 * 10**3, 20, 2, 1, True)
AapAAp = Enzyme('Aminoacyl Prolin, Alanyl Aminopeptidasen', '', 22 * 10**3, 20, 2, 1, True)
PPAp = Enzyme('Prolyl Aminopetidase, Aminopeptidase P', '', 109 * 10**3, 20, 2, 1, True)
subt = Enzyme('Subtilisin', '', 440 * 10**3, 20, 2, 1, True)
enzymes = [iPep, lpPep, fln, dpAp, LAMAp, AapAAp, PPAp, subt]
