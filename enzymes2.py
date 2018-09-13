import pandas as pd
import dataimport
import datainterpretation


class Enzyme(object):
    def __init__(self, name, category, pf_id_old, k_cat_k_m, max_enzyme_abundance, max_size, min_size, step_size, end_based):

            self.name = name
            self.category = category  # see types list
            self.pfidOld = pf_id_old  # pfid of reference enzyme used for abundance
            self.kCatKm = k_cat_k_m  # in 1/(M*s)
            self.line = 0  # line in proteome table
            self.abundance = []  # relative protein concentration in 2hpi resolution
            self.MAX_ENZYME_ABUNDANCE = max_enzyme_abundance  # in fmol
            self.maxSize = max_size  # maximum protein length for cleavage
            self.minSize = min_size  # minimum protein length for cleavage
            self.stepSize = step_size  # interval for cleavage sites
            self.endBased = end_based  # works only on the n terminus of the protein
            if category != -1:
                self.indices = datainterpretation.get_index_list(self)
            # list with start and end index for species wFpp followed by woFpp for ode solver


types = ['Initial peptidase', 'Long peptide peptidase', 'Dipeptidyl aminopeptidase', 'Dipeptid aminopeptidase']
# -1 is for other types that should not be processed by a function requiring types
plas1Concentration = 30 * 10**(-6)
iPepConcentration = 30 * 10**(-6)
fal2Concentration = 50 * 10**(-6)
hdpConcentration = 50 * 10**(-6)
lpPepConcentration = 80 * 10**(-6)
flnConcentration = 80 * 10**(-6)
dpapConcentration = 100 * 10**(-6)
apConcentration = 80 * 10**(-6)


dfProt = pd.read_excel('proteome_plas.xlsx', sep='\t', header=1)  # dataframe for protein abundance
plas1 = Enzyme('Plasmepsin I', -1, 'PF14_0076', 600 * 10**3, plas1Concentration, 0, 0, 0, False)
plas1.line = dataimport.det_data_frame_line(dfProt, plas1)     # 2 lines in dataframe only first is used
plas1.abundance = dataimport.get_abundance(dfProt, plas1)
iPep = Enzyme('Initial peptidases', 0, 'PF14_0077', 500 * 10**3, 2 * iPepConcentration, 146, 80, 20, False)
# plas1 used for abundance values (theory = initial concentration spikes of plas4 are for membrane degradation)
iPep.line = dataimport.det_data_frame_line(dfProt, iPep)
iPep.abundance = dataimport.get_abundance(dfProt, iPep)
fal2 = Enzyme('Falcipain II', 1, 'PF11_0165',  1511 * 10**3, fal2Concentration, 80, 20, 10, False)
fal2.line = dataimport.det_data_frame_line(dfProt, fal2)
fal2.abundance = dataimport.get_abundance(dfProt, fal2)
hdp = Enzyme('Heme Detoxification Protein', -1, 'PF14_0446', 4179 * 10**3, hdpConcentration, 0, 0, 0, False)
lpPep = Enzyme('Long Peptide Peptidases', 1, 'PF11_0161', 0.4 * 10**3, 2 * lpPepConcentration, 80, 30, 10, False)
lpPep.line = dataimport.det_data_frame_line(dfProt, lpPep)
lpPep.abundance = dataimport.get_abundance(dfProt, lpPep)
# hap is used for abundance values (no data on fal3)
fln = Enzyme('Falcilysin', 2, 'PF13_0322', 148 * 10**3, flnConcentration, 26, 8, 4, True)
fln.line = dataimport.det_data_frame_line(dfProt, fln)
fln.abundance = dataimport.get_abundance(dfProt, fln)
dpAp = Enzyme('Dipeptidyl aminopeptidase', 2, 'PF11_0174', 70 * 10**3, dpapConcentration, 20, 4, 2, True)
LAMAp = Enzyme('Leucyl, Aspartyl, Methionyl Aminopeptidasen', 3, '', 1.5 * 10**3, 3 * apConcentration, 20, 2, 1, True)
AapAAp = Enzyme('Aminoacyl Prolin, Alanyl Aminopeptidasen', 3, '', 22 * 10**3, 2 * apConcentration, 20, 2, 1, True)
PPAp = Enzyme('Prolyl Aminopetidase, Aminopeptidase P', 3, '', 109 * 10**3, 2 * apConcentration, 20, 2, 1, True)
subt = Enzyme('Subtilisin', 3, '', 440 * 10**3, apConcentration, 20, 2, 1, True)
enzymes = [iPep, fal2, lpPep, fln, dpAp, LAMAp, AapAAp, PPAp, subt]
