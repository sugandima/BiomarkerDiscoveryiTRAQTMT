################################################################################################
# Author: S.M. Vidanagamachchi
# Function: main()
#           Run several functions to find biomarker discovery of labeled (iTRAQ/TMT) experiments
#
# Output:
#         A set of proteins as biomarkers
##############################################################################################
source("trypsinDigestion.R")
source("massToCharge.R")
source("Quantification.R")
source("searching.R")
source("proteinIdentification_ratios.R")
source("proteinIdentificationSeperateRuns.R")
source("selectMarkers_wTest.R")

trypsinDigestion("RefSeqsequence.fasta")
massToCharge(2.7,"cystein","lysine")# 2 static modifications
#massToCharge(2.7,"","lysine")# 1 static modification
Quantification("data")
searching("out",22,20)#identify peptides and reports peptide ratios
proteinIdentification_ratios(22,20)#identify unique proteins from all fractions and reports 3 different protein ratios (for comparison)
proteinIdentificationSeperateRuns(22,20)#identify proteins in each fraction seperately and reports intensities (and ratios) that can be used in p-value calculations of unique proteins 
selectMarkers_wTest(22,20)#calculate p-values for each protein
