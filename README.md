# SPCovR
Sparse Principal Covariates Regression

This repository contains 1) scripts and data to reproduce the analyses described in the paper and 2) the MATLAB implementation of the sparse Principal Covariates regression procedure (Algorithm 1), including stability selection (Algorithm 2).

To reproduce the analysis of the influenza data discussed in the paper, follow the following steps.

  1. CREATE DATA
  Gene expression data: Run the R script R/RSCRIPT_TIV_RMA.R. This creates RMA pre-processed expression matrices for both seasons that are stored in DATA/TIVRMA.txt (2008 season) and DATA/TIVRMA_2007.txt (2007 season)
  Antibody titers: To recreate the data files DATA/TIVtiters.txt & DATA/TIVtiters_2007.txt take the following steps:
 download the GSE29614_series_matrix.txt from https://www.ncbi.nlm.nih.gov/geo/ and copy the lines 45-50 into the data file
 download GSE29617_series_matrix.txt from https://www.ncbi.nlm.nih.gov/geo/ and copy the lines 45-50 into the data file
 These are the antibody titers for three influenza strains and measured just before and 3 and 7 days after vaccination

2. PRE-PROCESS DATA
Further pre-process the expression data and titers using MatLab
 a. 2007 sample
 REVScript_makedata_TIV3vsD7_2007.m
 REVScriptHAI_TIV_2007.m
 b. 2008 sample
 REVScript_makedata_TIV3vsD7.m
 REVScriptHAI_TIV.m

4. analyze data
  1. SPLS (for comparison) using R
  Script_sgcca_spls.R
  2. Ordinary PCovR (for comparison) using Matlab
  Script_PCovR.m
  3. Sparse PCovR using Matlab
  Script_StabilitySelection.m


5. Annotate the probe-sets with non-zero weights resulting from SPCovR and SGCCA, SPLS
For SGCCA and SPLS: use GENEIDS_spls.txt and GENEIDS_sgcca.txt as input to AmiGO: http://amigo.geneontology.org/amigo
For SPCovR, retrieve the AFFYIDS and convert these to the official gene symbols using DAVID (as the annotation
file included in the data folder may be outdated): https://david.ncifcrf.gov/. 
Next submit the official gene symbols to AmiGO.
