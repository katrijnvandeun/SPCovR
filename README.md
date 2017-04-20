# SPCovR
**Sparse Principal Covariates Regression**

This repository contains 1) scripts and data to reproduce the analyses described in the paper and 2) the MATLAB implementation of the sparse Principal Covariates regression procedure (Algorithm 1), including stability selection (Algorithm 2).

Reference:
Van Deun, K., Crompvoets, E.A.V., & Ceulemans, E. (2017). Sparse principal covariates regression. Manuscript submitted for publication.

To reproduce the analysis of the influenza data discussed in the paper, take the following steps.

## 1. CREATE DATA
  
-Gene expression data: Run the R script R/RSCRIPT_TIV_RMA.R. This retrieves the expression data from the online database and creates RMA pre-processed expression matrices for both seasons. The resulting pre-processed gene expression data are stored in DATA/TIVRMA.txt (2008 season) and DATA/TIVRMA_2007.txt (2007 season).  
-Antibody titers: To recreate the data files DATA/TIVtiters.txt and DATA/TIVtiters_2007.txt take the following steps:  
  * download the GSE29614_series_matrix.txt from https://www.ncbi.nlm.nih.gov/geo/ and copy the lines 45-50 into the data file  
  * download GSE29617_series_matrix.txt from https://www.ncbi.nlm.nih.gov/geo/ and copy the lines 45-50 into the data file  
  
These files contain the antibody titers for three influenza strains and measured just before and 28 days after vaccination

Because we use the baseline-corrected gene expression data obtained three days after vaccination, the day 0 and day 3 data have to be matched on their subject identifiers and the difference scores have to be calculated. These difference scores are standardized (=centering + scaling to unit variance) per variable. To reproduce these steps, run the following MATLAB scripts.  

 a. 2007 sample  
  * Script_TIVD3vsD0_2007.m creates TIVD3_2007_rev.mat and TIVD3_2007_rev.txt
  * ScriptHAI_TIVD28vsD0_2007.m creates TIVtiter2007.m and TIVtiter2007.txt  
  
 b. 2008 sample  
  * MATLAB/Script_TIVD3vsD0_2008.m creates DATA/TIVD3_rev.mat and DATA/TIVD3_rev.txt
  * MATLAB/ScriptHAI_TIVD28vsD0_2008.m creates DATA/TIVtiter.m and DATA/TIVtiter.txt

## 2. ANALYZE DATA + POST-PROCESS RESULTS

1. SPLS (for comparison) using R and the packages SGCCA and SPLS; also includes trying out SPCR (gives out-of-memory)  
  *R/Script_sgcca_spls.R*
2. Ordinary PCovR analysis using MATLAB and creation of the Figures 1 & 2 in the paper  
  *MATLAB/plot_PCovR.m*: this script requires two external matlab functions: fig.m and exportfig.m
	Available from: http://www.mathworks.com/matlabcentral/fileexchange/30736 and
	https://nl.mathworks.com/matlabcentral/fileexchange/727-exportfig
3. Sparse PCovR using MATLAB  
  *MATLAB/Script_SPCovRanalysis.m*: Calls different function that implement Algorithm 1 and Algorithm 2.


## 3. ANNOTATION OF SELECTED PROBE SETS

Annotate the probe-sets with non-zero weights resulting from SPCovR and SGCCA, SPLS
  * For SGCCA and SPLS: use DATA/GENEIDS_spls.txt and DATA/GENEIDS_sgcca.txt as input to AmiGO: http://amigo.geneontology.org/amigo
  * For SPCovR, retrieve the AFFYIDS within MATLAB/Script_SPCovRanalysis.m and convert these to the official gene symbols using DAVID (as the annotation file included in the data folder may be outdated): https://david.ncifcrf.gov/. 
Next submit the official gene symbols to AmiGO.
