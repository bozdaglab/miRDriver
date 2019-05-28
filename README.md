# miRDriver
miRDriver identifies frequently aberrated regions using GISTIC and computes differentially expressed genes between the two patient groups as frequently aberrated group and not frequently aberrated group. Utilizing copy number aberration, DNA methylation, gene expression, transcription factor expression and miRNA expression datasets, miRDriver applies a LASSO-based method to select miRNAs-target gene interactions. We tested miRDriver on breast cancer and ovarian cancer data from the Cancer Genome Atlas (TCGA) database. 

The order of running the code are as given below

1. createTransGenePredictorFile.R (Helper file: createTransGenePredictorFileHelper.R)
2. runLassoParallelForTransGene.R (Helper file: runLassoParallelForTransGeneHelper.R)
   
Folder structure required is shown below. Also created in miRDriver.

OvarianLASSO, 
OvarianLASSO\DATAFRAMES,
OvarianLASSO\ERROR_FILES,
OvarianLASSO\gene_all_DEcis,
OvarianLASSO\gene_all_mirna,
OvarianLASSO\LassoMinCoeff,
OvarianLASSO\PipeTransGeneDEPredFile,

DATAFRAMES folder contains the following files
OV_All_methylation.Rda (Sample Data),
OV_filtered_mirna_0.01_30_gistic.Rda,
OV_LASSO_CNV.Rda (Sample Data),
RSeq_RPKM.Rda (Sample Data),
TF_in_Rseq.Rda,

gene_all_DEcis folder contains trans gene and cis gene dataframe files

gene_all_mirna folder contains trans gene and cis miRNA dataframe files

LassoMinCoeff, PipeTransGeneDEPredFile and ERROR_FILES are empty folders used in the code to store the files.
