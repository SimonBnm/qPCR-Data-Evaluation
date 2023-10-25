# qPCR-Data-Evaluation
R-based script based on the 2−ΔΔCt method for automated qPCR data evaluation and visualisation.


Allows for automated evaluation of the relative changes in gene expression data from real-time, quantitative PCR experiments using the 2−ΔΔCt method. 
The standard script provided is only suitable for quantification of 2 groups (Control vs Experimental = Treatment Group), but can be customised with simple adjustments.
Excel files can be used as input data. They must contain the following columns as shown in the example data: Target (=gene), Sample (=Group, Replicate), Cq (Ct values). 
To exclude specific files from the analysis folder, add a '#' before the file name.


Before running the script, certain parameters need to be set:

1) housekeeper <- Name of the houskeeper gene (indicated in "Target" column)
2) cont_name <- Name of the Control Group (Group of samples to use as a calibrator/reference when calculating the ∆∆Ct, indicated in "Sample" column)
3) foldername <- Name of folder containing input data (Excel files)
4) target_number <- Specify how many genes (targets) were analysed, including the Housekeeper.
5) con_name <- Name of Control Group (shown in Plot)
6) exp_name <- Name of Experimental Group (shown in Plot)
7) target_remove <- add Targets (genes) that should not be included in the evaluation
8) stat_method <- Name of the statistical test (Wilcoxon rank-sum test as standard)
9) rename <- Set new target names (alphabetical order)

Successful execution of the script will result in the generation of an evaluation folder containing an Excel file with all calculations of the 2−ΔΔCt method and the results of the statistical test, 
as well as a visualisation of the data as a bar plot (pdf, svg or tif format).

If the number of target genes specified does not match the number of unique target names, the code will stop and an error message will be displayed.

