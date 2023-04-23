# load libraries
library(parallel)
library(edgeR)

# load arguments
covariate_path <- snakemake@input[['covariate_path']]
deconvolution_path <- snakemake@input[['deconvolution_path']]
methylation_files <- snakemake@input[['methylation_files']]
control_group <- snakemake@params[['control']]
output_file <- snakemake@output[['path']]
comparison_group <- %function of output file name %

# 
