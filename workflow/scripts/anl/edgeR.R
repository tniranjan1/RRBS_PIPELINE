# load libraries
library(parallel)
library(edgeR)

# load arguments
covariate_path <- snakemake@input[['covariate_path']]
deconvolution_path <- snakemake@input[['deconvolution_path']]
methylation_files <- snakemake@input[['methylation_files']]
control_group <- snakemake@params[['control']]
comparison_group <- snakemake@params[['comparison']]
output_file <- snakemake@output[['path']]

# load enhancer data
enhancer_gene <- snakemake@input[['enhancer_hs']]
enhancer_all <- snakemake@input[['enhancer_sp']]
load(enhancer_all)
load(enhancer_gene)
