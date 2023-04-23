# load libraries
library(parallel)
library(edgeR)

# load arguments
covariate_path <- snakemake@input[['covariate_path']]
