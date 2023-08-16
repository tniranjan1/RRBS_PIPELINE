# load libraries
library(parallel)
library(edgeR)

# load arguments
covariate_path <- snakemake@input[['covariate_path']]
deconvolution_path <- snakemake@input[['deconvolution_path']]
methylation_files <- snakemake@input[['methylation_files']]
remove_samples <- snakemake@input[['remove_samples']]
tissue <- snakemake@wildcards[['tissue']]
control_group <- 'Ctrl'
comparison_group <- snakemake@wildcards[['comparison']]
do_deconvolution <- snakemake@wildcards[['deconvolution']] # with_deconvolution | without_deconvolution
output_file <- snakemake@output[['path']]

# load enhancer data
enhancer_gene <- snakemake@input[['enhancer_hs']]
enhancer_all <- snakemake@input[['enhancer_sp']]
load(enhancer_all)
load(enhancer_gene)

# long comparison name
comp_shortcut <- list(
  "A" = "FTLD-TDP43-A",
  "B" = "FTLD-TDP43-B",
  "C" = "FTLD-TDP43-C",
  "C9" = "FTLD-TDP43-C9-Positive",
  "F" = "FTLD-FUS",
  "G" = "FTLD-TDP43-GRN",
  "X" = c("FTLD-TDP43-A", "FTLD-TDP43-B", "FTLD-TDP43-C"),
  "Y" = c("FTLD-TDP43-A", "FTLD-TDP43-B", "FTLD-TDP43-C", "FTLD-TDP43-C9-Positive", "FTLD-TDP43-GRN")
)

# process covariate data
covariates <- read.table(covariate_path, stringsAsFactors=F, header=T)
covariates <- subset(covariates, Tissue==tissue)

# sample names
control_samples <- subset(covariates, SampleGroup=='Ctrl')$SampleID

# process deconvolution data
  # check do_deconvolution

# isolate methylkit file paths matching control

# isolate methylkit file paths matching comparison group
