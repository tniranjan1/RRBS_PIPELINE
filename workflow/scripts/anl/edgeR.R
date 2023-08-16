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

# get sample names to remove
removal <- scan(remove_samples, what=character(), sep="\n")
removal <- grep('^#', removal, value=T, invert=T)
removal <- as.vector(sapply(removal, function(x) { strsplit(x, "\t")[[1]][1] }))
removal <- gsub('s_', '', removal)
removal <- grep(tissue, removal, value=T)
for(n in names(comp_shortcut)) removal <- gsub(paste('^', n, '_', sep=''), '', removal)
removal <- gsub('^Ctrl_', '', removal)
removal <- gsub(paste('_', tissue, '$', sep=''), '', removal)

# process covariate data
covariates <- read.table(covariate_path, stringsAsFactors=F, header=T)
covariates <- subset(covariates, Tissue==tissue)
## remove outlier samples
for(sample in removal) covariates <- subset(covariates, SampleID != sample)

# sample names
control_samples <- subset(covariates, SampleGroup=='Ctrl')$SampleID
comp_samples <- sapply(comp_shortcut[[comparison_group]], function(x) {subset(covariates, SampleGroup==x)$SampleID})
comp_samples <- as.vector(unlist(comp_samples))

# split covariate data
rownames(covariates) <- covariates$SampleID
control_covariates <- covariates[control_samples,]
comparison_covariates <- covariates[comp_samples,]

# process deconvolution data
  if(do_deconvolution == 'with_deconvolution')
  {
    decon_df <- read.table(deconvolution_path, header=T, stringsAsFactors=F)
    decon_df <- subset(decon_df, grepl(tissue, group))
    decon_df <- subset(decon_df, method=='nnls' & references=='all' & CpGsites=='all')
    special_split <- function(x) { gsub(paste(decon_df[x,'group'], '.', sep=''), '', decon_df[x,'sample_name']) }
    cells <- c('Neuron', 'Oligodendrocyte', 'Other')
    rownames(decon_df) <- sapply(seq(nrow(decon_df)), special_split)
    control_covariates[control_samples, cells] <- decon_df[gsub('[_-]', '.', control_samples, perl=T), cells]
    comparison_covariates[comp_samples, cells] <- decon_df[gsub('[_-]', '.', comp_samples, perl=T), cells]
  }

# isolate methylkit file paths matching control and comparison group
control_files <- as.vector(sapply(control_samples, function(x) { grep(x, methylation_files, value=T) }))
comparison_files <- as.vector(sapply(comp_samples, function(x) { grep(x, methylation_files, value=T) }))
