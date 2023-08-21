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

# read methylation files
cat("Step 1: Reading files\n", file=output, append=TRUE)
DMR_input <- readBismark2DGE(c(control_files, comparison_files), sample.names=c(control_samples, comp_samples), readr=F)
# order the methylation dataset by chromosome and position
chrNames <- paste("chr", c(1:22,"X","Y"), sep="")
DMR_input$genes$Chr <- factor(DMR_input$genes$Chr, levels=chrNames)
o <- order(DMR_input$genes$Chr, DMR_input$genes$Locus)
DMR_input <- DMR_input[o,]

# reduce input (insufficient coverage or usable samples)
cat("Step 2: Reducing input\n", file=output, append=TRUE)
DMR_input <- DMR_input[!is.na(DMR_input$gene$Chr),]
methylation <- gl(2,1,ncol(DMR_input), labels=c("Me","Un"))
me <- DMR_input$counts[,methylation=="Me"]
un <- DMR_input$counts[,methylation=="Un"]
coverage <- me + un
case_v_control <- c(rep(0, times=length(control_samples)), rep(1, times=length(comp_samples)))
adequate_callrate_control <- rowSums(coverage[,case_v_control==0] > 4) >= 10
adequate_callrate_cases <- rowSums(coverage[,case_v_control==1] > 4) >= 10
adequate_callrate <- adequate_callrate_cases & adequate_callrate_control
hasBoth <- rowSums(me) > 0 & rowSums(un) > 0
diffme <- me/coverage
control_means <- rowMeans(diffme[,1:length(control_samples)], na.rm=TRUE)
comparison_means <- rowMeans(diffme[,(length(control_samples)+1):ncol(diffme)], na.rm=TRUE)
diffme <- abs(control_means - comparison_means)	> 0.025
DMR_input <- DMR_input[adequate_callrate & hasBoth & diffme,,keep.lib.sizes=FALSE]

# build experimental design
targets <- as.data.frame(cbind(Sample=c(control_samples, comp_samples),
                               Group =c(rep("Controls",times=length(control_samples)),
                                        rep("Cases", times=length(comp_samples)))
                                       ))
rownames(targets) <- targets$Sample
cov_names <- colnames(control_covariates)[min(grep('Covariate', colnames(control_covariates))):ncol(control_covariates)]
targets[control_samples, gsub('Covariate_', '', cov_names)] <- control_covariates[control_samples, cov_names]
targets[comp_samples, gsub('Covariate_', '', cov_names)] <- comparison_covariates[comp_samples, cov_names]
targets$Group = (targets$Group == "Cases") + 0
targets$Gender = (targets$Gender == "F") + 0
formula <- as.formula(paste('~0 +', paste(c(names(targets)[-c(1:2)], names(targets)[2]), collapse=' + ')))
designSL <- model.matrix(formula, data=targets)
design <- modelMatrixMeth(designSL)

# perform linear regression
empty_sites <- DMR_input$samples[,"lib.size"] > 0
y1 <- DMR_input[,empty_sites]
redesign <- design[empty_sites,]
if(sum(!(empty_sites)) > 0) redesign <- redesign[,-c(unique(round(which(DMR_input$samples[,"lib.size"] == 0)/2)))]
## estimate dispersion
y1 <- estimateDisp(y1, design=redesign, trend.method="none")
  ## fit model for LRT
    fit.LRT <- glmFit(y1, redesign)
  ## fit model for QLT
    fit.QLT <- glmQLFit(y1, redesign)
  ## run glm using likelihood ratio test (LRT)
    test.LRT <- glmLRT(fit.LRT)
  ## run glm using quasi-likelihood test (QLT)
    test.QLT <- glmQLFTest(fit.QLT)
