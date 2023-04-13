# load libraries
library(parallel)
library(e1071)
library(ggplot2)
library(nnls)

# load arguments
args <- commandArgs(trailingOnly=TRUE)
threads <- as.integer(args[1])
rrbs_gse_merge_path <- args[2]
neuN_path <- args[3]
output_path <- args[4]

# get gse column types
gse_data <- read.csv(rrbs_gse_merge_path, sep="\t", header=T, nrow=10)
keepNA <- max(grep('FCX|CER', colnames(gse_data), perl=T))
colClasses_vec <- rep('numeric', times=ncol(gse_data))
colClasses_vec[1:keepNA] <- 'NULL'
colClasses_vec[(keepNA+1)] <- 'character'
colClasses_vec[(keepNA+2)] <- 'integer'
colClasses_vec[(keepNA+3)] <- 'integer'

# load gse data
gse_data <- read.csv(rrbs_gse_merge_path, sep="\t", header=T, colClasses=colClasses_vec)
colnames(gse_data) <- gsub('\\.1$', '', colnames(gse_data), perl=T)

# get rrbs column types
rrbs_data <- read.csv(rrbs_gse_merge_path, sep="\t", header=T, nrow=10)
keepNA <- max(grep('FCX|CER', colnames(rrbs_data), perl=T))
colClasses_vec <- rep('numeric', times=ncol(rrbs_data))
colClasses_vec[1] <- 'character'
colClasses_vec[2] <- 'integer'
colClasses_vec[3] <- 'integer'
colClasses_vec[(keepNA+1):length(colClasses_vec)] <- 'NULL'

# load rrbs data
rrbs_data <- read.csv(rrbs_gse_merge_path, sep="\t", header=T, colClasses=colClasses_vec)

# remove duplicates
are_duplicates <- duplicated(paste(gse_data[,1], gse_data[,2], sep='_'))
gse_data <- gse_data[!are_duplicates,]
rrbs_data <- rrbs_data[!are_duplicates,]

# load neuN ref data
neuN_ref <- read.csv(neuN_path, sep="\t", header=F, skip=1)

# remove neuN duplicates
are_duplicates <- duplicated(paste(neuN_ref[,1], neuN_ref[,2], sep='_'))
neuN_ref <- neuN_ref[!are_duplicates,]

# find and restrict to common loci between rrbs/gse data and neuN_ref
common <- intersect(paste(gse_data[,1], gse_data[,2], sep='_'), paste(neuN_ref[,1], neuN_ref[,2], sep='_'))
present <- paste(gse_data[,1], gse_data[,2], sep='_') %in% common
gse_data <- gse_data[present,]
rrbs_data <- rrbs_data[present,]
present <- paste(neuN_ref[,1], neuN_ref[,2], sep='_') %in% common
neuN_ref <- neuN_ref[present,]

# fix column and row names in neuN ref
colnames(neuN_ref) <- c('chrom', 'start', 'end', paste(rep(c('x_neg', 'x_neuPos'), times=6), rep(1:6, each=2), sep='.'))
rownames(neuN_ref) <- rownames(gse_data)

# adjust methylation values to 0-1 range
neuN_ref[,4:15] <- neuN_ref[,4:15] / 100

# bind all reference data and remove redundant tables
data <- cbind(gse_data, neuN_ref[,4:15])
rm(gse_data, neuN_ref)
gc(verbose=F)

# identify categorical columns
  ## neuron columns
  neuron_cols <- grep('neu', colnames(data), ignore.case=T)
  non_neuron_cols <- grep('neu', colnames(data), ignore.case=T, invert=T)
  non_neuron_cols <- non_neuron_cols[4:length(non_neuron_cols)]

	## oligodendrocyte columns
	oligo_cols <- grep('oligo', colnames(data), ignore.case=T)
	non_oligo_cols <- grep('oligo', colnames(data), ignore.case=T, invert=T)
	non_oligo_cols <- non_oligo_cols[4:length(non_oligo_cols)]

	## non-neuron, non-oligodendrocyte columns
	other_groups <- colnames(data)[setdiff(non_neuron_cols, oligo_cols)]
	other_groups <- sapply(other_groups, function(o) { paste(strsplit(o, '_')[[1]][-c(1)], collapse='_') })
	other_groups <- sapply(other_groups, function(o) { paste(rev(rev(strsplit(o, '\\.')[[1]])[-c(1)]), collapse='.') })

# build signature methylation (meth_sign) table per cell type
neuron_mean <- apply(data[,neuron_cols], 1, mean, na.rm=T)
oligo_mean <- apply(data[,oligo_cols], 1, mean, na.rm=T)
meth_sign <- data.frame('Neuron'=neuron_mean, 'Oligodendrocyte'=oligo_mean)

# build reduced signature methylation table (non-neuron. non-oligo cell types collapsed to one 'Other' column)
reduced_meth_sign <- meth_sign
reduced_meth_sign[,'Other'] <- apply(data[,setdiff(non_neuron_cols, oligo_cols)], 1, mean, na.rm=T)

# build mean methylation signatures for each non-neuron/non-oligo cell to add to main meth_sign
#' Mean methylation signature for a cell type
#' @description for a given cell type, gather the mean methylation per locus for use as a signature
#' T Niranjan, tejasvi.niranjan@uantwerpen.vib.be
#' Apr 13, 2023
#'
#' Parallelizable function to get cell-type specific methylation means for all available loci
#' Assumes that unique cell type names have already been applied to the reference dataset (data) above.
#'
#' mean_meth <- other_mean(o)
#'
#' @param o (string) name of a cell type
#'
#' @return  (numeric) vector of mean methylations at each locus for a given cell type "o"
other_mean <- function(o) { apply(as.matrix(data[,names(other_groups)[other_groups==o]]), 1, mean, na.rm=T) }
add_to_meth <- mclapply(unique(other_groups), other_mean, mc.cores=threads)
names(add_to_meth) <- unique(other_groups)
for(o in unique(other_groups)) meth_sign[,o] <- add_to_meth[[o]]
rm(add_to_meth)

# isolate the 'most useful' methylation loci (greatest variance between neurons, oligo, and others)
most_useful <- rep(F, times=nrow(reduced_meth_sign))
for(l in list(neuron_cols, oligo_cols, setdiff(non_neuron_cols, oligo_cols)))
{
  maximum <- apply(as.matrix(data[,setdiff(seq(ncol(data)), l)[-c(1:3)]]), 1, max, na.rm=T)
	minimum <- apply(as.matrix(data[,setdiff(seq(ncol(data)), l)[-c(1:3)]]), 1, min, na.rm=T)
  top <- apply(as.matrix(data[,l]), 1, min, na.rm=T) > maximum
  bot <- apply(as.matrix(data[,l]), 1, max, na.rm=T) < minimum
  most_useful <- most_useful | top
  most_useful <- most_useful | bot
}
rm(maximum, minimum, top, bot)

#' NNLS Deconvolution
#' @description Use nnls to estimate the cell count percentage
#' T Niranjan, tejasvi.niranjan@uantwerpen.vib.be
#' Apr 13, 2023
#'
#' Applies non-negative least squares from r-nnls package to estimate cell type proportions.
#' Assumes that the methylation signature matrix and rrbs methylation data are already generated above.
#'
#' proportions_df <- nnls_cell_prop(i, useful)
#'
#' @param i      (integer) index value of column in rrbs data to analyze for cell proportion
#' @param useful (boolean) if analysis should be restricted to 'most_useful' methylation loci
#' @param ret    (boolean) if instead of performing analysis, return length of index column (default = F)
#'
#' @return a data frame with cell proportions for the indexed sample
nnls_cell_prop <- function(i, useful, ret = F)
{
	# hold proportion data
  prop_holder <- data.frame(sample_name=character(), group=character())
	# isolate sample name and sample group
	sample_name <- colnames(rrbs_data)[i]
	sample_group <- strsplit(sample_name, '\\.')[[1]]
	sample_group <- paste(sample_group[1:grep('FCX|CER', sample_group, perl=T)], collapse='.')
	# isolate sample methylation, removing NAs
	sample_meth <- if(useful) { rrbs_data[most_useful,i] / 100 } else { rrbs_data[,i] / 100 }
	sample_meth <- sample_meth[!is.na(sample_meth)]
	# isolate reference methylation signature
	sig_matrix <- if(useful) { meth_sign[most_useful,] } else { meth_sign }
	sig_matrix <- sig_matrix[!is.na(sample_meth),]
	# get intersection of sample methylation (sample_meth) and methylation signature matrix (sig_matrix)
	sample_meth <- sample_meth[rowSums(is.na(sig_matrix))==0]
	sig_matrix <- sig_matrix[rowSums(is.na(sig_matrix))==0,]
	# perform cell type deconvolution using nnls
	if(ret)
	{
		return(length(sample_meth))
	} else {
		prop <- nnls(as.matrix(sig_matrix), sample_meth)$x
		prop_holder[sample_name, c('sample_name', 'group')] <- c(sample_name, sample_group)
		prop_holder[sample_name,colnames(meth_sign)] <- prop
		rm(sample_meth, sig_matrix)
		gc(verbose = F)
		return(prop_holder)
	}
}

# get cell type proportions using nnls using all cell types and only most useful loci
proportions <- data.frame(sample_name=character(), group=character())
prop_list <- mclapply(seq(4, ncol(rrbs)), nnls_cell_prop, useful = T, mc.cores=threads)
for(p in prop_list) proportions[rownames(p), colnames(p)] <- p
prop.allCellTypes.usefulLoci <- proportions

# get cell type proportions using nnls using all cell types and all loci
proportions <- data.frame(sample_name=character(), group=character())
prop_list <- mclapply(seq(4, ncol(rrbs)), nnls_cell_prop, useful = F, mc.cores=threads)
for(p in prop_list) proportions[rownames(p), colnames(p)] <- p
prop.allCellTypes.allLoci <- proportions

# switch meth_sign table from full to reduced; store main meth_sign in large_meth_sign variable
large_meth_sign <- meth_sign
meth_sign <- reduced_meth_sign

# get cell type proportions using nnls using reduced cell types and only most useful loci
proportions <- data.frame(sample_name=character(), group=character())
prop_list <- mclapply(seq(4, ncol(rrbs)), nnls_cell_prop, useful = T, mc.cores=threads)
for(p in prop_list) proportions[rownames(p), colnames(p)] <- p
prop.reducedCellTypes.usefulLoci <- proportions

# get cell type proportions using nnls using reduced cell types and all loci
proportions <- data.frame(sample_name=character(), group=character())
prop_list <- mclapply(seq(4, ncol(rrbs)), nnls_cell_prop, useful = F, mc.cores=threads)
for(p in prop_list) proportions[rownames(p), colnames(p)] <- p
prop.reducedCellTypes.allLoci <- proportions

#' Support vector machine deconvolution
#' @description Use SVMDECONV to estimate the cell count percentage
#' David L Gibbs, dgibbs@systemsbiology.org
#' June 9, 2017
#' Adapted for multiple nu by T Niranjan, tejasvi.niranjan@uantwerpen.vib.be
#' Apr 13, 2023
#'
#' v-SVR is applied with a linear kernel to solve for f,
#' and the best result from multiple values of v = seq(19)/20
#' is saved, where ‘best’ is defined as the lowest root mean squared error
#' between m and the deconvolution result, f x B.
#'
#' Our current implementation executes v-SVR using the
#' ‘svm’ function in the R package, ‘e1071’.
#'
#' w2 <- SVMDECON(m, B)
#'
#' @param m  a matrix represenging the mixture (genes X 1 sample)
#' @param B  a matrix representing the references (genes X cells), m should be subset to match B
#'
#' @return A matrix with cell type estimates for each samples
SVMDECON <- function(m,B) {
  # multiple models are fit with different values of nu
  fits <- lapply(seq(19)/20, function(k) { e1071::svm(m~B, nu=k, kernel="linear", scale=T, type="nu-regression") })
  # these w's are the cell fractions
  ws <- lapply(fits, function(f) { weightNorm(t(f$coefs) %*% f$SV) })
  # return the model with the smallest mean sq error, adjust for NAs
  errs <- lapply(ws, function(w) { sqrt( sum( (m - B %*% t(w))^2 , na.rm=T)/nrow(m) ) + sum(is.na(w)) })
  resIdx <- which.min(unlist(errs))
  return(ws[[resIdx]])
}

#' SVMDECONV helper function
#' @description Use weightNorm to normalize the SVM weights.  Used for SVMDECONV
#'
#' w1 <- weightNorm(w)
#'
#' @param w  The weight vector from fitting an SVM, something like something like t(fit1$coefs) \%*\% fit1$SV, where fit comes from  <- svm(m~B, nu=0.25, kernel="linear"))
#' @return a weight vector
weightNorm <- function(w) { w[w<0] <- 0; return(w/sum(w)) }

#' SVR Deconvolution
#' @description Use svr to estimate the cell count percentage
#' T Niranjan, tejasvi.niranjan@uantwerpen.vib.be
#' Apr 13, 2023
#'
#' Applies support vector regression from SVMDECON package to estimate cell type proportions.
#' Assumes that the methylation signature matrix and rrbs methylation data are already generated above.
#' Assumes that the methylation signature matrix is already reduced to neuron, oligodendrocyte, and other types.
#' Assumes that only 'most_useful' loci will be used. The reduced signatures and useful loci improves speed.
#'
#' proportions_df <- nnls_cell_prop(i)
#'
#' @param i      (integer) index value of column in rrbs data to analyze for cell proportion
#'
#' @return a data frame with cell proportions for the indexed sample
svr_cell_prop <- function(i)
{
	# hold proportion data
  prop_holder <- data.frame(sample_name=character(), group=character())
	# isolate sample name and sample group
  sample_name <- colnames(rrbs_data)[i]
  sample_group <- strsplit(sample_name, '\\.')[[1]]
  sample_group <- paste(sample_group[1:grep('FCX|CER', sample_group, perl=T)], collapse='.')
  sample_meth <- rrbs_data[most_useful,i] / 100
  sig_matrix <- as.matrix(meth_sign[most_useful,])[!is.na(sample_meth),]
  sample_meth <- sample_meth[!is.na(sample_meth)]
  sample_meth <- as.matrix(sample_meth)
  ctp <- SVMDECON(sample_meth, sig_matrix)
  prop_holder[sample_name, c('sample_name', 'group')] <- c(sample_name, sample_group)
  prop_holder[sample_name,gsub('^B', '', colnames(ctp))] <- ctp
  rm(sample_meth, sig_matrix)
  gc()
  return(prop_holder)
}




sample_name <- colnames(rrbs_data)[i]
sample_group <- strsplit(sample_name, '\\.')[[1]]
sample_group <- paste(sample_group[1:grep('FCX|CER', sample_group, perl=T)], collapse='.')
# isolate sample methylation, removing NAs
sample_meth <- if(useful) { rrbs_data[most_useful,i] / 100 } else { rrbs_data[,i] / 100 }
sample_meth <- sample_meth[!is.na(sample_meth)]
# isolate reference methylation signature
sig_matrix <- if(useful) { meth_sign[most_useful,] } else { meth_sign }
sig_matrix <- sig_matrix[!is.na(sample_meth),]
# get intersection of sample methylation (sample_meth) and methylation signature matrix (sig_matrix)
sample_meth <- sample_meth[rowSums(is.na(sig_matrix))==0]
sig_matrix <- sig_matrix[rowSums(is.na(sig_matrix))==0,]
# perform cell type deconvolution using nnls





#
