input.rrbs = #
input.deco = #

# Load deconvSeq library
if(!("deconvSeq" %in% installed.packages()[,"Package"]))
{
  library(devtools)
  install_github("rosedu1/deconvSeq", dependencies=TRUE)
}
library(deconvSeq)

# build methyl table for deconvo reference
f <- input.deco[1]
data <- read.table(f, sep="\t", header=F, stringsAsFactors=F)
data <- data[!duplicated(data),]
rownames(data) <- paste(data[,1], data[,2], sep="_")
data <- data[,1:3]
colnames(data) <- c('chr', 'start', 'end')
data[,'strand'] <- "+"
deco_samples <- c()
deco_celltype <- c()
for(f in input.deco)
{
  insert <- read.table(f, sep="\t", header=F, stringsAsFactors=F)
  insert <- insert[!duplicated(insert),]
  rownames(insert) <- paste(insert[,1], insert[,2], sep="_")
  insert <- insert[rownames(data),]
  insert <- insert[,4:5]
  insert <- cbind(rowSums(insert), insert)
  insert <- insert / insert[,1] * 100
  for(x in c(2 , 3)) insert[,x] <- as.integer(insert[,x])
  insert[,1] <- insert[,2] + insert[,3]
  f_name <- rev(strsplit(f, '/')[[1]])[1]
  f_name <- strsplit(f_name, '\\.')[[1]]
  sample <- f_name[3]
  celltype <- rev(strsplit(grep('NeuN', f_name, value=T), '-')[[1]])[1]
  f_names <- paste(c('coverage', 'numCs', 'numTs'), (length(deco_samples)+1), sep="")
  deco_samples <- c(deco_samples, sample)
  deco_celltype <- c(deco_celltype, celltype)
  colnames(insert) <- f_names
  data[,f_names] <- insert
}
deconvo_table <- data
design.deco <- data.frame(neg=ifelse(deco_celltype=="neg",1,0),pos=ifelse(deco_celltype=="pos",1,0))

# compute projection matrix from deconvolution reference
set.seed(1234)
b0.deco = getb0.biseq(deconvo_table, design.deco, sigg=NULL)
# get top signature CpGs that will most contribute to deconvolution
# identify top CpGs using F-statistic Bonferroni-adjusted p-values of less than 0.05
top_x1 = getx1.biseq(NB0="top_bonferroni", b0.deco, deconvo_table, deco_samples, c("neg", "pos"))

# build methyl table for rrbs data
f <- input.rrbs[1]
data <- read.table(f, sep="\t", header=F, stringsAsFactors=F)
data <- data[!duplicated(data),]
rownames(data) <- paste(data[,1], data[,2], sep="_")
data <- data[,1:3]
colnames(data) <- c('chr', 'start', 'end')
data[,'strand'] <- "+"
rrbs_samples <- c()
rrbs_tissue <- c()
rrbs_group <- c()

for(f in input.rrbs)
{
  insert <- read.table(f, sep="\t", header=F, stringsAsFactors=F)
  insert <- insert[!duplicated(insert),]
  rownames(insert) <- paste(insert[,1], insert[,2], sep="_")
  insert <- insert[rownames(data),]
  insert <- insert[,4:5]
  insert <- cbind(rowSums(insert), insert)
  insert <- insert / insert[,1] * 100
  for(x in c(2 , 3)) insert[,x] <- as.integer(insert[,x])
  insert[,1] <- insert[,2] + insert[,3]
  f_name <- rev(strsplit(f, '/')[[1]])[1]
  f_name <- strsplit(f_name, '\\.')[[1]]
  sample <- f_name[3]
  tissue <- f_name[2]
  group <- f_name[1]
  group <- rev(strsplit(group, '_')[[1]])[1]
  f_names <- paste(c('coverage', 'numCs', 'numTs'), (length(rrbs_samples)+1), sep="")
  rrbs_samples <- c(rrbs_samples, sample)
  rrbs_tissue <- c(rrbs_tissue, tissue)
  rrbs_group <- c(rrbs_group, group)
  colnames(insert) <- f_names
  data[,f_names] <- insert
}
rrbs_table <- data
rrbs_full_names <- paste(rrbs_group, rrbs_samples, rrbs_tissue, sep="_")

# perform deconvolution
result_x1.rrbs <- getx1.biseq(NB0='top_bonferroni', b0.deco, rrbs_table, rrbs_full_names, c("neg", "pos"))
fixed_groups <- sapply(rownames(result_x1.rrbs$x1), function(r) { paste(strsplit(r, '_')[[1]][c(1,3)], collapse="_") })
fixed_ctrls <- unique(grep('^Ctrl', fixed_groups, value=T))
fixed_cases <- unique(grep('^Ctrl', fixed_groups, value=T, invert=T))
p_vals <- c()
for(f in fixed_cases)
{
  case_val <- result_x1.rrbs$x1[names(fixed_groups[fixed_groups == f]),'beta.neg']
  g <- paste('Ctrl', strsplit(f, '_')[[1]][2], sep='_')
  ctrl_val <- result_x1.rrbs$x1[names(fixed_groups[fixed_groups == g]),'beta.neg']
  p_vals[f] <- t.test(case_val, ctrl_val)$p.val
}

# add deconvolution to covariates
covariates_file <- '/home/tejasvi/rrbs_cristina_scratch/scripts/covariates.txt'
covariates <- read.table(covariates_file, header=T, sep="\t", stringsAsFactors=F, comment.char="")
colnames(covariates) <- gsub('^X.', '#', colnames(covariates))
rownames(covariates) <- sapply(covariates[,1], function(f) { strsplit(f, '_')[[1]][2] })
# CER
here <- grep('_CER$', rownames(result_x1.rrbs$x1))
covariates[as.vector(sapply(rownames(result_x1.rrbs$x1)[here], function(f) { strsplit(f, '_')[[1]][2] })),'CER'] <- result_x1.rrbs$x1[here,'beta.pos']
# FCX
here <- grep('_FCX$', rownames(result_x1.rrbs$x1))
covariates[as.vector(sapply(rownames(result_x1.rrbs$x1)[here], function(f) { strsplit(f, '_')[[1]][2] })),'FCX'] <- result_x1.rrbs$x1[here,'beta.pos']
covariates_file <- '/home/tejasvi/rrbs_cristina_scratch/scripts/covariates_and_deconv.txt'
write.table(covariates, file=covariates_file, sep="\t", row.names=F, col.names=T, quote=F)
