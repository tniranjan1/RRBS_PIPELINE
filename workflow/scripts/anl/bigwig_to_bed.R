# Load arguments
args <- commandArgs(trailingOnly=TRUE)
dir <- args[1]
threads <- as.integer(args[2])
rrbs <- args[3]

# Load libraries
library(parallel)

# Load file list
files <- list.files(path=dir, pattern='hg38.bigwig', full.names=T)

# function to convert bigwig to bedgraph format
convert <- function(f)
{
  new_name <- gsub('bigwig', 'bed', f)
  mid_name <- paste(f, 'tmp', sep='.')
  c1 <- paste("bigWigToBedGraph", f, mid_name)
  c2 <- paste('cat', mid_name, "| awk '($4 >= 0)' >", new_name)
  c3 <- paste('rm', mid_name)
  c4 <- paste('bgzip', new_name)
  system(c1)
  system(c2)
  system(c3)
  system(c4)
  new_name <- paste(new_name, 'gz', sep='.')
  return(new_name)
}

#new_names <- mclapply(files, convert, mc.cores=threads)

# files variable should be same as new_names
files <- list.files(path=dir, pattern='hg38.bed.gz', full.names=T)
# names lack parent directories and extension
names <- list.files(path=dir, pattern='hg38.bed.gz', full.names=F)
names <- sapply(names, function(n) { strsplit(n, '.hg38')[[1]][1] })

# get union of bed files
tool <- 'bedtools unionbedg -i'
params <- '-header -filler NA -names'
out <- paste(dir, 'merge.bed', sep='/')
command <- paste(tool, paste(files, collapse=' '), params, paste(names, collapse=' '), '>', out)
#system(command)

# compress union file
input <- out
awk_command <- 'awk \'BEGIN{FS=OFS="\t"}{ if($2 != "start"){$2 = $2 + 1}; print}\''
compress <- paste('bgzip -@', threads, '-c')
out <- paste(dir, 'merge.bed.gz', sep='/')
command <- paste('cat', input, '|', awk_command, '|', compress, '>', out)
#system(command)

# remove unnecessary files
extensions <- c('merge.bed', '*.bigwig', '*.beta', '*.hg38.bed.gz', 'GSE186458_RAW.tar', '*.pat.gz*')
for(e in extensions)
{
  command <- paste('rm ', dir, '/', e, sep='')
  #system(command)
}

# separate header and position sort, split by chromosome
header <- paste(dir, 'header.txt', sep='/')
command <- paste('zcat', out, '| head -n 1 >', header)
#system(command)
command <- paste('zcat', out, '| tail -n +2 | cut -f1 | uniq | sort -k 1,1 | uniq')
#chromosomes <- system(command, intern=T)

# function to position sort, one thread per chromosome for speed-up
sort_by_chr <- function(chr)
{
  c1 <- paste('zcat', out)
  c2 <- paste("grep -P '^", chr, "\t'", sep='')
  new_out <- paste(dir, '/merge.', chr, '.sorted.bed', sep='')
  c3 <- paste('sort -k 1,1 -k2,2n > ', new_out, sep='')
  command <- paste(c1, c2, c3, sep=' | ')
  system(command)
  return(new_out)
}
#out_by_chr <- unlist(mclapply(chromosomes, sort_by_chr, mc.cores=threads))

# bind chromosome-split sorted files
out <- paste(dir, 'merge.sorted.bed.gz', sep='/')
#command <- paste('cat', header, paste(out_by_chr, collapse=' '), '|', compress, '>', out)
#system(command)
command <- paste('rm', header)
#system(command)

# merge with rrbs methylation file
bedtools_command <- paste('bedtools intersect -sorted -wb -a stdin -b', out)
intersect <- paste(dir, 'merge.intersect.bed', sep='/')
command <- paste('zcat', rrbs, '| tail -n +2 | sort -k 1,1 -k2,2n |', bedtools_command, '>', intersect)
#system(command)

# build merged header
#command <- paste('zcat', rrbs, '| head -n 1')
#header <- strsplit(system(command, intern=T), '\t')[[1]]
#command <- paste('zcat', paste(dir, 'merge.bed.gz', sep='/'), '| head -n 1')
#header <- c(header, strsplit(system(command, intern=T), '\t')[[1]])
#header <- paste(header, collapse="\t")
#header <- paste(header, '\n', sep='')
#cat(header, file=paste(dir, 'header.txt', sep='/'))

# combined header and bed file, then compress
out <- paste(dir, 'merge.intersect.bed.gz', sep='/')
command <- paste('cat', paste(dir, 'header.txt', sep='/'), intersect, '|', compress, '>', out)
system(command)

# remove unnecessary files
command <- paste('rm', paste(dir, 'header.txt', sep='/'), intersect)
system(command)
