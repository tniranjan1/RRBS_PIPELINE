# load libraries
library(methylKit)

# load arguments
input_path <- snakemake@input[['path']]
output_path <- snakemake@input[['path']]
sample_name <- snakemake@wildcards[['sample']]
context <- snakemake@wildcards[['context']]

# read methylation data from input_path file
rawMethylObj <- methRead(location=input_path,
                         sample.id=sample_name,
                         assembly="hg38",
                         context=context,
                         mincov=1,
                         pipeline='bismarkCoverage')

# tile the methylation data
tiles <- tileMethylCounts(rawMethylObj, win.size=500, step.size=250, cov.bases=5)

# read the header line of the input_path file
gz_file <- gzfile(input_path, "r")
header <- readLines(gz_file, n = 1)
close(gz_file)

# organize the tiled methylation data into bismarkCoverage format
out <- data.frame(chr=tiles$chr,
                  start=tiles$start,
                  end=tiles$end,
                  meth=as.integer(tiles$numCs / tiles$coverage * 100),
                  C=tiles$numCs,
                  T=tiles$numTs)

# write the header and tiled data to file
gz_file <- gzfile(output_path, "w")
cat(header, file = gz_file, sep = "\n")
write.table(out, gz_file, sep = "\t", row.names=F, col.names=F, quote=F, append=T)
close(gz_file)
