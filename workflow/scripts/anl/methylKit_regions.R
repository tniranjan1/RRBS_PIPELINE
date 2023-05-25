# load arguments
input_path <- snakemake@input[['path']]
output_path <- snakemake@output[['path']]
sample_name <- snakemake@wildcards[['sample']]
context <- snakemake@wildcards[['context']]
log_path <- snakemake@log[[1]]
log <- file(log_path, open="w")
sink(log, type = "output")
sink(log, type = "message", append=TRUE)

# load libraries
library(methylKit)

# Log input parameters (sanity check)
sprintf("Input path:  %s\n", input_path)
sprintf("Output path: %s\n", output_path)
sprintf("Sample name: %s\n", sample_name)
sprintf("Context:     %s\n", context)
sprintf("Log path:    %s\n", log_path)

# read methylation data from input_path file
cat("Reading methylation data from input path.\n")
rawMethylObj <- methRead(location=input_path,
                         sample.id=sample_name,
                         assembly="hg38",
                         context=context,
                         mincov=1,
                         pipeline='bismarkCoverage')

# tile the methylation data
cat("Tiling methylation data.\n")
tiles <- tileMethylCounts(rawMethylObj, win.size=500, step.size=250, cov.bases=5)

# read the header line of the input_path file
cat("Reading header from input path.")
gz_file <- gzfile(input_path, "r")
header <- readLines(gz_file, n = 1)
close(gz_file)

# organize the tiled methylation data into bismarkCoverage format
cat("Organizing tiled methylation data.")
out <- data.frame(chr=tiles$chr,
                  start=tiles$start,
                  end=tiles$end,
                  meth=as.integer(tiles$numCs / tiles$coverage * 100),
                  C=tiles$numCs,
                  T=tiles$numTs)

# write the header and tiled data to file
cat("Writing tiled methylation data to file.\n")
gz_file <- gzfile(output_path, "w")
cat(header, file = gz_file, sep = "\n")
write.table(out, gz_file, sep = "\t", row.names=F, col.names=F, quote=F, append=T)
close(gz_file)

sink(type = "output")
sink(type = "message")
invisible()
