from snakemake import shell
import os
import re
import sys

## This script will extract methylation from bwameth-aligned BAM files.
##   Extraction will be done in all three contexts [ 'CpG', 'CHG', 'CHH' ].
##   Extracted files will then be bgzipped.

#----------------------------------------------------------------------------------------------------------------------#

# Get input arguments from snakemake object
input = snakemake.input
output = snakemake.output
params = snakemake.params
wildcards = snakemake.wildcards
threads = snakemake.threads
log = snakemake.log[0]

# Print standard output and error to log file
f = open(log, "w")
sys.stderr = sys.stdout = f

#----------------------------------------------------------------------------------------------------------------------#

mKit = "--methylKit" if wildcards.suffix == "methylKit" else ""
repeats = "" if wildcards.repeats == "with_repeats" else "-l " + input.ir
prefix = wildcards.path + "/methylation_calls/samples/MethylDackel_" + wildcards.sample + "." + wildcards.repeats
MDkl = "MethylDackel extract -q " + str(params.mapq) + " -d " + str(params.min_cov) + " --CHG --CHH"
print("Command submitted:")
print(f"  {MDkl} -@ {threads} {repeats} {mKit} -o {prefix} {input.ref} {input.bam} > {log} 2> {log}")
shell("{MDkl} -@ {threads} {repeats} {mKit} -o {prefix} {input.ref} {input.bam} > {log} 2> {log}")
if params.gzip:
    shell("bgzip -@ {threads} -c {prefix}_CpG.{wildcards.suffix} > {prefix}_CpG.{wildcards.suffix}.gz")
    shell("bgzip -@ {threads} -c {prefix}_CHG.{wildcards.suffix} > {prefix}_CHG.{wildcards.suffix}.gz")
    shell("bgzip -@ {threads} -c {prefix}_CHH.{wildcards.suffix} > {prefix}_CHH.{wildcards.suffix}.gz")

#----------------------------------------------------------------------------------------------------------------------#

f.close() # Close log file
