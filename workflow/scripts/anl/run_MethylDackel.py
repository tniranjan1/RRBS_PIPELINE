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
repeats = "-l " + input.ir if wildcards.repeats == "without_repeats" else ""
<<<<<<< HEAD
prefix = wildcards.path + "/methylation_calls/samples/MethylDackel_" + wildcards.sample + "." + wildcards.repeats
=======
prefix = wildcards.path + "/" + wildcards.sample + "." + wildcards.repeats
>>>>>>> 3a03b8e5e3705f5b325e67b7624b8d32cea8e3c6
MDkl = "MethylDackel extract -q 20 -d 5 --CHG --CHH"
print("Command submitted:")
print(f"  {MDkl} -@ {threads} {repeats} {mKit} -o {prefix} {input.ref} {input.bam} > {log} 2> {log}")
shell("{MDkl} -@ {threads} {repeats} {mKit} -o {prefix} {input.ref} {input.bam} > {log} 2> {log}")
shell("bgzip -@ {threads} -c {prefix}_CpG.{wildcards.suffix} > {prefix}_CpG.{wildcards.suffix}.gz 2>> {log}")
shell("bgzip -@ {threads} -c {prefix}_CHG.{wildcards.suffix} > {prefix}_CHG.{wildcards.suffix}.gz 2>> {log}")
shell("bgzip -@ {threads} -c {prefix}_CHH.{wildcards.suffix} > {prefix}_CHH.{wildcards.suffix}.gz 2>> {log}")

#----------------------------------------------------------------------------------------------------------------------#

f.close() # Close log file
