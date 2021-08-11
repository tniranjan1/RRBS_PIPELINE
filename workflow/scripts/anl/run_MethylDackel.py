from snakemake import shell
import os
import re
import sys

## This script will generate ...
##   ...

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
CHG = "--CHG" if len([ c for c in wildcards.context if c == "CHG" ]) == 1 else ""
CHH = "--CHH" if len([ c for c in wildcards.context if c == "CHH" ]) == 1 else ""
repeats = "-l " + input.ir if wildcards.repeats == "without_repeats" else ""
prefix = wildcards.path + "/methylation_calls/" + wildcards.sample + "." + wildcards.repeats
MDkl = "MethylDackel extract -q 20 -d 5"
shell("{MDkl} -@ {threads} {repeats} {CHG} {CHH} {mKit} -o {prefix} {input.ref} {input.bam} > {log} 2> {log}")

#----------------------------------------------------------------------------------------------------------------------#

f.close() # Close log file
