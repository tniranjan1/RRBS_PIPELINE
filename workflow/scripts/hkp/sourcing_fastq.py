from snakemake import shell
import os
import re
import sys

## This script will generate or link the raw fastq files needed for alignment by bwa-meth.
##   If the fastq file source is from an SRA repository, this script will download the repository temporarily
##   and dump the fastq files. If the fastq file source is from a BAM file, this script will read sort the
##   the bam file, and generate the fastq files from the bam. If the fastq file source is already fastq format,
##   the fastq files will be hard-linked to other original files.

#----------------------------------------------------------------------------------------------------------------------#

# Get input arguments from snakemake object
output = snakemake.output
params = snakemake.params
wildcards = snakemake.wildcards
threads = snakemake.threads
log = snakemake.log[0]

# Print standard output and error to log file
f = open(log, "w")
sys.stderr = sys.stdout = f

#----------------------------------------------------------------------------------------------------------------------#

if params.type == 'SRR':
    # Get path to local SRR download repository
    sra_repository_command = "vdb-config -o n -p | grep /repository/user/main/public/root | cut -d= -f2 | sed 's/ //g'"
    sra_loc_repo = os.popen(sra_repository_command).read()
    sra_loc_repo = sra_loc_repo.replace('"', '').replace("'", "").strip("\n")
    if len(sra_loc_repo) == 0:
        sra_loc_repo = "."
    # Download SRR from SRA repository
    shell("prefetch --max-size 100000000 {params.path}")
    # multi-thread fastq dump from local SRR
    shell("mkdir -p ../resources/sra_temp")
    pfd="parallel-fastq-dump"
    shell("{pfd} --sra-id {params.path} --threads {threads} --outdir {sra_loc_repo}/{params.path} --split-files --gzip")
    # move fastq dump files to desired output paths
    shell("mv {sra_loc_repo}/{params.path}/{params.path}_1.fastq.gz {output.fq1}")
    shell("mv {sra_loc_repo}/{params.path}/{params.path}_2.fastq.gz {output.fq2}")
    # delete local SRR download
    shell("rm -r {sra_loc_repo}/{params.path}")

#----------------------------------------------------------------------------------------------------------------------#

elif params.type == 'bam':
    # Read sort bam file
    tmp_fldr = "../resources/bam_temp"
    shell("mkdir -p {tmp_fldr}")
    tmp_bam = f"{tmp_fldr}/{wildcards.sample}.bam"
    shell("samtools sort -@ {threads} -n -O BAM -o {tmp_bam} {params.path}")
    # Generate fastq files from bam file. Samtools will auto-gzip the fastq files due to ".gz" extension in output names
    shell("samtools fastq -@ {threads} -n -1 {output.fq1} -2 {output.fq2} -0 /dev/null -s /dev/null {tmp_bam}")
    shell("rm {tmp_bam}")

#----------------------------------------------------------------------------------------------------------------------#

elif params.type == 'fq':
    # Get fastq paths and determine if already gzipped
    fastq_paths = params.path
    fastq_paths = fastq_paths.split(",")
    fq1_input_path = fastq_paths[0]
    fq2_input_path = fastq_paths[1]
    fq1_output_path = output.fq1
    fq2_output_path = output.fq2
    if os.path.splitext(fq1_input_path)[1] != '.gz':
        fq1_output_path = fq1_output_path.strip(".gz")
    if os.path.splitext(fq2_input_path)[1] != '.gz':
        fq2_output_path = fq2_output_path.strip(".gz")
    # Link input fastq files to rule-established destination
    shell("ln {fq1_input_path} {fq1_output_path}")
    shell("ln {fq2_input_path} {fq2_output_path}")
    # Gzip files if not already gzipped
    if os.path.splitext(fq1_input_path)[1] != '.gz':
        shell("bgzip -@ {threads} -c {fq1_output_path} > {output.fq1}")
        shell("rm {fq1_output_path}")
    if os.path.splitext(fq2_input_path)[1] != '.gz':
        shell("bgzip -@ {threads} -c {fq2_output_path} > {output.fq2}")
        shell("rm {fq2_output_path}")

#----------------------------------------------------------------------------------------------------------------------#

else:
    print('Error: When generating fastq source, param type ({}) is unknown. Workflow aborted.'.format(params.type))
    sys.exit()

#----------------------------------------------------------------------------------------------------------------------#

f.close() # Close log file
