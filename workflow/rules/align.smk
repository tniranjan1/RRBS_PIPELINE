# Generate raw fastq source by calling housekeeping "scripts/hkp/sourcing_fastq.py"
rule prep_fastq_from_source:
  output:
    fq1=temp("{path}/raw/{sample}.fq1.gz"),
    fq2=temp("{path}/raw/{sample}.fq2.gz")
  params:
    type=lambda wildcards: merged_sample_sheet['type'].loc[wildcards.sample],
    path=lambda wildcards: merged_sample_sheet['Path'].loc[wildcards.sample]
  conda: f"{workflow_dir}/envs/align.yaml"
  resources: disk_gb=lambda wildcards: get_disk_gb(merged_sample_sheet['type'].loc[wildcards.sample])
  threads: 8
#  log: workflow_dir + "/logs/align_rules/prep_fastq_from_source/{sample}.log"
  script: f"{workflow_dir}/scripts/hkp/sourcing_fastq.py"

#----------------------------------------------------------------------------------------------------------------------#

# Generate bgzipped file from input, specifically for fastq input, to indicate output is temporary
rule fastq_gzip:
  input: "{path}/{sample}.fq{num}"
  output: temp("{path}/{sample}.fq{num}.gz")
  conda: f"{workflow_dir}/envs/align.yaml"
  threads: 4
  log: f"{workflow_dir}/logs/align_rules/fastq_gzip/{{sample}}.fq{{num}}.log"
  shell: "bgzip -@ {threads} -c {input} > {output} 2> {log}"

#----------------------------------------------------------------------------------------------------------------------#

# Generate bwa-meth aligned bam files, not position sorted
rule bwa_meth_align:
  input:
    fq1="{path}/raw/{sample}.fq1.gz",
    fq2="{path}/raw/{sample}.fq2.gz",
    ref=reference_genome_path,
    bwamethidx=reference_genome_path + ".bwameth.c2t"
  output:
    bam=temp("{path}/alignments/{sample}.bam")
  threads: 4
  conda: f"{workflow_dir}/envs/align.yaml"
  log: f"{workflow_dir}/logs/align_rules/bwa_meth_align/{{sample}}.log"
  shell:
    """
    bwameth.py --reference {input.ref} -t {threads} {input.fq1} {input.fq2} | \
      samtools view -bh -o {output} 2> {log}
    """

#----------------------------------------------------------------------------------------------------------------------#

# Position sort a bam file
rule bam_position_sort:
  input: "{path}/{sample}.bam"
  output: "{path}/{sample}.POSsort.bam"
  threads: 4
  conda: f"{workflow_dir}/envs/align.yaml"
  log: f"{workflow_dir}/logs/align_rules/bam_position_sort/{{sample}}.log"
  shell: "samtools sort -@ {threads} -O BAM -o {output} {input} 2> {log}"

#----------------------------------------------------------------------------------------------------------------------#

# Generate bam index
rule bam_index:
  input: "{path}/{sample}.bam"
  output: "{path}/{sample}.bam.bai"
  threads: 4
  conda: f"{workflow_dir}/envs/align.yaml"
  log: f"{workflow_dir}/logs/align_rules/bam_index/{{sample}}.log"
