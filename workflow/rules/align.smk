# Generate raw fastq source by calling housekeeping "scripts/hkp/sourcing_fastq.py"
rule prep_fastq_from_source:
  output:
    fq1=temp("{path}/raw/{sample}.fq1.gz"),
    fq2=temp("{path}/raw/{sample}.fq2.gz")
  params:
    type=lambda wildcards: merged_sample_sheet['type'].loc[wildcards.sample],
    path=lambda wildcards: merged_sample_sheet['Path'].loc[wildcards.sample]
  log: "{path}/raw/.{sample}.rule-align.prep_fastq_from_source.log"
  resources: disk_gb=lambda wildcards: get_disk_gb(merged_sample_sheet['type'].loc[wildcards.sample])
  conda: f"{workflow_dir}/envs/align.yaml"
  threads: 8
  script: f"{workflow_dir}/scripts/hkp/sourcing_fastq.py"

#----------------------------------------------------------------------------------------------------------------------#

# Generate bwa-meth aligned bam files, not position sorted
rule bwa_meth_align:
  input:
    fq1="{top_path}/resources/{bot_path}/raw/{sample}.fq1.gz",
    fq2="{top_path}/resources/{bot_path}/raw/{sample}.fq2.gz",
    ref=reference_genome_path,
    bwamethidx=reference_genome_path + ".bwameth.c2t"
  output:
    bam=temp("{top_path}/results/{bot_path}/alignments/{sample}.bam")
  log: "{top_path}/results/{bot_path}/alignments/.{sample}.rule-align.bwa_meth_align.log"
  conda: f"{workflow_dir}/envs/align.yaml"
  threads: 4
  shell:
    """
    exec > {log}; exec 2> {log}
    bwameth.py --reference {input.ref} -t {threads} {input.fq1} {input.fq2} | \
      samtools view -bh -o {output}
    """

#----------------------------------------------------------------------------------------------------------------------#

# Position sort a bam file
rule bam_position_sort:
  input: "{path}/{sample}.bam"
  output: "{path}/{sample}.POSsort.bam"
  log: "{path}/.{sample}.POSsort.rule-align.bam_position_sort.log"
  threads: 4
  conda: f"{workflow_dir}/envs/align.yaml"
  shell: "samtools sort -@ {threads} -O BAM -o {output} {input} > {log} 2> {log}"

#----------------------------------------------------------------------------------------------------------------------#

# Generate bam index
rule bam_index:
  input: "{path}/{sample}.bam"
  output: "{path}/{sample}.bam.bai"
  log: "{path}/.{sample}.rule-align.bam_index.log"
  threads: 4
  conda: f"{workflow_dir}/envs/align.yaml"
  shell: "samtools index -@ {input} > {log} 2> {log}"
