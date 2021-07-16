# Generate raw fastq source
rule prep_fastq_from_source:
    input:
        type=lambda wildcards: merged_sample_sheet['type'].loc[wildcards.sample],
        path=lambda wildcards: merged_sample_sheet['Path'].loc[wildcards.sample]
    output:
        fq1=temp("{path}/raw/{sample}.fq1.gz"),
        fq2=temp("{path}/raw/{sample}.fq2.gz")
    shell:""
        # python code if input.type = 'SRR'
        # python code if input.type = 'bam'
        # python code if input.type = 'fq'

rule bwa_meth_align:
    input:
        fq1="{path}/raw/{sample}.fq1.gz",
        fq2="{path}/raw/{sample}.fq2.gz"
    output:
        bam=temp("{path}/alignments/{sample}.bam")
    threads: 8
#    conda:
#        "envs/align.yaml"
    shell:
        """

        """

rule bam_sort:
    input: "{sample}.bam"
    output: "{sample}.POSsort.bam"
    shell: ""
