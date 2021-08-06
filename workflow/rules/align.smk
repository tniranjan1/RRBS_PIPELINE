# Generate raw fastq source
rule prep_fastq_from_source:
    output:
        fq1=temp("{path}/raw/{sample}.fq1.gz"),
        fq2=temp("{path}/raw/{sample}.fq2.gz")
    params:
        type=lambda wildcards: merged_sample_sheet['type'].loc[wildcards.sample],
        path=lambda wildcards: merged_sample_sheet['Path'].loc[wildcards.sample]
    conda: "envs/align.yaml"
    resources: disk_gb=get_disk_gb
    threads: 8
    script:
        "scripts/hkp/sourcing_fastq.py"

rule fastq_gzip:
    input: "{sample}"
    output: temp("{sample}.gz")
    conda: "envs/align.yaml"
    threads: 4
    shell: "bgzip -@ {threads} -c {input} > {output}"

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
