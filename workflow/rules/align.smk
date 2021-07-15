# Generate empty source file indicating if source is SRA, BAM, or fastq
rule prep_source_type:
    output: "{path}/source/{type}.{sample}"
    shell: "touch {output}"

# Description (link to new bam and make fastq)
rule prep_from_bam:
    input:
        type="{path}/source/{type}.sample",
        bam_path=get_source
    output:
        fq1=temp("{path}/raw/{sample}.fq1.gz"),
        fq2=temp("{path}/raw/{sample}.fq2.gz")
    shell:
        """

        """

rule prep_from_fastq:
    input:
        type="{path}/source/{type}.sample",
        fq_path=get_source
    output:
        fq1=temp("{path}/raw/{sample}.fq1.gz"),
        fq2=temp("{path}/raw/{sample}.fq2.gz")
    shell:
        """

        """

rule prep_from_SRA:
    input:
        type="{path}/source/{type}.sample",
        SRR=get_source
    output:
        fq1=temp("{path}/raw/{sample}.fq1.gz"),
        fq2=temp("{path}/raw/{sample}.fq2.gz")
    shell:
        """

        """

rule bwa_meth_align:
    input:
        fq1="{path}/raw/{sample}.fq1.gz",
        fq2="{path}/raw/{sample}.fq2.gz"
    output:
        bam=temp("{path}/alignments/{sample}.bam")
    threads: 8
    conda:
        "envs/align.yaml"
    shell:
        """

        """
