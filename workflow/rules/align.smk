# Generate raw fastq source
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
    log: "{path}/raw/{sample}.fq.log"
    script:
        f"{workflow_dir}/scripts/hkp/sourcing_fastq.py"

rule fastq_gzip:
    input: "{sample}"
    output: temp("{sample}.gz")
    conda: f"{workflow_dir}/envs/align.yaml"
    threads: 4
    shell: "bgzip -@ {threads} -c {input} > {output}"

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
    shell:
        "bwameth.py --reference {input.ref} -t {threads} {input.fq1} {input.fq2} | samtools view -bh -o {output}"

rule bam_sort:
    input: "{sample}.bam"
    output: "{sample}.POSsort.bam"
    threads: 4
    conda: f"{workflow_dir}/envs/align.yaml"
    shell: "samtools sort -@ {threads} -O BAM -o {output} {input}"
