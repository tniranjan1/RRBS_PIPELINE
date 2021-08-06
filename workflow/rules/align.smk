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
    run:
        if params.type == 'SRR':
            # Get path to local SRR download repository
            sra_repository_command = "vdb-config -o n -p | grep /repository/user/main/public/root | cut -d= -f2 | sed 's/ //g'"
            sra_download_repository = os.popen(sra_repository_command).read()
            sra_download_repository = sra_download_repository.replace('"', '').replace("'", "").strip("\n") + "/sra"
            # Download SRR from SRA repository
            shell("prefetch --max-size 100000000 {params.path}")
            # multi-thread fastq dump from local SRR
            shell("mkdir -p ../resources/sra_temp")
            shell("parallel-fastq-dump --sra-id {params.path} --threads {threads} --outdir ../resources/sra_temp/{params.path} --split-files --gzip")
            # move fastq dump files to desired output paths
            shell("mv ../resources/sra_temp/{params.path}/{params.path}_1.fastq.gz {output.fq1}")
            shell("mv ../resources/sra_temp/{params.path}/{params.path}_2.fastq.gz {output.fq2}")
            # delete local SRR download
            shell("rm {sra_download_repository}/{params.path}.sra")
        elif params.type == 'bam':
            # Read sort bam file
            shell("mkdir -p ../resources/bam_temp")
            shell("samtools sort -@ {threads} -n -O BAM -o ../resources/bam_temp/{wildcards.sample}.bam {params.path}")
            # Generate fastq files (gzipped) from bam file
            shell("samtools fastq -@ {threads} -n -1 {output.fq1} -2 {output.fq2} -0 /dev/null -s /dev/null ../resources/bam_temp/{wildcards.sample}.bam")
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
            # Link input fastq files to standard destiantion
            shell("ln {fq1_input_path} {fq1_output_path}")
            shell("ln {fq2_input_path} {fq2_output_path}")
        else:
            print('Error: When generating fastq source, param type ({}) is unknown. Workflow aborted.'.format(param.type))
            sys.exit()

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
