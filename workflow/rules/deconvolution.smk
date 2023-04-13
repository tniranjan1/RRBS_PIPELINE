# Get additional reference
rule get_GSE186458:
  output: "{path1}/deconvo_ref_samples/GSE186458/merge.intersect.bed.gz"
  params:
    url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186458/suppl/GSE186458_RAW.tar",
    dir=lambda wc: wc.path1 + "/deconvo_ref_samples/GSE186458",
    script=f"{workflow_dir}/scripts/anl/bigwig_to_bed.R"
  threads: 12
  conda: f"{workflow_dir}/envs/deconvolution.yaml"
  log: "{path1}/deconvo_ref_samples/GSE186458/.get_GSE186458.log"
  shell:
    """
    wget -O {params.dir}/GSE186458_RAW.tar {params.url}
    tar -xvf {params.dir}/GSE186458_RAW.tar -C {params.dir}
    Rscript {params.script} {params.dir} {threads}
    """

#----------------------------------------------------------------------------------------------------------------------#

# Perform deconvolution
rule deconvolute:
  input:
    merge1="{path1}/deconvo_ref_samples/GSE186458/merge.intersect.bed.gz",
    neuN="{path1}/deconvo_ref_samples/methylation_calls/merged/merged_methylation.with_repeats.CpG.bedGraph.gz"
  output: "{path1}/deconvo_ref_samples/deconvolution/deconvo_values.txt"
  params:
    script=f"{workflow_dir}/scripts/anl/deconvolute.R"
  threads: 12
  conda: f"{workflow_dir}/envs/deconvolution.yaml"
  log: "{path1}/deconvo_ref_samples/deconvolution/.deconvo_values.log"
  shell:
    """
    exec > {log}; exec 2> {log}
    Rscript {threads} {input.merge1} {input.neuN} {output}
    """

#----------------------------------------------------------------------------------------------------------------------#

# rule reduce methylation data to deconvolution site list for deconvulation
rule reduce_methylation_for_deconvolution:
  input:
    sample="{path}/{file}.bedGraph.gz",
    bed=deconvo_site_list
  output: "{path}/{file}.deco.tmp"
  threads: 1
  conda: f"{workflow_dir}/envs/deconvolution.yaml"
  log: "{path}/.{file}.deco.tmp.rule-get_methylation.reduce_methylation_for_deconvolution.log"
  shell:
    "zcat {input.sample} | bedtools intersect -wb -a {input.bed} -b - | cut -f4-6,8,9 > {output}"

#----------------------------------------------------------------------------------------------------------------------#

# rule to perform deconvolution using deconvSeq in R
rule deconvSeq:
  input:
    rrbs=expand("{dir}/{analyte}/methylation_calls/samples/MethylDackel_{sample}.{repeats}_{context}.deco.tmp",
                dir = results_dir, analyte = "rrbs_samples", sample=rrbs_sample_names,
                repeats = "without_repeats", context = "CpG"),
    deco=expand("{dir}/{analyte}/methylation_calls/samples/MethylDackel_{sample}.{repeats}_{context}.deco.tmp",
                dir = results_dir, analyte = "deconvo_ref_samples", sample=deconvo_sample_names,
                repeats = "without_repeats", context = "CpG")
  output: f"{results_dir}/rrbs_samples/deconvolution/deconvolution.txt"
  threads: 1
  conda: f"{workflow_dir}/envs/deconvolution.yaml"
  log: f"{results_dir}/rrbs_samples/deconvolution/.deconvolution.log"
  shell: "echo 'done' > {output}"
