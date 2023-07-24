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
    neuN="{path1}/deconvo_ref_samples/methylation_calls/merged/merged_methylation.{repeats}.CpG.bedGraph.gz"
  output: "{path1}/deconvo_ref_samples/deconvolution/deconvo_values.{repeats}.significance.txt"
  params:
    script=f"{workflow_dir}/scripts/anl/deconvolute.R"
  threads: 12
  conda: f"{workflow_dir}/envs/deconvolution.yaml"
  resources:
    mem_gb=100
  log: "{path1}/deconvo_ref_samples/deconvolution/.deconvo_values.{repeats}.log"
  shell:
    """
    exec > {log}; exec 2> {log}
    Rscript {params.script} {threads} {input.merge1} {input.neuN} {output}
    """

#----------------------------------------------------------------------------------------------------------------------#
