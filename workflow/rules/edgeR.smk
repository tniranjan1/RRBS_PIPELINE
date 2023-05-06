# Merge methylation data in windows using methylKit
rule run_methylKit:
  input:
    path="{path}/methylation_calls/samples/MethylDackel_{sample}.{repeats}_{context}.bedGraph.gz"
  output:
    path="{path}/methylation_calls/samples/methylKit_{sample}.{repeats}_{context}.bedGraph.gz"
  conda: f"{workflow_dir}/envs/edgeR.yaml"
  threads: 1
  log: "{path}/methylation_calls/samples/.methylKit_{sample}.{repeats}_{context}.log"
  script: f"{workflow_dir}/scripts/anl/methylKit_regions.R"

#----------------------------------------------------------------------------------------------------------------------#

# Perform edgeR analysis of methylation
rule run_edgeR:
  input: expand("{{path}}/methylation_calls/samples/methylKit_{sample}.{repeats}_{context}.bedGraph.gz", sample=all_sample_names, repeats=repeats, context='CpG')
  output: "{path}/methylation_calls/samples/edgeR.tmp"
  threads: 12
  conda: f"{workflow_dir}/envs/edgeR.yaml"
  log: "{path}/methylation_calls/samples/.edgeR.tmp.log"
  shell: "exec > {log}; exec 2> {log}; echo run > {output}"
  #script: f"{workflow_dir}/scripts/anl/edgeR.R"

#----------------------------------------------------------------------------------------------------------------------#
