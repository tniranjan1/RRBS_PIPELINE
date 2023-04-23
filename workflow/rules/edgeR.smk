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
  input:
  output:
  params:
  threads: 12
  conda: f"{workflow_dir}/envs/edgeR.yaml"
  log:
  script: f"{workflow_dir}/scripts/anl/edgeR.R"

#----------------------------------------------------------------------------------------------------------------------#
