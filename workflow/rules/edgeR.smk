# Perform edgeR analysis of methylation
rule deconvolute:
  input:
  output:
  params:
  threads: 12
  conda: f"{workflow_dir}/envs/edgeR.yaml"
  log:
  script: f"{workflow_dir}/scripts/anl/edgeR.R"
