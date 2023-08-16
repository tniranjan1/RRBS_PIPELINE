# Merge methylation data in windows using methylKit
rule run_methylKit:
  input:
    path=results_dir + "/rrbs_samples/methylation_calls/samples/MethylDackel_{sample}.{repeats}_{context}.bedGraph.gz"
  output:
    path=results_dir + "/rrbs_samples/methylation_calls/samples/methylKit_{sample}.{repeats}_{context}.bedGraph.gz"
  conda: f"{workflow_dir}/envs/edgeR.yaml"
  threads: 1
  log: results_dir + "/rrbs_samples/methylation_calls/samples/.methylKit_{sample}.{repeats}_{context}.log"
  script: f"{workflow_dir}/scripts/anl/methylKit_regions.R"

#----------------------------------------------------------------------------------------------------------------------#

# Make temporary trackless file for input into edgeR
rule trackless:
    input: "{prefix}.bedGraph.gz"
    output: temp("{prefix}.detrack.gz")
    conda: f"{workflow_dir}/envs/edgeR.yaml"
    threads: 1
    shell: "exec > /dev/null; exec 2> /dev/null; zcat {input} | tail -n +2 | bgzip -c > {output}"

#----------------------------------------------------------------------------------------------------------------------#

# Perform edgeR analysis of methylation
rule run_edgeR:
  input: expand(results_dir+"/rrbs_samples/methylation_calls/samples/methylKit_{sample}.{repeats}_{context}.detrack.gz",
                sample=set_options['rrbs_samples'], repeats=repeats, context='CpG')
  output: results_dir + "/rrbs_samples/methylation_calls/samples/edgeR.tmp"
  threads: 12
  conda: f"{workflow_dir}/envs/edgeR.yaml"
  log: results_dir + "/rrbs_samples/methylation_calls/samples/.edgeR.tmp.log"
  shell: "exec > {log}; exec 2> {log}; exit 1"
  #script: f"{workflow_dir}/scripts/anl/edgeR.R"

#----------------------------------------------------------------------------------------------------------------------#
