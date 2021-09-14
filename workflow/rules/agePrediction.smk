# Retrieve infinium450k data
rule retrieve_infinium450k_hg38:
  output: resource_dir + "/infinium450k/infinium450k_hg38.tsv.gz"
  params:
    url=config['infinium450k']['url']
  log: resource_dir + "/infinium450k/.infinium450k_hg38.rule-agePrediction.retrieve_infinium450k_hg38.log"
  conda: f"{workflow_dir}/envs/agePrediction.yaml"
  threads: 1
  shell: "wget -o {log} -O {output} {params.url}"

#----------------------------------------------------------------------------------------------------------------------#

# Convert infinium450k data from tsv to bed format
rule infinium450k_tsvTObed:
  input: resource_dir + "/infinium450k/infinium450k_hg38.tsv.gz"
  output: resource_dir + "/infinium450k/infinium450k_hg38.bed"
  log: resource_dir + "/infinium450k/.infinium450k_hg38.rule-agePrediction.infinium450k_tsvTObed.log"
  conda: f"{workflow_dir}/envs/agePrediction.yaml"
  threads: 1
  shell:
    """
    exec > {log} 2> {log}
    echo -e "CHR\tStart\tStop\tProbeID" > {output}
    zcat {input} | grep '^chr' | awk 'BEGIN{{OFS="\t"}}{{$2=$2+1; print $0}}' | cut -f 1-3,5 >> {output}
    """

#----------------------------------------------------------------------------------------------------------------------#

# Restrict lrs methylation files to infinium bases
rule restrict_LRS_methyl_toInfinium:
  input:
    methyl=lambda wildcards: lrs_methyl_sample_sheet['Path'].loc[wildcards.sample],
    infinium=resource_dir + "/infinium450k/infinium450k_hg38.bed"
  output: results_dir + "/lrs-methyl/samples/{sample}.methylation.bedGraph.gz"
  log: results_dir + "/lrs-methyl/samples/.{sample}.rule-agePrediction.restrict_LRS_methyl_toInfinium.log"
  conda: f"{workflow_dir}/envs/agePrediction.yaml"
  threads: 1
  script: ""

#----------------------------------------------------------------------------------------------------------------------#
