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
    zcat {input} | grep '^chr' | awk 'BEGIN{{OFS="\t"}}{{$2=$2+1; print $0}}' | cut -f 1-3,5 > {output}
    """

#----------------------------------------------------------------------------------------------------------------------#

# Restrict lrs methylation files to infinium bases
rule restrict_LRS_methyl_toInfinium:
  input:
    methyl=lambda wildcards: lrs_methyl_sample_sheet['Path'].loc[wildcards.sample],
    infinium=resource_dir + "/infinium450k/infinium450k_hg38.bed"
  output: results_dir + "/lrs-methyl/samples/{lrs_sample}.methylationForEpiclock.bedGraph.gz"
  log: results_dir + "/lrs-methyl/samples/.{lrs_sample}.rule-agePrediction.restrict_LRS_methyl_toInfinium.log"
  conda: f"{workflow_dir}/envs/agePrediction.yaml"
  threads: 2
  shell:
    """
    exec > {log} 2> {log}

    zcat {input.methyl} | awk '($5>4)' | cut -f 1-3,7 |    \
      awk 'BEGIN{{OFS="\t"}}{{$2=$2+1;$3=$3+1;print $0}}' | \
      tail -n +2 > {output}.tmp

    bedtools unionbedg -i {output}.tmp {input.infinium} -filler 'remove' | grep -v 'remove' > {output}.uncomp
    rm {output}.tmp
    bgzip -@ {threads} -c {output}.uncomp > {output}
    rm {output}.uncomp
    """

#----------------------------------------------------------------------------------------------------------------------#

# Run epiclock on a sample bedgraph.gz
rule run_epiclock:
  input: "{path}/{sample}.methylationForEpiclock.bedGraph.gz"
  output: temp("{path}/{sample}.agePrediction.txt")
  log: "{path}/.{sample}.epiclockPrediction.rule-agePrediction.run_epiclock.log"
  conda: f"{workflow_dir}/envs/agePrediction.yaml"
  threads: 1
  shell: "touch {output}" # change output as non-temporary, and shift shell to Rscript

#----------------------------------------------------------------------------------------------------------------------#

# Merge epiclock results for lrs sample group
rule merge_and_markdown_epiclock:
  input: lambda wildcards: expand(f"{results_dir}/{{lrsORrrbs}}/samples/{{name}}.agePrediction.txt",
                           name = lrs_sample_names if wildcards.lrsORrrbs == 'lrs-methyl' else rrbs_sample_names)
  output: results_dir + "/{lrsORrrbs}/merged/merged.agePrediction.txt"
  log: results_dir + "/{lrsORrrbs}/merged/.merged.agePredction.rule-agePredction.merge_and_markdown_epiclock.log"
  conda: f"{workflow_dir}/envs/agePrediction.yaml"
  threads: 1
  shell: "touch {output}"
