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
  output: results_dir + "/lrs_methyl/samples/{sample}.methylationForEpiclock.bedGraph.gz"
  log: results_dir + "/lrs_methyl/samples/.{sample}.rule-agePrediction.restrict_LRS_methyl_toInfinium.log"
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
  output: "{path}/{sample}.agePrediction.txt"
  log: "{path}/.{sample}.epiclockPrediction.rule-agePrediction.run_epiclock.log"
  conda: f"{workflow_dir}/envs/agePrediction.yaml"
  threads: 12
  shell: f"Rscript {workflow_dir}/scripts/anl/run_epigeneticClock.R {{input}} {{threads}} > {{output}} 2> {{log}}"

#----------------------------------------------------------------------------------------------------------------------#

# Merge epiclock results for lrs sample group
rule merge_and_markdown_epiclock:
  input: lambda wildcards: distinguish_lrs_rrbs(middle_folder=wildcards.path, return_value="files")
  output: results_dir + "/{path}/merged/merged.agePrediction.txt"
  params: lambda wildcards: distinguish_lrs_rrbs(middle_folder=wildcards.path, return_value="samples")
  log: results_dir + "/{path}/merged/.merged.agePrediction.rule-agePrediction.merge_and_markdown_epiclock.log"
  threads: 1
  run:
    sample_names = params[0]
    output_df = lrs_methyl_sample_sheet[ [ 'SampleID', 'Covariate_Age' ] ]
    output_df['EpiToc'] = [ np.nan ] * len(output_df)
    output_df['Horvath'] = [ np.nan ] * len(output_df)
    output_df['Hannum'] = [ np.nan ] * len(output_df)
    for index in range(len(input)):
        filename = input[index]
        sample_name = sample_names[index]
        predictions = pd.read_table(filename, sep="=", names=[ 'Value' ], index_col=0)
        for index in predictions.index:
            output_df.at[sample_name,index] = predictions['Value'][index]
    output_df.to_csv(output[0], sep='\t', index=False)
