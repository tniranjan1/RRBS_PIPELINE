# Get condensed sites for deconvolution
rule condensed_sites_list:
  input:
    rrbs="{path1}/rrbs_samples/{path2}/merged_methylation.without_repeats.noNAs.CpG.bedGraph.gz",
    deconv="{path1}/deconvo_ref_samples/{path2}/merged_methylation.without_repeats.noNAs.CpG.bedGraph.gz"
  output: "{path1}/deconvo_ref_samples/{path2}/deconvolution_sites/sites.bed"
  threads: 1
  conda: f"{workflow_dir}/envs/deconvolution.yaml"
  log: "{path1}/deconvo_ref_samples/{path2}/deconvolution_sites/.sites.rule-get_methylation.condensed_sites_list.log"
  shell:
    """
    zcat {input.rrbs} | cut -f-3 > {input.rrbs}.fix
    zcat {input.deconv} | cut -f-3 > {input.deconv}.fix
    bedtools intersect -a {input.deconv}.fix -b {input.rrbs}.fix > {output}
    rm {input.rrbs}.fix {input.deconv}.fix
    """

#----------------------------------------------------------------------------------------------------------------------#

# rule reduce methylation data to deconvolution site list for deconvulation
rule reduce_methylation_for_deconvolution:
  input:
    sample="{path}/{file}.gz",
    bed=deconvo_site_list
  output: temp("{path}/{file}.deco.tmp")
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
