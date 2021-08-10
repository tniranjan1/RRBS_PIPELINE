# Bgzip methylation tables for storage space reduction
rule bgzip_table:
  input: "{path}/methylation_calls/{sample}.{suffix}"
  output: "{path}/methylation_calls/{sample}.{suffix}.gz"
  wildcard_constraints:
    suffix="\w+"
  threads: 8
  conda: f"{workflow_dir}/envs/methylation.yaml"
  log: f"{workflow_dir}/logs/get_methylation_rules/bgzip_table/{{sample}}.{{suffix}}.log"
  shell: "bgzip -@ {threads} -c {input} > {output}"

#----------------------------------------------------------------------------------------------------------------------#

# Extract methylation using MethylDackel
rule extract_methylation:
  input:
    bam="{path}/alignments/{sample}.POSsort.bam",
    bai="{path}/alignments/{sample}.POSsort.bam.bai",
    ref=reference_genome_path,
    ir=inverted_repeats
  output:
    CpG=temp("{path}/methylation_calls/{sample}_CpG.{repeats}.{suffix}") if context_truth['CpG'] else [],
    CHG=temp("{path}/methylation_calls/{sample}_CHG.{repeats}.{suffix}") if context_truth['CHG'] else [],
    CHH=temp("{path}/methylation_calls/{sample}_CHG.{repeats}.{suffix}") if context_truth['CHH'] else [],
  wildcard_constraints:
    suffix="bedGraph|methylKit"
  threads: 4
  conda: f"{workflow_dir}/envs/methylation.yaml"
  log: f"{workflow_dir}/logs/get_methylation_rules/extract_methylation/{{sample}}.log"
  run:
      mKit = "--methylKit" if wildcards.suffix == "methylKit" else ""
      CHG = "--CHG" if context_truth['CHG'] else ""
      CHH = "--CHH" if context_truth['CHH'] else ""
      repeat = "-l " + input.ir if wildcards.repeats == "without_repeats" else ""
      prefix = wildcards.path + "/methylation_calls/" + wildcards.sample
      MDackel = "MethylDackel extract"
      shell("{MDackel} -@ {threads} -q 20 -d 5 {repeat} {CHG} {CHH} {mKit} -o {prefix} {input.ref} {input.bam}")

#----------------------------------------------------------------------------------------------------------------------#

# Merge methylation calls for all samples in a sheet, subdivided by chromosome
rule merge_methylation_by_chr:
  input:
    orig: lambda wildcards: get_extracted_files(sample_sheet=rrbs_samples, wildcards=wildcards, gz=False),
    gzip: lambda wildcards: get_extracted_files(sample_sheet=rrbs_samples, wildcards=wildcards, gz=True )
  output: temp("{path}/methylation_calls/merged/merged_methylation.{repeats}.{context}.{chr}.bedGraph")
  params:
    sample_names=rrbs_samples.index
  log: f"{workflow_dir}/logs/get_methylation_rules/merge_methylation_by_chr/merged_methylation.{{repeats}}.{{context}}.{{chr}}.log"
  conda: f"{workflow_dir}/envs/methylation.yaml"
  shell:
      """
      chr_input=""
      for i in {input.orig}
      do
        grep '^{wildcards.chr}' $i | cut -f1-4 > $i.methy_tmp
        chr_input="$chr_input $i.methy_tmp"
      done
      bedtools unionbedg -filler NA -i $chr_input > {output}
      for i in {input.orig}; do rm $i.methy_tmp; done
      """

#----------------------------------------------------------------------------------------------------------------------#

# Merge coverage information for all samples in a sheet, subdivided by chromosome
rule merge_coverage_by_chr:
  input:
    orig: lambda wildcards: get_extracted_files(sample_sheet=rrbs_samples, wildcards=wildcards, gz=False),
    gzip: lambda wildcards: get_extracted_files(sample_sheet=rrbs_samples, wildcards=wildcards, gz=True )
  output: temp("{path}/methylation_calls/merged/merged_coverage.{repeats}.{context}.{chr}.bedGraph")
  params:
    sample_names=rrbs_samples.index
  log: f"{workflow_dir}/logs/get_methylation_rules/merge_coverage_by_chr/merged_methylation.{{repeats}}.{{context}}.{{chr}}.log"
  conda: f"{workflow_dir}/envs/methylation.yaml"
  shell:
      """
      chr_input=""
      for i in {input.orig}
      do
        grep '^{wildcards.chr}' $i | awk '{{print $1 "\t" $2 "\t" $3 "\t" $5+$6}}' > $i.cov_tmp
        chr_input="$chr_input $i.cov_tmp"
      done
      bedtools unionbedg -filler NA -i $chr_input > {output}
      for i in {input.orig}; do rm $i.cov_tmp; done
      """

#----------------------------------------------------------------------------------------------------------------------#

# Merge chromosome-separated tables
rule merge_table_from_chr:
  input:
    fai=reference_genome_path + ".fai",
    merge_list=lambda wildcards: get_merge_list(fai=reference_genome_path + ".fai", wildcards=wildcards)
  output: temp("{path}/methylation_calls/merged/merged_{meco}.{repeats}.{context}.bedGraph")
  params:
    sample_names=rrbs_samples.index
  conda: f"{workflow_dir}/envs/methylation.yaml"
  log: f"{workflow_dir}/logs/get_methylation_rules/merge_table_from_chr/merged_{{meco}}.{{repeats}}.{{context}}.log"
  shell:
        """
        echo "chrom" "start" "end" {params.sample_names} | sed 's/ /\t/g' > {output} ## Print header
        cat {input.merge_list} >> {output}
        """
