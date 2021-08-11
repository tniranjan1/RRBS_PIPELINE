# Bgzip methylation tables for storage space reduction
rule bgzip_table:
  input: "{path}/{file_name}"
  output: "{path}/{file_name}.gz"
  threads: 8
  conda: f"{workflow_dir}/envs/get_methylation.yaml"
  log: "{path}/.{file_name}.rule-get_methylation.bgzip_table.log"
  shell: "bgzip -@ {threads} -c {input} > {output} > {log} 2> {log}"

#----------------------------------------------------------------------------------------------------------------------#

# Extract methylation using MethylDackel
rule extract_methylation:
  input:
    bam="{path}/alignments/{sample}.POSsort.bam",
    bai="{path}/alignments/{sample}.POSsort.bam.bai",
    ref=reference_genome_path,
    ir=inverted_repeats
  output:
    temp(
        expand("{path}/methylation_calls/samples/{sample}.{repeats}_{context}.bedGraph",
            context=context_to_use, allow_missing=TRUE)
        )
  wildcard_constraints:
    suffix="bedGraph|methylKit"
  threads: 4
  conda: f"{workflow_dir}/envs/get_methylation.yaml"
  log: "{path}/methylation_calls/samples/.{sample}.{repeats}.rule-get_methylation.extract_methylation.log"
  run:
      mKit = "--methylKit" if wildcards.suffix == "methylKit" else ""
      CHG = "--CHG" if context_truth['CHG'] else ""
      CHH = "--CHH" if context_truth['CHH'] else ""
      repeat = "-l " + input.ir if wildcards.repeats == "without_repeats" else ""
      prefix = wildcards.path + "/methylation_calls/" + wildcards.sample + "." + wildcards.repeats
      MDkl = "MethylDackel extract -q 20 -d 5"
      shell("{MDkl} -@ {threads} {repeat} {CHG} {CHH} {mKit} -o {prefix} {input.ref} {input.bam} > {log} 2> {log}")

#----------------------------------------------------------------------------------------------------------------------#

# Merge methylation calls for all samples in a sheet, subdivided by chromosome
rule merge_methylation_by_chr:
  input:
    orig=lambda wildcards: get_extracted_files(sample_sheet=rrbs_samples, wildcards=wildcards, gz=False),
    gzip=lambda wildcards: get_extracted_files(sample_sheet=rrbs_samples, wildcards=wildcards, gz=True )
  output: temp("{path}/methylation_calls/merged/merged_methylation.{repeats}.{context}.{chr}.bedGraph")
  params:
    sample_names=rrbs_samples.index
  log: "{path}/methylation_calls/merged/" + \
       ".merged_methylation.{repeats}.{context}.{chr}.rule-get_methylation.merge_methylation_by_chr.log"
  conda: f"{workflow_dir}/envs/get_methylation.yaml"
  shell:
      """
      exec > {log}; exec 2> {log}
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
    orig=lambda wildcards: get_extracted_files(sample_sheet=rrbs_samples, wildcards=wildcards, gz=False),
    gzip=lambda wildcards: get_extracted_files(sample_sheet=rrbs_samples, wildcards=wildcards, gz=True )
  output: temp("{path}/methylation_calls/merged/merged_coverage.{repeats}.{context}.{chr}.bedGraph")
  params:
    sample_names=rrbs_samples.index
  log: "{path}/methylation_calls/merged/" + \
       ".merged_coverage.{repeats}.{context}.{chr}.rule-get_methylation.merge_coverage_by_chr.log"
  conda: f"{workflow_dir}/envs/get_methylation.yaml"
  shell:
      """
      exec > {log}; exec 2> {log}
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
  conda: f"{workflow_dir}/envs/get_methylation.yaml"
  log: "{path}/methylation_calls/merged/" + \
       ".merged_{meco}.{repeats}.{context}.rule-get_methylation.merge_table_from_chr.log"
  shell:
        """
        exec > {log}; exec 2> {log}
        echo "chrom" "start" "end" {params.sample_names} | sed 's/ /\t/g' > {output} ## Print header
        cat {input.merge_list} >> {output}
        """
