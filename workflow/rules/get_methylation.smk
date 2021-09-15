# Extract methylation using MethylDackel
rule extract_methylation:
  input:
    bam="{path}/alignments/{sample}.POSsort.bam",
    bai="{path}/alignments/{sample}.POSsort.bam.bai",
    ref=reference_genome_path,
    ir=inverted_repeats
  output:
    orig=temp(
         expand("{path}/methylation_calls/samples/MethylDackel_{sample}.{repeats}_{context}.{suffix}",
                context=[ 'CpG', 'CHG', 'CHH' ], allow_missing=True)
             ),
    gzip=expand("{path}/methylation_calls/samples/MethylDackel_{sample}.{repeats}_{context}.{suffix}.gz",
                context=[ 'CpG', 'CHG', 'CHH' ], allow_missing=True)
  wildcard_constraints:
    suffix="bedGraph|methylKit"
  threads: 4
  conda: f"{workflow_dir}/envs/get_methylation.yaml"
  log: "{path}/.MethylDackel_{sample}.{repeats}.{suffix}.rule-get_methylation.extract_methylation.log"
  script: f"{workflow_dir}/scripts/anl/run_MethylDackel.py"

#----------------------------------------------------------------------------------------------------------------------#

# Merge methylation calls for all samples in a sheet, subdivided by chromosome
rule merge_methylation_by_chr:
  input:
    orig=expand("{path}/samples/MethylDackel_{sample}.{repeats}_{context}.bedGraph",
                sample=rrbs_samples.index, allow_missing=True),
    gzip=expand("{path}/samples/MethylDackel_{sample}.{repeats}_{context}.bedGraph.gz",
                sample=rrbs_samples.index, allow_missing=True)
  output: temp("{path}/merged/merged_methylation.{repeats}.{context}.{chr}.bedGraph")
  log: "{path}/merged/.merged_methylation.{repeats}.{context}.{chr}.rule-get_methylation.merge_methylation_by_chr.log"
  conda: f"{workflow_dir}/envs/get_methylation.yaml"
  shell:
      """
      exec > {log}; exec 2> {log}
      chr_input=""
      for i in {input.orig}
      do
        grep '^{wildcards.chr}' $i | cut -f1-4 > $i.{wildcards.chr}.methy_tmp
        chr_input="$chr_input $i.{wildcards.chr}.methy_tmp"
      done
      bedtools unionbedg -filler NA -i $chr_input > {output}
      for i in {input.orig}; do rm $i.{wildcards.chr}.methy_tmp; done
      """

#----------------------------------------------------------------------------------------------------------------------#

# Merge coverage information for all samples in a sheet, subdivided by chromosome
rule merge_coverage_by_chr:
  input:
    orig=expand("{path}/samples/MethylDackel_{sample}.{repeats}_{context}.bedGraph",
                sample=rrbs_samples.index, allow_missing=True),
    gzip=expand("{path}/samples/MethylDackel_{sample}.{repeats}_{context}.bedGraph.gz",
                sample=rrbs_samples.index, allow_missing=True)
  output: temp("{path}/merged/merged_coverage.{repeats}.{context}.{chr}.bedGraph")
  log: "{path}/merged/.merged_coverage.{repeats}.{context}.{chr}.rule-get_methylation.merge_coverage_by_chr.log"
  conda: f"{workflow_dir}/envs/get_methylation.yaml"
  shell:
      """
      exec > {log}; exec 2> {log}
      chr_input=""
      for i in {input.orig}
      do
        grep '^{wildcards.chr}' $i | awk '{{print $1 "\t" $2 "\t" $3 "\t" $5+$6}}' > $i.{wildcards.chr}.cov_tmp
        chr_input="$chr_input $i.{wildcards.chr}.cov_tmp"
      done
      bedtools unionbedg -filler NA -i $chr_input > {output}
      for i in {input.orig}; do rm $i.{wildcards.chr}.cov_tmp; done
      """

#----------------------------------------------------------------------------------------------------------------------#

# Merge chromosome-separated tables
rule merge_table_from_chr:
  input:
    fai=reference_genome_path + ".fai",
    merge_list=expand("{path}/merged_{meco}.{repeats}.{context}.{chr}.bedGraph", chr=chromosomes, allow_missing=True)
  output: temp("{path}/merged_{meco}.{repeats}.{context}.bedGraph")
  params:
    sample_names = rrbs_sample_names
  conda: f"{workflow_dir}/envs/get_methylation.yaml"
  log: "{path}/.merged_{meco}.{repeats}.{context}.rule-get_methylation.merge_table_from_chr.log"
  shell:
        """
        exec > {log}; exec 2> {log}
        echo "chrom" "start" "end" {params.sample_names} | sed 's/ /\t/g' > {output} ## Print header
        cat {input.merge_list} >> {output}
        """

#----------------------------------------------------------------------------------------------------------------------#

# Bgzip methylation tables for storage space reduction
rule bgzip_merged_tables:
  input: "{path}/merged_{file_name}"
  output: "{path}/merged_{file_name}.gz"
  threads: 4
  conda: f"{workflow_dir}/envs/get_methylation.yaml"
  log: "{path}/.merged_{file_name}.rule-get_methylation.bgzip_table.log"
  shell: "bgzip -@ {threads} -c {input} > {output} 2> {log}"

#----------------------------------------------------------------------------------------------------------------------#

# Remove sites with <100% call rate (remove NAs) from merged methylation/coverage table.
rule remove_NA:
  input: "{path}/merged_{meco}.{repeats}.{context}.bedGraph.gz"
  output: temp("{path}/merged_{meco}.{repeats}.noNAs.{context}.bedGraph")
  threads: 2
  conda: f"{workflow_dir}/envs/get_methylation.yaml"
  log: "{path}/.merged_{meco}.{repeats}.noNAs.{context}.rule-get_methylation.remove_NA.log"
  shell: "zcat {input} | grep -vP '\tNA' > {output}"

#----------------------------------------------------------------------------------------------------------------------#
