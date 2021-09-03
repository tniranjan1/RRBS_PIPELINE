# Link source vcf to resource dir

rule link_vcf:
  input: vcf_source
  output: vcf_file
  shell: "ln {input} {output}"

#----------------------------------------------------------------------------------------------------------------------#

# Compress input vcf file genotypes.
#  This rule will take an input vcf file (from config file) containing genotypes for one or more samples in the study,
#  extract A <-> T transversion sites, and compress this information into a tab file for simpler access during
#  comparisons. This will be done for an individual chromosome to allow for parallelization.
rule compress_input_vcf_by_chr:
  input: vcf_file
  output: temp("{genotype_dir}/sample.compressed_genotypes.{chr}.txt")
  threads: 1
  conda: f"{workflow_dir}/envs/check_genotypes.yaml"
  log: "{genotype_dir}/.sample.compressed_genotypes.{chr}.rule-check_genotypes.compress_input_vcf_by_chr.log"
  script: f"{workflow_dir}/scripts/anl/compress_genotypes.py"

#----------------------------------------------------------------------------------------------------------------------#

# Rule to merge compressed genotypes across chromosomes.
rule merge_compressed_input_vcf:
  input: expand("{genotype_dir}/sample.compressed_genotypes.{chr}.txt", chr=chromosomes, allow_missing=True)
  output: "{genotype_dir}/sample.compressed_genotypes.txt"
  threads: 1
  conda: f"{workflow_dir}/envs/check_genotypes.yaml"
  log: "{genotype_dir}/.sample.compressed_genotypes.rule-check_genotypes.merge_compressed_input_vcf.log"
  shll:

#----------------------------------------------------------------------------------------------------------------------#
