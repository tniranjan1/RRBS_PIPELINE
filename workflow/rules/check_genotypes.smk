# Link source vcf to resource dir
rule link_vcf:
  input: vcf_source
  output: vcf_file
  threads: 1
  run:
    import shutil
    with open(input, 'rb') as test_f:
        if test_f.read(2) == b'\x1f\x8b':
            shutil.copy(input, output)
        else:
            test_f.close()
            import gzip
            with open(input, 'rb') as f_in, gzip.open(output, 'wb') as f_out:
                f_out.writelines(f_in)

#----------------------------------------------------------------------------------------------------------------------#

# Compress input vcf file genotypes.
#  This rule will take an input vcf file (from config file) containing genotypes for one or more samples in the study,
#  extract A <-> T transversion sites, and compress this information into a tab file for simpler access during
#  comparisons.
rule compress_input_vcf:
  input: vcf_file
  output: f"{genotype_dir}/sample.AT_genotypes.txt"
  threads: 1
  conda: f"{workflow_dir}/envs/check_genotypes.yaml"
  log: f"{genotype_dir}/.sample.AT_genotypes.rule-check_genotypes.compress_input_vcf.log"
  shell:
    """
    zcat {input} | \
        grep -v '#' | \
        awk '{ if( ($4 == "A" && $5 == "T") || ($4 == "T" && $5 =="A")) print $1 "\t" $2 "\t" $2+1}' > {output}
    """

#----------------------------------------------------------------------------------------------------------------------#

# Rule to get coverage at potential genotyping sites.
rule compress_adequate_cov_input_vcf:
  input:
    bam=expand(f"{results_dir}/rrbs_samples/alignments/{sample}.POSsort.bam", sample=rrbs_sample_names),
    bai=expand(f"{results_dir}/rrbs_samples/alignments/{sample}.POSsort.bam.bai", sample=rrbs_sample_names),
    bed=f"{genotype_dir}/sample.AT_genotypes.txt"
  output: f"{genotype_dir}/sample.potential_gt_cov.txt"
  threads: 1
  conda: f"{workflow_dir}/envs/check_genotypes.yaml"
  log: f"{genotype_dir}/.sample.potential_gt_cov.rule-check_genotypes.compress_adequate_cov_input_vcf.log"
  shell:
    """
    samtools mpileup -d 10 -l [input.bed] {input.bam} | \
      awk '{s=0; for(i=1;i<=NF;i++){if((i-1) % 3 == 0){if($i < 5)s++}};if(s==0) print $1 "\t" $2 "\t" $2+1}' > {output}
    """

#----------------------------------------------------------------------------------------------------------------------#

# get sample genotypes
rule get_sample_genotypes:
  input:
    bam=expand("{path}/alignments/{sample}.POSsort.bam", path=results_dir, sample=rrbs_sample_names),
    bai=expand("{path}/alignments/{sample}.POSsort.bam.bai", path=results_dir, sample=rrbs_sample_names),
    bed=f"{genotype_dir}/sample.potential_gt_cov.txt",
    ref=reference_genome_path
  output: f"{genotype_dir}/all_samples.gt.vcf"
  threads: 1
  conda: f"{workflow_dir}/envs/bcftools.yaml"
  log: f"{results_dir}/genotypes/.all_samples.gt.rule-check_genotypes.get_sample_genotypes.log"
  shell:
    """
    bcftools mpileup -R [input.bed] --fasta-ref {input.ref} {input.bam} | bcftools call -m - > {output}
    """
