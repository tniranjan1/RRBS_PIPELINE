# Generate reference dictionary
rule reference_dict:
  input: "{reference}.fa"
  output: "{reference}.dict"
  params:
    build=config['ref']['hg_build'],
    species=config['ref']['species']
  conda: f"{workflow_dir}/envs/index.yaml"
  shell: "samtools dict -a {params.build} -s {params.species} -o {output} {input}"

#----------------------------------------------------------------------------------------------------------------------#

# Generate reference index
rule reference_idx:
  input: "{reference}.fa"
  output: "{reference}.fa.fai"
  conda: f"{workflow_dir}/envs/index.yaml"
  shell: "samtools faidx {input}"

#----------------------------------------------------------------------------------------------------------------------#

# Generate bwa-meth index
rule bwa_meth_index:
  input: "{reference}.fa"
  output:
    main="{reference}.fa.bwameth.c2t",
    amb="{reference}.fa.bwameth.c2t.amb",
    ann="{reference}.fa.bwameth.c2t.ann",
    pac="{reference}.fa.bwameth.c2t.pac"
  conda: f"{workflow_dir}/envs/index.yaml"
  priority: 10
  shell: "bwameth.py index {input}"

#----------------------------------------------------------------------------------------------------------------------#

# Hard link reference repeats bed file to resource directory
rule link_repeats:
  input: abspath(config['ref']['repeats_path'])
  output: reference_repeats
  shell: "ln {input} {output}"

#----------------------------------------------------------------------------------------------------------------------#

# Invert repeats (isolate regions exclusive of repeat regions)
rule invert_repeat_regions:
  input:
    repeats="{path}/{prefix}.repeats{extra}bed"#,reference_repeats,
    fai=reference_genome_path + ".fai"
  output: "{path}/{prefix}.invert_repeats{extra}bed"#inverted_repeats
  conda: f"{workflow_dir}/envs/get_methylation.yaml"
  shell:
      """
      cat {input.repeats} | \
        awk '{{ split($1, a, "_"); gsub("chr", "", a[1]); print a[1] "\t" $1 "\t" $2 "\t" $3 }}' | \
        sed 's/^X/23/' | \
        sed 's/^Y/24/' | \
        sed 's/^M/25/' | \
        sed 's/^Un/26/' | \
        sort -k 1,1n -k 2,2 -k 3,3n | \
        cut -f2-4 > {input.repeats}.tmp

      cat {input.fai} | \
        awk '{{ split($1, a, "_"); gsub("chr", "", a[1]); print a[1] "\t" $1 "\t1\t" $2 }}' | \
        sed 's/^X/23/' | \
        sed 's/^Y/24/' | \
        sed 's/^M/25/' | \
        sed 's/^Un/26/' | \
        sort -k 1,1n -k 2,2 -k 3,3n | \
        cut -f2-4 > {input.repeats}.ref.bed

      bedtools subtract -a {input.repeats}.ref.bed -b {input.repeats}.tmp > {output}
      rm {input.repeats}.tmp {input.repeats}.ref.bed
      """
#----------------------------------------------------------------------------------------------------------------------#

# Bed telomere setions
rule telomere_bed:
  input: reference_repeats
  output: reference_repeats.replace('.bed', '.telomere.bed')
  conda: f"{workflow_dir}/envs/check_genotypes.yaml"
  log: log_file(reference_repeats.replace('.bed', '.telomere.log'))
  priority: 8
  shell:
    """
    exec > {log}; exec 2> {log}
    grep TAACCC {input} | bedtools sort -i stdin | awk '($1 !~ "_")' > {output}
    """
