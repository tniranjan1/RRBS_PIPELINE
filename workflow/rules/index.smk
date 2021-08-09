# Generate reference dictionary
rule reference_dict:
  input: "{reference}.fa"
  output: "{reference}.dict"
  params:
    build=config['ref']['hg_build'],
    species=config['ref']['species']
  conda: f"{workflow_dir}/envs/index.yaml"
  log: f"{workflow_dir}/logs/index_rules/reference_dict.log"
  shell: "samtools dict -a {params.build} -s {params.species} -o {output} {input} 2> {log}"

# Generate reference index
rule reference_idx:
  input: "{reference}.fa"
  output: "{reference}.fa.fai"
  conda: f"{workflow_dir}/envs/index.yaml"
  log: f"{workflow_dir}/logs/index_rules/reference_idx.log"
  shell: "samtools faidx {input} 2> {log}"

# Generate bwa-meth index
rule bwa_meth_index:
  input: "{reference}.fa"
  output:
    main="{reference}.fa.bwameth.c2t",
    amb="{reference}.fa.bwameth.c2t.amb",
    ann="{reference}.fa.bwameth.c2t.ann",
    pac="{reference}.fa.bwameth.c2t.pac"
  conda: f"{workflow_dir}/envs/index.yaml"
  log: f"{workflow_dir}/logs/index_rules/bwa_meth_index.log"
  shell: "bwameth.py index {input} 2> {log}"
