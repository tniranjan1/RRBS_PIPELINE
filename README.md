# Snakemake workflow: RRBS of Brain Tissue

[![Python](https://img.shields.io/badge/python-=3.8.10-brightgreen.svg)](https://docs.python.org/3.8/)
[![Snakemake](https://img.shields.io/badge/snakemake-=6.5.3-blueviolet.svg)](https://snakemake.github.io)
[![Conda](https://img.shields.io/badge/conda-=4.10.3-blue.svg)](https://docs.conda.io/en/latest/)
[![MIT License](https://img.shields.io/apm/l/atomic-design-ui.svg?)](https://github.com/tniranjan1/RRBS_PIPELINE/blob/master/LICENSE.md)

---
### Introduction
The methylome is the complement of DNA methylation, particularly at CpG dinucleotides, in cells/tissues that influence gene expression. Classically, hypermethylation at enhancers, gene promotors, and CpG islands locally promotes closed chromatin and reduced gene expression. Hypomethylation is associated with the opposite. Aberrant methylation has been associated with various brain diseases. For example, hypermethylation at a CGG repeat expansion near the gene *FMR1* results in decreased *FMR1* expression, resulting in Fragile X Syndrome. Loss-of-function mutations in *MeCP2*, which is a maintenance gene for DNA methylation, results in genome-wide methylation defects, causing Rett Syndrome.

Aberrant methylation is also associated with neurodegenerative disorders. For example, a hexanucleotide repeat expansion at *C9orf72* is the most common known genetic cause of ALS/FTLD (Amyotrophic Lateral Sclerosis / Frontotemporal Lobar Degeneration). Expansion of this repeat is associated with hypermethylation of the *C9orf72* promoter region, decreased *C9orf72* expression, and toxic aggregates (Jackson JL 2020). Los-of-function mutations in the Progranulin (*GRN*)gene are associated with a TDP-43 positive form of FTLD. The *GRN* promoter is more frequently hypermethylated in peripheral blood mononuclear cells of patients with sporadic FTLD, and this is associated with decreased gene expression, but this effect is not seen in patients with Alzheimer's disease (Banzhaf-Strathmann J 2013, Cooper YA 2018). *APOE* hypomethylation is associated with Lewy Body Dementia (LBD), and *DRD2* hypomethylation is associated with LBD and Parkinson Disease (PD).

To further characterize methylation patterns in FTLD, we are analyzing RRBS (reduced representation bisulfite sequencing) data generated from post-mortem brain tissue of >100 individuals with various forms of FTLD, as well as controls. This snakemake workflow and associate analytical scripts are made available for the purposes of data reproducibility and extension to other bisulfite-sequencing datasets in the brain.

---
### Datasets
FTLD RRBS datasets used for this project may be available after publications.
Access to reference and publicly available datasets are described in the [config](https://github.com/tniranjan1/RRBS_PIPELINE/blob/master/config/config.yaml) file.

---
### System Requirements
Worklfow and all analyses were run and tested on local server running [CentOS Linux release 7.9.2009](https://www.centos.org/).
Minimum CPUs used for most analyses is 4. Max memory usage is 12 GB, but may be higher for more CPU availability.

---
### Dependencies
All packages and dependencies are managed using conda.
See [conda installation](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) options before proceeding with this workflow.
The following packages will be managed by conda and do not need to be pre-installed.

     1. samtools [version='>=1.10']
     2. bedtools
     3. r [version='>=4.0.3']
     4. gatk4
     5. bowtie2
     6. bwa
     7. maven
     8. bcftools
     9. bwameth
    10. methyldackel
    11. bioconductor-methylkit
    12. bioconductor-edger
    13. bioconductor-annotationdbi
    14. bioconductor-org.hs.eg.db
    15. bioconductor-rnbeads
    16. git
    17. pip

### Installation
In your desired conda environment (e.g.: base) git and pip can be installed by:

    $ conda install -c anaconda git
    $ conda install -c anaconda pip
In your desired RRBS workflow directory, RRBS_PIPELINE can be installed by:

    $ pip install git+git://github.com/tniranjan1/RRBS_PIPELINE.git

### Workflow Configuration
Modify as appropriate snakemake arguements as described in the [profile](https://github.com/tniranjan1/RRBS_PIPELINE/blob/master/workflow/profile/config.yaml) file.

Modify as appropriate configuration parameters and sample sheet files as described in the [config](https://github.com/tniranjan1/RRBS_PIPELINE/blob/master/config/config.yaml) file.

### Citation
Please cite this GitHub link until such time as publication of our results. Then please cite the publication (DOI to be provided here).

### Contact
For advise, bug reports, or assistance please post to [RRBS_PIPELINE Issues](https://github.com/tniranjan1/RRBS_PIPELINE/issues).
