import pandas as pd
import snakemake
import os
import re
import sys

#----------------------------------------------------------------------------------------------------------------------#

# Distinguish between lrs and rrbs sample names
def distingish_lrs_rrbs(wildcards):
    if wildcards.path == 'lrs_methyl':
        use = lrs_sample_names
    else:
        use = rrbs_sample_names
    full_path = results_dir + "/" + wildcards.path + "/samples/{sample}.agePrediction.txt"
    return expand(full_path, sample = use)

#----------------------------------------------------------------------------------------------------------------------#
