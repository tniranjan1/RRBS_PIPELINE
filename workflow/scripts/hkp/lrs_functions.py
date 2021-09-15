import pandas as pd
import snakemake
import os
import re
import sys

#----------------------------------------------------------------------------------------------------------------------#

# Distinguish between lrs and rrbs sample names
def distinguish_lrs_rrbs(middle_folder, path_prefix, path_suffix, lrs_sample_names, rrbs_sample_names):
    if middle_folder == 'lrs_methyl':
        use = lrs_sample_names
    else:
        use = rrbs_sample_names
    expansion = []
    for u in sample: expansion.append(path_prefix + u + path_suffix)
    return expansion

#----------------------------------------------------------------------------------------------------------------------#
