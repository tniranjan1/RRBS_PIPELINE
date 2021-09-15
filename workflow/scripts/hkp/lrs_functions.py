import pandas as pd
import snakemake
import os
import re
import sys

#----------------------------------------------------------------------------------------------------------------------#

# Distinguish between lrs and rrbs sample names
def distinguish_lrs_rrbs(middle_folder, top_path, lrs, rrbs):
    path_prefix=f"{top_path}/{middle_folder}/samples/",
    path_suffix=".agePrediction.txt",
    if middle_folder == 'lrs_methyl':
        use = lrs
    else:
        use = rrbs
    expansion = []
    for u in sample: expansion.append(path_prefix + u + path_suffix)
    return expansion

#----------------------------------------------------------------------------------------------------------------------#
