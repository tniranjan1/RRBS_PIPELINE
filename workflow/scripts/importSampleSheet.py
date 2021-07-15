import pandas as pd
import snakemake
import re
import sys

def importSampleSheet(sample_path, schema_path):
    sample_sheet = pd.read_table(sample_path)
    try: snakemake.utils.validate(sample_sheet, schema=schema_path)
    except snakemake.exceptions.WorkflowError as e:
        print(e)
        print("Validation failed on ", sample_path, ". Workflow aborted.", sep="")
        sys.exit(-1)
    # Name each sample with the format: SampleID.SampleGroup.Tissue
    # Set sample name as index (row names) for rrbs_samples dataframe
    sample_names = []
    for i in range(0, len(sample_sheet)):
        new_name = sample_sheet[['SampleID', 'SampleGroup', 'Tissue']].iloc[0].str.cat(sep=".")
        slugify = re.sub(r'(?u)[^-\w.]', '', new_name)
        sample_names.append(new_name)
        if slugify != new_name:
            print("For SampleID.SampleGroup.Tissue = ", new_name, ", invalid character(s) detected.", sep="")
            print("Valid characters are Aa-Zz0-9.-_.")
            print("Remove invalid characters from ", sample_path, ". Workflow aborted.", sep="")
            sys.exit(-1)
    sample_sheet.index = sample_names
    return sample_sheet
