import pandas as pd
import snakemake
import os
import re
import sys

#----------------------------------------------------------------------------------------------------------------------#

# Short write for absolute path
def abspath(path):
    return os.path.abspath(path)

#----------------------------------------------------------------------------------------------------------------------#

# Subroutine to correctly import and validate a sample sheet for analysis.
def importSampleSheet(sample_path, schema_path):
    # Import sample sheet as pandas dataframe
    sample_sheet = pd.read_table(sample_path)

    # Convert appropriate columns to string, if not already in string format
    stringed_columns = ['SampleID','SampleGroup','CaseControl','Tissue','Path','Remove_Reason']
    for c in stringed_columns:
        sample_sheet[c] = sample_sheet[c].astype('str')

    try: snakemake.utils.validate(sample_sheet, schema=schema_path)
    except snakemake.exceptions.WorkflowError as e:
        print(e)
        print("Validation failed on ", sample_path, ". Workflow aborted.", sep="")
        sys.exit(-1)

    # Name each sample with the format: SampleGroup.Tissue.SampleID
    sample_names = []
    for i in range(0, len(sample_sheet)):
        new_name = sample_sheet[['SampleGroup', 'Tissue', 'SampleID']].iloc[i].str.cat(sep=".")
        # Ensure new sample name can be used as filename (no invalid characters)
        slugify = re.sub(r'(?u)[^-\w.]', '', new_name)
        sample_names.append(new_name)
        # Print error and abort if there are invalid characters
        if slugify != new_name:
            print("For SampleGroup.Tissue.SampleID = ", new_name, ", invalid character(s) detected.", sep="")
            print("Valid characters are Aa-Zz0-9.-_.")
            print("Remove invalid characters from ", sample_path, ". Workflow aborted.", sep="")
            sys.exit(-1)

    # Set sample name as index (row names) for rrbs_samples dataframe
    sample_sheet.index = sample_names
    return sample_sheet

#----------------------------------------------------------------------------------------------------------------------#

def mergeSampleSheet(sheetA, sheetB):
    mergedSheet = sheetA[['SampleID', 'Path']].append(sheetB[['SampleID', 'Path']])

    #get source type (SRR, BAM, fastq)
    patterns = {
        "fq1$" : "fq",
        "fq2$" : "fq",
        "fq1.gz$" : "fq",
        "fq2.gz$" : "fq",
        "fastq" : "fq",
        "fastq.gz" : "fq",
        "bam$" : "bam",
        "^SRR" : "SRR"
    }
    type_list = []
    for i in range(0, len(mergedSheet)):
        # match path to type
        path = mergedSheet['Path'].iloc[i]
        value = []
        for key in patterns:
            if re.search(key, path):
                value.append(patterns[key])
        if len(value) == 1:
            type_list.append(value)
        else: # Abort with error specifying invalid or too many path types.
            if len(value) == 0:
                print("Usable path type  not found for SampleGroup.Tissue.SampleID = {}.".format(mergedSheet.index[i]))
                print("Ensure path to sample reads is SRR***, ***.bam, or ***.fq1(.gz),***.fq2(.gz).")
            else:
                print("Too many  path types found for SampleGroup.Tissue.SampleID = {}.".format(mergedSheet.index[i]))
                print("Currently, the workflow can only accept one path type (SRR, BAM, fastq) per sample.")
            print("Workflow aborted.")
            sys.exit(-1)
    type_list = [ y for x in type_list for y in x ]
    mergedSheet['type'] = type_list
    # Return dataframe for each sample with path and path type.
    return mergedSheet

#----------------------------------------------------------------------------------------------------------------------#

# Subroutine to obtain the names of all inital output files (alignments) for a given sample sheet
# output name = "../resources/{sample_destinition from function input}/alignments/{sample name from sample sheet}.bam"
def get_initial_output(sample_sheet, sample_destination):
    output_files = []
    prefix = "../results/" + sample_destination + "/alignments/"
    suffix = ".POSsort.bam"
    for i in range(0, len(sample_sheet)):
        output_files.append(os.path.abspath(prefix + sample_sheet.index[i] + suffix))
    return output_files

#----------------------------------------------------------------------------------------------------------------------#

# Subroutine to obtain list of merged methylation or merged coverage files that have been separated by chr
def get_merge_list(fai, wildcards):
    path=wildcards.path
    meco=wildcards.meco
    repeats=wildcards.repeats
    context=wildcards.context
    if os.path.isfile(fai): # get names of chromosomes
        fai_table = pd.read_table(fai, header=None)
        chr_list = [ chr for chr in fai_table[0].tolist() if not ('_' in chr) ]
    else:
        chr_list = "chrM"
    file_template = f"{path}/methylation_calls/merged/merged_{meco}.{repeats}.{context}"
    file_list = [ f"{file_template}.{chr}.bedGraph" for chr in chr_list if ('c' in chr) ]
    return file_list

#----------------------------------------------------------------------------------------------------------------------#

# Subroutine to get max allowed disk usage
def get_disk_gb(source_type):
    if source_type == 'SRR':
        return 50
    elif source_type == 'bam':
        return 10
    elif source_type == 'fq':
        return 1
    return 1

#----------------------------------------------------------------------------------------------------------------------#

# Subroutine to restrict to full chromosomes in reference
def chromosome_constraint(reference_genome_path):
    fai=reference_genome_path + '.fai'
    if os.path.isfile(fai): # get names of chromosomes
        fai_table = pd.read_table(fai, header=None)
        chr_list = [ chr for chr in fai_table[0].tolist() if not ('_' in chr) ]
    else:
        chr_list = "chrM"
    return chr_list
