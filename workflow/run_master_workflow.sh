#!/bin/bash

# This is a master workflow initiator script. This is the only script that needs to be run. This script will check
# conda installations and then start the snakemake workflow.

std_log="./stdout.log"
err_log="./stderr.log"

echo -e  "Here you will find the standard output of the main workflow.\n\n---\n\n" > $std_log
echo -e "Here you will find the standard error of the main workflow.\n\n---\n\n" > $err_log

echo "Starting RRBS workflow." >> $std_log 2>> $err_log
## Is conda installed?
if ! type conda > /dev/null
then
  echo "Conda not installed or not aliased. Exiting." >> $std_log 2>> $err_log
  exit -1
else
  condaver=$(conda --version)
  echo "$condaver installed." | sed 's/conda/Conda/' >> $std_log 2>> $err_log
  CONDA_BASE="$(conda info --base)/etc/profile.d/conda.sh"
  source $CONDA_BASE
fi

## Is RRBS conda environment alrady created?
RRBS_made=$(conda env list | cut -f1 -d " " | grep 'RRBS' | wc -l)
if [ "$RRBS_made" -eq "0" ]
then
  ## yes, activate RRBS env
  conda create -n RRBS python=3.8 snakemake mamba
  echo "Creating RRBS conda environment." >> $std_log 2>> $err_log
  conda activate RRBS
else
  ## no, install and activate RRBS env
  conda activate RRBS
  echo "RRBS conda environment already created." >> $std_log 2>> $err_log
fi
echo "Activating RRBS conda environment." >> $std_log 2>> $err_log

echo "Starting Snakemake workflow. See ./stderr.log and ./stdout.log for logs."
snakemake --profile profile/ >> $std_log 2>> $err_log &
