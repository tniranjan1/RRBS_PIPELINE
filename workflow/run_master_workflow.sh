#!/bin/bash

## This is a master workflow initiator script.
## This is the only script that needs to be run.
## This script will check conda installations
## and then start the snakemake workflow.

echo "Starting RRBS workflow."
## Is conda installed?
if ! type conda > /dev/null
then
  echo "Conda not installed or not aliased. Exiting."
  exit -1
else
  condaver=$(conda --version)
  echo "$condaver installed." | sed 's/conda/Conda/'
  CONDA_BASE="$(conda info --base)/etc/profile.d/conda.sh"
  source $CONDA_BASE
fi

## Is RRBS conda environment alrady created?
RRBS_made=$(conda env list | cut -f1 -d " " | grep 'RRBS' | wc -l)
if [ "$RRBS_made" -eq "0" ]
then
  ## yes, activate RRBS env
  conda create -n RRBS python=3.8 snakemake
  echo "Creating RRBS conda environment."
  conda activate RRBS
else
  ## no, install and activate RRBS env
  conda activate RRBS
  echo "RRBS conda environment already created."
fi
echo "Activating RRBS conda environment."

echo "Starting Snakemake workflow. See ./stderr.log and ./stdout.log for logs."
snakemake --profile ./profile/config.yaml
