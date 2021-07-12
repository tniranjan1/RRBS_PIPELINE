#!/bin/bash

echo "Starting RRBS workflow"
## Is conda installed?
if ! type conda > /dev/null
  then
    echo "Conda not installed or not aliased. Exiting."
    exit -1
else
  condaver=$(conda --version)
  echo "$condaver installed" | sed 's/conda/Conda/'
fi

## Is RRBS conda environment alrady created?
  ## yes, activate RRBS env
  ## no, install and activate RRBS env
