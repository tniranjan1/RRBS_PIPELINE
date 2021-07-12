#!/bin/bash

echo "Starting RRBS workflow"
## Is conda installed?
{
  condaver=$(conda --version)
} || {
  "Error: Conda is not installed."
  set -e
}
echo "$condaver installed."

## Is RRBS conda environment alrady created?
  ## yes, activate RRBS env
  ## no, install and activate RRBS env
