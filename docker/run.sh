#!/bin/bash

export PATH="/opt/miniconda/bin:$PATH"
source activate snakemake
cd /opt/snakemake_rundir && snakemake -p --cores "$DCKR_THREADS" --resources ram=360

