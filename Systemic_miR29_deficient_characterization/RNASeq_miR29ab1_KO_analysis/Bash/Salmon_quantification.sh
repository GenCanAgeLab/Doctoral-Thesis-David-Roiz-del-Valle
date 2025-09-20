#!/bin/bash

# Data is already trimmed, so we can directly quantify it with salmon


# Make dirs

mkdir Salmon_quants

# Function to run all jobs 

myfun(){

# Prepare input

SampID=$(basename $1 | sed 's/_1.fq.gz//')
fq1=$1
fq2=${fq1/1.fq.gz/2.fq.gz}

# Run salmon

salmon quant -l ISR -p 15 \
-i /data/genomes/GRCm39_GENCODE/salmon_index_M32/ \
-1 ${fq1} -2 ${fq2} \
--gcBias \
-o data/Salmon_quants/${SampID} && echo "done"

}

export -f myfun

micromamba activate salmon


# Do all files

nice parallel -j 2 --progress myfun ::: /data/miR29/20180518_RNASeq_Liver_Heart_miR29AB_BC/data/*_1.fq.gz