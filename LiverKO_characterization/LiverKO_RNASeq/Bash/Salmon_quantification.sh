#!/bin/bash

# Make dirs

mkdir Salmon_quants

mkdir Trimming_logs

# Function to run all jobs 

myfun(){

# Prepare input

SampID=$(basename $1 | sed 's/_1.fastq.gz//')
fq1=$1
fq2=${fq1/1.fastq.gz/2.fastq.gz}

#Create pipes

mkfifo ${SampID}_s1 ${SampID}_s2

# Run cutadapt + salmon

echo "Doing file ${fq1}"

cutadapt -m 31 \
-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-o ${SampID}_s1 -p ${SampID}_s2 \
${fq1} ${fq2} > "Trimming_logs/${SampID}_cutadapt_stats.txt" &

salmon quant -l ISR -p 15 \
-i /data/genomes/GRCm39_GENCODE/salmon_index_M32/ \
-1 ${SampID}_s1 -2 ${SampID}_s2 \
--gcBias \
-o "Salmon_quants/${SampID}" && echo "done"

rm ${SampID}_s1 ${SampID}_s2
}

export -f myfun

micromamba activate salmon


# Do all files

nice parallel -j 2 --progress myfun ::: FASTQ/*_1.fastq.gz