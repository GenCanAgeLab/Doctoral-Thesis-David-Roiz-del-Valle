# RNASEQ LAKI

This is a large dataset of tissues from WT and LAKI mice. 
Sequencing was done by novogene using stranded lncRNA library (Ribo-depeletion) with 12 Gbs coverage
Clea prepared the samples and novogene did the RNA extraction, library prep and sequencing. 
Data was sent in a hard drive and then was backed up in CLONAS as a tar file in Fastq folder
Project code is: X204SC23073534-Z01-F001




### Inspection
First, I copy the folder to Chronos to run the analysis. Before doing anything, I do a md5sum check to verify that none of the fastq files are corrupt.
Inside the X204SC23073534-Z01-F001 foder, there is a md5.txt file with the md5 hash and the relative path to the files. 

```Bash
md5sum -c MD5.txt > ../md5sum_check.txt
```
All are OK, so I go ahead

## Test with sample file
Since this dataset is huge, I want to try to run the trimming of adapters and alignment without and intermediate step. I will try then to used named pipes to transfer data between commands. 

### First, I check the adapters are found:

The 3' adaptor provided in the report is: 5'-GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG-3'
This adapter must be found in the pair 1 read without any modification, but an A must be prepended to account for the T/A ligation

The 5' adaptor provided in the report is: 5'-AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT-3'

This adaptor must be found reverse complementary in the pair 2


```Bash
# Check 3' adaptor
gzip -dc X204SC23073534-Z01-F001/01.RawData/Tube_1/Tube_1_EKRO230003864-1A_HMGWVDSX7_L4_1.fq.gz | grep --color=auto  AGATCGGAAGAGCACACGTCTGAACTCCAGTCA

# Check 5' adaptor
gzip -dc X204SC23073534-Z01-F001/01.RawData/Tube_1/Tube_1_EKRO230003864-1A_HMGWVDSX7_L4_2.fq.gz | grep --color=auto AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

```

The 3' adaptor sequences are found as expected. If I search the full adaptor or truncated versions is not found. I have to trim until the index sequence is removed. then, I can see that there are common parts sorrounding the index. I checked several tubes and the index seems to be 11-nt in length, which does not match with the unaligned part of the provided adaptor (9nt). Anyway, it seems that I understood what I should use for trimming.

AGATCGGAAGAGCACACGTCTGAACTCCAGTCA  CGGATGACT  ATCTCGTATGCCGTCTTCTGCTTG 3' adaptor
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA CGCTCTTAGAC ATCTCGTATGCCGTCTTCTGCTTG pair1 tube 1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA CCATACAGACC ATCTCGTATGCCGTCTTCTGCTTG pair1 tube 2
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA CGAAGACCTAC ATCTCGTATGCCGTCTTCTGCTTG pair1 tube 3
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA CGCAGTAATCG ATCTCGTATGCCGTCTTCTGCTTG pair1 tube 4

The 5' adaptor does not behave as expected. It should be detected at the 3' end of the pair2 reads in reverse complement orientation. However, if I do this search I get no matches. However, if I use it without reverse complement, I get hits that match the adaptor. Therefore, it seems that they have not provided with the right adaptor sequence (is already reverse complemented). In addition, the part found in the reads after the index does not match completely the 3' adaptor. To sum up, there is something wrong with provided adaptor, but I know which sequence I should use for trimming.  

AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT  AGATCTCG  GTGGTCGCCGTATCATT                  3' adaptor
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT CAGGTGTATG GTGTAGATCTCGGTGGGCGCCGTATCATT      pair1 tube 1
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT TTCGCTAGCT GTGTAGATCTCGGTGGGCGCCGTATCCTTAAACA pair1 tube 2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT ATCTGACGAG GTGTAGATCTCGGTGGTCGCCGTATCA        pair1 tube 3
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT GATTCAGAGC GTGTAGATCTCGGTGGTCGCCGTATCATT      pair1 tube 4
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT AGCTGCTTAG GTGTAGATCTCGGTGGTCGCCGTATCATTA     pair1 tube 62


Using notepad, I tried to reconstruct the adaptor strategy. It's the classic Y-shaped adapters that have a small complementary region with an overhanging T for T/A ligation. I think they provided the wrong 5' adaptor: it must be its complement. Or it migth be that they provided the sequence that must be found instead of the real oligonucleotide sequences. 

If I do that everything fits. This would be the ligation:

  5' adap comprev NNNNNNNNNNNNNNNNNNNNNNNNNNNNA 3'adaptor w/o modification
  3' adap rev    ANNNNNNNNNNNNNNNNNNNNNNNNNNNN  5' adaptor complement

And after the PCR, the resulting is:

  3' adap comp    NNNNNNNNNNNNNNNNNNNNNNNNNNNNN  5' Adapt
  3' adap rev     NNNNNNNNNNNNNNNNNNNNNNNNNNNNN  5' adapt comp


Therefore, in pair 1 I should detect the 3' adaptor just as the sequence they provided. However, in pair2, I should 

### Check pipeline for one sample
Some samples have more than one file, so I need to somehow group them for the STAR alignment. To do that, the simplest way is to use R with the report file to concatenate them with comma, as required by STAR. 

First, I use cat to paste the content of file "02.Report_X204SC23073534-Z01-F001\src\tables\qc.summary.xls" into a file called report.tsv



```R
setwd('~/20231017_RNASeq_LAKI/')
# Load report
    mydf <- read.delim('Novogene_report.tsv')
# Create file name
  mydf$Fname <- with(mydf, paste(Sample,Library_Flowcell_Lane,'1.fq.gz',sep='_'))
# Create path
  mydf$Fpath <- with(mydf, file.path('X204SC23073534-Z01-F001/01.RawData',Sample,Fname))
  # Check all exist
  all(file.exists(mydf$Fpath))

# Concatenate names
  res <- with(mydf, aggregate(Fpath,list(Sample),paste0,collapse=' '))
# Save data frame
  write.table(res,'Aggregated_samples.tsv',sep='\t',row.names = F,col.names = F,quote=F)  
```

## Run all jobs
I make use of named pipes to connect cutadapt output with STAR input. A problem I have to dealed with is that some samples have more than one file. To solve that, I use R to group de file names and save it as an input file for the bash function. Since cutadapt does not support multiple files as input, I used zcat to unzip and concatenate the reads from multiple files. I used array slicing to automaticaly handle 1-file or multifile samples.


## Group sample files to feed bash script

Some samples have more than one file, so I need to somehow group them. The first program will be cutadapt. Cutadapt does not support multiple files as input. As an alternative, I decided to use zcat to decompress/concatenate the files and pipe them to cutadapt. The simplest way to group the files is to use R with the report file to concatenate them separated by spaces. Each row of the output file have the sample name followed by the corresponding files separated by spaces.

First, I use cat to paste the content of file "02.Report_X204SC23073534-Z01-F001\src\tables\qc.summary.xls" into a file called report.tsv

```r

# Load report
  mydf <- read.delim('Novogene_report.tsv')
# Create file name
  mydf$Fname <- with(mydf, paste(Sample,Library_Flowcell_Lane,'1.fq.gz',sep='_'))
# Create path
  mydf$Fpath <- with(mydf, file.path('X204SC23073534-Z01-F001/01.RawData',Sample,Fname))
  # Check all exist
  all(file.exists(mydf$Fpath))

# Concatenate names
  res <- with(mydf, aggregate(Fpath,list(Sample),paste0,collapse=' '))
# Save data frame
  write.table(res,'Aggregated_samples.tsv',sep='\t',row.names = F,col.names = F,quote=F)  
```

```Bash
# Create named pipes
mkfifo zcat1 zcat2 star1 star2

rm zcat1 zcat2 star1 star2

# Function to run all jobs 
myfun(){
mkfifo zcat1 zcat2 star1 star2
# Prepare input
line=($1)
SampID="${line[0]}"
fq1="${line[@]:1}"
fq2=${fq1//1.fq.gz/2.fq.gz}
# Start reading files
zcat ${fq1[@]}  > zcat1 &
zcat ${fq2[@]}  > zcat2 &
echo "Doing file(s) ${fq1[*]}" &
cutadapt -m 31 \
-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-o star1 -p star2 \
zcat1 zcat2 > "STAR_aligns/${SampID}_cutadapt_stats.txt" &
STAR --runThreadN 20 --genomeDir /data/genomes/GRCm39_GENCODE/STAR_index \
--genomeLoad LoadAndKeep --readFilesIn star1 star2 \
--outFileNamePrefix "STAR_aligns/${SampID}_" \
--outSAMtype BAM Unsorted \
--quantMode GeneCounts
rm zcat1 zcat2 star1 star2
}
export -f myfun
conda activate deepseq
# Test a 2-file sample
  #cat Aggregated_samples.tsv | grep 'Tube_10 ' | parallel myfun
# Do all files
cat Aggregated_samples.tsv | parallel -j 1 --progress myfun
cat Aggregated_subset.tsv | parallel -j 1 --progress myfun

# Remove index from memory
STAR --genomeLoad Remove --genomeDir /data/genomes/GRCm39_GENCODE/STAR_index
```


## Sort and Index BAM files to reduce space

```Bash

myfun(){
  samtools sort -@ 10 -o ${1/.bam/.sorted.bam} $1 && rm $1
}
export -f myfun
parallel -j 2 myfun ::: *.bam

# Index

parallel --progress -j 10 samtools index ::: *.sorted.bam
