##### Using salmon ########

I am concerned that the strict mode of function of STAR/htseq-count migth remove valuable information. I am going to try in parallel a multiple-alignment aware software, like salmon. Another possibility would be to run again STAR with --quantMode transcriptomeSAM, which projects the genomic alignments to the transcriptome. However, given that salmon is so fast, I think I will give a try to this option first. 

I install first the last version of salmon: v.1.10

Then, I build the index using the same files than for STAR:

```bash
# WD: chronos: /data/genomes/GRCm39_GENCODE
cat gencode.vM32.transcripts.fa.gz GRCm39.primary_assembly.genome.fa.gz > gentrome.fa.gz
gzip -dc GRCm39.primary_assembly.genome.fa.gz | grep '^>' | cut -d " " -f 1 | sed 's/>//g' > decoys.txt
salmon index -t gentrome.fa.gz -d decoys.txt -p 12 -i salmon_index --gencode
```

Now, I run salmon for all samples

```bash
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
salmon quant -l ISR -p 10 \
-i /data/genomes/GRCm39_GENCODE/salmon_index \
-1 star1 -2 star2 \
--gcBias \
-o "salmon_quants/${SampID}"
rm zcat1 zcat2 star1 star2
}
export -f myfun

# Test a 2-file sample
  #cat Aggregated_samples.tsv | grep 'Tube_10 ' | parallel myfun
# Do all files
cat Aggregated_samples.tsv | parallel -j 1 --progress myfun