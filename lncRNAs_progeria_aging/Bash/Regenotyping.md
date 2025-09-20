# Genotyping
I will call for variants in the LMNA gene to confirm the genotypes. I use BCFTOOLS mpileup and call. Before that, I have to index the reference genome fasta file. To do that, I have to first generate a bgzip compress fasta file and the run samtools faix

```bash
# Generate fasta index
gzip -dc GRCm39.primary_assembly.genome.fa.gz | bgzip -c > GRCm39.primary_assembly.genome.fa.bgz
samtools faidx GRCm39.primary_assembly.genome.fa.bgz

# Variant calling
ls -1 STAR_aligns/*.bam > bamlist.txt
bcftools mpileup -Ou -q 255 -f /data/genomes/GRCm39_GENCODE/GRCm39.primary_assembly.genome.fa.bgz -r chr3:88388455-88410642 -b bamlist.txt | bcftools call -mv -Ov -o RNASeq_LAKI_genotypes.vcf
```
Then, I R:
```R
require(VariantAnnotation)
setwd('~/20231017_RNASeq_LAKI/')
# Read VCF
myvcf <- readVcf('RNASeq_LAKI_genotypes.vcf')
# Extract genotypes for LAKI mutation
genos <- geno(myvcf)$GT['chr3:88389797_G/A',]
names(genos) <- gsub('_Aligned.*','',basename(names(genos)))
# Merge to sample info
sampinfo <- read.delim('R_analysis/Info/Sample_info_complete.tsv')
sampinfo$InferedGeno <- factor(genos[sampinfo$SampleID],
                               levels=c('0/0','0/1','1/1'),
                               labels=c('WT','HT','KO'))
colnames(sampinfo) <- c('SampleID','MouseNo','Genotype','Gender','Tissue',
                        'DAB','Notes','Infered_Genotype')
sampinfo[sampinfo$Genotype != sampinfo$Infered_Genotype,]
write.table(sampinfo,
            'R_analysis/Info/Sample_info_complete_genotiped.tsv',
            col.names =T, sep='\t',quote=F)

```

**There are two samples wrongly labeled**
SampleID MouseNo Genotype Gender Tissue DAB              Notes
Tube_90     913       WT   Male  Colon 150                   
Tube_101     991       KO Female  Heart 153  
