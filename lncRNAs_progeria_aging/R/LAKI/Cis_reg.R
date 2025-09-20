
# Load Gencode ----

mygtf <- import('/data/genomes/GRCm39_GENCODE/M32/gencode.vM32.primary_assembly.annotation.gtf')



# ncRNA subset ----

# Keep only ncRNA genes
mygtf <- mygtf[mygtf$type == 'gene']
ncgtf <- mygtf[mygtf$gene_type != 'protein_coding']



# protein coding subset ----

pcgtf <- mygtf[mygtf$gene_type == 'protein_coding']



# Get cis genes ----

# Get nearby genes and distance (from TSS to TSS)
maxgap = 1E6 # 1 Mb

# Calculate overlaps between ncRNAs and protein-coding
overs <- findOverlapPairs(resize(ncgtf,1,'start'), # this resize reduces the genes to 1, thus getting the TSS of the lncRNA
                      resize(pcgtf,1,'start'), # this gets the TSS of the protein coding genes
                      maxgap=maxgap, # max distance of 1 Mb
                      ignore.strand=T)

# Now we create a tibble with that info
cispcgenes <- tibble(
  ncRNA = first(overs)$gene_id,
  ncRNA_name = first(overs)$gene_name,
  ncRNA_TSS= as.character(first(overs)),
  target_id = second(overs)$gene_id,
  target_name = second(overs)$gene_name,
  target_TSS=as.character(second(overs)),
  distance = distance(first(overs),second(overs), ignore.strand=T),
  relDist = start(second(overs)) - start(first(overs)) #Distance with sign to indicate direction of cis gene
)


# Save
write_tsv(cispcgenes,'Out/LAKI/ncRNA_Cis_pcGenes.tsv.gz')
