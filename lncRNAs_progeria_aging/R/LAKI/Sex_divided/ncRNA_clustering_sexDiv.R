# Load libraries ----

# library(rtracklayer)
# library(GenomicRanges)

# Load Gencode ----

mygtf <- import(
  '/data/genomes/GRCm39_GENCODE/M32/gencode.vM32.primary_assembly.annotation.gtf'
)

# Keep only genes
mygtf <- mygtf[mygtf$type == 'gene']

# Load high significant ncRNAs ----

hsdf <- read_tsv('Out/LAKI/HighlySign_NcRNAs_cisProtGenes_sexDiv_LAKI.tsv')
ncdf <- read_tsv("Out/LAKI/ncRNA_cisGenes_TissueSex_LAKI.tsv.gz")

# Cluster lncRNAs ----

#Filter GTF to leave only highly significant ncRNAs
hsgtf <- mygtf[mygtf$gene_id %in% hsdf$gene_id, ]


## Cluster with resize (equivalent to bedtools cluster) ----

maxgap = 1E6 # Window to clusterize. We leave it in 1Mb

hsclusts <- reduce(
  resize(hsgtf, fix = 'center', width = width(hsgtf) + maxgap),
  ignore.strand = T
)


## Get ncRNAs for every cluster ----

#Use findOverlaps to match hsclusts and hsgtf
overs <- findOverlaps(hsclusts, hsgtf, ignore.strand = T)
clustdf <- tibble(
  cluster = as.character(hsclusts[from(overs)]), #Clusters for every overlap
  gene_id = hsgtf[to(overs)]$gene_id #ncRNA gene_id for every overlap
)


# Merge info ----

hsdf <- hsdf |>
  left_join(clustdf, by = join_by(gene_id))

ncdf <- ncdf |>
  left_join(clustdf, by = join_by(gene_id))

# Cluster stats ----

hsdf <- hsdf |>
  group_by(cluster, Sex) |>
  mutate(
    clusterNgenes = n_distinct(gene_id), # Calculate how many ncRNAs in cluster
    medianLFCnc_cluster = median(unique(medianLFC)), # Median of all tissues-median LFC of hs ncRNAs in cluster
    medianLFCpc_cluster = median(unique(medianLFC_target)) # Median of all-tissues-median LFC of pc genes in cluster
  )

ncdf <- ncdf |>
  group_by(cluster) |>
  mutate(
    clusterNgenes = n_distinct(gene_id), # Calculate how many ncRNAs in cluster
    medianLFCnc_cluster = median(unique(log2FoldChange_ncRNA)), # Median of all tissues-median LFC of hs ncRNAs in cluster
    medianLFCpc_cluster = median(unique(log2FoldChange_target)) # Median of all-tissues-median LFC of pc genes in cluster
  )

# Save
write_tsv(
  hsdf,
  'Out/LAKI/HighlySign_NcRNA_cisProtGenes_clusters_sexDiv_LAKI.tsv'
)
write_tsv(
  ncdf,
  file = "Out/LAKI/ncRNA_cisProtGenes_clusters_TissueSex_LAKI.tsv.gz"
)

#Remove environment objects but no functions
rm(list = setdiff(ls(), lsf.str()))
