# Load gene metadata ----

geneMetadata <- read_tsv("Data/geneMetadata.tsv")

# Load results and lncRNA annotations ----

resultsList <- readRDS("Out/Zhou/DESeq_results_Zhou_oldvsyoung.rds")
tissueList <- names(resultsList)

# ncRNA results ----

ncRNAList <- lapply(resultsList, function(x) {
  x <- as_tibble(x, rownames = "gene_id") |>
    left_join(geneMetadata, by = join_by(gene_id == gene_id2)) |>
    filter(gene_type != "protein_coding") |>
    # Also delete IG genes/pseudogenes
    filter(!str_detect(gene_type, "IG"))
})

pcGeneList <- lapply(resultsList, function(x) {
  x <- as_tibble(x, rownames = "gene_id") |>
    left_join(geneMetadata, by = join_by(gene_id == gene_id2)) |>
    filter(gene_type == "protein_coding")
})

# Data merge in a single data table and cis genes annotation
# Cis annotation

allNcRNA <- bind_rows(ncRNAList, .id = "Tissue")

allPcGene <- bind_rows(pcGeneList, .id = "Tissue")

write_tsv(allPcGene, file = "Out/Zhou/pcGene_Zhou.tsv.gz")

# Merge corr data and ncRNA DE data

corrData <- read_tsv("Out/Zhou/ncRNA_pcGene_Age_corr_Zhou.tsv.gz")

ncRNACorr <- corrData |>
  left_join(
    allNcRNA,
    by = join_by(Tissue, gene_id, gene_name, gene_type),
    relationship = "many-to-many"
  ) |>
  left_join(
    allPcGene |>
      dplyr::select(
        Tissue,
        target_id = gene_id,
        target_name = gene_name,
        log2FoldChange_target = log2FoldChange,
        padj_target = padj
      ),
    by = join_by(Tissue, target_id, target_name),
    relationship = "many-to-many"
  )

write_tsv(ncRNACorr, "Out/Zhou/ncRNA_pcGene_Age_corr_Zhou.tsv.gz")

## Significance per tissue ----

signNcRNA <- allNcRNA |>
  filter(padj < 0.05) |>
  dplyr::select(gene_id, gene_name, gene_type, Tissue) |>
  group_by(gene_id, gene_name, gene_type) |>
  summarize(
    signTissuesN = n(),
    signTissues = paste0(Tissue, collapse = ", ")
  )

# Now we want to add to our results the info about significance in other tissues

allNcRNA <- allNcRNA |>
  left_join(signNcRNA, by = join_by(gene_id, gene_type, gene_name))


write_tsv(allNcRNA, file = "Out/Zhou/ncRNA_Zhou.tsv.gz")

# Select top ncRNAs per tissue ----

LFCcutoff <- 2

UpregNcRNA <- allNcRNA |>
  filter(log2FoldChange > LFCcutoff) |>
  group_by(Tissue) |>
  slice_min(padj, n = 20)

DownregNcRNA <- allNcRNA |>
  filter(log2FoldChange < -LFCcutoff) |>
  group_by(Tissue) |>
  slice_min(padj, n = 20)

TopNcRNA <- bind_rows(UpregNcRNA, DownregNcRNA)

# Save it

write_delim(TopNcRNA, file = "Out/Zhou/Top_ncRNA_Zhou.tsv", delim = "\t")

# Remove environment
rm(list = setdiff(ls(), lsf.str()))
