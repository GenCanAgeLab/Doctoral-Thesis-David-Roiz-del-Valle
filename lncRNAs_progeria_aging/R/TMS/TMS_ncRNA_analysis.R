# Load gene metadata ----

geneMetadata <- read_tsv("Data/geneMetadata.tsv")

# Load results and lncRNA annotations ----

resultsList <- readRDS("Out/TMS/DESeq_results_TMS_oldvsyoung.rds")
# lncRNA_mouseDB <- readRDS("Data/resources/lncRNA_mouseDB.rds")

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


# # Annotated results ----
#
# annotatedResults <- lapply(ncRNAList,
#                            function(res) {
#                              res |>
#                                left_join(lncRNA_mouseDB,
#                                          by = join_by(gene_name == ncName))
#                            })
# names(annotatedResults) <- tissueList

# Data merge in a single data table and cis genes annotation
# Cis annotation

allNcRNA <- bind_rows(ncRNAList, .id = "Sex_Tissue") |>
  separate_wider_delim(Sex_Tissue, names = c("Sex", "Tissue"), delim = " ")

allPcGene <- bind_rows(pcGeneList, .id = "Sex_Tissue") |>
  separate_wider_delim(Sex_Tissue, names = c("Sex", "Tissue"), delim = " ")

write_tsv(allPcGene, file = "Out/TMS/pcGenes_TMS.tsv.gz")

# Merge corr data and ncRNA DE data

corrData <- read_tsv("Out/TMS/ncRNA_pcGene_Age_corr_TMS.tsv.gz")

ncRNACorr <- corrData |>
  left_join(
    allNcRNA,
    by = join_by(Tissue, Sex, gene_id, gene_name, gene_type),
    relationship = "many-to-many"
  ) |>
  left_join(
    allPcGene |>
      dplyr::select(
        Tissue,
        Sex,
        target_id = gene_id,
        target_name = gene_name,
        log2FoldChange_target = log2FoldChange,
        padj_target = padj
      ),
    by = join_by(Tissue, Sex, target_id, target_name),
    relationship = "many-to-many"
  )

write_tsv(ncRNACorr, "Out/TMS/ncRNA_DE_pcGene_Age_corr_TMS.tsv.gz")


# ## Significance per tissue ----
#
# signNcRNA <- allNcRNA |>
#   filter(padj < 0.05) |>
#   dplyr::select(gene_id, gene_name, gene_type, Tissue) |>
#   group_by(gene_id, gene_name, gene_type) |>
#   summarize(
#     signTissuesN = n(),
#     signTissues = paste0(Tissue, collapse = ", ")
#   )

# # Now we want to add to our results the info about significance in other tissues
#
# allNcRNA <- allNcRNA |>
#   left_join(signNcRNA,
#             by = join_by(gene_id, gene_type, gene_name))

write_tsv(allNcRNA, file = "Out/TMS/ncRNA_TMS.tsv.gz")

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

# Save the data

write_delim(TopNcRNA, file = "Out/TMS/Top_ncRNA_TMS.tsv", delim = "\t")

#Remove environment
rm(list = setdiff(ls(), lsf.str()))
