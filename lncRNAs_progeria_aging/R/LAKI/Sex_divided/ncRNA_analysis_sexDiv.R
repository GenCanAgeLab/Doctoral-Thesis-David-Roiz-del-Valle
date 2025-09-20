# Load gene metadata ----

geneMetadata <- read_delim("Data/geneMetadata.tsv", delim = "\t")

# Load results and lncRNA annotations ----

resultsList <- readRDS("Out/LAKI/DESeq_results_sexDiv_LAKI.rds")
# lncRNA_mouseDB <- readRDS("Data/resources/lncRNA_mouseDB.rds")

# ncRNA results ----

ncRNAList <- lapply(resultsList, function(x) {
  x <- as_tibble(x, rownames = "gene_id") |>
    filter(gene_type != "protein_coding") |>
    # Also delete IG genes/pseudogenes
    filter(!str_detect(gene_type, "IG"))
})

protCoddf <- lapply(resultsList, function(x) {
  x <- as_tibble(x, rownames = "gene_id") |>
    filter(gene_type == "protein_coding")
}) |>
  bind_rows(.id = "Sex_Tissue") |>
  separate_wider_delim(Sex_Tissue, names = c("Sex", "Tissue"), delim = " ")

write_tsv(protCoddf, file = "Out/LAKI/pcGenes_LAKI.tsv")


## Annotate results ----

# trans interactions lncRNAs
# annotatedResults <- lapply(ncRNAList,
#                            function(res) {
#                              res |>
#                                left_join(lncRNA_mouseDB,
#                                          by = join_by(gene_name == ncName))
#                            })

# Data merge in a single data table and cis genes annotation
allNcRNA <- bind_rows(ncRNAList, .id = "Sex_Tissue") |>
  separate_wider_delim(Sex_Tissue, names = c("Sex", "Tissue"), delim = " ")

write_tsv(allNcRNA, file = "Out/LAKI/ncRNA_LAKI.tsv")

# Tissue-specific ncRNA-pcGene pairs

ncPcdf <- allNcRNA |>
  left_join(
    read_tsv('Out/LAKI/ncRNA_Cis_pcGenes.tsv.gz'), # Load cis genes info
    by = join_by(gene_id == ncRNA),
    relationship = "many-to-many"
  ) |>
  left_join(
    protCoddf,
    by = join_by(target_id == gene_id, target_name == gene_name, Tissue, Sex),
    suffix = c("_ncRNA", "_target")
  ) |>
  dplyr::select(
    Tissue,
    Sex,
    gene_id,
    gene_name,
    gene_type_ncRNA,
    log2FoldChange_ncRNA,
    padj_ncRNA,
    target_id,
    target_name,
    log2FoldChange_target,
    padj_target
  )

write_tsv(ncPcdf, file = "Out/LAKI/ncRNA_cisGenes_TissueSex_LAKI.tsv.gz")


# Highly significant lncRNAs ----

# First, we subset those genes that change in at least one tissue
# and get their info from every tissue
hsignNcRNA <- allNcRNA |>
  filter(padj < 0.05 & abs(log2FoldChange) > 1) |> #Select significant changes
  distinct(gene_id) |> # Get gene_id of changing ncRNAs
  left_join(
    allNcRNA, # Now re-join so we get info for these highly sign ncRNA from all tissues
    by = join_by(gene_id)
  ) |>

  # We will also add info about their direction of change in every tissue
  mutate(
    # We are going to codify changes in ncRNA  genes
    # 1 for sign UP and -1 for sign DOWN. 0 for NS change
    ncRNA_change = case_when(
      log2FoldChange > 0 & padj < 0.05 ~ 1,
      log2FoldChange < 0 & padj < 0.05 ~ -1,
      .default = 0
    ),
  )

# Now we want stats about those genes: in which and how many tissues they change
hsignNcRNAstats <- hsignNcRNA |>
  group_by(Sex, gene_id, gene_name, gene_type, ncRNA_change) |>
  summarise(
    signTissuesN = n(),
    signTissues = paste0(Tissue, collapse = ", "),
    medianLFC = median(log2FoldChange, na.rm = T)
  ) |>
  filter(ncRNA_change != 0) |>
  arrange(desc(signTissuesN), desc(medianLFC))


### Cis interactions ----

# We first merge the stats info into the hsign df
hsignNcRNA_cisGenes <- hsignNcRNA |>
  left_join(hsignNcRNAstats) |>

  # Then, we are adding info about cis ProtCod genes
  left_join(
    read_tsv('Out/LAKI/ncRNA_Cis_pcGenes.tsv.gz'), # Load cis genes info
    by = join_by(gene_id == ncRNA),
    relationship = "many-to-many"
  ) |>

  # We add prot cod expression info
  left_join(
    protCoddf |>
      dplyr::select(Sex, Tissue, gene_id, log2FoldChange, padj),
    by = join_by(target_id == gene_id, Tissue, Sex),
    suffix = c("", "_target"),
    relationship = "many-to-many"
  ) |>

  mutate(
    # We are going to codify changes in cis protCod genes
    # 1 for sign UP and -1 for sign DOWN. 0 for NS change

    protCod_change = case_when(
      log2FoldChange_target > 0 & padj_target < 0.05 ~ 1,
      log2FoldChange_target < 0 & padj_target < 0.05 ~ -1,
      .default = 0
    ),

    # Then, we are going to create a change direction variable
    # 1 for UP-UP/DOWN-DOWN and -1 to UP-DOWN/DOWN-UP
    changeCorr = ncRNA_change * protCod_change
  ) |>

  group_by(Sex, gene_id, gene_name, target_id, target_name) |>
  mutate(corrPower = sum(changeCorr)) |>
  filter(changeCorr != 0) |>
  group_by(
    Sex,
    gene_id,
    gene_name,
    target_id,
    target_name,
    corrPower,
    changeCorr,
    signTissues,
    signTissuesN,
    medianLFC
  ) |>
  summarise(
    corrTissuesN = n(),
    corrTissues = paste0(Tissue, collapse = ", "),
    medianLFC_target = median(log2FoldChange_target, na.rm = T)
  ) |>
  arrange(desc(corrPower), desc(medianLFC_target))


# Save data ----

write_tsv(
  hsignNcRNA_cisGenes,
  file = "Out/LAKI/HighlySign_NcRNAs_cisProtGenes_sexDiv_LAKI.tsv"
)

#Remove environment
rm(list = setdiff(ls(), lsf.str()))
