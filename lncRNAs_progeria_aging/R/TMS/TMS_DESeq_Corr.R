Cores <- 6

# Load gene metadata ----

geneMetadata <- read_tsv("Data/geneMetadata.tsv")

# Load gse objects ----

gse <- readRDS(
  "/data/jandrulas/Aging_datasets/TabulaMS_Bulk/Reanalysis/Objects/TMS_reanalysis_se.rds"
)

tissueList <- unique(colData(gse)$Tissue)

# Quality check ----

# mitoReads <- assay(gse) |>
#   as_tibble(rownames = "gene_id") |>
#   left_join(geneMetadata, by = join_by(gene_id == gene_id2)) |>
#   filter(Chr == "chrM") |>
#   summarise(across(contains("SRR"), sum)) |>
#   pivot_longer(everything(), names_to = "Sample", values_to = "totalCounts")

# Number of mapped reads was already computed and included in the SE object
# ggplot(colData(gse), aes(Run, num_mapped)) +
#   geom_col() +
#   labs(title = "Total reads")

# ggplot(mitoReads, aes(Sample, totalCounts)) +
#   geom_col() +
#   labs(title = "Mitochondrial reads")

# Sample filtering ----

# We are going to exclude those samples with less than 1M reads

# Histogram of mapped reads (in M reads)
ggplot(colData(gse), aes(x = num_mapped / 1E6)) +
  geom_histogram() +
  theme_classic() +
  labs(
    title = "Density of mapped reads",
    x = "Mapped reads(in M reads)",
    y = "Number of samples"
  )
ggsave("Out/TMS/Mapped_reads_histogram.pdf", width = 4, height = 3)
# The quantile distribution shows that the median is 11M reads
# 80% of the samples have more than 6.5 M reads
quantile(colData(gse)$num_mapped / 1E6, seq(0, 1, 0.1))

# Thus, we filter those samples with > 5M reads

gse <- gse[, colData(gse)$num_mapped > 5E6]

# Explore the per-sex timepoint distribution

ggplot(colData(gse), aes(x = factor(AGE), col = sex, group = sex)) +
  stat_count(geom = "line") +
  stat_count(geom = "point") +
  facet_wrap(vars(Tissue), scale = "free_x") +
  theme_classic() +
  labs(x = "Timepoint", y = "Number of samples")

ggsave("Out/TMS/Samples_per_timepoint.pdf", width = 10, height = 6)

# Females don't have samples beyond 21 months, and in subcutaneous tissue there are no samples
# older than 15 month. In males, there are some tissues with just 1 or 2 mice per point (e.g. sub_adipose, heart, spleen)
# The low number of samples in males can be corrected by merging nearby ages (e.g 3+6 vs 24+27)
# But for females, we would have to do 3+6 against 18+21 and exclude sub_adipose

# DESeq2 ----

# DESeq old vs young age

# Recode age
gse$Age_recoded <- factor(
  gse$AGE > 12,
  levels = c(T, F),
  labels = c("Old", "Young")
)
gse$Age_recoded <- relevel(gse$Age_recoded, ref = "Young")

# Plot
# ggplot(colData(gse), aes(x=Age_recoded,col=sex,group=sex))+
#   stat_count(geom='bar')+
#   stat_count(geom='point')+
#   facet_wrap(~Tissue)

ggplot(colData(gse), aes(x = Age_recoded, fill = sex, group = sex)) +
  geom_bar(position = position_dodge()) +
  facet_wrap(vars(Tissue))


# Keep 3+6 vs. 24+27 for males. Maybe leave it in 18+21 to simplify?
# Keep 3+6 vs. 18+21 for females

gseSex <- list(
  "Male" = gse[
    ,
    colData(gse)$sex == "male" &
      colData(gse)$AGE %in% c(3, 6, 24, 27)
  ],
  "Female" = gse[
    ,
    colData(gse)$sex == "female" &
      colData(gse)$AGE %in% c(3, 6, 18, 21)
  ]
)

# stats <- bind_rows(
#   as_tibble(colData(gseSex[["Male"]])) |>
#     group_by(sex,Tissue,Age_recoded) |>
#     summarize(n=n()),
#   as_tibble(colData(gseSex[["Female"]])) |>
#     group_by(sex,Tissue,Age_recoded) |>
#     summarize(n=n())
# )

# DESeq

# First, we discard subcutaneous adipose of the analysis
tissueList <- tissueList[tissueList != "Subcutaneus_adipose"]


analysisDESeq <- function(gse, Tissue) {
  tmpGse <- gse[, colData(gse)[, "Tissue"] == Tissue] # First, select the samples from our target tissue

  tmpDds <- DESeqDataSet(tmpGse, design = ~Age_recoded)

  keepRow <- rowSums(counts(tmpDds) >= 10) >= (ncol(tmpDds) / 2)
  tmpDds <- tmpDds[keepRow, ]

  DESeq(tmpDds)
}

deseqList <- mapply(
  analysisDESeq,
  # Repeat each element of gseSex 16 times (once for each tissue)
  rep(gseSex, each = length(tissueList)),
  # Repeat the tissueList for both sexes
  rep(tissueList, times = length(gseSex)),
  SIMPLIFY = FALSE
)

names(deseqList) <- paste(
  rep(names(gseSex), each = length(tissueList)),
  rep(tissueList, times = length(gseSex))
)

resultsList <- lapply(deseqList, results)
names(resultsList) <- names(deseqList)

saveRDS(resultsList, "Out/TMS/DESeq_results_TMS_oldvsyoung.rds")

## Vsd and PCA

vsdList <- lapply(deseqList, vst, blind = TRUE)
# pcaList <- lapply(vsdList, plotPCA, intgroup = "Age_recoded")

# Correlations ----

## ncRNA-Age correlation ----

calc_cor <- function(vsd) {
  expr_matrix <- assay(vsd)
  age_vector <- vsd$AGE
  n_samples <- length(age_vector)

  # Vectorized correlation calculation (much faster than apply)
  estimates <- cor(t(expr_matrix), age_vector, method = "spearman")

  # Vectorized p-value calculation using t-distribution approximation
  t_stat <- estimates * sqrt((n_samples - 2) / (1 - estimates^2))
  p_values <- 2 * pt(-abs(t_stat), df = n_samples - 2)

  tibble(
    gene_id = rownames(estimates),
    p_value = p_values[, 1],
    estimate = estimates[, 1]
  )
}

corList <- mclapply(vsdList, calc_cor, mc.cores = Cores)
names(corList) <- names(resultsList)


corListdf <- corList |>
  bind_rows(.id = "Sex_Tissue") |>
  left_join(geneMetadata, by = join_by(gene_id == gene_id2)) |>
  filter(gene_type != "protein_coding") |>
  # Also delete IG genes/pseudogenes
  filter(!str_detect(gene_type, "IG")) |>
  separate_wider_delim(Sex_Tissue, delim = " ", names = c("Sex", "Tissue"))

write_tsv(corListdf, "Out/TMS/Correlation_results_TMS_all.tsv.gz")

# Age corr stats

corTissueSex <- corListdf |>
  filter(p_value < 0.05) |>
  group_by(Sex, gene_id, gene_name, gene_type, CorrDir = sign(estimate)) |>
  summarise(
    corrTissueN = n(),
    corrTissues = paste0(Tissue, collapse = ", "),
    medianEstimate = median(estimate)
  )


write_tsv(corTissueSex, "Out/TMS/Correlation_results_TMS.tsv")


## ncRNA - cis pcGene correlation ----

ncRNApc_combinations <- read_tsv("Out/LAKI/ncRNA_Cis_pcGenes.tsv.gz") |>
  left_join(
    geneMetadata |>
      dplyr::select(gene_id, gene_id2, gene_name),
    by = join_by(ncRNA == gene_id)
  ) |>
  dplyr::select(gene_id = gene_id2, target_id, gene_name) |>
  left_join(
    geneMetadata |>
      dplyr::select(gene_id, gene_id2, gene_name),
    by = join_by(target_id == gene_id)
  ) |>
  dplyr::select(
    gene_id,
    target_id = gene_id2,
    gene_name = gene_name.x,
    target_name = gene_name.y
  )


calculate_all_ncp_correlations <- function(vsd) {
  # Filter for pairs where both genes are present in the expression data
  relevant_pairs <- ncRNApc_combinations |>
    filter(gene_id %in% rownames(vsd) & target_id %in% rownames(vsd))

  # Get unique ncRNA and pcGene lists to create expression sub-matrices
  nc_genes <- unique(relevant_pairs$gene_id)
  pc_genes <- unique(relevant_pairs$target_id)

  vsd_assay <- assay(vsd)
  nc_expr <- vsd_assay[nc_genes, , drop = FALSE]
  pc_expr <- vsd_assay[pc_genes, , drop = FALSE]

  # Calculate all-vs-all Spearman correlations in a single matrix operation
  corr_matrix <- cor(t(nc_expr), t(pc_expr), method = "spearman")

  # Vectorized p-value calculation for the entire matrix
  n_samples <- ncol(vsd_assay)

  t_stat_matrix <- corr_matrix * sqrt((n_samples - 2) / (1 - corr_matrix^2))
  pval_matrix <- 2 * pt(-abs(t_stat_matrix), df = n_samples - 2)

  # Create a matrix of indices to look up results for our specific pairs
  pair_indices <- cbind(
    match(relevant_pairs$gene_id, rownames(corr_matrix)),
    match(relevant_pairs$target_id, colnames(corr_matrix))
  )

  # Efficiently extract the correlation and p-value for each pair
  relevant_pairs$corrEstimate <- corr_matrix[pair_indices]
  relevant_pairs$corrPval <- pval_matrix[pair_indices]

  return(relevant_pairs)
}


ncRNApcGene_corr <- mclapply(vsdList, calculate_all_ncp_correlations, mc.cores = Cores)

names(ncRNApcGene_corr) <- names(deseqList)


ncRNApcGene_corrdf <- bind_rows(ncRNApcGene_corr, .id = "Sex_Tissue") |>
  separate_wider_delim(Sex_Tissue, delim = " ", names = c("Sex", "Tissue")) |>
  dplyr::rename(ncPcCorrPval = corrPval, ncPcEstimate = corrEstimate)


# Merge Age corr and nc-pc corr
ncRNApcGene_Agecorrdf <- ncRNApcGene_corrdf |>
  left_join(
    corListdf |>
      dplyr::select(
        Tissue,
        Sex,
        gene_id,
        gene_type,
        AgePval = p_value,
        AgeEstimate = estimate
      ),
    by = join_by(Sex, Tissue, gene_id)
  )

write_tsv(ncRNApcGene_Agecorrdf, "Out/TMS/ncRNA_pcGene_Age_corr_TMS.tsv.gz")


# Remove environment
rm(list = setdiff(ls(), lsf.str()))
