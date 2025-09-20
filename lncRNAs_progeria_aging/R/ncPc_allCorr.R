Cores <- 6

# Import results ----

LAKIres <- readRDS("Out/LAKI/DESeq_obj_LAKI.rds")
TMSres <- readRDS("Out/TMS/DESeq_obj_TMS_oldvsyoung.rds")
Zhoures <- readRDS("Out/Zhou/DESeq_obj_Zhou_oldvsyoung.rds")
geneMetadata <- read_tsv("Data/geneMetadata.tsv")

## Generate vsd ----

LAKIcounts <- lapply(LAKIres, counts, normalized = TRUE)
LAKIcounts <- lapply(LAKIcounts, as_tibble, rownames = "gene_id")
LAKIcounts <- Reduce(
  function(df1, df2) full_join(df1, df2, by = "gene_id"),
  LAKIcounts
)

TMScounts <- lapply(TMSres, counts, normalized = TRUE)
TMScounts <- lapply(TMScounts, as_tibble, rownames = "gene_id")
TMScounts <- Reduce(
  function(df1, df2) full_join(df1, df2, by = "gene_id"),
  TMScounts
)

Zhoucounts <- lapply(Zhoures, counts, normalized = TRUE)
Zhoucounts <- lapply(Zhoucounts, as_tibble, rownames = "gene_id")
Zhoucounts <- Reduce(
  function(df1, df2) full_join(df1, df2, by = "gene_id"),
  Zhoucounts
)

allCounts <- full_join(TMScounts, Zhoucounts, by = join_by(gene_id)) |>
  left_join(
    geneMetadata |>
      dplyr::select(gene_id, gene_id2),
    by = join_by(gene_id == gene_id2)
  ) |>
  dplyr::select(-gene_id) |>
  dplyr::rename(gene_id = gene_id.y) |>
  full_join(LAKIcounts)

## ncRNA - cis pcGene correlation ----

ncRNApc_combinations <- read_tsv("Out/LAKI/ncRNA_Cis_pcGenes.tsv.gz") |>
  dplyr::select(
    gene_id = ncRNA,
    target_id,
    gene_name = ncRNA_name,
    target_name = target_name
  )


# --- Robust calcNcPcCorr (Data Frame version) ---
calcNcPcCorr <- function(ncRNA_id, pcGene_id, counts_df) {
  gene_id_col_index <- which(names(counts_df) == "gene_id")
  if (length(gene_id_col_index) == 0) {
    stop("Column 'gene_id' not found in counts_df")
  }

  # Find rows using logical indexing
  ncRNA_rows <- counts_df$gene_id == ncRNA_id
  pcGene_rows <- counts_df$gene_id == pcGene_id

  # Check if genes were found (sum logical vector)
  if (sum(ncRNA_rows) == 0 || sum(pcGene_rows) == 0) {
    # Return NA if either gene not found
    return(c(p.value = NA_real_, estimate = NA_real_))
  }
  # Use first match if duplicates exist
  idx_ncRNA <- which(ncRNA_rows)[1]
  idx_pcGene <- which(pcGene_rows)[1]

  # Extract numeric sample values and convert
  ncRNA_counts <- try(as.numeric(counts_df[idx_ncRNA, -gene_id_col_index]), silent = TRUE)
  pcGene_counts <- try(as.numeric(counts_df[idx_pcGene, -gene_id_col_index]), silent = TRUE)

  # Check if extraction worked
  if (inherits(ncRNA_counts, "try-error") || !is.numeric(ncRNA_counts) || inherits(pcGene_counts, "try-error") || !is.numeric(pcGene_counts)) {
    return(c(p.value = NA_real_, estimate = NA_real_))
  }

  # Check for sufficient variance BEFORE calling cor.test
  sd_ncRNA <- sd(ncRNA_counts, na.rm = TRUE)
  sd_pcGene <- sd(pcGene_counts, na.rm = TRUE)
  if (is.na(sd_ncRNA) || sd_ncRNA == 0 || is.na(sd_pcGene) || sd_pcGene == 0) {
    return(c(p.value = NA_real_, estimate = NA_real_))
  }

  # Perform correlation test within a try() block
  result <- try(cor.test(
    ncRNA_counts,
    pcGene_counts,
    method = "spearman",
    use = "complete.obs", # Handles missing values
    exact = FALSE
  ), silent = TRUE)

  # Check if cor.test itself failed (e.g., "not enough finite observations")
  if (inherits(result, "try-error")) {
    return(c(p.value = NA_real_, estimate = NA_real_))
  }

  # Return p-value and estimate (rho)
  unlist(result[c("p.value", "estimate")])
}


# --- Robust applyNcPcCorr ---
applyNcPcCorr <- function(counts_df) {
  # --- Input Validation ---
  if (!is.data.frame(counts_df) && !tibble::is_tibble(counts_df)) {
    stop("'counts_df' must be a data frame or tibble.")
  }
  if (!"gene_id" %in% names(counts_df)) {
    stop("'counts_df' must have a column named 'gene_id'.")
  }
  if (ncol(counts_df) <= 1) {
    stop("'counts_df' must have sample columns.")
  }
  counts_df$gene_id <- as.character(counts_df$gene_id) # Ensure character IDs

  # --- Filtering Combinations ---
  valid_gene_ids <- unique(counts_df$gene_id)
  df <- ncRNApc_combinations |>
    # Ensure combinations also use character IDs for matching
    mutate(gene_id = as.character(gene_id), target_id = as.character(target_id)) |>
    filter(gene_id %in% valid_gene_ids & target_id %in% valid_gene_ids)

  # Check if any pairs remain
  if (nrow(df) == 0) {
    warning("No valid ncRNA-PCG pairs found with both genes present in the counts data frame.")
    df$corrPval <- numeric(0)
    df$corrEstimate <- numeric(0)
    return(df)
  } else {
    message(paste("Found", nrow(df), "valid ncRNA-PCG pairs to test."))
  }

  # --- Parallel Correlation Calculation ---
  num_cores_mapply <- max(1, Cores)
  message(paste("Calculating correlations using mcmapply with", num_cores_mapply, "cores..."))

  # mcmapply will call the ROBUST calcNcPcCorr
  corResults <- try(mcmapply(
    FUN = calcNcPcCorr,
    ncRNA_id = df$gene_id,
    pcGene_id = df$target_id,
    MoreArgs = list(counts_df = counts_df), # Pass df explicitly
    mc.cores = num_cores_mapply,
    SIMPLIFY = TRUE # Request matrix output, will simplify if possible
  ), silent = TRUE)

  # --- Process Results ---
  if (inherits(corResults, "try-error")) {
    stop("mcmapply failed. Error: ", conditionMessage(attr(corResults, "condition")))
  }

  # Check if mcmapply returned the expected matrix format
  # The robust calcNcPcCorr ensures returning a 2-element vector or NA vector,
  # so mcmapply should reliably return a 2-row matrix here.
  if (!is.matrix(corResults) || nrow(corResults) != 2) {
    message("mcmapply did not return a 2-row matrix. Inspecting corResults:")
    # Use str() from utils package explicitly if needed
    print(utils::str(corResults))
    stop("Unexpected format returned by mcmapply. This shouldn't happen with the robust calcNcPcCorr.")
  }

  # Assign results
  df$corrPval <- corResults[1, ]
  df$corrEstimate <- corResults[2, ]

  message("Finished correlation calculations.")
  return(df)
}

# --- Execute ---
# Assume 'allCounts' is your prepared data frame
# Make sure prerequisite checks on allCounts structure are done if needed

ncRNApcGene_corr <- applyNcPcCorr(allCounts)

# Write results ----

write_tsv(ncRNApcGene_corr, file = "Out/allDb_ncPcCorr.tsv")
