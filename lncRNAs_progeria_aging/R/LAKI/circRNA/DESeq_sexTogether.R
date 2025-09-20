library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(DESeq2)

library(parallel)

# Import Ciri2 results ----

ciri2Res <- read_tsv("Data/circRNA_LAKI/Ciri2_circRNA_counts.tsv") |>
  column_to_rownames("ID")

circExp2Res <- read_tsv("Data/circRNA_LAKI/circExplorer2_circRNA_counts.tsv") |>
  column_to_rownames("ID")

findCircRes <- read_tsv("Data/circRNA_LAKI/findcirc_circRNA_counts.tsv") |>
  column_to_rownames("ID")

# terraceRes <- read_tsv("Data/circRNA_LAKI/terrace_circRNA_counts.tsv") |>
#   column_to_rownames("ID")

starchipRes <- read_tsv("Data/circRNA_LAKI/starchip_circRNA_counts.tsv") |>
  column_to_rownames("ID")


circExp2Res <- circExp2Res |>
  rename_with(
    # The function to apply: extract "Tube_NUMBER" part
    ~ sub("^(Tube_\\d+)_.*", "\\1", .),
    # Select columns to apply the function to: those starting with "Tube_" followed by a digit and an underscore
    .cols = matches("^Tube_\\d+_")
  )

# terraceRes <- terraceRes |>
#   rename_with(
#     # The function to apply: extract "Tube_NUMBER" part
#     ~ sub("^(Tube_\\d+)_.*", "\\1", .),
#     # Select columns to apply the function to: those starting with "Tube_" followed by a digit and an underscore
#     .cols = matches("^Tube_\\d+_")
#   )

circRes <- list(ciri2Res, circExp2Res, findCircRes, starchipRes)
names(circRes) <- c("CIRI2", "circExplorer2", "find_circ", "starchip")


circRes <- lapply(circRes, function(x) {
  as.matrix(x)
  # as.matrix(x) + 1 # Add a pseudo-count value of 1
})



# Load sample info ----

sampdf <- read.table(
  "Data/Info/Sample_info_complete_genotiped.tsv",
  header = T,
  sep = "\t"
) |>
  filter(Notes != "Ratones de Face-NR") |>
  # Fix wrong genotypes and remove HT animals (see Bash/Regenotyping)
  dplyr::rename(Regenotype = Infered_Genotype) |>
  filter(Regenotype != "HT") |>
  mutate(Regenotype = factor(Regenotype, levels = c("WT", "KO")))


rownames(sampdf) <- sampdf$SampleID


sampdf <- sampdf |>
  dplyr::select(
    Gender,
    MouseNo,
    Genotype = Regenotype,
    Gender,
    Tissue,
    DAB
  )

# Rearrange counts

circRes <- lapply(circRes, function(x) {
  x[, colnames(x) != "Tube_101"] # Delete the HT sample
})


# This is done bellow, isn't it?
# circRes <- lapply(circRes, function(x){
#   intersect <- intersect(rownames(sampdf), colnames(x))
#   x[, intersect]
# })

# DESeqDataSet

circDds <- lapply(circRes, function(x) {
  intersect <- intersect(rownames(sampdf), colnames(x))
  DESeqDataSetFromMatrix(
    countData = x,
    colData = sampdf[intersect, ],
    design = ~ Genotype + Gender
  )
})

# Use library sizes from all gene Salmon quantification ----

gse_Allgenes <- readRDS("Out/LAKI/gse_allSamples_LAKI.rds")

dds_Allgenes <- DESeqDataSet(gse_Allgenes, design = ~Genotype)

dds_Allgenes <- estimateSizeFactors(dds_Allgenes)

# Size factors are roughly mean cols of norm factors
# https://support.bioconductor.org/p/117781/
libSize <- colMeans(normalizationFactors(dds_Allgenes))

# This is not used anywhere in the code
# totalLib <- colSums(assay(gse_Allgenes))


# DESeq analysis per tissue

analysisDESeq <- function(Tissue, dds) {
  tmpDds <- dds[, colData(dds)[, "Tissue"] == Tissue] # First, select the samples from our target tissue

  keep <- rowSums(counts(tmpDds) >= 5) >= 5 # At least 5 junction reads in 5 samples
  tmpDds <- tmpDds[keep, ]

  # DESeq Wald test
  # sizeFactors(tmpDds) <- libSize[colnames(tmpDds)]
  # tmpDds <- estimateDispersionsGeneEst(tmpDds)
  # dispersions(tmpDds) <- mcols(tmpDds)$dispGeneEst
  # nbinomWaldTest(tmpDds)
  DESeq(tmpDds)
}

tissueList <- unique(sampdf$Tissue)


circTisDds <- mclapply(circDds,
  function(dds) {
    lapply(tissueList, analysisDESeq, dds = dds)
  },
  mc.cores = 5
)

circTisDds <- unlist(circTisDds)
names(circTisDds) <- paste(rep(names(circRes), each = length(tissueList)),
  rep(tissueList, times = length(circDds)),
  sep = "."
)

# Results

resultsList <- lapply(
  circTisDds,
  results,
  contrast = c("Genotype", "KO", "WT")
)


# Merge res

resultsList <- lapply(resultsList, function(x) {
  as_tibble(x, rownames = "circRNA_ID")
})

resAll <- bind_rows(resultsList, .id = "Tool.Tissue") |>
  separate_wider_delim(Tool.Tissue, names = c("Tool", "Tissue"), delim = ".")


# Save analysis

saveRDS(circTisDds, file = "Out/LAKI/circRNA/dds_circRNA_sexTogether.rds")
saveRDS(resultsList, file = "Out/LAKI/circRNA/results_circRNA_sexTogether.rds")

write_tsv(resAll, file = "Out/LAKI/circRNA/circRNA_results_all_sexTogether.tsv")

# Remove environment objects but no functions
rm(list = setdiff(ls(), lsf.str()))
