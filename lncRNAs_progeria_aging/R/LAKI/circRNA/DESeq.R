library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(DESeq2)
library(ggplot2)

library(parallel)

# Import Ciri2 results ----

ciri2Res <- read_tsv("/home/jandrulas/20231017_RNASeq_LAKI/circRNA_analysis/Results/Ciri2_circRNA_counts.tsv") |>
  column_to_rownames("ID")

circExp2Res <- read_tsv("/home/jandrulas/20231017_RNASeq_LAKI/circRNA_analysis/Results/circExplorer2_circRNA_counts.tsv") |>
  column_to_rownames("ID")

findCircRes <- read_tsv("/home/jandrulas/20231017_RNASeq_LAKI/circRNA_analysis/Results/findcirc_circRNA_counts.tsv") |>
  column_to_rownames("ID")

# terraceRes <- read_tsv("/home/jandrulas/20231017_RNASeq_LAKI/circRNA_analysis/Results/terrace_circRNA_counts.tsv") |>
#   column_to_rownames("ID")

starchipRes <- read_tsv("/home/jandrulas/20231017_RNASeq_LAKI/circRNA_analysis/Results/starchip_circRNA_counts.tsv") |>
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
  mutate(Regenotype = factor(Regenotype, levels = c("WT", "KO"))) |> 
  #Tube_100 and Tube_20 are incorrectly classified as males,
  #and Tube_91 and Tube_11 incorrectly as females. So we change them
  mutate(Gender = case_when(
    SampleID %in% c("Tube_100","Tube_20") ~ "Female",
    SampleID %in% c("Tube_91","Tube_11") ~ "Male",
    .default = Gender
  ))


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

# circRes <- lapply(circRes, function(x){
#   intersect <- intersect(rownames(sampdf), colnames(x))
#   x[, intersect]
# })

# DESeqDataSet

circDds <- lapply(circRes, function(x) {
  intersect <- intersect(rownames(sampdf), colnames(x))
  DESeqDataSetFromMatrix(
    countData = x[, intersect],
    colData = sampdf[intersect, ],
    design = ~Genotype
  )
})

# Use library sizes from all gene Salmon quantification ----

gse_Allgenes <- readRDS("Out/LAKI/gse_allSamples_LAKI.rds")

dds_Allgenes <- DESeqDataSet(gse_Allgenes, design = ~Genotype)

dds_Allgenes <- estimateSizeFactors(dds_Allgenes)

# Size factors are roughly mean cols of norm factors
# https://support.bioconductor.org/p/117781/
libSize <- colMeans(normalizationFactors(dds_Allgenes))

totalLib <- colSums(assay(gse_Allgenes))


# Percentage of BSJ reads over total reads ----

BSJperc <- lapply(circRes, function(x) {
  as.data.frame(colSums(x) / totalLib * 100)
})

BSJ_all <- bind_cols(BSJperc)
colnames(BSJ_all) <- names(circRes)

BSJ_all <- merge(BSJ_all, sampdf, by = "row.names") |>
  mutate(Genotype = factor(Genotype, levels = c("WT", "KO"))) |>
  pivot_longer(
    cols = c("CIRI2", "circExplorer2", "find_circ", "starchip"),
    names_to = "Tool",
    values_to = "PercBSJ"
  )

ggplot(BSJ_all, aes(x = Tissue, y = PercBSJ, fill = Genotype)) +
  stat_summary(fun = "median", geom = "bar", position = "dodge") +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.2),
    shape = 16,
    alpha = 0.5
  ) +
  scale_fill_manual(values = c("#0072B2", "#ff3b30")) +
  theme_classic() +
  labs(y = "BSJ (%)") +
  facet_grid(
    rows = vars(Tool),
    cols = vars(Gender)
  )

ggsave("Out/LAKI/circRNA/BSJ_percentage.pdf", width = 8, height = 6)

# Divide into sex

circSexDiv <- lapply(circDds, function(x) {
  list(
    "Male" = x[
      ,
      colData(x)$Gender == "Male"
    ],
    "Female" = x[
      ,
      colData(x)$Gender == "Female"
    ]
  )
})

# DESeq analysis per tissue

analysisDESeq <- function(dds, Tissue) {
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


circSexTisDds <- mcmapply(
  function(gseList, Tissue) {
    print(Tissue)
    list(
      "Male" = analysisDESeq(gseList[[1]], Tissue),
      "Female" = analysisDESeq(gseList[[2]], Tissue)
    )
  },
  # Repeat each element of gseSex 6 times (once for each tissue)
  rep(circSexDiv, each = length(tissueList)),
  # Repeat the tissueList for all circRNA methods
  rep(tissueList, times = length(circSexDiv)),
  SIMPLIFY = FALSE,
  mc.cores = 5
)

names(circSexTisDds) <- paste(
  rep(names(circRes), each = length(tissueList)),
  rep(tissueList, times = length(circSexDiv))
)


circSexTisDds <- unlist(circSexTisDds)

# Results

resultsList <- lapply(
  circSexTisDds,
  results,
  contrast = c("Genotype", "KO", "WT")
)


# Merge res

resultsList <- lapply(resultsList, function(x) {
  as_tibble(x, rownames = "circRNA_ID")
})

resAll <- bind_rows(resultsList, .id = "Tool.Sex") |>
  separate_wider_delim(Tool.Sex, names = c("Tool_Tissue", "Sex"), delim = ".") |>
  separate_wider_delim(Tool_Tissue, names = c("Tool", "Tissue"), delim = " ")


# Save analysis

saveRDS(circSexTisDds, file = "Out/LAKI/circRNA/dds_circRNA.rds")
saveRDS(resultsList, file = "Out/LAKI/circRNA/results_circRNA.rds")

write_tsv(resAll, file = "Out/LAKI/circRNA/circRNA_results_all.tsv")

# Remove environment objects but no functions
rm(list = setdiff(ls(), lsf.str()))
