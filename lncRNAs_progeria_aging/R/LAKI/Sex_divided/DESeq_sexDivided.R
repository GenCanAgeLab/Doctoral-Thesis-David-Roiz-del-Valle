# Load gene metadata ----

geneMetadata <- read_tsv("Data/geneMetadata.tsv") |>
  column_to_rownames("gene_id")

# Load and arrange sample info ----

## Some samples are from a different experiment, so they will be removed.

# Load Sample info
sampdf <- read.table(
  "Data/Info/Sample_info_complete_genotiped.tsv",
  header = T,
  sep = "\t"
) |> # Cannot use readr as tibbles do not support rownames
  as_tibble() |> # We convert it to tibble as tximeta has problems with rownames

  # Add a column with file name
  mutate(
    files = file.path("Data/Salmon_quants", SampleID, "quant.sf")
  ) |>
  # Remove some animals from another experiment
  filter(Notes != "Ratones de Face-NR") |>
  # Fix wrong genotypes and remove HT animals (see Bash/Regenotyping)
  dplyr::rename(Regenotype = Infered_Genotype) |>
  filter(Regenotype != "HT")

stats <- sampdf |> 
  group_by(Gender,Regenotype,Tissue) |> 
  summarize(n=n())

# Load data with tximeta ----

# Check file paths
# all(file.exists(sampdf$files))

# Create names column
sampdf <- sampdf |>
  dplyr::rename(names = SampleID) |>
  dplyr::select(
    names,
    files,
    MouseNo,
    Genotype = Regenotype,
    Gender,
    Tissue,
    DAB
  )

# Using tximeta
txi <- tximeta(sampdf)
gse <- summarizeToGene(txi)

# Add rowData

rowData(gse) <- geneMetadata[rownames(gse), ]

# First, dds of all samples to check sex ----

## Check sex info ----

#check that sex metainfo is correct plotting gene expression of Xist

ddsAll <- DESeqDataSet(gse, design = ~Genotype + Tissue + Gender)
ddsAll <- DESeq(ddsAll)

Xist_cts <- plotCounts(ddsAll, 'ENSMUSG00000086503.5', intgroup = 'Gender', returnData = T)

ggplot(Xist_cts, aes(x = Gender, y = log(count), label = rownames(Xist_cts))) +
  geom_label(position = position_jitter(width = 0.3)) +
  theme_classic()

#It looks that Tube_100 and Tube_20 are incorrectly classified as males,
#and Tube_91 and Tube_11 incorrectly as females. So we change them

sampdf <- sampdf |> 
  #Tube_100 and Tube_20 are incorrectly classified as males,
  #and Tube_91 and Tube_11 incorrectly as females. So we change them
  mutate(Gender = case_when(
    names %in% c("Tube_100","Tube_20") ~ "Female",
    names %in% c("Tube_91","Tube_11") ~ "Male",
    .default = Gender
  ))

# Using tximeta
txi <- tximeta(sampdf)
gse <- summarizeToGene(txi)

# Add rowData

rowData(gse) <- geneMetadata[rownames(gse), ]

saveRDS(gse, file = "Out/LAKI/gse_allSamples_LAKI.rds")

# Quality check ----

sampleReads <- assay(gse) |>
  as_tibble(rownames = "gene_id")

mitoReads <- sampleReads |>
  left_join(
    geneMetadata |>
      rownames_to_column("gene_id"),
    by = join_by(gene_id)
  ) |>
  filter(Chr == "chrM") |>
  summarise(across(contains("Tube_"), sum)) |>
  pivot_longer(everything(), names_to = "Sample", values_to = "totalCounts")

totalReads <- sampleReads |>
  summarise(across(contains("Tube_"), sum)) |>
  pivot_longer(everything(), names_to = "Sample", values_to = "totalCounts")

# ggplot(totalReads, aes(Sample, totalCounts)) +
#   geom_col() +
#   labs(title = "Total reads")

# ggplot(mitoReads, aes(Sample, totalCounts)) +
#   geom_col() +
#   labs(title = "Mitochondrial reads")

ggplot(totalReads, aes(x = totalCounts / 1E6)) +
  geom_histogram() +
  theme_classic() +
  labs(
    title = "Density of mapped reads",
    x = "Mapped reads(in M reads)",
    y = "Number of samples"
  )
ggsave("Out/LAKI/Mapped_reads_histogram.pdf", width = 4, height = 3)

# Divide in male/female

gseSex <- list(
  "Male" = gse[
    ,
    colData(gse)$Gender == "Male"
  ],
  "Female" = gse[
    ,
    colData(gse)$Gender == "Female"
  ]
)


# DESeq2 ----

analysisDESeq <- function(gse, Tissue) {
  tmpGse <- gse[, colData(gse)[, "Tissue"] == Tissue] # First, select the samples from our target tissue

  tmpDds <- DESeqDataSet(tmpGse, design = ~Genotype)

  # smallestGroupSize <- as_tibble(colData(tmpDds)) |>
  #   group_by(Genotype, Gender) |>
  #   summarize(n = n()) |>
  #   pull(n) |>
  #   min()

  # Keep it simple: 116 / 2 genotypes / 6 tissues = almost 5
  smallestGroupSize <- 5

  keep <- rowSums(counts(tmpDds) >= 10) >= smallestGroupSize

  tmpDds <- tmpDds[keep, ]

  # DESeq Wald test
  DESeq(tmpDds)
}

tissueList <- unique(colData(gse)$Tissue)

ddsList <- mapply(
  analysisDESeq,
  # Repeat each element of gseSex 6 times (once for each tissue)
  rep(gseSex, each = length(tissueList)),
  # Repeat the tissueList for both sexes
  rep(tissueList, times = length(gseSex)),
  SIMPLIFY = FALSE
)

names(ddsList) <- paste(
  rep(names(gseSex), each = length(tissueList)),
  rep(tissueList, times = length(gseSex))
)

saveRDS(ddsList, "Out/LAKI/DESeq_obj_LAKI.rds")

# PCA plots ----

vsdList <- lapply(ddsList, vst, blind = TRUE)

pcaList <- lapply(vsdList, plotPCA, intgroup = "Genotype")


# Full results ----

resultsList <- lapply(
  ddsList,
  results,
  contrast = c("Genotype", "KO", "WT"),
  saveCols = c("gene_type", "gene_name")
)
names(resultsList) <- names(ddsList)

# When we use name = scaleAge, the log2FoldChange is the increase/decrease
# in expression per unit of scaleAge.

# Save results ----

saveRDS(resultsList, file = "Out/LAKI/DESeq_results_sexDiv_LAKI.rds")

# Remove environment objects but no functions
rm(list = setdiff(ls(), lsf.str()))
