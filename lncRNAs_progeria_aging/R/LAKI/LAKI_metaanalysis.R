LAKI_ncRNA <- read_tsv("Out/LAKI/ncRNA_LAKI.tsv")
LAKI_pcRNA <- read_tsv("Out/LAKI/pcGenes_LAKI.tsv")
hSignLaki <- read_tsv(
  "Out/LAKI/HighlySign_NcRNA_cisProtGenes_clusters_sexDiv_LAKI.tsv"
)

resLAKI <- readRDS("Out/LAKI/DESeq_results_sexDiv_LAKI.rds")
geneMetadata <- read_tsv("Data/geneMetadata.tsv")

# Tissue stats ----

# ncRNA
DEstats <- LAKI_ncRNA |>
  group_by(Tissue, Sex) |>
  summarize(
    Upreg. = sum(padj < 0.05 & log2FoldChange > 0, na.rm = TRUE),
    Downreg. = sum(padj < 0.05 & log2FoldChange < 0, na.rm = TRUE)
  ) |>
  pivot_longer(cols = c(Upreg., Downreg.), names_to = "Change", values_to = "Num") |>
  mutate(Change = fct_relevel(Change, "Upreg."))

ggplot(DEstats, aes(x = Tissue, y = Num, fill = Change)) +
  geom_col(position = position_dodge(width = 0.8, preserve = "single"), width = 0.8) +
  theme_classic() +
  labs(y = "Differentially expressed genes") +
  scale_fill_manual(values = c("#ff3b30", "#4cd964")) +
  facet_wrap(vars(Sex))

ggsave("Out/LAKI/DE_tissue.pdf", width = 8, height = 3)

# pcGene
DEstats <- LAKI_pcRNA |>
  group_by(Tissue, Sex) |>
  summarize(
    Upreg. = sum(padj < 0.05 & log2FoldChange > 0, na.rm = TRUE),
    Downreg. = sum(padj < 0.05 & log2FoldChange < 0, na.rm = TRUE)
  ) |>
  pivot_longer(cols = c(Upreg., Downreg.), names_to = "Change", values_to = "Num") |>
  mutate(Change = fct_relevel(Change, "Upreg."))

ggplot(DEstats, aes(x = Tissue, y = Num, fill = Change)) +
  geom_col(position = position_dodge(width = 0.8, preserve = "single"), width = 0.8) +
  theme_classic() +
  labs(y = "Differentially expressed genes") +
  scale_fill_manual(values = c("#ff3b30", "#4cd964")) +
  facet_wrap(vars(Sex))

ggsave("Out/LAKI/DE_tissue_pcGenes.pdf", width = 8, height = 3)

resList <- lapply(resLAKI, function(x) {
  x <- as_tibble(x, rownames = "gene_id")
})

resAll <- bind_rows(resList, .id = "Tissue")

DE_type_stats <- resAll |>
  mutate(Tissue = str_replace(Tissue, " ", "\n ")) |>
  group_by(Tissue) |>
  mutate(totalDE = sum(padj < 0.05, na.rm = TRUE)) |>
  group_by(Tissue, gene_type, totalDE) |>
  summarise(
    DE = sum(padj < 0.05, na.rm = TRUE)
  ) |>
  ungroup() |>
  mutate(percDE = DE / totalDE * 100) |>
  filter(percDE > 2)

ggplot(DE_type_stats, aes(x = Tissue, y = log10(DE), fill = gene_type)) +
  geom_col() +
  labs(y = "log10(number of DE genes)") +
  theme_classic()

ggsave("Out/LAKI/DE_gene_type.pdf", width = 10, height = 5)

# Volcano plot ncRNAs ----

ggplot(
  LAKI_ncRNA |>
    mutate(
      padj = if_else(padj < 1e-10, 1e-10, padj),
      log2FoldChange = case_when(
        log2FoldChange > 5 ~ 5,
        log2FoldChange < -5 ~ -5,
        .default = log2FoldChange
      )
    ),
  aes(x = log2FoldChange, y = -log10(padj), color = padj < 0.05)
) +
  geom_point(shape = 16) +
  facet_wrap(vars(paste(Tissue, Sex))) +
  theme_classic(base_size = 15) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey", "red"))

ggsave("Out/LAKI/volcano.pdf", width = 10, height = 10)

# Highly sign ncRNAs ----

hsNcRNA <- LAKI_ncRNA |>
  filter(padj < 0.05) |>
  group_by(Sex, gene_id, gene_name, gene_type, sign(log2FoldChange)) |>
  summarize(
    nTissues = n(),
    signTissues = paste(Tissue, collapse = ", "),
    medianLFC = median(log2FoldChange)
  )

write_tsv(
  hsNcRNA |>
    filter(Sex == "Male"),
  file = "Out/LAKI/hsNcRNA_male.tsv"
)

write_tsv(
  hsNcRNA |>
    filter(Sex == "Female"),
  file = "Out/LAKI/hsNcRNA_female.tsv"
)

hSignLaki_pc <- hSignLaki |> 
  filter(corrTissuesN >= 4) |> 
  mutate(changeCorr = if_else(changeCorr == 1, "DIR.","INV.")) |> 
  mutate(corrTissues = str_replace_all(corrTissues,c("Ileon" = "IL",
                                         "Muscle" = "MU",
                                         "Kidney" = "KI",
                                         "Liver" = "LI",
                                         "Colon" = "CO",
                                         "Heart" = "HE"))) |> 
  arrange(desc(Sex), desc(corrTissuesN),changeCorr,desc(medianLFC)) |> 
  dplyr::select(Sex, gene_name, target_name, changeCorr, corrTissues, medianLFC, medianLFC_target, cluster)

write_tsv(hSignLaki_pc,
          "Out/LAKI/LAKI_hsNcRNA_pcGenes.tsv")

# ncRNA stats ----

# Summarize and calculate percentage
ncStats <- hsNcRNA |>
  group_by(nTissues, Sex) |>
  summarize(n = n()) |>
  group_by(Sex) |>
  mutate(perc = n / sum(n) * 100)

# Plot with percentages
ggplot(ncStats, aes(x = "", y = n, fill = as.factor(nTissues))) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  # geom_text(aes(label = label),
  #           position = position_stack(vjust = 0.5)) +
  labs(fill = "Sign. tissues") +
  theme_void() +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(vars(Sex), scales = "free")

ggsave("Out/LAKI/ncRNA_tissue_specif.pdf", width = 5, height = 3)

# Tissue comparative ----

SexSign <- function(sex) {
  # UP
  UPNcRNA <- LAKI_ncRNA |>
    filter(Sex == sex) |>
    filter(padj < 0.05 & log2FoldChange > 0) |>
    dplyr::select(Tissue, gene_id) |>
    group_by(Tissue) |>
    nest() |>
    pull(data, name = Tissue)

  UPNcRNA <- lapply(UPNcRNA, function(x) {
    unlist(as.vector(x))
  })

  # DN
  DNNcRNA <- LAKI_ncRNA |>
    filter(Sex == sex) |>
    filter(padj < 0.05 & log2FoldChange < 0) |>
    dplyr::select(Tissue, gene_id) |>
    group_by(Tissue) |>
    nest() |>
    pull(data, name = Tissue)

  DNNcRNA <- lapply(DNNcRNA, function(x) {
    unlist(as.vector(x))
  })

  return(list(UPNcRNA, DNNcRNA))
}

signSex <- lapply(c("Male", "Female"), SexSign)


pdf("Out/LAKI/Upset.pdf", width = 10, height = 5)

# Male

print(upset(fromList(signSex[[1]][[1]]), nsets = 6))
grid::grid.text("Male UP", just = "top")


print(upset(fromList(signSex[[1]][[2]]), nsets = 6))
grid::grid.text("Male DN", just = "top")

# Female

print(upset(fromList(signSex[[2]][[1]]), nsets = 6))
grid::grid.text("Female UP", just = "top")


print(upset(fromList(signSex[[2]][[2]]), nsets = 6))
grid::grid.text("Female DN", just = "top")

dev.off()

# Sex comparative ----

TissueSign <- function(tis) {
  # Male
  MaleNcRNA <- LAKI_ncRNA |>
    filter(Tissue == tis) |>
    filter(Sex == "Male") |>
    filter(padj < 0.05) |>
    mutate(Change = factor(log2FoldChange > 0, labels = c("UP", "DOWN"))) |>
    dplyr::select(Change, gene_id) |>
    group_by(Change) |>
    nest() |>
    pull(data, name = Change)

  MaleNcRNA <- lapply(MaleNcRNA, function(x) {
    unlist(as.vector(x))
  })

  # Female
  FemaleNcRNA <- LAKI_ncRNA |>
    filter(Tissue == tis) |>
    filter(Sex == "Female") |>
    filter(padj < 0.05) |>
    mutate(Change = factor(log2FoldChange > 0, labels = c("UP", "DOWN"))) |>
    dplyr::select(Change, gene_id) |>
    group_by(Change) |>
    nest() |>
    pull(data, name = Change)

  FemaleNcRNA <- lapply(FemaleNcRNA, function(x) {
    unlist(as.vector(x))
  })

  return(list(
    "Male UP" = MaleNcRNA$UP,
    "Male DOWN" = MaleNcRNA$DOWN,
    "Female UP" = FemaleNcRNA$UP,
    "Female DOWN" = FemaleNcRNA$DOWN
  ))
}

signTiss <- lapply(unique(LAKI_ncRNA$Tissue), TissueSign)

names(signTiss) <- unique(LAKI_ncRNA$Tissue)

# Out file for nVenn2 ----

signTiss <- unlist(signTiss, recursive = FALSE)

# Function to pad lists
pad_list <- function(lst) {
  max_len <- max(sapply(lst, length))
  lapply(lst, function(x) {
    length(x) <- max_len
    x
  })
}

signTiss <- pad_list(signTiss)
signTiss <- as_tibble(signTiss) # This will create columns with NA for shorter lists

write_tsv(signTiss, file = "Out/LAKI/nVenn2_Male_Female_Tissue.tsv")


# Correlation ----

LAKI_corr <- LAKI_ncRNA |>
  mutate(Group = paste(Sex, Tissue)) |>
  dplyr::select(gene_id, Group, log2FoldChange) |>
  pivot_wider(
    names_from = Group,
    values_from = log2FoldChange
  ) |>
  dplyr::select(-gene_id)


# Compute the pairwise Spearman correlation between samples
# Using "pairwise.complete.obs" ensures that missing values (NA) are handled appropriately
cor_LAKI <- cor(LAKI_corr, use = "pairwise.complete.obs", method = "spearman")

library(ggcorrplot)

ggcorrplot(cor_LAKI, hc.order = TRUE)
ggsave("Out/LAKI/Heatmap_TisSex.pdf", width = 5, height = 5)


# NcRNA-pcGene plot ----

# Think about LFC values to choose.... for both target and ncGene
cisGenesSelect <- hSignLaki |>
  filter(gene_name == "Terc") |>
  arrange(desc(corrPower), desc(medianLFC_target)) |>
  distinct(target_name, .keep_all = TRUE) |>
  dplyr::select(gene_id = target_id, medianLFC = medianLFC_target) |>
  left_join(geneMetadata, by = join_by(gene_id)) |>
  arrange(start) |>
  mutate(yPos = rep(c(1, 1.1, 1.2), length.out = n()))

ncRNASelect <- hSignLaki |>
  filter(gene_name == "Terc") |>
  arrange(desc(corrPower), desc(medianLFC)) |>
  slice_head(n = 1) |>
  dplyr::select(gene_id, medianLFC) |>
  left_join(geneMetadata, by = join_by(gene_id)) |>
  arrange(start) |>
  mutate(yPos = 1.1)

cisGenesSelect <- bind_rows(cisGenesSelect, ncRNASelect)


cisGenePlot <- function(ncRNA) {
  # Select cis Genes and ncRNAs
  cisGenesSelect <- hSignLaki |>
    filter(gene_name == ncRNA) |>
    arrange(desc(corrPower), desc(medianLFC_target)) |>
    distinct(target_name, .keep_all = TRUE) |>
    dplyr::select(gene_id = target_id, medianLFC = medianLFC_target) |>
    left_join(geneMetadata, by = join_by(gene_id)) |>
    arrange(start) |>
    mutate(yPos = rep(c(1, 1.1, 1.2), length.out = n()))

  ncRNASelect <- hSignLaki |>
    filter(gene_name == ncRNA) |>
    arrange(desc(corrPower), desc(medianLFC)) |>
    slice_head(n = 1) |>
    dplyr::select(gene_id, medianLFC) |>
    left_join(geneMetadata, by = join_by(gene_id)) |>
    arrange(start) |>
    mutate(yPos = 1.1)


  # Bind cisGenes and ncRNA
  cisGenesSelect <- bind_rows(cisGenesSelect, ncRNASelect)

  # Plot

  ggplot(cisGenesSelect, aes(xmin = start, xmax = end, ymin = 1, ymax = yPos + 0.8, fill = medianLFC)) +
    theme_classic() +
    geom_label_repel(
      # Data is inherited from ggplot() call
      # Map aesthetics, including conditional ones for styling:
      aes(
        x = (start + end) / 2,
        y = yPos,
        label = gene_name,
        # Conditionally set size based on gene_name:
        size = ifelse(gene_name == ncRNA, 4, 2),
        # Conditionally set fontface based on gene_name:
        fontface = ifelse(gene_name == ncRNA, "bold", "plain") # "plain" is default normal face
      ),
      # Constant arguments applying to all labels in this layer:
      color = "black", # Text color
      label.padding = unit(0.15, "lines"),
      segment.colour = NA,
      label.r = unit(0.1, "lines"), # Corner rounding
      label.size = 0, # No border around the label box itself
      box.padding = 0.5, # Padding for repulsion calculation (adjust as needed)
      min.segment.length = 0.1, # Shorten minimum segment length to draw more lines
      max.overlaps = Inf # Allow overlaps if necessary, or set a limit
      # The 'fill' aesthetic is automatically inherited from the global aes()
    ) +
    scale_fill_gradient2(
      low = "green", mid = "white", high = "red",
      midpoint = 0, na.value = "grey50", # Color for NA LFC values
      name = "log2FoldChange"
    ) +
    scale_size_identity() +
    ylim(0.9, 1.25) +
    labs(y = NULL, x = "Genomic coordinates (bp)")
}


cisGenePlot("Terc")
ggsave("Out/LAKI/cisTert.pdf", width = 10, height = 3)

cisGenePlot("Reno1")
ggsave("Out/LAKI/cisReno1.pdf", width = 10, height = 3)

library(ggrepel)

ggplot(cisGenesSelect, aes(xmin = start, xmax = end, ymin = 1, ymax = yPos + 0.8, fill = medianLFC)) +
  theme_classic() +
  geom_label_repel(
    # Data is inherited from ggplot() call
    # Map aesthetics, including conditional ones for styling:
    aes(
      x = (start + end) / 2,
      y = yPos,
      label = gene_name,
      # Conditionally set size based on gene_name:
      size = ifelse(gene_name == "Terc", 4, 2),
      # Conditionally set fontface based on gene_name:
      fontface = ifelse(gene_name == "Terc", "bold", "plain") # "plain" is default normal face
    ),
    # Constant arguments applying to all labels in this layer:
    color = "black", # Text color
    label.padding = unit(0.15, "lines"),
    segment.colour = NA,
    label.r = unit(0.1, "lines"), # Corner rounding
    label.size = 0, # No border around the label box itself
    box.padding = 0.5, # Padding for repulsion calculation (adjust as needed)
    min.segment.length = 0.1, # Shorten minimum segment length to draw more lines
    max.overlaps = Inf # Allow overlaps if necessary, or set a limit
    # The 'fill' aesthetic is automatically inherited from the global aes()
  ) +
  scale_fill_gradient2(
    low = "green", mid = "white", high = "red",
    midpoint = 0, na.value = "grey50", # Color for NA LFC values
    name = "log2FoldChange"
  ) +
  scale_size_identity() +
  ylim(0.9, 1.25) +
  labs(y = NULL, x = "Genomic coordinates (bp)")

ggsave("Out/LAKI/cisTert.pdf", width = 10, height = 3)


# Cis-pcGenes heatmap plot

require(rtracklayer)
grs <- import("/data/genomes/GRCm39_GENCODE/M32/gencode.vM32.primary_assembly.annotation.gtf")
grs <- subset(grs, type == "gene")

query_g <- "ENSMUSG00000118406.2"
lakids <- readRDS("Out/LAKI/DESeq_results_sexDiv_LAKI.rds")
lakiexp <- cbind(
  do.call("rbind", lakids),
  "names" = rep(names(lakids), times = sapply(lakids, nrow)),
  "gene_id" = unlist(lapply(lakids, row.names))
)


prepare_cisexp <- function(query_id, genegtf, deseqres) {
  # Get nearby genes and order (by gene start and ignoring strand)
  query_gr <- subset(genegtf, gene_id == query_g)
  stopifnot(length(query_gr) == 1)
  target_g <- subsetByOverlaps(genegtf, query_gr, maxgap = 1E6)
  target_g <- sort(unstrand(resize(target_g, 1)))
  # target_g$distance <- distance(query_gr,target_g)/1000
  target_g$distance <- (start(target_g) - start(query_gr)) / 1000
  target_g <- as.data.frame(target_g)[, c("gene_id", "gene_name", "distance")]

  # Get expression
  target_exp <- deseqres[deseqres$gene_id %in% target_g$gene_id, ]
  # Fill empty info
  target_exp <- merge(
    expand_grid("names" = unique(deseqres$names), target_g),
    target_exp,
    by = c("gene_id", "names"),
    all = T
  )
  metadata(target_exp) <- as.list(target_g)
  target_exp$gene_id <- factor(target_exp$gene_id, levels = target_g$gene_id)
  return(target_exp)
}

plotdt <- prepare_cisexp("ENSMUSG00000118406.2", grs, lakiexp)
plotdt$gender <- gsub("\\s+.*", "", plotdt$names)
plotdt$tissue <- gsub(".*\\s+", "", plotdt$names)


# Plot
require(ggplot2)
ggplot(plotdt, aes(x = gene_id, y = tissue, fill = log2FoldChange)) +
  geom_tile(col = "black") +
  scale_fill_gradient2(low = "blue4", high = "darkred", mid = "white", na.value = "grey70") +
  facet_wrap(~gender, ncol = 1) +
  scale_x_discrete(
    breaks = metadata(plotdt)$gene_id,
    labels = sprintf(
      "%s (%.0f)",
      metadata(plotdt)$gene_name,
      metadata(plotdt)$distance
    )
  ) +
  coord_fixed() +
  theme_bw(8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
