TMS_ncRNA <- read_tsv("Out/TMS/ncRNA_TMS.tsv.gz")
TMS_pcGene <- read_tsv("Out/TMS/pcGenes_TMS.tsv.gz")

ncPc_TMS <- read_tsv("Out/TMS/ncRNA_pcGene_Age_corr_TMS.tsv.gz")

TMS_ageCorr <- read_tsv("Out/TMS/Correlation_results_TMS_all.tsv.gz")
hsTMS_ageCorr <- read_tsv("Out/TMS/Correlation_results_TMS.tsv")

resultsTMS <- readRDS("Out/TMS/DESeq_results_TMS_oldvsyoung.rds")
geneMetadata <- read_tsv("Data/geneMetadata.tsv")

# Tissue stats ----

# ncRNA
DEstats <- TMS_ncRNA |>
  mutate(
    Tissue = case_when(
      Tissue == "Mesenteric_adipose" ~ "Mes. adipose",
      Tissue == "Gonadal_adipose" ~ "Gon. adipose",
      Tissue == "Brown_adipose" ~ "Brown adipose",
      Tissue == "Small_intestine" ~ "Small intestine",
      Tissue == "Limb_muscle" ~ "Limb muscle",
      .default = Tissue
    )
  ) |>
  mutate(Tissue = str_replace(Tissue, " ", "\n ")) |>
  mutate(Tissue = str_replace(Tissue, "_", "\n ")) |>
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
  facet_wrap(vars(Sex), ncol = 1, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Out/TMS/DE_tissue.pdf", width = 8, height = 5)

# pcGenes
DEstats <- TMS_pcGene |>
  mutate(
    Tissue = case_when(
      Tissue == "Mesenteric_adipose" ~ "Mes. adipose",
      Tissue == "Gonadal_adipose" ~ "Gon. adipose",
      Tissue == "Brown_adipose" ~ "Brown adipose",
      Tissue == "Small_intestine" ~ "Small intestine",
      Tissue == "Limb_muscle" ~ "Limb muscle",
      .default = Tissue
    )
  ) |>
  mutate(Tissue = str_replace(Tissue, " ", "\n ")) |>
  mutate(Tissue = str_replace(Tissue, "_", "\n ")) |>
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
  facet_wrap(vars(Sex), ncol = 1, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Out/TMS/DE_tissue_pcGenes.pdf", width = 8, height = 5)


Corstats <- TMS_ageCorr |>
  group_by(Tissue, Sex) |>
  summarise(
    DIR_COR = sum(p_value < 0.05 & estimate > 0, na.rm = TRUE),
    INV_COR = sum(p_value < 0.05 & estimate < 0, na.rm = TRUE)
  )

resList <- lapply(resultsTMS, function(x) {
  x <- as_tibble(x, rownames = "gene_id") |>
    left_join(geneMetadata, by = join_by(gene_id == gene_id2))
})

resAll <- bind_rows(resList, .id = "Tissue")

DE_type_stats <- resAll |>
  mutate(
    Tissue = case_when(
      Tissue == "Mesenteric_adipose" ~ "Mes. adipose",
      Tissue == "Gonadal_adipose" ~ "Gon. adipose",
      Tissue == "Brown_adipose" ~ "Brown adipose",
      Tissue == "Small_intestine" ~ "Small intestine",
      Tissue == "Limb_muscle" ~ "Limb muscle",
      .default = Tissue
    )
  ) |>
  mutate(Tissue = str_replace(Tissue, " ", "\n ")) |>
  mutate(Tissue = str_replace(Tissue, "_", "\n ")) |>
  group_by(Tissue) |>
  mutate(totalDE = sum(padj < 0.05, na.rm = TRUE)) |>
  group_by(Tissue, gene_type, totalDE) |>
  summarise(
    DE = sum(padj < 0.05, na.rm = TRUE)
  ) |>
  ungroup() |>
  mutate(percDE = DE / totalDE * 100) |>
  filter(percDE > 5)

ggplot(DE_type_stats, aes(x = Tissue, y = log10(DE), fill = gene_type)) +
  geom_col() +
  scale_fill_brewer(palette = "Set2") +
  labs(y = "log10(number of DE genes)") +
  theme_classic()

ggsave("Out/TMS/DE_gene_type.pdf", width = 20, height = 5)

# MA plot ncRNAs ----

ggplot(
  TMS_ncRNA |>
    mutate(
      padj = if_else(padj < 1e-5, 1e-5, padj),
      log2FoldChange = case_when(
        log2FoldChange > 5 ~ 5,
        log2FoldChange < -5 ~ -5,
        .default = log2FoldChange
      )
    ) |>
    mutate(
      Tissue = case_when(
        Tissue == "Mesenteric_adipose" ~ "Mes. adipose",
        Tissue == "Gonadal_adipose" ~ "Gon. adipose",
        Tissue == "Brown_adipose" ~ "Brown adipose",
        Tissue == "Small_intestine" ~ "Small intestine",
        Tissue == "Limb_muscle" ~ "Limb muscle",
        .default = Tissue
      )
    ),
  aes(x = log2FoldChange, y = -log10(padj), color = padj < 0.05)
) +
  geom_point(shape = 16) +
  facet_wrap(vars(paste(Tissue, Sex))) +
  theme_classic(base_size = 15) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey", "red"))

ggsave("Out/TMS/volcano.pdf", width = 12, height = 10)

# Highly sign ncRNAs ----

hsNcRNA <- TMS_ncRNA |>
  filter(padj < 0.05) |>
  group_by(Sex, gene_id, gene_name, Change = sign(log2FoldChange)) |>
  summarize(
    nDETis = n(),
    TisDE = paste(Tissue, collapse = ", "),
    medianLFC = median(log2FoldChange)
  )

hsNcRNA_DE_Cor <- full_join(
  hsNcRNA,
  hsTMS_ageCorr,
  by = join_by(Sex, gene_id, gene_name),
  relationship = "many-to-many"
) |>
  mutate(Change = if_else(Change == 1, "UP", "DOWN"))

TophsNc <- hsNcRNA_DE_Cor |>
  filter(nDETis >= 3 | corrTissueN >= 5) |>
  # We keep only the Top relationships for the table
  arrange(desc(nDETis), desc(corrTissueN)) |>
  distinct(Sex, gene_id, gene_name, .keep_all = TRUE) |>
  drop_na()

write_tsv(TophsNc, "Out/TMS/Cons_changed_ncRNA_TMS.tsv")

# nc-pc interactions of hsNcRNAs ----

Top_ncPc <- ncPc_TMS |>
  full_join(
    TMS_ncRNA |>
      dplyr::select(Sex, Tissue, gene_id, gene_name, log2FoldChange),
    by = join_by(Sex, Tissue, gene_id, gene_name)
  ) |>
  full_join(
    TMS_pcGene |>
      dplyr::select(Sex, Tissue, gene_id, gene_name, log2FoldChange_target = log2FoldChange),
    by = join_by(
      Sex, Tissue, target_id == gene_id,
      target_name == gene_name
    )
  ) |>
  filter(gene_id %in% TophsNc$gene_id) |>
  filter(ncPcCorrPval < 0.05) |>
  group_by(
    Sex,
    gene_id,
    target_id,
    gene_name,
    target_name,
    corrDir = sign(ncPcEstimate)
  ) |>
  summarize(
    nCorPc = n(),
    TisCor = paste(Tissue, collapse = ", "),
    medianEstimate = median(ncPcEstimate),
    medianLFC_nc = median(log2FoldChange),
    medianLFC_pc = median(log2FoldChange_target)
  ) |>
  filter(nCorPc >= 7)

write_tsv(Top_ncPc, "Out/TMS/Top_ncPC_TMS.tsv")

# ncRNA stats ----

# Summarize and calculate percentage
ncStats <- hsNcRNA |>
  group_by(nDETis, Sex) |>
  summarize(n = n()) |>
  group_by(Sex) |>
  mutate(perc = n / sum(n) * 100)

# Plot with percentages
ggplot(ncStats, aes(x = "", y = n, fill = as.factor(nDETis))) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  # geom_text(aes(label = label),
  #           position = position_stack(vjust = 0.5)) +
  labs(fill = "Sign. tissues") +
  scale_fill_brewer(palette = "Set1") +
  theme_void() +
  facet_wrap(vars(Sex), scales = "free")

ggsave("Out/TMS/ncRNA_tissue_specif.pdf", width = 5, height = 3)

# Tissue comparative ----

SexSign <- function(sex) {
  # UP
  UPNcRNA <- TMS_ncRNA |>
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
  DNNcRNA <- TMS_ncRNA |>
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


pdf("Out/TMS/Upset.pdf")

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
  MaleNcRNA <- TMS_ncRNA |>
    filter(Tissue == tis) |>
    filter(Sex == "Male") |>
    filter(padj < 0.05) |>
    mutate(
      Change = case_when(
        log2FoldChange > 0 ~ "UP",
        log2FoldChange < 0 ~ "DOWN"
      )
    ) |>
    dplyr::select(Change, gene_id) |>
    group_by(Change) |>
    nest() |>
    pull(data, name = Change)

  MaleNcRNA <- lapply(MaleNcRNA, function(x) {
    unlist(as.vector(x))
  })

  # Female
  FemaleNcRNA <- TMS_ncRNA |>
    filter(Tissue == tis) |>
    filter(Sex == "Female") |>
    filter(padj < 0.05) |>
    mutate(
      Change = case_when(
        log2FoldChange > 0 ~ "UP",
        log2FoldChange < 0 ~ "DOWN"
      )
    ) |>
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

signTiss <- lapply(unique(TMS_ncRNA$Tissue), TissueSign)

names(signTiss) <- unique(TMS_ncRNA$Tissue)

outVenn <- paste0("Out/TMS/nVenn/", unique(TMS_ncRNA$Tissue), ".svg")


## First, delete empty elements inside each list
signTiss <- map(signTiss, ~ keep(.x, ~ length(.) > 0))

## Now, tissue_sex that are empty
signTiss <- compact(signTiss)


library(nVennR)

for (i in 1:length(signTiss)) {
  plotVenn(signTiss[[i]], outFile = outVenn[i])
}

# Correlation ----

TMS_corr <- TMS_ncRNA |>
  mutate(Group = paste(Sex, Tissue)) |>
  dplyr::select(gene_id, Group, log2FoldChange) |>
  pivot_wider(
    names_from = Group,
    values_from = log2FoldChange
  ) |>
  dplyr::select(-gene_id)


# Compute the pairwise Spearman correlation between samples
# Using "pairwise.complete.obs" ensures that missing values (NA) are handled appropriately
cor_TMS <- cor(TMS_corr, use = "pairwise.complete.obs", method = "spearman")

library(ggcorrplot)

ggcorrplot(cor_TMS, hc.order = TRUE)
ggsave("Out/TMS/Heatmap_TisSex.pdf", width = 9, height = 9)


# Remove environment
rm(list = setdiff(ls(), lsf.str()))
