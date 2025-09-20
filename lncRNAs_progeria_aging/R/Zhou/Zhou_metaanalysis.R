Zhou_ncRNA <- read_tsv("Out/Zhou/ncRNA_Zhou.tsv.gz")
Zhou_pcGenes <- read_tsv("Out/Zhou/pcGene_Zhou.tsv.gz")


hSignZhou <- read_tsv("Out/Zhou/ncRNA_pcGene_Age_corr_Zhou.tsv.gz")
Zhou_ageCorr <- read_tsv("Out/Zhou/Correlation_results_Zhou.tsv")

resultsZhou <- readRDS("Out/Zhou/DESeq_results_Zhou_oldvsyoung.rds")
geneMetadata <- read_tsv("Data/geneMetadata.tsv")

# Tissue stats ----

# ncRNA
DEstats <- Zhou_ncRNA |>
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
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Out/Zhou/DE_tissue.pdf", width = 6, height = 3)

# pcGenes
DEstats <- Zhou_pcGenes |>
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
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Out/Zhou/DE_tissue_pcGenes.pdf", width = 6, height = 3)

Corstats <- hSignZhou |>
  select(!contains("target")) |>
  select(!contains("ncPc")) |>
  drop_na() |>
  distinct(Tissue, gene_id, gene_name, .keep_all = TRUE) |>
  group_by(Tissue) |>
  summarize(
    DIR_COR = sum(AgePval < 0.05 & AgeEstimate > 0, na.rm = TRUE),
    INV_COR = sum(AgePval < 0.05 & AgeEstimate < 0, na.rm = TRUE)
  )

resList <- lapply(resultsZhou, function(x) {
  x <- as_tibble(x, rownames = "gene_id") |>
    left_join(geneMetadata, by = join_by(gene_id == gene_id2))
})

resAll <- bind_rows(resList, .id = "Tissue")

DE_type_stats <- resAll |>
  mutate(
    Tissue = case_when(
      Tissue == "Bone_marrow" ~ "Bone\n marrow",
      Tissue == "Brown_Adipose" ~ "Brown\n adipose",
      Tissue == "Gastric_muscle" ~ "Gastric\n muscle",
      Tissue == "Inginal_WAT" ~ "Inguinal\n WAT",
      Tissue == "Epididymal_WAT" ~ "Epididymal\n WAT",
      .default = Tissue
    )
  ) |>
  group_by(Tissue) |>
  mutate(totalDE = sum(padj < 0.05, na.rm = TRUE)) |>
  group_by(Tissue, gene_type, totalDE) |>
  summarise(
    DE = sum(padj < 0.05, na.rm = TRUE)
  ) |>
  ungroup() |>
  mutate(percDE = DE / totalDE * 100) |>
  filter(percDE > 1)

ggplot(DE_type_stats, aes(x = Tissue, y = log10(DE), fill = gene_type)) +
  geom_col() +
  scale_fill_brewer(palette = "Set2") +
  labs(y = "log10(number of DE genes)") +
  theme_classic()

ggsave("Out/Zhou/DE_gene_type.pdf", width = 10, height = 5)

# Volcano plot ncRNAs ----

ggplot(
  Zhou_ncRNA |>
    mutate(
      padj = if_else(padj < 1e-10, 1e-10, padj),
      log2FoldChange = case_when(
        log2FoldChange > 2 ~ 2,
        log2FoldChange < -2 ~ -2,
        .default = log2FoldChange
      )
    ) |>
    mutate(
      Tissue = case_when(
        Tissue == "Bone_marrow" ~ "Bone marrow",
        Tissue == "Brown_Adipose" ~ "Brown adipose",
        Tissue == "Gastric_muscle" ~ "Gastric muscle",
        Tissue == "Inginal_WAT" ~ "Inguinal WAT",
        Tissue == "Epididymal_WAT" ~ "Epididymal WAT",
        .default = Tissue
      )
    ),
  aes(x = log2FoldChange, y = -log10(padj), color = padj < 0.05)
) +
  geom_point(shape = 16) +
  facet_wrap(vars(Tissue)) +
  theme_classic(base_size = 15) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey", "red"))

ggsave("Out/Zhou/volcano.pdf", width = 10, height = 10)

# Consistently changed sign ncRNAs ----

hsNcRNA <- Zhou_ncRNA |>
  filter(padj < 0.05) |>
  group_by(gene_id, gene_name, Change = sign(log2FoldChange)) |>
  summarize(
    nDETis = n(),
    TisDE = paste(Tissue, collapse = ", "),
    medianLFC = median(log2FoldChange)
  )

hsNcRNA_DE_Cor <- full_join(
  hsNcRNA,
  Zhou_ageCorr,
  by = join_by(gene_id, gene_name),
  relationship = "many-to-many"
) |>
  mutate(Change = if_else(Change == 1, "UP", "DOWN"))

TophsNc <- hsNcRNA_DE_Cor |>
  filter(nDETis >= 5 | corrTissueN >= 7) |>
  # We keep only the Top relationships for the table
  arrange(desc(nDETis), desc(corrTissueN)) |>
  distinct(gene_id, gene_name, .keep_all = TRUE) |>
  drop_na()


write_tsv(TophsNc, "Out/Zhou/Cons_changed_ncRNA_Zhou.tsv")

# nc-pc interactions of consistently changed ncRNAs ----

Top_ncPc <- hSignZhou |>
  filter(gene_id %in% TophsNc$gene_id) |>
  filter(ncPcCorrPval < 0.05) |>
  group_by(
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
  filter(nCorPc >= 5)

write_tsv(Top_ncPc, "Out/Zhou/Cons_changed_ncRNA_pcGenes_Zhou.tsv")

# ncRNA stats ----

# Summarize and calculate percentage
ncStats <- hsNcRNA |>
  mutate(nDETis = if_else(nDETis > 5, 6, nDETis)) |>
  mutate(
    nDETis = as.character(nDETis),
    nDETis = if_else(nDETis == "6", "5+", nDETis)
  ) |>
  group_by(nDETis) |>
  summarize(n = n()) |>
  ungroup() |>
  mutate(perc = n / sum(n) * 100)

# Plot with percentages
ggplot(ncStats, aes(x = "", y = n, fill = as.factor(nDETis))) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  # geom_text(aes(label = label),
  #           position = position_stack(vjust = 0.5)) +
  scale_fill_brewer(palette = "Set1") +
  labs(fill = "Sign. tissues") +
  theme_void()

ggsave("Out/Zhou/ncRNA_tissue_specif.pdf", width = 5, height = 3)

# Tissue comparative ----

SexSign <- function() {
  # UP
  UPNcRNA <- Zhou_ncRNA |>
    filter(padj < 0.05 & log2FoldChange > 0) |>
    dplyr::select(Tissue, gene_id) |>
    group_by(Tissue) |>
    nest() |>
    pull(data, name = Tissue)

  UPNcRNA <- lapply(UPNcRNA, function(x) {
    unlist(as.vector(x))
  })

  # DN
  DNNcRNA <- Zhou_ncRNA |>
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

signSex <- SexSign()


pdf("Out/Zhou/Upset.pdf")

# Male

upset(fromList(signSex[[1]]), nsets = 16)
grid::grid.text("Male UP", just = "top")


upset(fromList(signSex[[2]]), nsets = 16)
grid::grid.text("Male DN", just = "top")


dev.off()


# Correlation ----

Zhou_corr <- Zhou_ncRNA |>
  dplyr::select(gene_id, Tissue, log2FoldChange) |>
  pivot_wider(
    names_from = Tissue,
    values_from = log2FoldChange
  ) |>
  dplyr::select(-gene_id)


# Compute the pairwise Spearman correlation between samples
# Using "pairwise.complete.obs" ensures that missing values (NA) are handled appropriately
cor_Zhou <- cor(Zhou_corr, use = "pairwise.complete.obs", method = "spearman")

library(ggcorrplot)

ggcorrplot(cor_Zhou, hc.order = TRUE)
ggsave("Out/Zhou/Heatmap_Tis.pdf", width = 5, height = 5)

# Remove environment
rm(list = setdiff(ls(), lsf.str()))
