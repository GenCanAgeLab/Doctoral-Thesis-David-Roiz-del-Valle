geneMetadata <- read_tsv("Data/geneMetadata.tsv")

# Load PhysioAging ----

TMS_results <- readRDS("Out/TMS/DESeq_results_TMS_oldvsyoung.rds") |>
  lapply(as_tibble, rownames = "gene_id") |>
  bind_rows(.id = "Sex_Tissue") |>
  separate_wider_delim(Sex_Tissue, names = c("Sex", "Tissue"), delim = " ")

Zhou_results <- readRDS("Out/Zhou/DESeq_results_Zhou_oldvsyoung.rds") |>
  lapply(as_tibble, rownames = "gene_id") |>
  bind_rows(.id = "Tissue") |>
  mutate(Sex = "Male")

Srivastava_results <- readRDS("Out/Srivastava/DESeq_results_Srivastava_oldvsyoung.rds") |>
  as_tibble(rownames = "gene_id") |>
  mutate(Sex = "Male", Tissue = "Liver")

# Change Tissue names in Zhou to make it comparable to TMS
Zhou_results <- Zhou_results |>
  mutate(
    Tissue = case_when(
      Tissue == "Bone_marrow" ~ "Marrow",
      Tissue == "Brown_Adipose" ~ "Brown adipose",
      Tissue == "Inginal_WAT" ~ "Subcutaneus adipose",
      Tissue == "Epididymal_WAT" ~ "Gonadal adipose",
      .default = Tissue
    )
  )

TMS_results <- TMS_results |>
  mutate(
    Tissue = case_when(
      Tissue == "Bone_marrow" ~ "Marrow",
      Tissue == "Brown_adipose" ~ "Brown adipose",
      Tissue == "Subcutaneus_adipose" ~ "Subcutaneus adipose",
      Tissue == "Gonadal_adipose" ~ "Gonadal adipose",
      .default = Tissue
    )
  )

# Merge PhysioAging ----

PhysioAging <- full_join(
  TMS_results,
  Zhou_results,
  by = join_by(gene_id, Sex, Tissue),
  suffix = c("_TMS", "_Zhou")
) |>
  # dplyr::select only those ncRNA with data on both datasets
  drop_na(log2FoldChange_TMS, log2FoldChange_Zhou) |>
  # Merge gene metadata info

  left_join(
    geneMetadata |>
      dplyr::select(gene_id2, gene_name, gene_type),
    by = join_by(gene_id == gene_id2)
  )

# DE stat plot

# DE_get <- function(res) {
#   res |>
#     group_by(Tissue, Sex) |>
#     summarize(
#       Upreg. = sum(padj < 0.05 & log2FoldChange > 0, na.rm = TRUE),
#       Downreg. = sum(padj < 0.05 & log2FoldChange < 0, na.rm = TRUE)
#     ) |>
#     pivot_longer(cols = c(Upreg., Downreg.), names_to = "Change", values_to = "Num") |>
#     mutate(Change = fct_relevel(Change, "Upreg."))
# }
#
# DE_stats <- lapply(list("TMS" = TMS_results,
#                         "Zhou" = Zhou_results),
#                    DE_get)
#
#
# DE_stats <- bind_rows(DE_stats, .id = "Dataset")
#
# common_tissues <- DE_stats |>
#   filter(Change == "Upreg.") |>
#   group_by(Tissue) |>
#   filter(n_distinct(Dataset) == 2) |>
#   distinct(Tissue) |>
#   pull(Tissue)



ggplot(
  PhysioAging |>
    filter(padj_TMS < 0.05 & padj_Zhou < 0.05) |>
    filter(sign(log2FoldChange_TMS) == sign(log2FoldChange_Zhou)) |>
    group_by(Tissue, Sex) |>
    summarize(
      Upreg. = sum(log2FoldChange_TMS > 0, na.rm = TRUE),
      Downreg. = sum(log2FoldChange_TMS < 0, na.rm = TRUE)
    ) |>
    pivot_longer(cols = c(Upreg., Downreg.), names_to = "Change", values_to = "Num") |>
    mutate(Change = fct_relevel(Change, "Upreg.")),
  aes(x = Tissue, y = Num, fill = Change)
) +
  geom_col(position = position_dodge(width = 0.8, preserve = "single"), width = 0.8) +
  theme_classic() +
  labs(y = "Differentially expressed genes") +
  scale_fill_manual(values = c("#CC79A7", "#0072B2")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Out/PhysioAging/DE_tissue.pdf", width = 5, height = 3)


# Plot

ggplot(
  PhysioAging |>
    filter(gene_type == "protein_coding") |>
    filter(padj_TMS < 0.05 & padj_Zhou < 0.05),
  aes(x = log2FoldChange_TMS, y = log2FoldChange_Zhou)
) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(vars(Tissue), scales = "free") +
  theme_bw()

# Calculate correlation

# TMS_Zhou_corr <- PhysioAging |>
#   filter(gene_type == "protein_coding") |>
#   filter(padj_TMS < 0.05 | padj_Zhou < 0.05) |>
#   group_by(Tissue) |>
#   summarise(
#     pval = cor.test(x = log2FoldChange_TMS,
#                     y = log2FoldChange_Zhou,
#                     method = "spearman")[["p.value"]],
#     rho = cor.test(x = log2FoldChange_TMS,
#                     y = log2FoldChange_Zhou,
#                     method = "spearman")[["estimate"]]
#   )

TMS_Zhou_corr_pc <- PhysioAging |>
  filter(gene_type == "protein_coding") |>
  filter(padj_TMS < 0.05 | padj_Zhou < 0.05) |>
  dplyr::select(
    Tissue,
    gene_id,
    log2FoldChange_TMS,
    log2FoldChange_Zhou
  ) |>
  pivot_longer(
    cols = starts_with("log2"),
    names_prefix = "log2FoldChange_",
    names_to = "Experiment",
    values_to = "log2FoldChange"
  ) |>
  pivot_wider(
    names_from = c(Tissue, Experiment),
    values_from = log2FoldChange
  ) |>
  dplyr::select(-gene_id)

TMS_Zhou_corr_nc <- PhysioAging |>
  filter(gene_type != "protein_coding") |>
  filter(!str_detect(gene_type, "IG")) |> # Also delete IG genes/pseudogenes
  filter(padj_TMS < 0.05 | padj_Zhou < 0.05) |>
  dplyr::select(
    Tissue,
    gene_id,
    log2FoldChange_TMS,
    log2FoldChange_Zhou
  ) |>
  pivot_longer(
    cols = starts_with("log2"),
    names_prefix = "log2FoldChange_",
    names_to = "Experiment",
    values_to = "log2FoldChange"
  ) |>
  pivot_wider(
    names_from = c(Tissue, Experiment),
    values_from = log2FoldChange
  ) |>
  dplyr::select(-gene_id)


# Compute the pairwise Spearman correlation between samples
# Using "pairwise.complete.obs" ensures that missing values (NA) are handled appropriately
cor_TMS_Zhou_pc <- cor(
  TMS_Zhou_corr_pc,
  use = "pairwise.complete.obs",
  method = "spearman"
)

cor_TMS_Zhou_nc <- cor(
  TMS_Zhou_corr_nc,
  use = "pairwise.complete.obs",
  method = "spearman"
)

ggcorrplot(cor_TMS_Zhou_pc, hc.order = TRUE)
ggsave("Out/PhysioAging/Heatmap_PhysioAging_pc.pdf", width = 8, height = 8)

ggcorrplot(cor_TMS_Zhou_nc, hc.order = TRUE)
ggsave("Out/PhysioAging/Heatmap_PhysioAging_nc.pdf", width = 8, height = 8)

# HSign ncRNA comparative ----

TopHsNcRNA_TMS <- read_tsv("Out/TMS/hsNcRNA_DE_Cor_TMS.tsv")
TopHsNcRNA_Zhou <- read_tsv("Out/Zhou/Cons_changed_ncRNA_Zhou.tsv")


DE_ncRNA_TMS <- TopHsNcRNA_TMS |>
  filter(Sex == "Male") |>
  filter(nDETis >= 3) |>
  dplyr::distinct(Sex, Change, gene_name) |>
  nest(data = gene_name) |>
  mutate(Change = if_else(Change == "UP", "TMS UP", "TMS DOWN")) |>
  pull(data, name = Change)

Cor_ncRNA_TMS <- TopHsNcRNA_TMS |>
  filter(Sex == "Male") |>
  filter(corrTissueN >= 5) |>
  dplyr::distinct(Sex, CorrDir, gene_name) |>
  nest(data = gene_name) |>
  mutate(CorrDir = if_else(CorrDir == 1, "TMS DirCorr", "TMS InvCorr")) |>
  pull(data, name = CorrDir)


DE_ncRNA_Zhou <- TopHsNcRNA_Zhou |>
  filter(nDETis >= 5) |>
  dplyr::distinct(Change, gene_name) |>
  nest(data = gene_name) |>
  mutate(Change = if_else(Change == "UP", "ZHOU UP", "ZHOU DOWN")) |>
  pull(data, name = Change)

Cor_ncRNA_Zhou <- TopHsNcRNA_Zhou |>
  filter(corrTissueN >= 7) |>
  dplyr::distinct(CorrDir, gene_name) |>
  nest(data = gene_name) |>
  mutate(CorrDir = if_else(CorrDir == 1, "Zhou DirCorr", "Zhou InvCorr")) |>
  pull(data, name = CorrDir)

DE_PhysioAging <- c(DE_ncRNA_TMS, DE_ncRNA_Zhou)
Corr_PhysioAging <- c(Cor_ncRNA_TMS, Cor_ncRNA_Zhou)

DE_PhysioAging <- lapply(DE_PhysioAging, function(x) {
  unlist(as.vector(x))
})
Corr_PhysioAging <- lapply(Corr_PhysioAging, function(x) {
  unlist(as.vector(x))
})


# Function to pad lists
pad_list <- function(lst) {
  max_len <- max(sapply(lst, length))
  lapply(lst, function(x) {
    length(x) <- max_len
    x
  })
}

DE_PhysioAging <- pad_list(DE_PhysioAging) |>
  as_tibble() # This will create columns with NA for shorter lists

Corr_PhysioAging <- pad_list(Corr_PhysioAging) |>
  as_tibble() # This will create columns with NA for shorter lists


write_tsv(DE_PhysioAging, file = "Out/nVenn_DE_TMS_Zhou.tsv")
write_tsv(Corr_PhysioAging, file = "Out/nVenn_Corr_TMS_Zhou.tsv")

UP_UP <- intersect(DE_PhysioAging$`TMS UP`, DE_PhysioAging$`ZHOU UP`)

DN_DN <- intersect(DE_PhysioAging$`TMS DOWN`, DE_PhysioAging$`ZHOU DOWN`)

DN_UP <- intersect(DE_PhysioAging$`TMS DOWN`, DE_PhysioAging$`ZHOU UP`)

DIR_DIR <- intersect(
  Corr_PhysioAging$`TMS DirCorr`,
  Corr_PhysioAging$`Zhou DirCorr`
)

INV_DIR <- intersect(
  Corr_PhysioAging$`TMS InvCorr`,
  Corr_PhysioAging$`Zhou DirCorr`
)
