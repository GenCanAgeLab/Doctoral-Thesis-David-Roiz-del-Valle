source("R/init.R")

# GSEA  ---------------------------

## Load data ----

gene_res <- readRDS("Data/results.rds")

## Load GMTs  ---------------------

mouseHallmarks <- gmtPathways(
  "/data/genomes/MsigDB/mh.all.v2024.1.Mm.symbols.gmt"
)
mouseM2 <- gmtPathways(
  "/data/genomes/MsigDB/m2.all.v2024.1.Mm.symbols.gmt"
)

mouseM3 <- gmtPathways(
  "/data/genomes/MsigDB/m3.all.v2024.1.Mm.symbols.gmt"
)

mouseM5 <- gmtPathways(
  "/data/genomes/MsigDB/m5.all.v2024.1.Mm.symbols.gmt"
)

mouseM8 <- gmtPathways(
  "/data/genomes/MsigDB/m8.all.v2024.1.Mm.symbols.gmt"
)

Progeria_datasets <- c(
  Ortega_Molina_progeria <- read_delim("/data/droiva/20240509_Guille_RNASeq/20240509_Guille_RNASeq_analysis/data/Ortega_Molina_2024/Signatures_OrtegaMolina_progeria.csv"),
  gmtPathways("/data/droiva/20240509_Guille_RNASeq/20240509_Guille_RNASeq_analysis/data/Varela_Zmpste24/VARELA_ZMPSTE24_TARGETS_UP.v2024.1.Hs.gmt"),
  gmtPathways("/data/droiva/20240509_Guille_RNASeq/20240509_Guille_RNASeq_analysis/data/Varela_Zmpste24/VARELA_ZMPSTE24_TARGETS_DN.v2024.1.Hs.gmt"),
  readRDS("/data/droiva/20240509_Guille_RNASeq/20240509_Guille_RNASeq_analysis/data/Osorio_LAKI_geneset.rds"),
  readRDS("/data/droiva/20240509_Guille_RNASeq/20240509_Guille_RNASeq_analysis/data/CS_Pluijm_geneset.rds"),
  readRDS("/data/droiva/20240509_Guille_RNASeq/20240509_Guille_RNASeq_analysis/data/ERCC_Schumacher_geneset.rds"),
  readRDS("/data/droiva/20240509_Guille_RNASeq/20240509_Guille_RNASeq_analysis/data/ERCC_Niedernhofer_geneset.rds"),
  readRDS("/data/droiva/20240509_Guille_RNASeq/20240509_Guille_RNASeq_analysis/data/SIRT6_Strand_geneset.rds"),
  readRDS("/data/droiva/20240509_Guille_RNASeq/20240509_Guille_RNASeq_analysis/data/TERC_Ferrara_geneset.rds"),
  readRDS("/data/droiva/20240509_Guille_RNASeq/20240509_Guille_RNASeq_analysis/data/POLG_Hauser_geneset.rds"),
  readRDS("/data/droiva/20240509_Guille_RNASeq/20240509_Guille_RNASeq_analysis/data/POLG_Hiona_geneset.rds"),
  readRDS("/data/droiva/20240509_Guille_RNASeq/20240509_Guille_RNASeq_analysis/data/SIRT6_Masri_geneset.rds"),
  readRDS("/data/droiva/20240509_Guille_RNASeq/20240509_Guille_RNASeq_analysis/data/SIRT1_Masri_geneset.rds")
)
## Rank  ---------------------

Rank <- gene_res |>
  drop_na(pvalue, log2FoldChange, gene_name) |> #  CollapsePathways fails if NA gene names are present
  arrange(desc(stat)) |> # Using stat as reference for ranking
  distinct(gene_name, .keep_all = T) |>
  pull(stat, name = gene_name)

## GSEA and plotting  ---------------------

set.seed(42) # Fix seed!

GSEA_H <- fgsea(mouseHallmarks, Rank, nproc = 1) |>
  filter(padj < 0.05) |>
  mutate(pathway = str_remove(pathway, "HALLMARK_")) |>
  rowwise() |>
  mutate(leadingEdge = paste(unlist(leadingEdge),
    collapse = ", "
  ))

GSEA_M2 <- fgsea(mouseM2, Rank, nproc = 1) |>
  filter(padj < 0.05) |>
  rowwise() |>
  mutate(leadingEdge = paste(unlist(leadingEdge),
    collapse = ", "
  )) |>
  ungroup()

GSEA_M3 <- fgsea(mouseM3, Rank, nproc = 1) |>
  filter(padj < 0.05) |>
  rowwise() |>
  mutate(leadingEdge = paste(unlist(leadingEdge),
    collapse = ", "
  )) |>
  ungroup()

GSEA_M5 <- fgsea(mouseM5, Rank, nproc = 1) |>
  filter(padj < 0.01) |>
  rowwise() |>
  mutate(leadingEdge = paste(unlist(leadingEdge),
    collapse = ", "
  )) |>
  ungroup()

collapsedM2 <- collapsePathways(
  as.data.table(GSEA_M2)[order(pval)],
  mouseM2, Rank
)

collapsedM5 <- collapsePathways(
  as.data.table(GSEA_M5)[order(pval)],
  mouseM5, Rank
)

GSEA_Progeria <- fgsea(Progeria_datasets, Rank)

# Save results in table format ----

write_tsv(GSEA_H, file = "Out/GSEA/GSEA_H.tsv")

write_tsv(
  GSEA_M2 |>
    filter(pathway %in% collapsedM2$mainPathways),
  file = "Out/GSEA/GSEA_M2_databases.tsv"
)

write_tsv(GSEA_M3, file = "Out/GSEA/GSEA_M3.tsv")

write_tsv(
  GSEA_M5 |>
    filter(pathway %in% collapsedM5$mainPathways),
  file = "Out/GSEA/GSEA_M5_GO.tsv"
)


# GSEA_all <- bind_rows(
#   list("Hallmarks" = GSEA_H,
#        "Curated_GS" = GSEA_M2 |>
#          filter(pathway %in% collapsedM2$mainPathways)),
#   .id = "Dataset"
# )

ggplot(
  GSEA_H,
  aes(
    x = NES,
    y = fct_reorder(pathway, NES),
    size = -log10(padj)
  )
) +
  geom_point() +
  theme_classic() +
  xlim(-3,NA) +
  labs(
    x = "Normalized enrichment score",
    y = "Pathways",
    fill = "Statistical\n significance"
  )

ggsave(filename = "Out/GSEA_H.pdf", width = 6, height = 2.5)

ggplot(
  GSEA_M2 |>
    filter(pathway %in% collapsedM2$mainPathways),
  aes(
    x = NES,
    y = fct_reorder(pathway, NES),
    size = -log10(padj)
  )
) +
  geom_point() +
  theme_classic() +
  labs(
    x = "Normalized enrichment score",
    y = "Pathways",
    fill = "Statistical\n significance"
  )

ggsave(filename = "Out/GSEA_M2.pdf", width = 10, height = 10)

ggplot(
  GSEA_M5 |>
    filter(pathway %in% collapsedM5$mainPathways),
  aes(
    x = NES,
    y = fct_reorder(pathway, NES),
    size = -log10(padj)
  )
) +
  geom_point() +
  theme_classic() +
  labs(
    x = "Normalized enrichment score",
    y = "Pathways",
    fill = "Statistical\n significance"
  )

ggsave(filename = "Out/GSEA_M5.pdf", width = 10, height = 15)

# ORA miR-29 UP targets ----

UP_targets <- gene_res |>
  filter(padj < 0.05 & log2FoldChange > 0) |>
  filter(isTarget) |>
  pull(gene_name)

allORA <- list(
  ORA_H = fora(mouseHallmarks, UP_targets, universe = unique(gene_res$gene_name)) |>
    rowwise() |>
    mutate(overlapGenes = paste(unlist(overlapGenes),
      collapse = ", "
    )) |>
    ungroup(),
  ORA_M2 = fora(mouseM2, UP_targets, universe = unique(gene_res$gene_name)) |>
    rowwise() |>
    mutate(overlapGenes = paste(unlist(overlapGenes),
      collapse = ", "
    )) |>
    ungroup(),
  ORA_M3 = fora(mouseM3, UP_targets, universe = unique(gene_res$gene_name)) |>
    rowwise() |>
    mutate(overlapGenes = paste(unlist(overlapGenes),
      collapse = ", "
    )) |>
    ungroup(),
  ORA_M5 = fora(mouseM5, UP_targets, universe = unique(gene_res$gene_name)) |>
    rowwise() |>
    mutate(overlapGenes = paste(unlist(overlapGenes),
      collapse = ", "
    )) |>
    ungroup()
)

allORA <- bind_rows(allORA, .id = "Database")

write_tsv(allORA, "Out/GSEA/ORA_targets.tsv")

# Remove environment
rm(list = ls())
