# Source scripts ----

source("R/init.R")
source("R/annotation.R")
source("R/data_import.R")
source("R/Liver_analysis.R")
source("R/Heart_analysis.R")
source("R/drug_metabolism_liver.R")

# Import data ----

## Load RDS files  ---------------------

res_Liver_MA <- readRDS(file = "data/R_data/liver_dds_shrinked.rds")
# filter(gene_type == "protein_coding" & is.na(tag))

res_Heart_MA <- readRDS(file = "data/R_data/heart_dds_shrinked.rds")
# filter(gene_type == "protein_coding" & is.na(tag))

res_Liver <- readRDS(file = "data/R_data/liver_results.rds")
res_Heart <- readRDS(file = "data/R_data/heart_results.rds")

vsd_Liver <- readRDS(file = "data/R_data/liver_vsd.rds")
vsd_Heart <- readRDS(file = "data/R_data/heart_vsd.rds")

# RNASeq stats

totalGenes <- length(res_Liver$gene_id)
res_Liver |>
  filter(padj < 0.05) |>
  group_by(sign(log2FoldChange)) |>
  summarize(
    n = n(),
    perc = n / totalGenes * 100
  )

# MA plot  ---------------------

plotMA <- function(shrinkedRes) {
  title <- deparse(substitute(shrinkedRes))
  ggplot(
    shrinkedRes,
    aes(
      x = log2(baseMean),
      y = log2FoldChange,
      color = padj < 0.05,
      label = gene_name
    )
  ) +
    geom_point(shape = 16, alpha = 0.7, size = 0.8) +
    geom_label_repel(
      data = shrinkedRes |>
        filter(gene_type == "protein_coding" & is.na(tag)) |>
        filter(padj < 0.05) |>
        slice_max(abs(log2FoldChange), n = 10) |>
        ungroup(),
      max.overlaps = 200,
      point.size = NA,
      size = 3
    ) +
    scale_color_manual(values = c("grey", "black")) +
    labs(x = " log2 mean reads", y = "log2FC") +
    theme_bw()
}

pdf(file = "Out/MA_plots.pdf", width = 5, height = 3)

lapply(list(res_Liver_MA, res_Heart_MA), plotMA)

dev.off()

# GSEA  ---------------------------

## Load GMTs  ---------------------

mouseHallmarks <- gmtPathways(
  "/data/genomes/MsigDB/mh.all.v2023.2.Mm.symbols.gmt"
)
mouseWP <- gmtPathways(
  "/data/genomes/MsigDB/m2.cp.wikipathways.v2024.1.Mm.symbols.gmt"
)

mouseGO <- gmtPathways(
  "/data/genomes/MsigDB/c5.go.v2023.2.Hs.symbols.gmt"
)

## Rank  ---------------------

Rank_Liver <- res_Liver |>
  drop_na(pvalue, log2FoldChange, gene_name) |> #  CollapsePathways fails if NA gene names are present
  arrange(desc(stat)) |>
  distinct(gene_name, .keep_all = T) |>
  pull(stat, name = gene_name)

Rank_Heart <- res_Heart |>
  drop_na(pvalue, log2FoldChange, gene_name) |> #  CollapsePathways fails if NA gene names are present
  arrange(desc(stat)) |>
  distinct(gene_name, .keep_all = T) |>
  pull(stat, name = gene_name)

set.seed(42) # Fix seed!

GSEA_H_Liver <- fgsea(mouseHallmarks, Rank_Liver, nproc = 1) |>
  filter(padj < 0.05) |>
  mutate(pathway = str_remove(pathway, "HALLMARK_"))

GSEA_H_Heart <- fgsea(mouseHallmarks, Rank_Heart, nproc = 1) |>
  filter(padj < 0.05) |>
  mutate(pathway = str_remove(pathway, "HALLMARK_"))

GSEA_WP_Liver <- fgsea(mouseWP, Rank_Liver, nproc = 1) |>
  filter(padj < 0.05) |>
  mutate(pathway = str_remove(pathway, "WP_"))

GSEA_GO_Liver <- fgsea(mouseGO, Rank_Liver, nproc = 1) |>
  filter(padj < 0.05)

GSEA_H <- bind_rows(
  list("Liver" = GSEA_H_Liver, "Heart" = GSEA_H_Heart),
  .id = "Tissue"
)

ggplot(
  GSEA_H,
  aes(
    x = NES,
    y = fct_reorder(pathway, NES),
    size = -log10(padj)
  )
) +
  geom_point() +
  facet_wrap(vars(Tissue)) +
  theme_classic() +
  labs(
    x = "Normalized enrichment score",
    y = "Pathways",
    fill = "Statistical\n significance"
  )

ggsave(filename = "Out/GSEA.pdf", width = 10, height = 7)

# Liver-only

ggplot(
  GSEA_H_Liver,
  aes(
    x = NES,
    y = fct_reorder(pathway, NES),
    size = -log10(padj)
  )
) +
  geom_point() +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey") +
  xlim(-3, 3) +
  theme_classic() +
  labs(
    x = "Normalized enrichment score",
    y = "Pathways",
    fill = "Statistical\n significance"
  )

ggsave(filename = "Out/GSEA_liver.pdf", width = 5, height = 4)

ggplot(
  GSEA_WP_Liver,
  aes(
    x = NES,
    y = fct_reorder(pathway, NES),
    size = -log10(padj)
  )
) +
  geom_point() +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey") +
  xlim(-3, 3) +
  theme_classic() +
  labs(
    x = "Normalized enrichment score",
    y = "Pathways",
    fill = "Statistical\n significance"
  )

ggsave(filename = "Out/GSEA_WP_liver.pdf", width = 7, height = 4)

########################################################################################

# Merge results file

## Merge res  ---------------------

resAll <- res_Liver |>
  left_join(
    res_Heart,
    suffix = c("_Liver", "_Heart"),
    by = join_by(
      gene_id,
      gene_name,
      gene_type,
      tag,
      TargetScan_nseed,
      miRDB_score,
      microT_miRNA,
      microT_score,
      miRDB_miRNA,
      isTarget
    )
  )

write_delim(resAll, file = "Out/results.tsv", delim = "\t")

# Top genes heatmap  ---------------------

## Gene list  ---------------------

TopGenes_Liver <- resAll |>
  filter(gene_type == "protein_coding" & is.na(tag)) |>
  filter(padj_Liver < 0.05) |>
  slice_max(abs(log2FoldChange_Liver), n = 30) |>
  pull(gene_id) |>
  unique()

TopGenes_Heart <- resAll |>
  filter(gene_type == "protein_coding" & is.na(tag)) |>
  filter(padj_Heart < 0.05) |>
  slice_max(abs(log2FoldChange_Heart), n = 30) |>
  pull(gene_id) |>
  unique()

createHeatmap <- function(vsd, geneList, title) {
  # Vsd object from genelist

  vsdHeatmap <- assay(vsd) |>
    as_tibble(rownames = "gene_id") |>
    mutate(gene_name = rowData(vsd)$gene_name) |>
    filter(gene_id %in% geneList) |>
    dplyr::select(-c(gene_id)) |>
    tibble::column_to_rownames(var = "gene_name") |>
    drop_na()

  # Order samples and annotation col

  df <- colData(vsd) |>
    as.data.frame() |>
    dplyr::select(Genotype)

  orderSamples <- colData(vsd) |>
    as.data.frame() |>
    arrange(Genotype) |>
    pull(names)

  vsdHeatmap <- vsdHeatmap[, match(orderSamples, colnames(vsdHeatmap))]

  # Heatmap

  pheatmap(
    vsdHeatmap,
    cluster_rows = T,
    show_rownames = T,
    cluster_cols = F,
    annotation_col = df,
    color = colorRampPalette(c("blue", "white", "red"))(12),
    scale = "row",
    main = title
  )
}

## Plot heatmaps ----------

dev.off()

pdf(file = "Out/Top_genes_heatmaps.pdf")

createHeatmap(vsd_Liver, TopGenes_Liver, "Top genes Liver")
createHeatmap(vsd_Heart, TopGenes_Heart, "Top genes Heart")

dev.off()
