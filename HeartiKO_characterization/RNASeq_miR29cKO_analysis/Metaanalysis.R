# Source R scripts ------

source("R/init.R") # Load libraries and multicore parameters
source("R/mmus_annotation.R")
source("R/hsap_annotation.R")
source("R/DE_Bulk_analysis.R")
source("R/DE_CMspec_analysis.R")
source("R/TMS_correlation.R")

# Source python script ----

system2(
  "/home/droiva/.local/bin/micromamba",
  args = c("run", "-n", "scanpy", "Python/TMS_single_cell_expression.py")
)

# Import data  ---------------------

## Load RDS files  ---------------------

res_Bulk_MA <- readRDS(file = "Data/dds_shrinked_Bulk.rds") |>
  filter(gene_name != "Esr1") |>
  filter(gene_type == "protein_coding" & is.na(tag))

res_CMspec_MA <- readRDS(file = "Data/dds_shrinked_CMspec.rds") |>
  filter(gene_name != "Esr1") |>
  filter(gene_type == "protein_coding" & is.na(tag))

res_iPSCs_MA <- readRDS(file = "Data/dds_shrinked_iPSCs.rds") |>
  filter(gene_type == "protein_coding" & is.na(tag))

res_Bulk <- readRDS(file = "Data/Results_Bulk.rds")
res_CMspec <- readRDS(file = "Data/Results_CMspec.rds")
res_iPSCs <- readRDS(file = "Data/Results_iPSCs.rds")

vsd_Bulk <- readRDS(file = "Data/vsd_Bulk.rds")
vsd_CMspec <- readRDS(file = "Data/vsd_CMspec.rds")
vsd_iPSCs <- readRDS(file = "Data/vsd_iPSCs.rds")

# Sstats

totalGenes <- length(res_CMspec$gene_id)
res_CMspec |>
  filter(padj < 0.05) |>
  group_by(sign(log2FoldChange)) |>
  summarize(
    n = n(),
    perc = n / totalGenes * 100
  )

# MA plot  ---------------------

plotMA <- function(shrinkedRes) {
  shrinkedRes <- shrinkedRes |>
    mutate(signif = if_else(padj < 0.05, TRUE, FALSE, missing = FALSE))

  ggplot(
    shrinkedRes,
    aes(
      x = log2(baseMean),
      y = log2FoldChange,
      color = signif,
      label = gene_name
    )
  ) +
    geom_point(shape = 16) +
    geom_label_repel(
      data = shrinkedRes |>
        filter(padj < 0.05) |>
        slice_min(padj, n = 10) |>
        ungroup(),
      max.overlaps = Inf
    ) +
    scale_color_manual(values = c("grey", "black")) +
    labs(x = " log2 mean expression", y = "log2FC") +
    theme_bw()
}

pdf(file = "Out/MA_plots.pdf")

lapply(list(res_Bulk_MA, res_CMspec_MA, res_iPSCs_MA), plotMA)

dev.off()

pdf(file = "Out/MA_plots_print.pdf", width = 5, height = 3)

lapply(list(res_Bulk_MA, res_CMspec_MA, res_iPSCs_MA), plotMA)

dev.off()

# GSEA  ---------------------------

## Load GMTs  ---------------------

mouseHallmarks <- gmtPathways(
  "/data/genomes/MsigDB/mh.all.v2023.2.Mm.symbols.gmt"
)
humanHallmarks <- gmtPathways(
  "/data/genomes/MsigDB/h.all.v2023.2.Hs.symbols.gmt"
)
mouseWP <- gmtPathways(
  "/data/genomes/MsigDB/m2.cp.wikipathways.v2024.1.Mm.symbols.gmt"
)

## Rank  ---------------------

Rank_Bulk <- res_Bulk |>
  filter(gene_name != "Esr1") |> #Delete Estrogen receptor, as it is an artifact
  drop_na(pvalue, log2FoldChange, gene_name) |> #  CollapsePathways fails if NA gene names are present
  arrange(desc(stat)) |>
  distinct(gene_name, .keep_all = T) |>
  pull(stat, name = gene_name)

Rank_CMspec <- res_CMspec |>
  filter(gene_name != "Esr1") |> #Delete Estrogen receptor, as it is an artifact
  drop_na(pvalue, log2FoldChange, gene_name) |> #  CollapsePathways fails if NA gene names are present
  arrange(desc(stat)) |>
  distinct(gene_name, .keep_all = T) |>
  pull(stat, name = gene_name)

Rank_iPSCs <- res_iPSCs |>
  drop_na(pvalue, log2FoldChange, gene_name) |> #  CollapsePathways fails if NA gene names are present
  arrange(desc(stat)) |>
  distinct(gene_name, .keep_all = T) |>
  pull(stat, name = gene_name)

set.seed(42) # Fix seed!

GSEA_H_Bulk <- fgsea(mouseHallmarks, Rank_Bulk, nproc = 1) |>
  filter(padj < 0.05) |>
  mutate(pathway = str_remove(pathway, "HALLMARK_")) |> 
  mutate(pathway = str_replace_all(pathway, "_"," "))

GSEA_H_CMspec <- fgsea(mouseHallmarks, Rank_CMspec, nproc = 1) |>
  filter(padj < 0.05) |>
  mutate(pathway = str_remove(pathway, "HALLMARK_")) |> 
  mutate(pathway = str_replace_all(pathway, "_"," "))

# GSEA_H_iPSCs <- fgsea(humanHallmarks, Rank_iPSCs, nproc = 1) |>
#   filter(padj < 0.05) |>
#   mutate(pathway = str_remove(pathway, "HALLMARK_"))

GSEA_WP_Bulk <- fgsea(mouseWP, Rank_Bulk, nproc = 1) |>
  filter(padj < 0.05) |>
  mutate(pathway = str_remove(pathway, "WP_")) |>
  filter(str_detect(
    pathway,
    "HYPERTROPHY|ELECTRON|RIBOSOMAL|OXIDATIVE_PHOSPHORYLATION"
  )) |> 
  mutate(pathway = str_replace_all(pathway, "_"," "))

GSEA_WP_CMspec <- fgsea(mouseWP, Rank_CMspec, nproc = 1) |>
  filter(padj < 0.05) |>
  mutate(pathway = str_remove(pathway, "WP_")) |>
  filter(str_detect(
    pathway,
    "HYPERTROPHY|ELECTRON|RIBOSOMAL|OXIDATIVE_PHOSPHORYLATION"
  )) |> 
  mutate(pathway = str_replace_all(pathway, "_"," "))

GSEA_H <- bind_rows(
  list(
    "Bulk" = GSEA_H_Bulk,
    "CM specific" = GSEA_H_CMspec
  ),
  .id = "Model"
)

GSEA_WP <- bind_rows(
  list(
    "Bulk" = GSEA_WP_Bulk,
    "CM specific" = GSEA_WP_CMspec
  ),
  .id = "Model"
)

GSEA <- bind_rows(
  list(
  "MSigDB" = GSEA_H, 
  "WP" = GSEA_WP),
  .id = "Database") |> 
  # Only first 30 leading edge genes
  mutate(leadingEdge = map_chr(leadingEdge, ~ paste(head(.x, 30), collapse = ", ")))

write_tsv(GSEA,"Out/GSEA_df.tsv")

# ggplot(
#   GSEA_H,
#   aes(
#     x = NES,
#     y = fct_reorder(pathway, NES),
#     size = -log10(padj),
#     color = -log10(padj)
#   )
# ) +
#   geom_point() +
#   labs(
#     x = "Normalized enrichment score",
#     y = "Pathways",
#     fill = "Statistical\n significance"
#   ) +
#   geom_vline(xintercept = 0, linewidth = 0.4, color = "grey") +
#   xlim(-3.5, 3) +
#   theme_classic() +
#   facet_wrap(vars(Model)) +
#   scale_colour_steps2()
# 
# ggsave(filename = "Out/GSEA_print.pdf", width = 7, height = 5)
# 
# ggplot(
#   GSEA_WP,
#   aes(
#     x = NES,
#     y = fct_reorder(pathway, NES),
#     size = -log10(padj),
#     color = -log10(padj)
#   )
# ) +
#   geom_point() +
#   labs(
#     x = "Normalized enrichment score",
#     y = "Pathways",
#     fill = "Statistical\n significance"
#   ) +
#   geom_vline(xintercept = 0, linewidth = 0.4, color = "grey") +
#   theme_classic() +
#   xlim(-3.5, 2) +
#   facet_wrap(vars(Model)) +
#   scale_colour_steps2()

ggsave(filename = "Out/GSEA_WP_print.pdf", width = 7, height = 3)

ggplot(
  GSEA,
  aes(
    x = NES,
    y = fct_reorder(pathway, NES),
    size = -log10(padj)
  )
) +
  geom_point() +
  labs(
    x = "Normalized enrichment score",
    y = "Pathways",
    fill = "Statistical\n significance"
  ) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey") +
  theme_classic() +
  xlim(-3.5, 2) +
  facet_grid(rows = vars(Database),
             cols = vars(Model),
             scales = "free_y",
             space = "free_y")

ggsave(filename = "Out/GSEA_print.pdf", width = 7, height = 6)


########################################################################################

# Merge results file

## Merge res  ---------------------

resAll <- res_Bulk |>
  left_join(
    res_CMspec |> 
      select(-miRDB_miRNA),
    suffix = c("_Bulk", "_CMspec"),
    by = join_by(
      gene_id,
      gene_name,
      gene_type,
      tag,
      TargetScan_nseed,
      miRDB_score,
      microT_miRNA,
      microT_score,
      isTarget
    )
  )
  # left_join(
  #   res_iPSCs |>
  #     dplyr::select(-c(gene_id, gene_name, gene_type, tag)),
  #   suffix = c("_Bulk", "_iPSCs"),
  #   by = join_by(gene_name == Mouse_gene_name)
  # )

write_delim(resAll, file = "Out/results.tsv", delim = "\t")

# Individual files

write_tsv(res_Bulk, file = "Out/results_Bulk.tsv")
write_tsv(res_CMspec, file = "Out/results_CMspec.tsv")

# Top genes heatmap  ---------------------

## Gene list  ---------------------

TopGenes <- resAll |>
  filter(padj_Bulk < 0.05 | padj_CMspec < 0.05) |>
  filter(abs(log2FoldChange_Bulk) > 2 | abs(log2FoldChange_CMspec) > 2) |>
  pull(gene_id) |>
  unique()

## createHeatmaps function  ---------------------

createHeatmap <- function(vsd, geneList, title) {
  # Vsd object from genelist

  vsdHeatmap <- assay(vsd) |>
    as_tibble(rownames = "gene_id") |>
    mutate(
      gene_name = rowData(vsd)$gene_name,
      gene_type = rowData(vsd)$gene_type,
      tag = rowData(vsd)$tag
    ) |>
    filter(gene_name != "Esr1") |> #Delete Estrogen receptor, as it is an artifact
    filter(gene_type == "protein_coding" & is.na(tag)) |>
    filter(gene_id %in% TopGenes) |>
    dplyr::select(-c(gene_id, gene_type, tag)) |>
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

## Plot heatmaps --------------------

dev.off()
pdf(file = "Out/Top_genes_heatmaps.pdf")

createHeatmap(vsd_Bulk, TopGenes, "Top genes Bulk RNASeq")
createHeatmap(vsd_CMspec, TopGenes, "Top genes CMspec RNASeq")

dev.off()

# Venn diagrams ----

## Creating gene subsets ----

UP_Bulk <- resAll |>
  filter(padj_Bulk < 0.05 & log2FoldChange_Bulk > 1) |>
  pull(gene_name)

DOWN_Bulk <- resAll |>
  filter(padj_Bulk < 0.05 & log2FoldChange_Bulk < -1) |>
  pull(gene_name)

UP_CMspec <- resAll |>
  filter(padj_CMspec < 0.05 & log2FoldChange_CMspec > 1) |>
  pull(gene_name)

DOWN_CMspec <- resAll |>
  filter(padj_CMspec < 0.05 & log2FoldChange_CMspec < -1) |>
  pull(gene_name)

Venn_plot <- plotVenn(
  list(
    "Upregulated in Bulk" = UP_Bulk,
    "Upregulated in CM specific" = UP_CMspec,
    "Downregulated in Bulk" = DOWN_Bulk,
    "Downregulated in CM specific" = DOWN_CMspec
  ),
  outFile = "Out/Venn.svg"
)

DOWN_Bulk_only <- setdiff(DOWN_Bulk, DOWN_CMspec)
UP_Bulk_only <- setdiff(UP_Bulk, UP_CMspec)

UP_both <- intersect(UP_Bulk, UP_CMspec)
DOWN_both <- intersect(DOWN_Bulk, DOWN_CMspec)

# miR-29 target list ----

UP_targets <- resAll |> 
  filter(isTarget) |> 
  filter(padj_CMspec < 0.05 | padj_Bulk < 0.05) |> 
  filter(log2FoldChange_CMspec > 0.5 | log2FoldChange_Bulk > 0.5)

write_tsv(UP_targets, "Out/miR_29_targets.tsv")

