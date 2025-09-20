source("R/init.R")

gene_res <- readRDS("Data/results.rds")
shrinkedRes <- readRDS("Data/res_shrinked.rds")

# Summary of DE ----

totalGenes <- length(gene_res$gene_id)

DE_stats <- gene_res |>
  filter(padj < 0.05) |>
  mutate(Change = if_else(log2FoldChange > 0,
    "UP",
    "DOWN"
  )) |>
  group_by(Change) |>
  summarise(
    DE = n(),
    Perc = n() / totalGenes * 100
  )


# MA plot ----

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
    labs(x = " log2 mean reads", y = "log2FC") +
    theme_classic()
}

plotMA(shrinkedRes)


ggsave(file = "Out/MA_plot.pdf", width = 5, height = 3)

# miR-29 target list ----

miR_29_targets <- gene_res |> 
  filter(isTarget) |> 
  filter(padj < 0.05) |> 
  filter(log2FoldChange > 0.5)

write_tsv(miR_29_targets, file = "Out/miR_29_targets.tsv")

# Remove environment
rm(list = ls())
