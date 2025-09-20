## Load RDS file  ---------------------

gene_res <- readRDS("Data/results.rds")
gene_vsd <- readRDS("Data/vsd.rds")

# Select drug metabolism enzymes ----

circGenes <- c("Ror",
               "Nr1d",
               "Arntl",
               "Nfil3",
               "Dbp", "Tef", "Hlf",
               "Clock",
               "Per", "Cry")

circGeneData <- lapply(circGenes,
                    function(x) {
                      gene_res |>
                        filter(str_detect(gene_name, x))
                    })

circGenesdf <- bind_rows(circGeneData) |>
  arrange(log2FoldChange) |>
  mutate(gene_name = factor(gene_name,
                            levels = gene_name)) |>
  ungroup()


ggplot(circGenesdf, aes(x = gene_name,
                        y = gene_name,
                        fill = log2FoldChange)) +
  geom_tile() +
  geom_text(aes(label = gene_name),
            color = "black",
            size = 2) +
  labs(
    x = NULL, 
    y = NULL, 
    fill = "log2FoldChange") +
  theme_void() +
  scale_fill_gradient2(
    low = "red",
    mid = "white",
    high = "green",
    midpoint = 0)

ggsave("Out/circadian_genes_all.pdf",
       width = 10,
       height = 10)

circGenesdf <- bind_rows(circGeneData) |>
  filter(padj < 0.05) |>
  arrange(log2FoldChange) |>
  mutate(gene_name = factor(gene_name,
                            levels = gene_name)) |>
  ungroup()


ggplot(circGenesdf, aes(x = gene_name,
                        y = gene_name,
                        fill = log2FoldChange)) +
  geom_tile() +
  geom_text(aes(label = gene_name),
            color = "black",
            size = 2) +
  labs(
    x = NULL, 
    y = NULL, 
    fill = "log2FoldChange") +
  theme_void() +
  scale_fill_gradient2(
    low = "red",
    mid = "white",
    high = "green",
    midpoint = 0)

ggsave("Out/circadian_genes_sign.pdf",
       width = 10,
       height = 10)

# Heatmap

design <- read_tsv("Design.tsv")

circGeneVsd <- lapply(circGenes,
                    function(x) {
                      tmp <- gene_vsd[str_detect(rowData(gene_vsd)$gene_name, x)]
                      rownames(tmp) <- rowData(tmp)$gene_name
                      colnames(tmp) <- paste(design$ID,design$Genotype,sep = "_")
                      as.data.frame(assay(tmp))
                    })

circGeneVsddf <- bind_rows(circGeneVsd)

dev.off()
pdf("Out/heatmap_circadian.pdf")

pheatmap(circGeneVsddf, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, scale = "row")

dev.off()

#Remove environment
rm(list = ls())