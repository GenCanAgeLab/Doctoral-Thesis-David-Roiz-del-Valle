## Load RDS file  ---------------------

gene_res <- readRDS("Data/results.rds")
gene_vsd <- readRDS("Data/vsd.rds")

# Select drug metabolism enzymes ----

enzymeFamilies <- c("Cyp", "Sult", "Ugt", "Ces", "Gst", "Slc", "Abc",
                    "Alas1","Nr1i3","Ppara","Pparg","Rxr","Nr1i2", "Acot",
                    "Por", "Fmo", "Paps")

xenoEnzymes <- lapply(enzymeFamilies,
                      function(x) {
                        gene_res |>
                          filter(str_detect(gene_name, x))
                      })

xenoEnzymesdf <- bind_rows(xenoEnzymes,
                           .id = "Enzyme_family") |>
  filter(gene_type == "protein_coding") |>
  filter(padj < 0.05) |>
  group_by(Enzyme_family) |>
  arrange(log2FoldChange,
          .by_group = TRUE) |>
  mutate(gene_name = factor(gene_name,
                            levels = gene_name)) |>
  ungroup()


ggplot(xenoEnzymesdf, aes(x = Enzyme_family,
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

ggsave("Out/drug_metabolism.pdf",
       width = 10,
       height = 30)

# Heatmap

design <- read_tsv("Design.tsv")

xenoEnzymes <- lapply(enzymeFamilies,
                      function(x) {
                        tmp <- gene_vsd[str_detect(rowData(gene_vsd)$gene_name, x)]
                        rownames(tmp) <- rowData(tmp)$gene_name
                        colnames(tmp) <- paste(design$ID,design$Genotype,sep = "_")
                        as.data.frame(assay(tmp))
                      })

xenoEnzymesdf <- bind_rows(xenoEnzymes)

dev.off()
pdf("Out/heatmap_drug_metabolism.pdf",
    height = 50)

pheatmap(xenoEnzymesdf, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, scale = "row")

dev.off()




#Remove environment
rm(list = ls())
