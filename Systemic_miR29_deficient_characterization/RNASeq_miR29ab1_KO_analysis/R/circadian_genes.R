## Load RDS file  ---------------------

res_Liver <- readRDS(file = "data/R_data/liver_results.rds")

# Select drug metabolism enzymes ----

circGenes <- c("Ror",
               "Nr1d",
               "Arntl",
               "Nfil3",
               "Dbp", "Tef", "Hlf",
               "Clock",
               "Per", "Cry")

circGenes <- lapply(circGenes,
                      function(x) {
                        res_Liver |>
                          filter(str_detect(gene_name, x))
                      })

circGenesdf <- bind_rows(circGenes) |>
  filter(gene_type == "protein_coding" & is.na(tag)) |>
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

ggsave("Out/circadian_genes.pdf",
       width = 10,
       height = 10)

#Remove environment
rm(list = ls())