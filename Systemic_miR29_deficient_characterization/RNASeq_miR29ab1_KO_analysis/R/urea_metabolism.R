## Load RDS file  ---------------------

res_Liver <- readRDS(file = "data/R_data/liver_results.rds")

# Select drug metabolism enzymes ----

ureaGenes <- c("Cps1",
               "Otc",
               "Ass1",
               "Asl",
               "Arg1")

ureaGenes <- lapply(ureaGenes,
                      function(x) {
                        res_Liver |>
                          filter(str_detect(gene_name, x))
                      })

ureaGenesdf <- bind_rows(ureaGenes) |>
  filter(gene_type == "protein_coding") |>
  # filter(padj < 0.05) |>
  arrange(log2FoldChange,
          .by_group = TRUE) |>
  mutate(gene_name = factor(gene_name,
                            levels = gene_name)) |>
  ungroup()


ggplot(ureaGenesdf, aes(x = gene_name,
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

ggsave("Out/urea_metabolism.pdf",
       width = 5,
       height = 5)

#Remove environment
rm(list = ls())
