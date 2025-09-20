## Load RDS file  ---------------------

res_Liver <- readRDS(file = "data/R_data/liver_results.rds")

# Select drug metabolism enzymes ----

enzymeFamilies <- c("Cyp", "Sult", "Ugt", "Ces", "Gst", "Slc", "Abc",
                    "Alas1","Nr1i3","Ppara","Pparg","Rxr","Nr1i2", "Acot",
                    "Por", "Fmo", "Paps")

xenoEnzymes <- lapply(enzymeFamilies,
                      function(x) {
                        res_Liver |>
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

#Remove environment
rm(list = ls())
