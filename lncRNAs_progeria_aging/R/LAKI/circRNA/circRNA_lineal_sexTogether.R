
source("R/init.R")

# Read circRNA data ----

circ_res <- read_tsv("Out/LAKI/circRNA/circRNA_results_annotated_sexTogether.tsv")

consensusCirc <- read_tsv("Out/LAKI/circRNA/circRNA_results_consensus_sexTogether.tsv")

consistent_circDE <- read_tsv("Out/LAKI/circRNA/circRNA_results_consistent_sexTogether.tsv")

# Load lineal transcript data ----

geneLAKIList <- readRDS("Out/LAKI/DESeq_results_LAKI.rds")

geneLAKI <- lapply(geneLAKIList, function(x) {
  x <- as_tibble(x, rownames = "gene_id")
})

geneLAKIdf <- bind_rows(geneLAKI, .id = "Tissue") |>
  dplyr::select(Tissue, gene_id, log2FoldChange, padj, gene_type, gene_name)

# Merge info ----

# circ_res_lineal <- circ_res |> 
#   left_join(geneLAKIdf,
#             by = join_by(Tissue, gene_id, gene_name, gene_type),
#             suffix = c("","_lineal"))

consensusCirc_lineal <- consensusCirc |> 
  left_join(geneLAKIdf |> 
              dplyr::rename(LFC_lineal = log2FoldChange,
                            padj_lineal = padj),
            by = join_by(Tissue, gene_id, gene_name, gene_type)) |> 
  mutate(changeLineal = case_when(
    # Use 0.1 as cutoff for lineal genes
    padj_lineal < 0.1 & LFC_lineal > 0 ~ "UP",
    padj_lineal < 0.1 & LFC_lineal < 0 ~ "DOWN",
    .default = "NS"
  )) |> 
  left_join(consistent_circDE |> 
              dplyr::select(circID, NtissuesDE, Tissues),
            by = join_by(circID))


# Correlation log2FC circ_lineal ----

# LM fit

circ_lineal <- lm(LFC_lineal ~ medianLFC,
   data = consensusCirc_lineal)

summary(circ_lineal)

# circ_lineal_all <- lm(log2FoldChange_lineal ~ log2FoldChange,
#                   data = circ_res_lineal)
# 
# summary(circ_lineal_all)


# Plot

ggplot(consensusCirc_lineal, aes(x = medianLFC, y = LFC_lineal)) +
  geom_point() +
  geom_smooth(method = "lm")




# Candidate search ----

# circRNAs in NS or opposed change lineal genes
# Also, info about changes in other tissues for each candidate

candidateCirc <- consensusCirc_lineal |> 
  filter(Change != "NS") |> 
  filter(Change != changeLineal)

geneCountAll <- candidateCirc |> 
  group_by(gene_name) |> 
  summarize(Gene_count = n())


# Kidney 

kidneyCirc <- candidateCirc |> 
  filter(Sex == "Male" & Tissue == "Kidney")

geneKidney <- kidneyCirc |> 
  group_by(gene_name) |> 
  summarize(Gene_count = n())

# Save results ----

# write_tsv(circ_res_lineal, file = "Out/LAKI/circRNA/circRNA_results_annotated_lineal.tsv")

write_tsv(consensusCirc_lineal, file = "Out/LAKI/circRNA/circRNA_results_consensus_lineal_sexTogether.tsv")
write_tsv(candidateCirc, file = "Out/LAKI/circRNA/circRNA_candidates_sexTogether.tsv")




