source("R/init.R")

# Read circRNA data ----

circ_res <- read_tsv("Out/LAKI/circRNA/circRNA_results_annotated.tsv")

consensusCirc <- read_tsv("Out/LAKI/circRNA/circRNA_results_consensus.tsv")

consistent_circDE <- read_tsv("Out/LAKI/circRNA/circRNA_results_consistent.tsv")

# Load lineal transcript data ----

geneLAKIList <- readRDS("Out/LAKI/DESeq_results_sexDiv_LAKI.rds")

geneLAKI <- lapply(geneLAKIList, function(x) {
  x <- as_tibble(x, rownames = "gene_id")
})

geneLAKIdf <- bind_rows(geneLAKI, .id = "Sex_Tissue") |>
  separate_wider_delim(Sex_Tissue, names = c("Sex", "Tissue"), delim = " ") |>
  dplyr::select(Sex, Tissue, gene_id, log2FoldChange, padj, gene_type, gene_name)

# Merge info ----

# circ_res_lineal <- circ_res |>
#   left_join(geneLAKIdf,
#             by = join_by(Tissue, Sex, gene_id, gene_name, gene_type),
#             suffix = c("","_lineal"))

consensusCirc_lineal <- consensusCirc |>
  left_join(
    geneLAKIdf |>
      dplyr::rename(
        LFC_lineal = log2FoldChange,
        padj_lineal = padj
      ),
    by = join_by(Tissue, Sex, gene_id, gene_name, gene_type)
  ) |>
  mutate(changeLineal = case_when(
    # Use 0.1 as cutoff for lineal genes
    padj_lineal < 0.1 & LFC_lineal > 0 ~ "UP",
    padj_lineal < 0.1 & LFC_lineal < 0 ~ "DOWN",
    .default = "NS"
  )) |>
  left_join(
    consistent_circDE |>
      dplyr::select(circID, Sex, NtissuesDE, Tissues),
    by = join_by(Sex, circID)
  )


# Correlation log2FC circ_lineal ----

# LM fit

circ_lineal <- lm(LFC_lineal ~ medianLFC,
  data = consensusCirc_lineal
)

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

# kidneyCirc <- candidateCirc |>
#   filter(Sex == "Male" & Tissue == "Kidney")

# geneKidney <- kidneyCirc |>
#   group_by(gene_name) |>
#   summarize(Gene_count = n())

# Save results ----

# write_tsv(circ_res_lineal, file = "Out/LAKI/circRNA/circRNA_results_annotated_lineal.tsv")

write_tsv(consensusCirc_lineal, file = "Out/LAKI/circRNA/circRNA_results_consensus_lineal.tsv")
write_tsv(candidateCirc, file = "Out/LAKI/circRNA/circRNA_candidates.tsv")

# Crosscheck with ORF info ----

library(readxl)
circRNA_ORF <- read_excel("Data/resources/Fan_etal_2022_circRNA_ORF.xlsx",
  sheet = "circPeptide"
) |>
  distinct(Proteins, `gene name`)
Human_Mouse <- read_tsv("/data/genomes/MsigDB/Human_Gene_Symbol_with_Remapping_Mouse_Orthologs_MSigDB.v2024.1.Mm.tsv") |>
  dplyr::select(
    Human = `Probe Set ID`,
    Mouse = `Gene Symbol`
  )

circRNA_ORF <- left_join(circRNA_ORF, Human_Mouse,
  by = join_by(`gene name` == Human)
)

candidateCirc_ORF <- left_join(candidateCirc, circRNA_ORF,
  by = join_by(gene_name == Mouse),
  relationship = "many-to-many"
)

write_tsv(candidateCirc_ORF, file = "Out/LAKI/circRNA/circRNA_candidates_Fan_etal_ORF.tsv")
