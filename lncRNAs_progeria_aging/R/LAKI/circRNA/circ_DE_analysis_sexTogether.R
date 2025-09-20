
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)

# Read files

circSexTisDds <- readRDS("Out/LAKI/circRNA/dds_circRNA_sexTogether.rds")
resultsList <- readRDS("Out/LAKI/circRNA/results_circRNA_sexTogether.rds")

resAll <- read_tsv("Out/LAKI/circRNA/circRNA_results_all_sexTogether.tsv")

# All circ annotation ----

circ_annotation <- read_tsv("/home/jandrulas/20231017_RNASeq_LAKI/circRNA_analysis/Results/circRNA_ID_across_tools.tsv")
circ_gene_annotation <- read_tsv("/home/jandrulas/20231017_RNASeq_LAKI/circRNA_analysis/Results/All_circRNAs_annotated_cross_with_circDBs.tsv")

geneMetadata <- read_tsv("Data/geneMetadata.tsv")

circ_annotation <- circ_annotation |> 
  pivot_longer(cols = c(CIRI2, Terrace, circExplorer2, find_circ, starchip),
               names_to = "Tool",
               values_to = "circID_tool") |> 
  left_join(circ_gene_annotation,
            by = join_by(circID, seqnames, start, end, width)) |> 
  left_join(geneMetadata |> 
              dplyr::select(gene_id, gene_type, gene_name),
            by = join_by(gene_id))

# Merge annotation with results ----

annotated_res <- resAll |> 
  left_join(circ_annotation,
            by = join_by(Tool, circRNA_ID == circID_tool))

# DE analysis per tool per tissue/sex ----

circ_DE_tis <- annotated_res |> 
  filter(padj < 0.05) |> 
  group_by(Tissue, Tool) |> 
  summarize(
    UP = sum(log2FoldChange > 0, na.rm = TRUE),
    DOWN = sum(log2FoldChange < 0, na.rm = TRUE)
  ) |> 
  pivot_longer(cols = c(UP, DOWN),
               names_to = "Change",
               values_to = "DE_genes") |> 
  mutate(Change = factor(Change, levels = c("UP", "DOWN")))

ggplot(circ_DE_tis, aes(x = Tissue, y = DE_genes, fill = Tool)) +
  geom_bar(position = "dodge", stat = "identity") +
  facet_wrap(vars(Change)) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic()

ggsave("Out/LAKI/circRNA/DE_tissue_sex_sexTogether.pdf",
       width = 8,
       height = 5)

# Consensus circRNAs ----

consensus_circDE <- annotated_res |> 
  drop_na(padj) |> 
  mutate(Change = factor(sign(log2FoldChange),
                         levels = c(1,-1),
                         labels = c("UP", "DOWN"))) |> 
  filter(padj < 0.05) |> 
  group_by(Tissue,
           circID, Cons_coords, type, gene_id, gene_type, gene_name,
           Change, circAtlas, circBase, circPedia) |> 
  summarize(
    NtoolsDE = n(),
    medianBMean = median(baseMean),
    medianLFC = median(log2FoldChange),
    medianPadj = median(padj),
    medianStat = median(stat)
  ) |> 
  filter(NtoolsDE > 2) # DE detected by more than 2 tools

consensus_circDE_tis <- consensus_circDE |> 
  filter(NtoolsDE > 2) |> # DE detected by more than 2 tools
  group_by(Tissue, Change, type) |> 
  summarize(
    DE_genes = n()
  )

ggplot(consensus_circDE_tis, aes(x = Tissue, y = DE_genes, fill = Change)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic()

ggsave("Out/LAKI/circRNA/DE_tissue_sex_consensus_sexTogether.pdf",
       width = 6,
       height = 3)

ggplot(consensus_circDE_tis, aes(x = Tissue, y = DE_genes, fill = type)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(vars(Change)) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic()

ggsave("Out/LAKI/circRNA/DE_tissue_sex_consensus_type_sexTogether.pdf",
       width = 8,
       height = 4)

# Consistently-changed circRNAs ----

consistent_circDE <- consensus_circDE |> 
  filter(NtoolsDE > 2) |> # DE detected by more than 2 tools
  group_by(circID, Cons_coords, type, gene_id, gene_type, gene_name,
           Change, circAtlas, circBase, circPedia) |> 
  summarize(
    NtissuesDE = n_distinct(Tissue),
    Tissues = paste(unique(Tissue),
                    collapse = ", "),
    medianLFC = median(medianLFC),
    medianPadj = median(medianPadj),
    medianBMean = median(medianBMean)
  )


# Save data ----

write_tsv(annotated_res,
          file = "Out/LAKI/circRNA/circRNA_results_annotated_sexTogether.tsv")

write_tsv(consensus_circDE,
          file = "Out/LAKI/circRNA/circRNA_results_consensus_sexTogether.tsv")

write_tsv(consistent_circDE,
          file = "Out/LAKI/circRNA/circRNA_results_consistent_sexTogether.tsv")




